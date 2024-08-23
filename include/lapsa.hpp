#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <queue>
#include <thread>
#include <vector>

#include "aviize.hpp"
#include "iestade.hpp"
#include "rododendrs.hpp"

namespace lapsa {

// tools
struct Settings {
    size_t n_states            = 1000000;
    double init_p_acceptance   = 0.99;
    size_t init_t_log_len      = 100;
    double cooling_rate        = 0.95;
    size_t cooling_round_len   = 1;
    size_t e_sma_fast_len      = 50;
    size_t e_sma_slow_len      = 200;
    size_t e_sma_update_period = 100;

    size_t progress_update_period = 1;
    std::string log_filename      = "log.csv";
    std::string stats_filename    = "stats.txt";

    size_t n_reports = 0;

    Settings() {}

    // clang-format off
    Settings(const std::string& config_filepath,
             const std::string& key_path_prefix) :
        n_states            (iestade::size_t_from_json(config_filepath, key_path_prefix + "/n_states")),
        init_p_acceptance   (iestade::double_from_json(config_filepath, key_path_prefix + "/init_p_acceptance")),
        init_t_log_len      (iestade::size_t_from_json(config_filepath, key_path_prefix + "/init_t_log_len")),
        cooling_rate        (iestade::double_from_json(config_filepath, key_path_prefix + "/cooling_rate")),
        cooling_round_len   (iestade::size_t_from_json(config_filepath, key_path_prefix + "/cooling_round_len")),
        e_sma_fast_len      (iestade::size_t_from_json(config_filepath, key_path_prefix + "/e_sma_fast_len")),
        e_sma_slow_len      (iestade::size_t_from_json(config_filepath, key_path_prefix + "/e_sma_slow_len")),
        e_sma_update_period (iestade::size_t_from_json(config_filepath, key_path_prefix + "/e_sma_update_period")),
        log_filename        (iestade::string_from_json(config_filepath, key_path_prefix + "/log_filename")),
        stats_filename      (iestade::string_from_json(config_filepath, key_path_prefix + "/stats_filename")),
        n_reports           (iestade::size_t_from_json(config_filepath, key_path_prefix + "/n_reports", true, 0))
    {
    }
    // clang-format on
};

class State {
protected:
    Settings _settings;
    bool _evaluated;
    double _energy;
    double _value;

    void reset_evaluation()
    {
        _evaluated = false;
        _energy    = -1;
        _value     = -1;
    }

public:
    explicit State(Settings &in_settings) :
        _settings(in_settings)
    {
        reset_evaluation();
    }

    // virtual destructor is required if virtual methods are used
    virtual ~State() {}

    virtual double get_energy()
    {
        std::cout << "error: get_energy method not implemented" << std::endl;
        return -1.0;
    }

    virtual double get_value()
    {
        std::cout << "error: get_value method not implemented" << std::endl;
        return -1.0;
    }

    virtual void randomize()
    {
        reset_evaluation();
        std::cout << "error: randomize method not implemented" << std::endl;
    }

    virtual void change()
    {
        reset_evaluation();

        // this method is very individual-specific, so to not overthink it
        // I leave it virtual
        std::cout << "error: change method not implemented" << std::endl;
    }
};

template <typename TState>
class Context {
public:
    Settings settings;

    double temperature = 0;
    bool do_cool       = false;
    size_t cooling_i   = 0;
    TState state;
    TState proposed_state;

    std::atomic<bool> q_pressed    = false;
    std::atomic<bool> stop_threads = false;

    std::vector<double> init_t_log;
    bool init_done = false;

    // important! states begin with 1st, not 0th
    double t_max   = 0;
    size_t state_i = 1;
    std::deque<double> e_log;
    size_t e_log_len;
    bool run_done = false;

    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> stop_time;
    double cycle_time_us = 0;

    aviize::Progress run_progress;
    std::ofstream log_f;

    std::queue<size_t> report_states_queue;
    bool do_report = false;

    explicit Context(Settings &s) :
        settings(s),
        state(settings),
        proposed_state(settings),
        e_log_len(settings.e_sma_slow_len)
    {
    }

    std::string get_stats()
    {
        const double runtime_s =
                std::chrono::duration_cast<std::chrono::seconds>(stop_time -
                                                                 start_time)
                        .count();
        const double state_s         = state_i / runtime_s;
        const size_t first_col_width = 22;
        std::stringstream ss{};
        // standard parameters
        ss << std::left << std::setw(first_col_width) << "states" << state_i
           << std::endl;
        // runtime stats
        ss << std::left << std::setw(first_col_width) << "runtime"
           << aviize::seconds_to_hhmmss_string(runtime_s) << std::endl;
        ss << std::left << std::setw(first_col_width) << "state/s" << state_s
           << std::endl;
        ss << std::left << std::setw(first_col_width) << "result energy"
           << state.get_energy() << std::endl;
        ss << std::left << std::setw(first_col_width) << "result value"
           << state.get_value() << std::endl;
        return ss.str();
    }
};

template <typename TState>
class StateMachine {
private:
    typedef std::function<void(Context<TState> &)> state_function_t;
    Context<TState> _context;

    void _input_handler()
    {
        while (!_context.stop_threads.load()) {
            if (std::cin.get() == 'q') {
                _context.q_pressed.store(true);
                return;
            }
        }
    }

public:
    std::vector<state_function_t> init_functions{};
    std::vector<state_function_t> init_loop_functions{};
    std::vector<state_function_t> run_loop_functions{};
    std::vector<state_function_t> finalize_functions{};

    explicit StateMachine(Settings &s) :
        _context(s)
    {
    }

    void run()
    {
        aviize::out_stream << "[enter q to interrupt]" << std::endl;
        std::thread input_thread =
                std::thread(&StateMachine::_input_handler, this);

        _context.start_time = std::chrono::steady_clock::now();
        for (const state_function_t &f : init_functions) {
            f(_context);
        }
        auto cycle_begin_time = std::chrono::steady_clock::now();
        while (init_loop_functions.size() > 0 && !_context.init_done) {
            for (const state_function_t &f : init_loop_functions) {
                if (_context.init_done) {
                    break;
                }
                f(_context);
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            _context.cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        while (run_loop_functions.size() > 0 && !_context.run_done) {
            for (const state_function_t &f : run_loop_functions) {
                if (_context.run_done) {
                    break;
                }
                f(_context);
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            _context.cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        _context.stop_time = std::chrono::steady_clock::now();
        for (const state_function_t &f : finalize_functions) {
            f(_context);
        }
        if (!_context.q_pressed.load()) {
            aviize::out_stream << "[press enter to quit]" << std::flush;
        }
        _context.stop_threads.store(true);
        input_thread.join();
    }
};

template <typename TState>
void log_init(Context<TState> &c)
{
    assert(!c.log_f.is_open());
    if (c.settings.log_filename.empty()) {
        return;
    }

    c.log_f.open(c.settings.log_filename);
    c.log_f << "state_i,temperature,energy,value" << std::endl;
}

template <typename TState>
void log_update(Context<TState> &c)
{
    // write the log to file
    assert(!c.settings.log_filename.empty());
    if (!c.log_f.is_open()) {
        return;
    }

    c.log_f << c.state_i;
    c.log_f << "," << c.temperature;
    c.log_f << "," << c.state.get_energy();
    c.log_f << "," << c.state.get_value();
    c.log_f << std::endl;
    c.log_f << std::flush;
}

template <typename TState>
void run_progress_init(Context<TState> &c)
{
    if (c.temperature > 0) {
        c.run_progress.n_min         = 1;
        c.run_progress.n_max         = c.settings.n_states;
        c.run_progress.update_period = c.settings.progress_update_period;
    }
}

template <typename TState>
void run_progress_text_reset(Context<TState> &c)
{
    c.run_progress.text.clear();
}

template <typename TState>
void run_progress_text_add_stats(Context<TState> &c)
{
    const double state_s = 1 / c.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " n/s " << state_s;

    c.run_progress.text += std::string(ss.str());
}

template <typename TState>
void run_progress_text_add_total(Context<TState> &c)
{
    if (!c.run_progress.text.empty()) {
        c.run_progress.text += " ";
    }
    c.run_progress.text += c.run_progress.str_total(c.state_i);
}

template <typename TState>
void run_progress_text_add_pct(Context<TState> &c)
{
    if (!c.run_progress.text.empty()) {
        c.run_progress.text += " ";
    }
    c.run_progress.text += c.run_progress.str_pct(c.state_i);
}

template <typename TState>
void run_progress_text_add_eta(Context<TState> &c)
{
    if (!c.run_progress.text.empty()) {
        c.run_progress.text += " ";
    }
    c.run_progress.text += c.run_progress.str_eta(c.state_i);
}

template <typename TState>
void run_progress_text_add_e(Context<TState> &c)
{
    if (!c.run_progress.text.empty()) {
        c.run_progress.text += " ";
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << c.state.get_energy();
    c.run_progress.text += "e " + oss.str();
}

template <typename TState>
void run_progress_text_add_t(Context<TState> &c)
{
    if (!c.run_progress.text.empty()) {
        c.run_progress.text += " ";
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << c.temperature;
    c.run_progress.text += "t " + oss.str();
}

template <typename TState>
void run_progress_text_add_v(Context<TState> &c)
{
    if (!c.run_progress.text.empty()) {
        c.run_progress.text += " ";
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << c.state.get_value();
    c.run_progress.text += "v " + oss.str();
}

template <typename TState>
void run_progress_text_add_freq(Context<TState> &c)
{
    if (!c.run_progress.text.empty()) {
        c.run_progress.text += " ";
    }
    const double state_s = 1 / c.cycle_time_us * 1000000;
    c.run_progress.text += "freq " + std::to_string(state_s);
}

template <typename TState>
void run_progress_print(Context<TState> &c)
{
    c.run_progress.print();
}

template <typename TState>
void run_progress_clear(Context<TState> &c)
{
    (void)c;
    aviize::erase_line();
}

template <typename TState>
void stats_print(Context<TState> &c)
{
    std::cout << c.get_stats();
}

template <typename TState>
void stats_create_file(Context<TState> &c)
{
    if (c.settings.stats_filename.empty()) {
        return;
    }

    std::ofstream f(c.settings.stats_filename);
    f << c.get_stats();
}

template <typename TState>
void state_randomize(Context<TState> &c)
{
    assert(c.state_i == 1);
    c.state.randomize();
}

template <typename TState>
void state_propose_new(Context<TState> &c)
{
    c.proposed_state = c.state;
    c.proposed_state.change();
}

template <typename TState>
void temperature_init_record(Context<TState> &c)
{
    // calculate intermediate init temperatures
    // - T = -(E_proposed - E) / ln(P)
    // - from Acceptance probability
    //   - if E_new < E_current: P = 1
    //   - else: P = exp(-(E_proposed - E) / T)
    //   - ref: https://en.wikipedia.org/wiki/Simulated_annealing
    // - only for cases where new state is worse and we need
    //   to use the p_acceptance
    // - smaller energy = better
    assert(c.init_t_log.size() < c.settings.init_t_log_len);

    const double dE = c.proposed_state.get_energy() - c.state.get_energy();
    if (dE > 0) {
        const double t = -dE / std::log(c.settings.init_p_acceptance);
        assert(t > 0);
        c.init_t_log.push_back(t);
    }
}

template <typename TState>
void temperature_init_select_as_max(Context<TState> &c)
{
    if (c.init_t_log.size() < c.settings.init_t_log_len) {
        return;
    }

    assert(c.init_t_log.size() == c.settings.init_t_log_len);
    assert(c.temperature == 0);
    assert(c.t_max == 0);

    // t_max = max(init_t at init_p_acceptance)
    c.t_max       = *max_element(c.init_t_log.begin(), c.init_t_log.end());
    c.temperature = c.t_max;
    assert(c.temperature > 0);
}

template <typename TState>
void init_done_decide(Context<TState> &c)
{
    if (c.init_t_log.size() == c.settings.init_t_log_len ||
        c.q_pressed.load()) {
        c.init_done = true;
    }
}

template <typename TState>
void state_update(Context<TState> &c)
{
    c.state_i++;
    const double dE = c.proposed_state.get_energy() - c.state.get_energy();

    if (dE < 0) {
        c.state = c.proposed_state;
        return;
    }

    // dE >= 0
    const double p_acceptance = std::exp(-dE / c.temperature);
    assert(p_acceptance <= 1);
    if (rododendrs::rnd01() <= p_acceptance) {
        c.state = c.proposed_state;
        return;
    }
}

template <typename TState>
void do_cool_always(Context<TState> &c)
{
    assert(!c.do_cool);
    c.do_cool = true;
}

template <typename TState>
void cool_at_rate(Context<TState> &c)
{
    assert(c.settings.cooling_rate > 0);
    assert(c.settings.cooling_rate <= 1);
    assert(c.temperature <= c.t_max);
    assert(c.temperature >= 0);
    if (c.do_cool) {
        // t = T0 * a^(i/R)
        // ref:
        // https://www.cicirello.org/publications/CP2007-Autonomous-Search-Workshop.pdf
        c.temperature =
                c.t_max * std::pow(c.settings.cooling_rate,
                                   std::floor(c.cooling_i /
                                              c.settings.cooling_round_len));
        c.cooling_i++;

        if (c.temperature < 0) {
            c.temperature = 0;
        }
    }
    c.do_cool = false;
}

template <typename TState>
void log_energy(Context<TState> &c)
{
    const double dE = c.proposed_state.get_energy() - c.state.get_energy();
    if (dE >= 0) {
        return;
    }

    assert(c.e_log_len > 0);
    c.e_log.push_front(c.state.get_energy());
    if (c.e_log.size() > c.e_log_len) {
        c.e_log.pop_back();
    }
    assert(c.e_log.size() <= c.e_log_len);
}

template <typename TState>
void do_cool_decide_sma(Context<TState> &c)
{
    assert(c.settings.e_sma_fast_len < c.settings.e_sma_slow_len);
    assert(!c.do_cool);
    if (c.e_log.size() < c.e_log_len) {
        return;
    }

    assert(c.settings.e_sma_update_period > 0);
    if (c.state_i % c.settings.e_sma_update_period) {
        return;
    }

    assert(c.temperature <= c.t_max);
    assert(c.temperature >= 0);

    // calculate avg(energy) at two intervals in the past
    double sum_e_sma_fast = 0;
    for (size_t i = 0; i < c.settings.e_sma_fast_len; i++) {
        sum_e_sma_fast += c.e_log[i];
    }
    const double e_sma_fast = sum_e_sma_fast / c.settings.e_sma_fast_len;

    double sum_e_sma_slow = 0;
    for (size_t i = 0; i < c.settings.e_sma_slow_len; i++) {
        sum_e_sma_slow += c.e_log[i];
    }
    const double e_sma_slow = sum_e_sma_slow / c.settings.e_sma_slow_len;

    if (e_sma_fast > e_sma_slow) {
        c.do_cool = true;
    }
}

template <typename TState>
void run_done_decide(Context<TState> &c)
{
    assert(c.temperature <= c.t_max);
    assert(c.temperature >= 0);
    assert(!c.run_done);
    if (c.state_i == c.settings.n_states || c.q_pressed.load()) {
        c.run_done = true;
    }
}

template <typename TState>
void report_linear_init(Context<TState> &c)
{
    // report state queue holds the
    // numbers of all states that must produce a report
    // - always starts with 1
    // - always end n_states
    assert(c.report_states_queue.empty());
    c.report_states_queue.push(1);
    assert(c.settings.n_reports > 0);
    const double report_step =
            c.settings.n_states / (double)(c.settings.n_reports - 1);
    for (size_t report_i = 1; report_i < c.settings.n_reports - 1;
         report_i++) {
        c.report_states_queue.push(
                static_cast<size_t>(report_step * report_i));
    }
    c.report_states_queue.push(c.settings.n_states);
    assert(c.report_states_queue.size() == c.settings.n_reports);
}

template <typename TState>
void report_at_rate_init(Context<TState> &c)
{
    // report state queue holds the
    // numbers of all states that must produce a report
    // - always starts with 1 and ends with n_states
    assert(c.report_states_queue.empty());

    // from
    // a k^0 = 1
    // a k^(n_reports-1) = n_states,
    // where a = 1, reports = [0, n_reports-1]
    // k ^ (n_reports-1) = n_states/a
    // k = (n_states/a) ^ 1/(n_reports-1)
    // k = n_states ^ 1/(n_reports-1)
    assert(c.settings.n_reports > 0);
    const double float_err_compensation = 1e-5;
    const double report_rate            = std::pow(
            c.settings.n_states, 1 / (double)(c.settings.n_reports - 1));

    double state_i = 1;
    for (size_t report_i = 1; report_i <= c.settings.n_reports; report_i++) {
#ifndef NDEBUG
        if (report_i == 1) {
            assert(static_cast<size_t>(state_i) == 1);
        }
        if (report_i == c.settings.n_reports) {
            assert(static_cast<size_t>(state_i + float_err_compensation) ==
                   c.settings.n_states);
        }
#endif
        c.report_states_queue.push(
                static_cast<size_t>(state_i + float_err_compensation));
        state_i *= report_rate;
    }
    assert(c.report_states_queue.size() == c.settings.n_reports);
}

template <typename TState>
void do_report_decide(Context<TState> &c)
{
    if (!c.report_states_queue.empty() &&
        c.state_i >= c.report_states_queue.front()) {
        c.report_states_queue.pop();
        c.do_report = true;
        return;
    }

    c.do_report = false;
}

}  // namespace lapsa
