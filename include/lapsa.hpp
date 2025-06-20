// #define LAPSA_DEBUG

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
#include <vector>

#include "aviize.hpp"
#include "iestaade.hpp"
#include "rododendrs.hpp"

#ifdef LAPSA_DEBUG
#define debug_lapsa(x)                                                  \
    std::cout << "debug (" << __FILE__ << ":" << __LINE__ << "): " << x \
              << std::endl;
#else
#define debug_lapsa(x) while (0)
#endif

namespace lapsa {

// utility
bool last_line_empty(const std::string &text)
{
    return (text.empty() || text.back() == '\n' || text.back() == '\r');
}

// tools
class Schedule {
private:
    struct Period {
        size_t i_start;
        size_t period;
    };

    std::vector<Period> _periods;

    void _init_periods(const std::string &config_path,
                       const std::string &key_path,
                       std::vector<Period> &periods)
    {
        const boost::json::array periods_json =
                iestaade::value_from_json(config_path, key_path).as_array();
        assert(periods_json.size() > 0);
        assert(periods.empty());

        for (const boost::json::value &period_json : periods_json) {
            Period new_period;
            new_period.i_start = period_json.as_array()[0].to_number<size_t>();
            new_period.period  = period_json.as_array()[1].to_number<size_t>();
            periods.push_back(new_period);
        }
    }

    size_t _i_period(size_t i, const std::vector<Period> &periods) const
    {
        assert(!periods.empty());

        // check i
        if (i < periods.front().i_start) {
            return 0;
        }

        // go over all periods
        for (size_t i_period = 0; i_period + 1 < periods.size(); i_period++) {
            // find first that fits
            if (i >= periods[i_period].i_start &&
                i < periods[i_period + 1].i_start) {
                return i_period;
            }
        }

        // only last period remains
        if (i >= periods.back().i_start) {
            return periods.size() - 1;
        }
        else {
            assert(false);
            return 0;
        }
    }

public:
    Schedule() {}

    Schedule(const std::string &config_path, const std::string &key_path)
    {
        _init_periods(config_path, key_path, _periods);
    }

    bool is_time(size_t i) const
    {
        const size_t i_period = _i_period(i, _periods);
        return i % _periods[i_period].period == 0;
    }
};

struct Settings {
    size_t n_states            = 1000000;
    double init_p_acceptance   = 0.99;
    size_t init_t_candidates_n = 100;
    size_t init_t_max_attempts = 100;
    double cooling_rate        = 0.95;
    size_t cooling_round_len   = 1;
    size_t e_decision_period   = 100;
    size_t e_sma_fast_len      = 50;
    size_t e_sma_slow_len      = 200;
    size_t e_window            = 100;
    size_t e_shift             = 100;
    double e_min_az_overlap    = 0.99;
    size_t e_history_len;

    size_t progress_update_period = 1;
    std::string log_file_name     = "log.csv";
    std::string stats_file_name   = "stats.txt";

    size_t n_records = 0;
    Schedule rec_schedule;

    Settings() {}

    // clang-format off
    Settings(const std::string& config_filepath,
             const std::string& key_path_prefix) :
        n_states            (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/n_states")),
        init_p_acceptance   (iestaade::double_from_json(config_filepath, key_path_prefix + "/init_p_acceptance")),
        init_t_candidates_n (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/init_t_candidates_n")),
        init_t_max_attempts (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/init_t_max_attempts")),
        cooling_rate        (iestaade::double_from_json(config_filepath, key_path_prefix + "/cooling_rate")),
        cooling_round_len   (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/cooling_round_len")),
        e_decision_period   (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/e_decision_period")),
        e_sma_fast_len      (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/e_sma_fast_len")),
        e_sma_slow_len      (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/e_sma_slow_len")),
        e_window            (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/e_window")),
        e_shift             (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/e_shift")),
        e_min_az_overlap    (iestaade::double_from_json(config_filepath, key_path_prefix + "/e_min_az_overlap")),
        e_history_len       (std::max(e_sma_slow_len, e_shift + e_window)),
        log_file_name       (iestaade::string_from_json(config_filepath, key_path_prefix + "/log_file_name")),
        stats_file_name     (iestaade::string_from_json(config_filepath, key_path_prefix + "/stats_file_name")),
        n_records           (iestaade::size_t_from_json(config_filepath, key_path_prefix + "/n_records", true, 0)),
        rec_schedule(config_filepath, key_path_prefix + "/rec_periods")
    {
        assert(init_t_candidates_n <= init_t_max_attempts);
    }
    // clang-format on
};

class State {
protected:
    bool _evaluated;
    double _energy;

public:
    void reset_evaluation()
    {
        _evaluated = false;
    }

    explicit State() :
        _energy(-1)
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

    virtual void randomize()
    {
        std::cout << "error: randomize method not implemented" << std::endl;
    }

    virtual void change()
    {
        // this method is very individual-specific, so to not overthink it
        // I leave it virtual
        std::cout << "error: change method not implemented" << std::endl;
    }
};

template <typename TState>
struct StateMachine {
    // settings
    typedef std::function<void(StateMachine<TState> &)> state_function_t;
    const Settings *p_settings;

    // init
    bool init_done = false;

    // runtime
    // important! states begin with 1st, not 0th
    double t         = 0;
    double t_max     = 0;
    bool do_cool     = false;
    size_t cooling_i = 0;
    size_t run_i     = 1;
    bool run_done    = false;
    bool do_rec      = false;

    // performance
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> stop_time;
    double cycle_time_us = 0;

    // history
    std::vector<double> init_t_candidates;
    std::deque<double> e_history;
    std::queue<size_t> rec_states_queue;

    // data display and logging
    aviize::Progress progress;
    std::ofstream log_f;

    // state
    TState state;
    TState proposed_state;

    // state machine functions
    std::vector<state_function_t> init_functions{};
    std::vector<state_function_t> init_loop_functions{};
    std::vector<state_function_t> run_functions{};
    std::vector<state_function_t> run_loop_functions{};
    std::vector<state_function_t> finalize_functions{};

    explicit StateMachine(const Settings *s) :
        p_settings(s)
    {
    }

    void run()
    {
        start_time = std::chrono::steady_clock::now();
        for (const state_function_t &f : init_functions) {
            f(*this);
        }
        auto cycle_begin_time = std::chrono::steady_clock::now();
        while (init_loop_functions.size() > 0 && !init_done) {
            for (const state_function_t &f : init_loop_functions) {
                if (init_done) {
                    break;
                }
                f(*this);
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        for (const state_function_t &f : run_functions) {
            f(*this);
        }
        while (run_loop_functions.size() > 0 && !run_done) {
            for (const state_function_t &f : run_loop_functions) {
                if (run_done) {
                    break;
                }
                f(*this);
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        stop_time = std::chrono::steady_clock::now();
        for (const state_function_t &f : finalize_functions) {
            f(*this);
        }
    }

    std::string get_stats()
    {
        const double runtime_s =
                std::chrono::duration_cast<std::chrono::seconds>(stop_time -
                                                                 start_time)
                        .count();
        const double run_s           = run_i / runtime_s;
        const size_t first_col_width = 16;
        std::stringstream ss{};
        // standard parameters
        ss << std::left << std::setw(first_col_width) << "states" << run_i
           << std::endl;
        // runtime stats
        ss << std::left << std::setw(first_col_width) << "runtime"
           << aviize::seconds_to_hhmmss_string(runtime_s) << std::endl;
        ss << std::left << std::setw(first_col_width) << "state/s" << run_s
           << std::endl;
        ss << std::left << std::setw(first_col_width) << "result energy"
           << state.get_energy() << std::endl;
        return ss.str();
    }
};

template <typename TState>
void init_log(StateMachine<TState> &sm)
{
    assert(!sm.log_f.is_open());
    if (sm.p_settings->log_file_name.empty()) {
        return;
    }

    sm.log_f.open(sm.p_settings->log_file_name);
    sm.log_f << "run_i,";
    sm.log_f << "t,";
    sm.log_f << "e";
    sm.log_f << std::endl;
}

template <typename TState>
void update_log(StateMachine<TState> &sm)
{
    // write the log to file
    assert(!sm.p_settings->log_file_name.empty());
    if (!sm.log_f.is_open()) {
        return;
    }

    sm.log_f << sm.run_i;
    sm.log_f << "," << sm.t;
    sm.log_f << "," << sm.state.get_energy();
    sm.log_f << std::endl;
}

template <typename TState>
void progress_init_init_loop(lapsa::StateMachine<TState> &sm)
{
    sm.progress.n_min         = 1;
    sm.progress.n_max         = sm.p_settings->init_t_max_attempts;
    sm.progress.update_period = sm.p_settings->progress_update_period;
}

template <typename TState>
void progress_init_run_loop(lapsa::StateMachine<TState> &sm)
{
    sm.progress.n_min         = 1;
    sm.progress.n_max         = sm.p_settings->n_states;
    sm.progress.update_period = sm.p_settings->progress_update_period;
}

template <typename TState>
void progress_text_reset(StateMachine<TState> &sm)
{
    sm.progress.reset();
}

template <typename TState>
void progress_text_add_nl(StateMachine<TState> &sm)
{
    sm.progress.text += "\n";
}

template <typename TState>
void progress_text_add_title_init_loop_nl(StateMachine<TState> &sm)
{
    std::stringstream ss;
    ss << "init loop:" << std::endl;
    sm.progress.text += ss.str();
}

template <typename TState>
void progress_text_add_title_run_loop_nl(StateMachine<TState> &sm)
{
    std::stringstream ss;
    ss << "run loop:" << std::endl;
    sm.progress.text += ss.str();
}

template <typename TState>
void progress_text_add_stats(StateMachine<TState> &sm)
{
    const double run_s = 1 / sm.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " n/s " << run_s;

    sm.progress.text += std::string(ss.str());
}

template <typename TState>
void progress_text_add_total(StateMachine<TState> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    sm.progress.text += sm.progress.str_total(sm.run_i);
}

template <typename TState>
void progress_text_add_pct(StateMachine<TState> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    sm.progress.text += sm.progress.str_pct(sm.run_i);
}

template <typename TState>
void progress_text_add_eta(StateMachine<TState> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    sm.progress.text += sm.progress.str_eta(sm.run_i);
}

template <typename TState>
void progress_text_add_e(StateMachine<TState> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << sm.state.get_energy();
    sm.progress.text += "e " + oss.str();
}

template <typename TState>
void progress_text_add_t(StateMachine<TState> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << sm.t;
    sm.progress.text += "t " + oss.str();
}

template <typename TState>
void progress_text_add_init_candidates(StateMachine<TState> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    std::ostringstream oss;
    oss << sm.init_t_candidates.size();
    sm.progress.text += "init candidates  " + oss.str();
}

template <typename TState>
void progress_text_add_freq(StateMachine<TState> &sm)
{
    if (!last_line_empty(sm.progress.text)) {
        sm.progress.text += " ";
    }
    const double run_s = 1 / sm.cycle_time_us * 1000000;
    sm.progress.text += "freq " + std::to_string(run_s);
}

template <typename TState>
void progress_print(StateMachine<TState> &sm)
{
    sm.progress.print();
}

template <typename TState>
void progress_clear_line(StateMachine<TState> &sm)
{
    aviize::erase_line(sm.progress.text);
    sm.progress.print();
}

template <typename TState>
void print_stats(StateMachine<TState> &sm)
{
    std::stringstream ss{};
    ss << sm.get_stats();
    aviize::print(ss);
}

template <typename TState>
void create_stats_file(StateMachine<TState> &sm)
{
    if (sm.p_settings->stats_file_name.empty()) {
        return;
    }

    std::ofstream f(sm.p_settings->stats_file_name);
    f << sm.get_stats();
}

template <typename TState>
void randomize_state(StateMachine<TState> &sm)
{
    assert(sm.run_i == 1);
    sm.state.reset_evaluation();
    sm.state.randomize();
}

template <typename TState>
void propose_new_state(StateMachine<TState> &sm)
{
    sm.proposed_state = sm.state;
    sm.proposed_state.reset_evaluation();
    sm.proposed_state.change();
}

template <typename TState>
void generate_init_t_candidates(StateMachine<TState> &sm)
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
    assert(sm.init_t_candidates.size() < sm.p_settings->init_t_candidates_n);

    const double dE = sm.proposed_state.get_energy() - sm.state.get_energy();
    if (dE > 0) {
        const double t = -dE / std::log(sm.p_settings->init_p_acceptance);
        assert(t > 0);
        sm.init_t_candidates.push_back(t);
    }
}

template <typename TState>
void select_init_t_as_max(StateMachine<TState> &sm)
{
    if (sm.init_t_candidates.size() < sm.p_settings->init_t_candidates_n) {
        return;
    }

    assert(sm.init_t_candidates.size() == sm.p_settings->init_t_candidates_n);
    assert(sm.t == 0);
    assert(sm.t_max == 0);

    // temperature_max = max(init_t at init_p_acceptance)
    sm.t_max = *max_element(sm.init_t_candidates.begin(),
                            sm.init_t_candidates.end());
    sm.t     = sm.t_max;
    assert(sm.t > 0);
}

template <typename TState>
void decide_init_done(StateMachine<TState> &sm)
{
    if (sm.run_i >= sm.p_settings->init_t_max_attempts) {
        std::cerr << "init max attempts reached: "
                  << sm.p_settings->init_t_max_attempts << std::endl;
        std::cerr << "valid candidates: " << sm.init_t_candidates.size()
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    if (sm.init_t_candidates.size() == sm.p_settings->init_t_candidates_n) {
        sm.init_done = true;
        sm.run_i     = 1;
    }
}

template <typename TState>
void run_i_zeroize(StateMachine<TState> &sm)
{
    sm.run_i = 0;
}

template <typename TState>
void run_i_inc(StateMachine<TState> &sm)
{
    sm.run_i++;
}

template <typename TState>
void update_state(StateMachine<TState> &sm)
{
    const double dE = sm.proposed_state.get_energy() - sm.state.get_energy();

    if (dE < 0) {
        sm.state = sm.proposed_state;
        return;
    }

    // dE >= 0
    const double p_acceptance = std::exp(-dE / sm.t);
    if (rododendrs::rnd01() <= p_acceptance) {
        sm.state = sm.proposed_state;
        return;
    }
}

template <typename TState>
void reset_state_evaluation(StateMachine<TState> &sm)
{
    sm.state.reset_evaluation();
}

template <typename TState>
void do_cool_set(StateMachine<TState> &sm)
{
    assert(!sm.do_cool);
    sm.do_cool = true;
}

template <typename TState>
void cool_at_rate(StateMachine<TState> &sm)
{
    assert(sm.p_settings->cooling_rate > 0);
    assert(sm.p_settings->cooling_rate <= 1);
    assert(sm.t <= sm.t_max);
    assert(sm.t >= 0);
    if (sm.do_cool) {
        // t = T0 * a^(i/R)
        // ref:
        // https://www.cicirello.org/publications/CP2007-Autonomous-Search-Workshop.pdf
        sm.t = sm.t_max *
               std::pow(sm.p_settings->cooling_rate,
                        std::floor(sm.cooling_i /
                                   sm.p_settings->cooling_round_len));
        sm.cooling_i++;

        if (sm.t < 0) {
            sm.t = 0;
        }
    }
    sm.do_cool = false;
}

template <typename TState>
void e_history_update(StateMachine<TState> &sm)
{
    const double dE = sm.proposed_state.get_energy() - sm.state.get_energy();
    if (dE >= 0) {
        return;
    }

    assert(sm.p_settings->e_history_len > 0);
    sm.e_history.push_front(sm.state.get_energy());
    if (sm.e_history.size() > sm.p_settings->e_history_len) {
        sm.e_history.pop_back();
    }
    assert(sm.e_history.size() <= sm.p_settings->e_history_len);
}

template <typename TState>
void decide_cool_sma(StateMachine<TState> &sm)
{
    assert(sm.p_settings->e_sma_fast_len < sm.p_settings->e_sma_slow_len);
    assert(!sm.do_cool);
    if (sm.e_history.size() < sm.p_settings->e_history_len) {
        return;
    }

    assert(sm.p_settings->e_decision_period > 0);
    if (sm.run_i % sm.p_settings->e_decision_period) {
        return;
    }

    assert(sm.t <= sm.t_max);
    assert(sm.t >= 0);

    // calculate avg(energy) at two intervals in the past
    double sum_e_sma_fast = 0;
    for (size_t i = 0; i < sm.p_settings->e_sma_fast_len; i++) {
        sum_e_sma_fast += sm.e_history[i];
    }
    const double e_sma_fast = sum_e_sma_fast / sm.p_settings->e_sma_fast_len;

    double sum_e_sma_slow = 0;
    for (size_t i = 0; i < sm.p_settings->e_sma_slow_len; i++) {
        sum_e_sma_slow += sm.e_history[i];
    }
    const double e_sma_slow = sum_e_sma_slow / sm.p_settings->e_sma_slow_len;

    if (e_sma_fast > e_sma_slow) {
        sm.do_cool = true;
    }
}

template <typename TState>
void decide_cool_min_sd(StateMachine<TState> &sm)
{
    assert(sm.p_settings->e_window > 0);
    assert(sm.p_settings->e_shift > 0);
    const size_t req_e_history_len =
            sm.p_settings->e_shift + sm.p_settings->e_window;
    assert(sm.p_settings->e_history_len == req_e_history_len);

    assert(!sm.do_cool);
    if (sm.e_history.size() < req_e_history_len) {
        return;
    }

    assert(sm.p_settings->e_decision_period > 0);
    if (sm.run_i % sm.p_settings->e_decision_period) {
        return;
    }

    assert(sm.t <= sm.t_max);
    assert(sm.t >= 0);

    // calculate energy at two intervals in the past
    const auto front_begin = sm.e_history.begin();
    const auto front_end   = sm.e_history.begin() + sm.p_settings->e_window;
    const std::vector<double> front_window(front_begin, front_end);
    assert(front_window.size() == sm.p_settings->e_window);
    const double front_sd = rododendrs::sd<std::vector>(front_window);

    const auto back_begin = sm.e_history.begin() + sm.p_settings->e_shift;
    const auto back_end   = sm.e_history.begin() + sm.p_settings->sm.e_window +
                          sm.p_settings->e_shift;
    const std::vector<double> back_window(back_begin, back_end);
    assert(back_window.size() == sm.p_settings->e_window);
    const double back_mean = rododendrs::mean<std::vector>(back_window);
    const double back_sd   = rododendrs::sd<std::vector>(back_window);

    if (std::abs(back_mean - sm.e_history[0]) < std::min(front_sd, back_sd)) {
        sm.do_cool = true;
    }
}

template <typename TState>
void decide_cool_az(StateMachine<TState> &sm)
{
    assert(sm.p_settings->e_window > 0);
    assert(sm.p_settings->e_shift > 0);
    const size_t req_e_history_len =
            sm.p_settings->e_shift + sm.p_settings->e_window;
    assert(sm.p_settings->e_history_len == req_e_history_len);

    assert(!sm.do_cool);
    if (sm.e_history.size() < req_e_history_len) {
        return;
    }

    assert(sm.p_settings->e_decision_period > 0);
    if (sm.run_i % sm.p_settings->e_decision_period) {
        return;
    }

    assert(sm.t <= sm.t_max);
    assert(sm.t >= 0);

    // calculate energy at two intervals in the past
    const auto front_begin = sm.e_history.begin();
    const auto front_end   = sm.e_history.begin() + sm.p_settings->e_window;
    const std::vector<double> front_window(front_begin, front_end);
    assert(front_window.size() == sm.p_settings->e_window);
    const double front_mean = rododendrs::mean<std::vector>(front_window);
    const double front_sd   = rododendrs::sd<std::vector>(front_window);

    const auto back_begin = sm.e_history.begin() + sm.p_settings->e_shift;
    const auto back_end   = sm.e_history.begin() + sm.p_settings->e_window +
                          sm.p_settings->e_shift;
    const std::vector<double> back_window(back_begin, back_end);
    assert(back_window.size() == sm.p_settings->e_window);
    const double back_mean = rododendrs::mean<std::vector>(back_window);
    const double back_sd   = rododendrs::sd<std::vector>(back_window);

    const double az_overlap =
            rododendrs::az_pdf_overlap<double>(front_mean,
                                               front_sd,
                                               sm.p_settings->e_window,
                                               back_mean,
                                               back_sd,
                                               sm.p_settings->e_window);
    if (az_overlap >= sm.p_settings->e_min_az_overlap) {
        sm.do_cool = true;
    }
}

template <typename TState>
void decide_run_done(StateMachine<TState> &sm)
{
    assert(sm.t <= sm.t_max);
    assert(sm.t >= 0);
    assert(!sm.run_done);
    if (sm.run_i == sm.p_settings->n_states) {
        sm.run_done = true;
    }
}

template <typename TState>
void init_rec_linear(StateMachine<TState> &sm)
{
    // record state queue holds the
    // numbers of all states that must produce a record
    // - always starts with 1
    // - always end n_states
    assert(sm.rec_states_queue.empty());
    sm.rec_states_queue.push(1);
    assert(sm.p_settings->n_records > 0);
    const double rec_step =
            sm.p_settings->n_states / (double)(sm.p_settings->n_records - 1);
    for (size_t rec_i = 1; rec_i < sm.p_settings->n_records - 1; rec_i++) {
        sm.rec_states_queue.push(static_cast<size_t>(rec_step * rec_i));
    }
    sm.rec_states_queue.push(sm.p_settings->n_states);
    assert(sm.rec_states_queue.size() == sm.p_settings->n_records);
}

template <typename TState>
void init_rec_at_rate(StateMachine<TState> &sm)
{
    // record state queue holds the
    // numbers of all states that must produce a record
    // - always starts with 1 and ends with n_states
    assert(sm.rec_states_queue.empty());

    // from
    // a k^0 = 1
    // a k^(n_records-1) = n_states,
    // where a = 1, records = [0, n_records-1]
    // k ^ (n_records-1) = n_states/a
    // k = (n_states/a) ^ 1/(n_records-1)
    // k = n_states ^ 1/(n_records-1)
    assert(sm.p_settings->n_records > 0);
    const double float_err_compensation = 1e-5;
    const double rec_rate =
            std::pow(sm.p_settings->n_states,
                     1 / (double)(sm.p_settings->n_records - 1));

    double run_i = 1;
    for (size_t rec_i = 1; rec_i <= sm.p_settings->n_records; rec_i++) {
#ifndef NDEBUG
        if (rec_i == 1) {
            assert(static_cast<size_t>(run_i) == 1);
        }
        if (rec_i == sm.p_settings->n_records) {
            assert(static_cast<size_t>(run_i + float_err_compensation) ==
                   sm.p_settings->n_states);
        }
#endif
        sm.rec_states_queue.push(
                static_cast<size_t>(run_i + float_err_compensation));
        run_i *= rec_rate;
    }
    assert(sm.rec_states_queue.size() == sm.p_settings->n_records);
}

template <typename TState>
void init_rec_periods(StateMachine<TState> &sm)
{
    // record state queue holds the
    // numbers of all states that must produce a record
    // - always starts with 1
    // - always end n_states
    assert(sm.rec_states_queue.empty());
    for (size_t i = 1; i <= sm.p_settings->n_states; i++) {
        if (sm.p_settings->rec_schedule.is_time(i)) {
            sm.rec_states_queue.push(i);
        }
    }
}

template <typename TState>
void decide_rec(StateMachine<TState> &sm)
{
    if (!sm.rec_states_queue.empty() &&
        sm.run_i >= sm.rec_states_queue.front()) {
        sm.rec_states_queue.pop();
        sm.do_rec = true;
        return;
    }

    sm.do_rec = false;
}

}  // namespace lapsa
