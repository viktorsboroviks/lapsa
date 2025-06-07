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
    }
    // clang-format on
};

class State {
protected:
    bool _evaluated;
    double _energy;
    double _value;

public:
    void reset_evaluation()
    {
        _evaluated = false;
    }

    explicit State() :
        _energy(-1),
        _value(-1)
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
    typedef std::function<void(void)> state_function_t;
    const Settings *p_settings;

    // init
#ifndef NDEBUG
    std::size_t init_t_candidates_failed_n = 0;
#endif
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
            f();
        }
        auto cycle_begin_time = std::chrono::steady_clock::now();
        while (init_loop_functions.size() > 0 && !init_done) {
            for (const state_function_t &f : init_loop_functions) {
                if (init_done) {
                    break;
                }
                f();
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        while (run_loop_functions.size() > 0 && !run_done) {
            for (const state_function_t &f : run_loop_functions) {
                if (run_done) {
                    break;
                }
                f();
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
            f();
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
        ss << std::left << std::setw(first_col_width) << "result value"
           << state.get_value() << std::endl;
        return ss.str();
    }

    void init_log()
    {
        assert(!log_f.is_open());
        if (p_settings->log_file_name.empty()) {
            return;
        }

        log_f.open(p_settings->log_file_name);
        log_f << "run_i,t,e,v" << std::endl;
    }

    void update_log()
    {
        // write the log to file
        assert(!p_settings->log_file_name.empty());
        if (!log_f.is_open()) {
            return;
        }

        log_f << run_i;
        log_f << "," << t;
        log_f << "," << state.get_energy();
        log_f << "," << state.get_value();
        log_f << std::endl;
    }

    void init_progress()
    {
        if (t > 0) {
            progress.n_min         = 1;
            progress.n_max         = p_settings->n_states;
            progress.update_period = p_settings->progress_update_period;
        }
    }

    void progress_text_reset()
    {
        progress.reset();
    }

    void progress_text_add_nl()
    {
        progress.text += "\n";
    }

    void progress_text_add_stats()
    {
        const double run_s = 1 / cycle_time_us * 1000000;

        std::stringstream ss;
        ss << " n/s " << run_s;

        progress.text += std::string(ss.str());
    }

    void progress_text_add_total()
    {
        if (!last_line_empty(progress.text)) {
            progress.text += " ";
        }
        progress.text += progress.str_total(run_i);
    }

    void progress_text_add_pct()
    {
        if (!last_line_empty(progress.text)) {
            progress.text += " ";
        }
        progress.text += progress.str_pct(run_i);
    }

    void progress_text_add_eta()
    {
        if (!last_line_empty(progress.text)) {
            progress.text += " ";
        }
        progress.text += progress.str_eta(run_i);
    }

    void progress_text_add_e()
    {
        if (!last_line_empty(progress.text)) {
            progress.text += " ";
        }
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(3) << state.get_energy();
        progress.text += "e " + oss.str();
    }

    void progress_text_add_t()
    {
        if (!last_line_empty(progress.text)) {
            progress.text += " ";
        }
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(3) << t;
        progress.text += "t " + oss.str();
    }

    void progress_text_add_v()
    {
        if (!last_line_empty(progress.text)) {
            progress.text += " ";
        }
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(3) << state.get_value();
        progress.text += "v " + oss.str();
    }

    void progress_text_add_freq()
    {
        if (!last_line_empty(progress.text)) {
            progress.text += " ";
        }
        const double run_s = 1 / cycle_time_us * 1000000;
        progress.text += "freq " + std::to_string(run_s);
    }

    void progress_print()
    {
        progress.print();
    }

    void progress_clear_line()
    {
        aviize::erase_line(progress.text);
        progress.print();
    }

    void stats_print()
    {
        std::stringstream ss{};
        ss << get_stats();
        aviize::print(ss);
    }

    void stats_create_file()
    {
        if (p_settings->stats_file_name.empty()) {
            return;
        }

        std::ofstream f(p_settings->stats_file_name);
        f << get_stats();
    }

    void randomize_state()
    {
        assert(run_i == 1);
        state.reset_evaluation();
        state.randomize();
    }

    void propose_new_state()
    {
        proposed_state = state;
        proposed_state.reset_evaluation();
        proposed_state.change();
    }

    void generate_init_t_candidates()
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
        assert(init_t_candidates.size() < p_settings->init_t_candidates_n);

        const double dE = proposed_state.get_energy() - state.get_energy();
        if (dE > 0) {
            const double t = -dE / std::log(p_settings->init_p_acceptance);
            assert(t > 0);
            init_t_candidates.push_back(t);
#ifndef NDEBUG
            init_t_candidates_failed_n = 0;
#endif
        }
#ifndef NDEBUG
        else {
            init_t_candidates_failed_n++;
        }
#endif
    }

    void select_init_t_as_max()
    {
        if (init_t_candidates.size() < p_settings->init_t_candidates_n) {
            return;
        }

        assert(init_t_candidates.size() == p_settings->init_t_candidates_n);
        assert(t == 0);
        assert(t_max == 0);

        // temperature_max = max(init_t at init_p_acceptance)
        t_max = *max_element(init_t_candidates.begin(),
                             init_t_candidates.end());
        t     = t_max;
        assert(t > 0);
    }

    void decide_init_done()
    {
        if (init_t_candidates.size() == p_settings->init_t_candidates_n) {
            init_done = true;
        }
#ifndef NDEBUG
        if (init_t_candidates_failed_n > p_settings->init_t_candidates_n) {
            std::cerr << "init_t_candidates_failed_n: "
                      << init_t_candidates_failed_n << std::endl;
            std::cerr << "init_t_candidates_n: "
                      << p_settings->init_t_candidates_n << std::endl;
            assert(false);
        }
#endif
    }

    void update_state()
    {
        run_i++;
        const double dE = proposed_state.get_energy() - state.get_energy();

        if (dE < 0) {
            state = proposed_state;
            return;
        }

        // dE >= 0
        const double p_acceptance = std::exp(-dE / t);
        if (rododendrs::rnd01() <= p_acceptance) {
            state = proposed_state;
            return;
        }
    }

    void reset_state_evaluation()
    {
        state.reset_evaluation();
    }

    void do_cool_set()
    {
        assert(!do_cool);
        do_cool = true;
    }

    void cool_at_rate()
    {
        assert(p_settings->cooling_rate > 0);
        assert(p_settings->cooling_rate <= 1);
        assert(t <= t_max);
        assert(t >= 0);
        if (do_cool) {
            // t = T0 * a^(i/R)
            // ref:
            // https://www.cicirello.org/publications/CP2007-Autonomous-Search-Workshop.pdf
            t = t_max * std::pow(p_settings->cooling_rate,
                                 std::floor(cooling_i /
                                            p_settings->cooling_round_len));
            cooling_i++;

            if (t < 0) {
                t = 0;
            }
        }
        do_cool = false;
    }

    void e_history_update()
    {
        const double dE = proposed_state.get_energy() - state.get_energy();
        if (dE >= 0) {
            return;
        }

        assert(p_settings->e_history_len > 0);
        e_history.push_front(state.get_energy());
        if (e_history.size() > p_settings->e_history_len) {
            e_history.pop_back();
        }
        assert(e_history.size() <= p_settings->e_history_len);
    }

    void decide_cool_sma()
    {
        assert(p_settings->e_sma_fast_len < p_settings->e_sma_slow_len);
        assert(!do_cool);
        if (e_history.size() < p_settings->e_history_len) {
            return;
        }

        assert(p_settings->e_decision_period > 0);
        if (run_i % p_settings->e_decision_period) {
            return;
        }

        assert(t <= t_max);
        assert(t >= 0);

        // calculate avg(energy) at two intervals in the past
        double sum_e_sma_fast = 0;
        for (size_t i = 0; i < p_settings->e_sma_fast_len; i++) {
            sum_e_sma_fast += e_history[i];
        }
        const double e_sma_fast = sum_e_sma_fast / p_settings->e_sma_fast_len;

        double sum_e_sma_slow = 0;
        for (size_t i = 0; i < p_settings->e_sma_slow_len; i++) {
            sum_e_sma_slow += e_history[i];
        }
        const double e_sma_slow = sum_e_sma_slow / p_settings->e_sma_slow_len;

        if (e_sma_fast > e_sma_slow) {
            do_cool = true;
        }
    }

    void decide_cool_min_sd()
    {
        assert(p_settings->e_window > 0);
        assert(p_settings->e_shift > 0);
        const size_t req_e_history_len =
                p_settings->e_shift + p_settings->e_window;
        assert(p_settings->e_history_len == req_e_history_len);

        assert(!do_cool);
        if (e_history.size() < req_e_history_len) {
            return;
        }

        assert(p_settings->e_decision_period > 0);
        if (run_i % p_settings->e_decision_period) {
            return;
        }

        assert(t <= t_max);
        assert(t >= 0);

        // calculate energy at two intervals in the past
        const auto front_begin = e_history.begin();
        const auto front_end   = e_history.begin() + p_settings->e_window;
        const std::vector<double> front_window(front_begin, front_end);
        assert(front_window.size() == p_settings->e_window);
        const double front_sd = rododendrs::sd<std::vector>(front_window);

        const auto back_begin = e_history.begin() + p_settings->e_shift;
        const auto back_end =
                e_history.begin() + p_settings->e_window + p_settings->e_shift;
        const std::vector<double> back_window(back_begin, back_end);
        assert(back_window.size() == p_settings->e_window);
        const double back_mean = rododendrs::mean<std::vector>(back_window);
        const double back_sd   = rododendrs::sd<std::vector>(back_window);

        if (std::abs(back_mean - e_history[0]) < std::min(front_sd, back_sd)) {
            do_cool = true;
        }
    }

    void decide_cool_az()
    {
        assert(p_settings->e_window > 0);
        assert(p_settings->e_shift > 0);
        const size_t req_e_history_len =
                p_settings->e_shift + p_settings->e_window;
        assert(p_settings->e_history_len == req_e_history_len);

        assert(!do_cool);
        if (e_history.size() < req_e_history_len) {
            return;
        }

        assert(p_settings->e_decision_period > 0);
        if (run_i % p_settings->e_decision_period) {
            return;
        }

        assert(t <= t_max);
        assert(t >= 0);

        // calculate energy at two intervals in the past
        const auto front_begin = e_history.begin();
        const auto front_end   = e_history.begin() + p_settings->e_window;
        const std::vector<double> front_window(front_begin, front_end);
        assert(front_window.size() == p_settings->e_window);
        const double front_mean = rododendrs::mean<std::vector>(front_window);
        const double front_sd   = rododendrs::sd<std::vector>(front_window);

        const auto back_begin = e_history.begin() + p_settings->e_shift;
        const auto back_end =
                e_history.begin() + p_settings->e_window + p_settings->e_shift;
        const std::vector<double> back_window(back_begin, back_end);
        assert(back_window.size() == p_settings->e_window);
        const double back_mean = rododendrs::mean<std::vector>(back_window);
        const double back_sd   = rododendrs::sd<std::vector>(back_window);

        const double az_overlap =
                rododendrs::az_pdf_overlap<double>(front_mean,
                                                   front_sd,
                                                   p_settings->e_window,
                                                   back_mean,
                                                   back_sd,
                                                   p_settings->e_window);
        if (az_overlap >= p_settings->e_min_az_overlap) {
            do_cool = true;
        }
    }

    void decide_run_done()
    {
        assert(t <= t_max);
        assert(t >= 0);
        assert(!run_done);
        if (run_i == p_settings->n_states) {
            run_done = true;
        }
    }

    void init_rec_linear()
    {
        // record state queue holds the
        // numbers of all states that must produce a record
        // - always starts with 1
        // - always end n_states
        assert(rec_states_queue.empty());
        rec_states_queue.push(1);
        assert(p_settings->n_records > 0);
        const double rec_step =
                p_settings->n_states / (double)(p_settings->n_records - 1);
        for (size_t rec_i = 1; rec_i < p_settings->n_records - 1; rec_i++) {
            rec_states_queue.push(static_cast<size_t>(rec_step * rec_i));
        }
        rec_states_queue.push(p_settings->n_states);
        assert(rec_states_queue.size() == p_settings->n_records);
    }

    void init_rec_at_rate()
    {
        // record state queue holds the
        // numbers of all states that must produce a record
        // - always starts with 1 and ends with n_states
        assert(rec_states_queue.empty());

        // from
        // a k^0 = 1
        // a k^(n_records-1) = n_states,
        // where a = 1, records = [0, n_records-1]
        // k ^ (n_records-1) = n_states/a
        // k = (n_states/a) ^ 1/(n_records-1)
        // k = n_states ^ 1/(n_records-1)
        assert(p_settings->n_records > 0);
        const double float_err_compensation = 1e-5;
        const double rec_rate               = std::pow(
                p_settings->n_states, 1 / (double)(p_settings->n_records - 1));

        double run_i = 1;
        for (size_t rec_i = 1; rec_i <= p_settings->n_records; rec_i++) {
#ifndef NDEBUG
            if (rec_i == 1) {
                assert(static_cast<size_t>(run_i) == 1);
            }
            if (rec_i == p_settings->n_records) {
                assert(static_cast<size_t>(run_i + float_err_compensation) ==
                       p_settings->n_states);
            }
#endif
            rec_states_queue.push(
                    static_cast<size_t>(run_i + float_err_compensation));
            run_i *= rec_rate;
        }
        assert(rec_states_queue.size() == p_settings->n_records);
    }

    void init_rec_periods()
    {
        // record state queue holds the
        // numbers of all states that must produce a record
        // - always starts with 1
        // - always end n_states
        assert(rec_states_queue.empty());
        for (size_t i = 1; i <= p_settings->n_states; i++) {
            if (p_settings->rec_schedule.is_time(i)) {
                rec_states_queue.push(i);
            }
        }
    }

    void decide_rec()
    {
        if (!rec_states_queue.empty() && run_i >= rec_states_queue.front()) {
            rec_states_queue.pop();
            do_rec = true;
            return;
        }

        do_rec = false;
    }
};

}  // namespace lapsa
