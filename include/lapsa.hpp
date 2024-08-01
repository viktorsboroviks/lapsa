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
#include <vector>

#include "rododendrs.hpp"

namespace lapsa {

// tools
std::string seconds_to_hhmmss_string(const double seconds)
{
    std::stringstream ss{};
    ss << std::setw(2) << std::setfill('0') << ((int)seconds / 60 / 60) % 60;
    ss << ":" << std::setw(2) << std::setfill('0') << ((int)seconds / 60) % 60;
    ss << ":" << std::setw(2) << std::setfill('0') << (int)seconds % 60;
    return ss.str();
}

class Progress {
private:
    std::chrono::time_point<std::chrono::steady_clock> last_update_time =
            std::chrono::steady_clock::now();

public:
    std::ostream &os = std::cerr;
    char c_opening_bracket = '[';
    char c_closing_bracket = ']';
    char c_fill = '.';
    char c_no_fill = ' ';
    size_t c_bar_len = 10;
    size_t update_period;
    size_t n_min;
    size_t n_max;

    Progress(size_t in_n_min, size_t in_n_max, size_t in_update_period = 1) :
        update_period(in_update_period),
        n_min(in_n_min),
        n_max(in_n_max)
    {
        assert(update_period > 0);
    }

    Progress() :
        Progress(0, 0)
    {
    }

    void update(double n)
    {
        update(n, "");
    }

    void update(size_t n, std::string text)
    {
        assert(n >= n_min);
        assert(n <= n_max);

        // if update_period specified:
        // to not overload the console - update only at update_period
        if (update_period > 1) {
            static size_t next_update_in = 0;

            if (next_update_in == 0) {
                next_update_in = update_period;
            }
            else {
                next_update_in--;
                return;
            }
        }

        // generate a string first and then write the whole string to `os`
        // to prevent blinking cursor from jumping all over the place
        std::stringstream ss;

        ss << c_opening_bracket;
        size_t n_fill = n / (double)(n_max - n_min) * c_bar_len;
        for (size_t i = 0; i < c_bar_len; i++) {
            if (i < n_fill) {
                ss << c_fill;
            }
            else {
                ss << c_no_fill;
            }
        }

        const size_t n_max_strlen = std::to_string(n_max).length();

        ss << c_closing_bracket;
        ss << " " << std::setfill('0') << std::setw(n_max_strlen) << n;
        ss << "/" << n_max;
        ss << " " << std::fixed << std::setprecision(1)
           << (double)n / n_max * 100 << "%";
        double eta_s = get_eta_s(n, update_period);
        ss << " ETA " << seconds_to_hhmmss_string(eta_s);
        ss << text;
        // overwrite remalining command line with ' '
        const size_t n_chars = 5;
        for (size_t i = 0; i < n_chars; i++) {
            ss << " ";
        }
        ss << "\r";
        os << ss.str();
    }

    double get_eta_s(size_t n, const size_t update_period)
    {
        const std::chrono::time_point<std::chrono::steady_clock> now =
                std::chrono::steady_clock::now();
        const double us_per_update_period =
                std::chrono::duration_cast<std::chrono::microseconds>(
                        now - last_update_time)
                        .count();
        last_update_time = now;
        const size_t remaining_n = n_max - n;
        const double eta_s =
                us_per_update_period / update_period * remaining_n / 1000000;
        return eta_s;
    }

    void os_clear_line()
    {
        // move the cursor to the beginning of the current line and clear it
        os << "\r\033[K";
        os << std::flush;
    }
};

// base classes:
// - Settings
// - Context
//   - holds all calculation data
// - State
// - StateMachine(Context)
//   - init
//   - init loop
//   - run loop
//   - closure

struct Settings {
    size_t n_states = 1000000;
    double init_p_acceptance = 0.99;
    size_t init_t_log_len = 100;
    double cooling_rate = 0.95;
    size_t cooling_round_len = 1;
    size_t e_sma_fast_len = 50;
    size_t e_sma_slow_len = 200;

    size_t progress_update_period = 1;
    std::string log_filename{""};
    std::string stats_filename{"stats.txt"};
};

class State {
protected:
    Settings _settings;
    bool _energy_calculated;
    double _energy;

    void reset_energy()
    {
        _energy_calculated = false;
        _energy = -1;
    }

public:
    State(Settings &in_settings) :
        _settings(in_settings)
    {
        reset_energy();
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
        reset_energy();
        std::cout << "error: randomize method not implemented" << std::endl;
    }

    virtual void change()
    {
        reset_energy();

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
    bool cool = false;
    size_t cooling_i = 0;
    TState state;
    TState proposed_state;

    std::vector<double> init_t_log;
    bool init_done = false;

    // important! states begin with 1st, not 0th
    double t_max = 0;
    size_t state_i = 1;
    std::deque<double> e_log;
    size_t e_log_len;
    bool run_done = false;

    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> stop_time;
    double cycle_time_us = 0;

    Progress run_progress;
    std::ofstream log_f;

    Context(Settings &s) :
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
        const double state_s = state_i / runtime_s;
        const size_t first_col_width = 22;
        std::stringstream ss{};
        // standard parameters
        ss << std::left << std::setw(first_col_width) << "states" << state_i
           << std::endl;
        ss << std::endl;
        // runtime stats
        ss << std::left << std::setw(first_col_width) << "runtime"
           << seconds_to_hhmmss_string(runtime_s) << std::endl;
        ss << std::left << std::setw(first_col_width) << "state/s" << state_s
           << std::endl;
        ss << std::left << std::setw(first_col_width) << "result energy"
           << state.get_energy() << std::endl;
        return ss.str();
    }
};

template <typename TState>
class StateMachine {
private:
    typedef std::function<void(Context<TState> &)> state_function_t;
    Context<TState> context;

public:
    std::vector<state_function_t> init_functions{};
    std::vector<state_function_t> init_loop_functions{};
    std::vector<state_function_t> run_loop_functions{};
    std::vector<state_function_t> finalize_functions{};

    StateMachine(Settings &s) :
        context(s)
    {
    }

    void run()
    {
        context.start_time = std::chrono::steady_clock::now();
        for (state_function_t &f : init_functions) {
            f(context);
        }
        auto cycle_begin_time = std::chrono::steady_clock::now();
        while (init_loop_functions.size() > 0 && !context.init_done) {
            for (state_function_t &f : init_loop_functions) {
                if (context.init_done) {
                    break;
                }
                f(context);
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            context.cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        while (run_loop_functions.size() > 0 && !context.run_done) {
            for (state_function_t &f : run_loop_functions) {
                if (context.run_done) {
                    break;
                }
                f(context);
            }
            auto cycle_end_time = std::chrono::steady_clock::now();
            context.cycle_time_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                            cycle_end_time - cycle_begin_time)
                            .count();
            cycle_begin_time = cycle_end_time;
        }
        context.stop_time = std::chrono::steady_clock::now();
        for (state_function_t &f : finalize_functions) {
            f(context);
        }
    }
};

template <typename TState>
void init_log(Context<TState> &c)
{
    assert(!c.log_f.is_open());
    if (c.settings.log_filename.empty()) {
        return;
    }

    c.log_f.open(c.settings.log_filename);
    c.log_f.is_open();
    c.log_f << "state_i,temperature,energy" << std::endl;
}

template <typename TState>
void update_log(Context<TState> &c)
{
    // write the log to file
    assert(!c.settings.log_filename.empty());
    if (!c.log_f.is_open()) {
        return;
    }

    c.log_f << c.state_i;
    c.log_f << "," << c.temperature;
    c.log_f << "," << c.state.get_energy();
    c.log_f << std::endl;
    c.log_f << std::flush;
}

template <typename TState>
void init_run_progress(Context<TState> &c)
{
    if (c.temperature > 0) {
        c.run_progress.n_min = 1;
        c.run_progress.n_max = c.settings.n_states;
        c.run_progress.update_period = c.settings.progress_update_period;
    }
}

template <typename TState>
void print_run_progress(Context<TState> &c)
{
    const double state_s = 1 / c.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " e " << c.state.get_energy();
    ss << " t " << c.temperature;
    ss << " n/s " << state_s;

    c.run_progress.update(c.state_i, std::string(ss.str()));
}

template <typename TState>
void clear_run_progress(Context<TState> &c)
{
    c.run_progress.os_clear_line();
}

template <typename TState>
void print_stats(Context<TState> &c)
{
    std::cout << c.get_stats();
}

template <typename TState>
void create_stats_file(Context<TState> &c)
{
    if (c.settings.stats_filename.empty()) {
        return;
    }

    std::ofstream f(c.settings.stats_filename);
    f.is_open();
    f << c.get_stats();
}

template <typename TState>
void randomize_state(Context<TState> &c)
{
    assert(c.state_i == 1);
    c.state.randomize();
}

template <typename TState>
void propose_new_state(Context<TState> &c)
{
    c.proposed_state = c.state;
    c.proposed_state.change();
}

template <typename TState>
void record_init_temperature(Context<TState> &c)
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
    if (dE >= 0) {
        const double t = -dE / std::log(c.settings.init_p_acceptance);
        c.init_t_log.push_back(t);
    }
}

template <typename TState>
void select_init_temperature_as_max(Context<TState> &c)
{
    if (c.init_t_log.size() < c.settings.init_t_log_len) {
        return;
    }

    assert(c.init_t_log.size() == c.settings.init_t_log_len);
    assert(c.temperature == 0);
    assert(c.t_max == 0);

    // t_max = max(init_t at init_p_acceptance)
    c.t_max = *max_element(c.init_t_log.begin(), c.init_t_log.end());
    c.temperature = c.t_max;
}

template <typename TState>
void check_init_done(Context<TState> &c)
{
    if (c.temperature > 0) {
        assert(c.init_t_log.size() == c.settings.init_t_log_len);
        c.init_done = true;
    }
}

template <typename TState>
void update_state(Context<TState> &c)
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
void decide_to_cool(Context<TState> &c)
{
    assert(!c.cool);
    c.cool = true;
}

template <typename TState>
void cool_at_rate(Context<TState> &c)
{
    assert(c.settings.cooling_rate > 0);
    assert(c.settings.cooling_rate < 1);
    assert(c.temperature <= c.t_max);
    assert(c.temperature >= 0);
    if (c.cool) {
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
    c.cool = false;
}

template <typename TState>
void record_energy(Context<TState> &c)
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
void decide_to_cool_sma(Context<TState> &c)
{
    assert(!c.cool);
    if (c.e_log.size() < c.e_log_len) {
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
        c.cool = true;
    }
}

template <typename TState>
void check_run_done(Context<TState> &c)
{
    assert(c.temperature <= c.t_max);
    assert(c.temperature >= 0);
    assert(!c.run_done);
    if (c.state_i == c.settings.n_states) {
        c.run_done = true;
    }
}

}  // namespace lapsa
