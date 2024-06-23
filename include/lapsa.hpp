#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <ostream>
#include <random>
#include <vector>

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

class Random {
    // thread-safe implementation of rnd01() borrowed from
    // https://github.com/Arash-codedev/openGA/blob/master/README.md
    // assuming those people knew what they were doing
private:
    std::mutex mtx_rand;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> unif_dist;

public:
    Random()
    {
        // initialize the random number generator with time-dependent seed
        uint64_t timeSeed = std::chrono::high_resolution_clock::now()
                                    .time_since_epoch()
                                    .count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff),
                         uint32_t(timeSeed >> 32)};
        rng.seed(ss);
        std::uniform_real_distribution<double> unif(0, 1);
    }

    Random(const Random &other)
    {
        // construct a new object on copy
        (void)other;
        Random();
    }

    double rnd01(void)
    {
        // prevent data race between threads
        std::lock_guard<std::mutex> lock(mtx_rand);
        return unif_dist(rng);
    }
};

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
    double n_min;
    double n_max;

    Progress(double in_n_min, double in_n_max, size_t in_update_period = 0) :
        update_period(in_update_period),
        n_min(in_n_min),
        n_max(in_n_max)
    {
    }

    ~Progress()
    {
        os_clean();
    }

    void update(double n)
    {
        update(n, "");
    }

    void update(double n, std::string text)
    {
        assert(n >= n_min);
        assert(n <= n_max);

        // if update_period specified:
        // to not overload the console - update only at update_period
        if (update_period) {
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
        size_t n_fill = n / (n_max - n_min) * c_bar_len;
        for (size_t i = 0; i < c_bar_len; i++) {
            if (i < n_fill) {
                ss << c_fill;
            }
            else {
                ss << c_no_fill;
            }
        }

        ss << c_closing_bracket;
        // TODO: add ETA
        ss << text;
        // overwrite remalining command line with ' '
        const size_t n_chars = 5;
        for (size_t i = 0; i < n_chars; i++) {
            ss << " ";
        }
        ss << "\r";
        os << ss.str();

        last_update_time = std::chrono::steady_clock::now();
    }

    void os_clean(size_t n_chars = 200)
    {
        // overwrite command line with n_chars ' '
        for (size_t i = 0; i < n_chars; i++) {
            os << " ";
        }
        os << "\r" << std::flush;
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
    size_t progress_update_period = 0;
    std::string log_filename{""};
    std::string stats_filename{"stats.txt"};
};

template <typename TData>
class State {
protected:
    bool _energy_calculated = false;
    double _energy = -1;
    Settings _settings;

public:
    TData data;

    State(Settings &in_settings) :
        _settings(in_settings)
    {
    }

    // virtual destructor is required if virtual methods are used
    virtual ~State() {}

    virtual double get_energy()
    {
        std::cout << "error: get_energy method not implemented" << std::endl;
        return -1.0;
    }

    virtual void randomize(const std::function<double(void)> &rnd01)
    {
        (void)rnd01();
        std::cout << "error: randomize method not implemented" << std::endl;
    }

    virtual void change(const std::function<double(void)> &rnd01)
    {
        // this method is very individual-specific, so to not overthink it
        // I leave it virtual
        (void)rnd01;
        std::cout << "error: change method not implemented" << std::endl;
    }
};

template <typename TState>
class Context {
public:
    Settings settings;
    Random random{};

    double temperature = 0;
    TState state;
    TState proposed_state;

    // important! states begin with 1st, not 0th
    double init_p_acceptance = 0.97;
    size_t init_state_i = 1;
    size_t init_state_imax = 100;
    std::vector<double> init_temperatures;
    bool init_done = false;

    // important! states begin with 1st, not 0th
    size_t state_i = 1;
    size_t state_imax = 1000000;
    bool run_done = false;

    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> stop_time;
    double cycle_time_us = 0;

    Progress progress;
    std::ofstream log_f;

    Context(Settings &s) :
        settings(s),
        // TODO: consider setting min/max temperature from real data
        state(settings),
        proposed_state(settings),
        progress(0, 100, s.progress_update_period)
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
    c.log_f << "state_i,energy" << std::endl;
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
void print_run_progress(Context<TState> &c)
{
    const double state_s = 1 / c.cycle_time_us * 1000000;

    std::stringstream ss;
    ss << " t " << c.temperature;
    ss << " e " << c.state.get_energy();
    ss << " state/s " << state_s;

    c.progress.update(c.temperature, std::string(ss.str()));
}

template <typename TState>
void update_state(Context<TState> &c)
{
    // TODO: review and change this function
    if (c.state_i > c.state_imax) {
        c.run_done = true;
        return;
    }

    c.state_i++;
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
    c.state.randomize([&c]() { return c.random.rnd01(); });
}

template <typename TState>
void propose_changed_state(Context<TState> &c)
{
    c.proposed_state = c.state;
    c.proposed_state.change([&c]() { return c.random.rnd01(); });
}

template <typename TState>
void init_temperature(Context<TState> &c)
{
    if (c.init_state_i <= c.init_state_imax) {
        // calculate intermediate init temperatures
        // - T = -(E_proposed - E) / ln(P)
        // - from Acceptance probability
        //   - if E_new < E_current: P = 1
        //   - else: P = exp(-(E_proposed - E) / T)
        //   - ref: https://en.wikipedia.org/wiki/Simulated_annealing
        // - only for cases where new state is worse and we need
        //   to use the p_acceptance
        // - smaller energy = better
        const double dE = c.proposed_state.get_energy() - c.state.get_energy();
        if (dE > 0) {
            const double t = -dE / std::log(c.init_p_acceptance);
            c.init_temperatures.push_back(t);
            c.init_state_i++;
        }
    }
    else {
        // calculate resulting init temperature
        assert(c.init_temperatures.size() == c.init_state_imax);
        c.temperature = *max_element(c.init_temperatures.begin(),
                                     c.init_temperatures.end());
        c.init_done = true;
    }
}

}  // namespace lapsa
