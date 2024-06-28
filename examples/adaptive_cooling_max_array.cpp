// #include <sstream>
#include <string>
#include <vector>

#include "lapsa.hpp"

size_t g_n_states = 1000000;
size_t g_progress_update_period = 100;
double g_init_p_acceptance = 0.97;
size_t g_init_t_log_len = 100;
double g_t_geom_k = (1 - 1e-4);
size_t g_e_sma_fast_len = 50;
size_t g_e_sma_slow_len = 100;

const size_t g_state_data_size = 100;
const std::string g_log_filename = "max_double_array_log.csv";

class MyState : lapsa::State<std::vector<double>> {
public:
    MyState(lapsa::Settings &in_settings) :
        State(in_settings)
    {
        data.resize(g_state_data_size);
    }

    double get_energy()
    {
        // if energy not calculated, do it now and store the result
        if (!_energy_calculated) {
            assert(data.size() != 0);
            _energy = 0;
            for (auto &d : data) {
                _energy -= d;
            }
        }
        _energy_calculated = true;
        return _energy;
    }

    void randomize(const std::function<double(void)> &rnd01)
    {
        assert(data.size() != 0);
        for (auto &d : data) {
            d = rnd01();
        }

        reset_energy();
    }

    void change(const std::function<double(void)> &rnd01)
    {
        assert(data.size() != 0);
        size_t changed_i = rnd01() * data.size();
        data[changed_i] = rnd01();

        reset_energy();
    }
};

int main()
{
    lapsa::Settings s{};
    s.n_states = g_n_states;
    s.progress_update_period = g_progress_update_period;
    s.init_p_acceptance = g_init_p_acceptance;
    s.init_t_log_len = g_init_t_log_len;
    s.t_geom_k = g_t_geom_k;
    s.log_filename = g_log_filename;
    s.e_sma_fast_len = g_e_sma_fast_len;
    s.e_sma_slow_len = g_e_sma_slow_len;

    lapsa::StateMachine<MyState> sm{s};
    sm.init_functions = {
            lapsa::init_log<MyState>,
            lapsa::init_state<MyState>,
    };
    sm.init_loop_functions = {
            lapsa::propose_new_state<MyState>,
            lapsa::record_init_temperature<MyState>,
            lapsa::select_init_temperature_as_max<MyState>,
            lapsa::init_run_progress<MyState>,
            lapsa::check_init_done<MyState>,
    };
    sm.run_loop_functions = {
            lapsa::propose_new_state<MyState>,
            lapsa::record_energy<MyState>,
            lapsa::update_temperature_with_geom_adaptive_cooling<MyState>,
            lapsa::update_state<MyState>,
            lapsa::check_run_done<MyState>,
            lapsa::update_log<MyState>,
            lapsa::print_run_progress<MyState>,
    };
    sm.finalize_functions = {
            lapsa::clear_run_progress<MyState>,
            lapsa::print_stats<MyState>,
            lapsa::create_stats_file<MyState>,
    };
    sm.run();
    return 0;
}
