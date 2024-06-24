// #include <sstream>
#include <string>
#include <vector>

#include "lapsa.hpp"

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
                _energy += d;
            }
        }

        return _energy;
    }

    void randomize(const std::function<double(void)> &rnd01)
    {
        assert(data.size() != 0);
        for (auto &d : data) {
            d = rnd01();
        }
    }

    void change(const std::function<double(void)> &rnd01)
    {
        assert(data.size() != 0);
        size_t changed_i = rnd01() * data.size();
        data[changed_i] = rnd01();
    }
};

int main()
{
    lapsa::Settings s{};
    s.log_filename = g_log_filename;

    lapsa::StateMachine<MyState> sm{s};
    sm.init_functions = {lapsa::init_log<MyState>,
                         lapsa::randomize_state<MyState>};
    sm.init_loop_functions = {
            lapsa::propose_changed_state<MyState>,
            lapsa::init_temperature<MyState>,
            // TODO: print_init_progress
    };
    sm.run_loop_functions = {lapsa::update_state<MyState>,
                             lapsa::propose_changed_state<MyState>,
                             // TODO: accept_state
                             // TODO: update_temperature (adaptive cooling)
// - calc sma(dE)
// - if < min -> dec(t)
                             // - if T==0 -> set return flag
                             // TODO: update_records
                             lapsa::update_log<MyState>,
                             lapsa::print_run_progress<MyState>};
    sm.finalize_functions = {lapsa::print_stats<MyState>,
                             lapsa::create_stats_file<MyState>};
    sm.run();
    return 0;
}
