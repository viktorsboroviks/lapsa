#include <string>
#include <vector>

#include "iestade.hpp"
#include "lapsa.hpp"
#include "rododendrs.hpp"

const std::string CONFIG_PATH =
        "examples/max_array_adaptive_cooling_config.json";

class MyState : lapsa::State {
private:
    std::vector<double> _data;

public:
    explicit MyState(lapsa::Settings& in_settings) :
        State(in_settings),
        _data(iestade::size_t_from_json(CONFIG_PATH, "state/data_size"))
    {
    }

    double get_energy() override
    {
        // if energy not calculated, do it now and store the result
        if (!_energy_calculated) {
            assert(_data.size() != 0);
            _energy = std::accumulate(_data.begin(), _data.end(), 0.0);
        }
        _energy_calculated = true;
        return 1.0 / _energy;
    }

    void randomize() override
    {
        assert(_data.size() != 0);
        std::generate(_data.begin(), _data.end(), []() {
            return rododendrs::rnd01();
        });

        reset_energy();
    }

    void change() override
    {
        assert(_data.size() != 0);
        size_t changed_i = rododendrs::rnd01() * _data.size();
        _data[changed_i] = rododendrs::rnd01();

        reset_energy();
    }
};

int main()
{
    lapsa::Settings s(CONFIG_PATH, "lapsa");
    lapsa::StateMachine<MyState> sm{s};
    sm.init_functions = {
            lapsa::log_init<MyState>,
            lapsa::state_randomize<MyState>,
    };
    sm.init_loop_functions = {
            lapsa::state_propose_new<MyState>,
            lapsa::temperature_init_record<MyState>,
            lapsa::temperature_init_select_as_max<MyState>,
            lapsa::run_progress_init<MyState>,
            lapsa::init_done_decide<MyState>,
    };
    sm.run_loop_functions = {
            // log and report
            lapsa::log_update<MyState>,
            lapsa::run_progress_text_reset<MyState>,
            lapsa::run_progress_text_add_total<MyState>,
            lapsa::run_progress_text_add_pct<MyState>,
            lapsa::run_progress_text_add_eta<MyState>,
            lapsa::run_progress_text_add_t<MyState>,
            lapsa::run_progress_text_add_e<MyState>,
            lapsa::run_progress_print<MyState>,
            // decide to proceed
            lapsa::run_done_decide<MyState>,
            // proceed
            lapsa::state_propose_new<MyState>,
            lapsa::log_energy<MyState>,
            lapsa::do_cool_decide_sma<MyState>,
            lapsa::cool_at_rate<MyState>,
            lapsa::state_update<MyState>,
    };
    sm.finalize_functions = {
            lapsa::run_progress_clear<MyState>,
            lapsa::stats_print<MyState>,
            lapsa::stats_create_file<MyState>,
    };
    sm.run();
    return 0;
}
