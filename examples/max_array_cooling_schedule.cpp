#include <string>
#include <vector>

#include "iestaade.hpp"
#include "lapsa.hpp"
#include "rododendrs.hpp"

const std::string CONFIG_PATH =
        "examples/max_array_cooling_schedule_config.json";

// cppcheck-suppress ctuOneDefinitionRuleViolation
class MyState : public lapsa::State {
private:
    std::vector<double> _data;

    void _evaluate()
    {
        // if energy not calculated, do it now and store the result
        if (!_evaluated) {
            assert(_data.size() != 0);
            _value  = std::accumulate(_data.begin(), _data.end(), 0.0);
            _energy = 1.0 / _value;
        }
        _evaluated = true;
    }

public:
    explicit MyState(lapsa::Settings &in_settings) :
        State(in_settings),
        _data(iestaade::size_t_from_json(CONFIG_PATH, "state/data_size"))
    {
    }

    double get_energy() override
    {
        _evaluate();
        return _energy;
    }

    double get_value() override
    {
        _evaluate();
        return _value;
    }

    void randomize() override
    {
        assert(_data.size() != 0);
        std::generate(_data.begin(), _data.end(), []() {
            return rododendrs::rnd01();
        });
    }

    void change() override
    {
        assert(_data.size() != 0);
        size_t changed_i = rododendrs::rnd01() * _data.size();
        _data[changed_i] = rododendrs::rnd01();
    }
};

int main()
{
    lapsa::Settings s(CONFIG_PATH, "lapsa");
    lapsa::StateMachine<MyState> sm{s};
    sm.init_functions = {
            lapsa::init_log<MyState>,
            lapsa::randomize_state<MyState>,
    };
    sm.init_loop_functions = {
            // update state
            lapsa::propose_new_state<MyState>,
            lapsa::init_t_history<MyState>,
            lapsa::init_t_select_max<MyState>,
            lapsa::init_progress<MyState>,
            // decide to proceed
            lapsa::decide_init_done<MyState>,
    };
    sm.run_loop_functions = {
            // log
            lapsa::update_log<MyState>,
            // progress bar
            lapsa::progress_text_reset<MyState>,
            lapsa::progress_text_add_total<MyState>,
            lapsa::progress_text_add_pct<MyState>,
            lapsa::progress_text_add_eta<MyState>,
            lapsa::progress_text_add_t<MyState>,
            lapsa::progress_text_add_e<MyState>,
            lapsa::progress_text_add_v<MyState>,
            lapsa::progress_print<MyState>,
            // decide to proceed and begin next run
            lapsa::decide_run_done<MyState>,
            // update state
            lapsa::propose_new_state<MyState>,
            lapsa::do_cool_set<MyState>,
            lapsa::cool_at_rate<MyState>,
            lapsa::update_state<MyState>,
    };
    sm.finalize_functions = {
            lapsa::progress_text_reset<MyState>,
            lapsa::stats_print<MyState>,
            lapsa::stats_create_file<MyState>,
    };
    sm.run();
    return 0;
}
