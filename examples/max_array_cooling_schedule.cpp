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
            _energy = -std::accumulate(_data.begin(), _data.end(), 0.0);
        }
        _evaluated = true;
    }

public:
    MyState() :
        _data(iestaade::from_json<size_t>(CONFIG_PATH, "state/data_size"))
    {
    }

    double get_energy() override
    {
        _evaluate();
        return _energy;
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
    lapsa::StateMachine<MyState> sm(CONFIG_PATH, "lapsa");
    // clang-format off
    sm.init_functions = {
            lapsa::init_log                     <MyState>,
            lapsa::randomize_state              <MyState>,
    };
    sm.init_loop_functions = {
            // update state
            lapsa::run_i_inc                    <MyState>,
            lapsa::propose_new_state            <MyState>,
            lapsa::generate_init_t_candidates   <MyState>,
            lapsa::select_init_t_as_max         <MyState>,
            // decide to proceed
            lapsa::decide_init_done             <MyState>,
    };
    sm.run_functions = {
            lapsa::progress_init_run_loop       <MyState>,
    };
    sm.run_loop_functions = {
            // log
            lapsa::update_log                   <MyState>,
            // progress bar
            lapsa::progress_text_reset          <MyState>,
            lapsa::progress_text_add_total      <MyState>,
            lapsa::progress_text_add_pct        <MyState>,
            lapsa::progress_text_add_eta        <MyState>,
            lapsa::progress_text_add_t          <MyState>,
            lapsa::progress_text_add_e          <MyState>,
            lapsa::progress_print               <MyState>,
            // decide to proceed and begin next run
            lapsa::decide_run_done              <MyState>,
            // update state
            lapsa::propose_new_state            <MyState>,
            lapsa::do_cool_set                  <MyState>,
            lapsa::cool_at_rate                 <MyState>,
            lapsa::run_i_inc                    <MyState>,
            lapsa::update_state                 <MyState>,
            lapsa::reset_state_evaluation       <MyState>,
    };
    sm.finalize_functions = {
            lapsa::progress_text_reset          <MyState>,
            lapsa::print_stats                  <MyState>,
            lapsa::create_stats_file            <MyState>,
    };
    // clang-format on
    sm.run();
    return 0;
}
