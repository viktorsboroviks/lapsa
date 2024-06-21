//#include <sstream>
#include <vector>
#include <string>

#include "lapsa.hpp"

const std::string g_log_filename = "max_double_array_log.csv";

class MyState : lapsa::State<std::vector<double> >
{

};

int main()
{
    lapsa::Settings s{};
    s.log_filename = g_log_filename;

    lapsa::StateMachine<MyState> sm{s};
    sm.init_functions = {lapsa::init_log<MyState>};
    sm.run_loop_functions = {
            lapsa::update_state<MyState>,
            lapsa::update_log<MyState>,
            lapsa::print_progress<MyState>};
    sm.finalize_functions = {lapsa::print_stats<MyState>,
                            lapsa::create_stats_file<MyState>};
    sm.run();
    return 0;
}
