.PHONY: \
	all \
	examples \
	plot \
	format \
	clang-format \
	jq-format \
	black-format \
	clean \
	distclean

all: examples

iestade:
	git clone git@github.com:viktorsboroviks/iestade.git
	cd iestade; git checkout v2.3

rododendrs:
	git clone git@github.com:viktorsboroviks/rododendrs.git
	cd rododendrs; git checkout v1.1

examples: \
	max_array_cooling_schedule.o \
	max_array_adaptive_cooling.o

max_array_cooling_schedule.o: \
		iestade \
		rododendrs \
		examples/max_array_cooling_schedule.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./iestade/include \
		-I./rododendrs/include \
		examples/max_array_cooling_schedule.cpp -o $@

max_array_adaptive_cooling.o: \
		iestade \
		rododendrs \
		examples/max_array_adaptive_cooling.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./iestade/include \
		-I./rododendrs/include \
		examples/max_array_adaptive_cooling.cpp -o $@

plot_examples: examples
	PYTHONPATH=${PYTHONPATH}:python python3 \
		examples/max_array_adaptive_cooling_plot.py \
		--config examples/max_array_adaptive_cooling_config.json

format: clang-format jq-format black-format

clang-format: \
		include/lapsa.hpp \
		examples/max_array_cooling_schedule.cpp \
		examples/max_array_adaptive_cooling.cpp
	clang-format -i $^

jq-format: \
		config.json \
		examples/max_array_cooling_schedule_config.json \
		examples/max_array_adaptive_cooling_config.json
	jq . config.json | \
		sponge config.json
	jq . examples/max_array_cooling_schedule_config.json | \
		sponge examples/max_array_cooling_schedule_config.json
	jq . examples/max_array_adaptive_cooling_config.json | \
		sponge examples/max_array_adaptive_cooling_config.json

black-format: \
		python/lapsa.py \
		examples/max_array_adaptive_cooling_plot.py
	black $^

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.csv"`
	rm -rf `find . -name "*.txt"`

distclean: clean
	m -rf iestade
	m -rf rododendrs
