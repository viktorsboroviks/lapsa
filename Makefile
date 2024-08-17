.PHONY: \
	all \
	examples \
	plot \
	format \
	clang-format \
	jq-format \
	black-format \
	lint \
	cppcheck-lint \
	clean \
	distclean

all: examples

aviize:
	git clone git@github.com:viktorsboroviks/aviize.git
	cd iestade; git checkout v1.3

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
		aviize \
		iestade \
		rododendrs \
		examples/max_array_cooling_schedule.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./aviize/include \
		-I./iestade/include \
		-I./rododendrs/include \
		examples/max_array_cooling_schedule.cpp -o $@

max_array_adaptive_cooling.o: \
		aviize \
		iestade \
		rododendrs \
		examples/max_array_adaptive_cooling.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./aviize/include \
		-I./iestade/include \
		-I./rododendrs/include \
		examples/max_array_adaptive_cooling.cpp -o $@

plot: examples
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

lint: cppcheck-lint

cppcheck-lint: \
		include/lapsa.hpp \
		examples/max_array_cooling_schedule.cpp \
		examples/max_array_adaptive_cooling.cpp
	cppcheck \
		--enable=warning,portability,performance \
		--enable=style,information \
		--enable=missingInclude \
		--inconclusive \
		--library=std,posix,gnu \
		--platform=unix64 \
		--language=c++ \
		--std=c++20 \
		--inline-suppr \
		--check-level=exhaustive \
		--suppress=missingIncludeSystem \
		--suppress=checkersReport \
		--checkers-report=cppcheck_report.txt \
		-I./include \
		-I./aviize/include \
		-I./iestade/include \
		-I./rododendrs/include \
		$^

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.csv"`
	rm -rf `find . -name "*.txt"`

distclean: clean
	m -rf aviize
	m -rf iestade
	m -rf rododendrs
