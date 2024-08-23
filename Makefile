.PHONY: \
	all \
	examples \
	plot \
	format \
	format-cpp \
	format-json \
	format-python \
	lint \
	lint-cpp \
	lint-python \
	clean \
	distclean

all: examples

aviize:
	git clone git@github.com:viktorsboroviks/aviize.git
	cd iestade; git checkout v1.5

iestade:
	git clone git@github.com:viktorsboroviks/iestade.git
	cd iestade; git checkout v2.5

rododendrs:
	git clone git@github.com:viktorsboroviks/rododendrs.git
	cd rododendrs; git checkout v1.4

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
	PYTHONPATH=${PYTHONPATH}:python python3 \
		scripts/plot_log.py \
		--config examples/max_array_adaptive_cooling_config.json

format: format-cpp format-json format-python

format-cpp: \
		include/lapsa.hpp \
		examples/max_array_cooling_schedule.cpp \
		examples/max_array_adaptive_cooling.cpp
	clang-format -i $^

format-json: \
		config.json \
		scripts/plot_log_config.json \
		examples/max_array_cooling_schedule_config.json \
		examples/max_array_adaptive_cooling_config.json
	for json_file in $^; do \
		jq . $$json_file | sponge $$json_file; \
    done

format-python: \
		python/lapsa.py \
		scripts/plot_log.py \
		examples/max_array_adaptive_cooling_plot.py
	black $^

lint: lint-cpp lint-python

lint-cpp: \
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

lint-python: \
		python/lapsa.py \
		scripts/plot_log.py \
		examples/max_array_adaptive_cooling_plot.py
	pylint $^
	flake8 $^

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.csv"`
	rm -rf `find . -name "*.txt"`
	rm -rf `find . -name "*.svg"`
	rm -rf `find . -name "*.html"`

distclean: clean
	m -rf aviize
	m -rf iestade
	m -rf rododendrs
