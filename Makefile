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
	cd aviize; git checkout v2.3

iestaade:
	git clone git@github.com:viktorsboroviks/iestaade.git
	cd iestaade; git checkout v2.5

rododendrs:
	git clone git@github.com:viktorsboroviks/rododendrs.git
	cd rododendrs; git checkout v1.6

vasarniica:
	git clone git@github.com:viktorsboroviks/vasarniica.git
	cd vasarniica; git checkout v1.2

examples: \
	max_array_cooling_schedule.o \
	max_array_adaptive_cooling.o

max_array_cooling_schedule.o: \
		aviize \
		iestaade \
		rododendrs \
		examples/max_array_cooling_schedule.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./aviize/include \
		-I./iestaade/include \
		-I./rododendrs/include \
		examples/max_array_cooling_schedule.cpp -o $@

max_array_adaptive_cooling.o: \
		aviize \
		iestaade \
		rododendrs \
		examples/max_array_adaptive_cooling.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./aviize/include \
		-I./iestaade/include \
		-I./rododendrs/include \
		examples/max_array_adaptive_cooling.cpp -o $@

plot: examples vasarniica
	PYTHONPATH=$(PYTHONPATH):vasarniica/python:python python3 \
		examples/max_array_adaptive_cooling_plot.py \
		--config examples/max_array_adaptive_cooling_config.json
	PYTHONPATH=$(PYTHONPATH):vasarniica/python:python python3 \
		scripts/plot_lapsa.py \
		--config examples/max_array_adaptive_cooling_config.json

saves:
	mkdir saves

saves_dir := saves/${shell date -Iseconds}
save: saves
	mkdir -p "${saves_dir}"
	cp stats.txt "${saves_dir}"
	cp output.html "${saves_dir}"
	cp examples/max_array_adaptive_cooling_config.json "${saves_dir}"
	cp examples/max_array_adaptive_cooling.cpp "${saves_dir}"
	cp include/lapsa.hpp "${saves_dir}"

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
		-I./iestaade/include \
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
	m -rf iestaade
	m -rf rododendrs
	m -rf vasarniica
