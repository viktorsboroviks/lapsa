.PHONY: all examples format clean distclean

all: examples

rododendrs:
	git clone git@github.com:viktorsboroviks/rododendrs.git
	cd rododendrs; git checkout v1.1

examples: \
	cooling_schedule_max_array.o \
	adaptive_cooling_max_array.o

cooling_schedule_max_array.o: rododendrs examples/cooling_schedule_max_array.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./rododendrs/include \
		examples/cooling_schedule_max_array.cpp -o $@

adaptive_cooling_max_array.o: rododendrs examples/adaptive_cooling_max_array.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		-I./rododendrs/include \
		examples/adaptive_cooling_max_array.cpp -o $@

format: \
		include/lapsa.hpp \
		examples/cooling_schedule_max_array.cpp \
		examples/adaptive_cooling_max_array.cpp
	clang-format -i $^

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.csv"`
	rm -rf `find . -name "*.txt"`

distclean: clean
	rm -rf lapsa
	rm -rf grafiins
	rm -rf rododendrs

