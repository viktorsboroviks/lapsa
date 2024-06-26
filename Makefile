.PHONY: all examples format clean

all: examples

examples: cooling_schedule_max_array.o

cooling_schedule_max_array.o: examples/cooling_schedule_max_array.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		examples/cooling_schedule_max_array.cpp -o cooling_schedule_max_array.o

format: include/lapsa.hpp examples/cooling_schedule_max_array.cpp
	clang-format -i $^

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.csv"`
	rm -rf `find . -name "*.txt"`
