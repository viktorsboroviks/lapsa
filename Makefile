.PHONY: all examples format clean

all: examples

examples: max_double_array.o

max_double_array.o: examples/max_double_array.cpp
	g++ -Wall -Wextra -Werror -Wpedantic \
		-std=c++20 -O3 \
		-I./include \
		examples/max_double_array.cpp -o max_double_array.o

format: include/lapsa.hpp examples/max_double_array.cpp
	clang-format -i $^

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.csv"`
	rm -rf `find . -name "*.txt"`
