# https://devhints.io/makefile
CC = g++
CCFLAGS = -std=c++17 -Wall -Wextra -Wno-sign-compare -pedantic -fPIC -O2 # run mode
# CCFLAGS = -std=c++17 -Wall -Wextra -Wno-sign-compare -pedantic -O0 -g # debug mode
GINAC = -lginac -lcln -lgtest -lgtest_main
BUILD_DIR = build
SRC_DIR = src
TEST_DIR = test

dependencies := $(wildcard src/*.cpp)
build_dependencies := $(dependencies:src/%.cpp=build/%.o)

all: main

clean:
	rm -rf $(BUILD_DIR)

prepare:
	mkdir -p $(BUILD_DIR)

reformat:
	clang-format -i -style=Google main.cpp $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h $(TEST_DIR)/*.cpp

main: clean reformat prepare run_main # for debug mode add "gdb" on the beginning
	$(BUILD_DIR)/main.out 

test: clean reformat prepare run_test # for debug mode add "gdb" on the beginning
	$(BUILD_DIR)/test.out

run_main: build/main.o $(build_dependencies)
	$(CC) $(CCFLAGS) $^ -I$(SRC_DIR) -o $(BUILD_DIR)/main.out $(GINAC)

run_test: build/test.o $(build_dependencies)
	$(CC) $(CCFLAGS) $^ -I$(SRC_DIR) -o $(BUILD_DIR)/test.out $(GINAC)

build/%.o: %.cpp
	$(CC) $(CCFLAGS) -I$(SRC_DIR) -c $^ -o $@

build/%.o: $(TEST_DIR)/%.cpp $(build_dependencies)
	$(CC) -I$(SRC_DIR) $(CCFLAGS) -c $< -o $@

build/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CCFLAGS) -c $< -o $@
