INC := -I/usr/local/include/eigen3/
INC += -I./include/
CC := g++
FLAGS := -g -Wall -std=c++11 -fno-omit-frame-pointer -fsanitize=address
LFLAGS := -lfftw3

SOURCE := $(wildcard ./src/*.cc)
OBJS := $(patsubst %.cc, %.o, $(SOURCE))
TEST_SRC := $(wildcard unit_test/*.cc)
TEST_PREFIX := $(patsubst %.cc, % , $(TEST_SRC))
TEST_EXE := $(patsubst %.cc, %.out, $(TEST_SRC))

EXAMPLE_SRC := $(wildcard examples/*.cc)
EXAMPLE_EXE := $(patsubst %.cc, %.out, $(EXAMPLE_SRC))

.PHONY: clean test all
.SECONDARY: $(OBJS)

all: test example

test: $(TEST_EXE)

example: $(EXAMPLE_EXE)

./examples/%.out: ./examples/%.cc $(OBJS)
	$(CC) -o $@ $< \
		$(OBJS) \
		$(INC) \
		$(FLAGS) \
		$(LFLAGS) 

./unit_test/%.out: ./unit_test/%.cc $(OBJS)
	$(CC) -o $@ $< \
		$(OBJS) \
		$(INC) \
		$(FLAGS) \
		$(LFLAGS) 

./src/%.o: ./src/%.cc ./include/*.hpp
	$(CC) -o $@ -c $< \
		$(INC) $(FLAGS)

clean:
	@rm -rf $(OBJS) unit_test/*.out* examples/*.out
