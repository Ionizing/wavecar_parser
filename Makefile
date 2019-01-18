INC := -I./include/
CC := c++
FLAGS := -Wextra -Wall -std=c++11 -O2 #-g -fno-omit-frame-pointer -fsanitize=address
LFLAGS := -lfftw3

SOURCE := $(wildcard ./src/*.cc)
OBJS := $(patsubst %.cc, %.o, $(SOURCE))

.PHONY: help clean test all veryclean
.SECONDARY: $(OBJS)

all: plotWave.out

test: plotWave.out
	./plotWave.out --poscar=./tests/MoS2/POSCAR --wavecar=./tests/MoS2/WAVECAR --prefix=MoS2
	./plotWave.out --poscar=./tests/fccSidos/POSCAR --wavecar=./tests/fccSidos/WAVECAR --prefix=Si

plotWave.out: $(OBJS)
	$(CC) -o $@ $(OBJS) $(FLAGS) $(LFLAGS)

./src/%.o: ./src/%.cc ./include/*.hpp
	$(CC) -o $@ -c $< \
		$(INC) $(FLAGS)

clean:
	@rm -rf $(OBJS) *.out

veryclean: clean
	@rm *.vasp
