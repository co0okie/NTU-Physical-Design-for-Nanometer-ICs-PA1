CC=g++
CFLAGS=-c -Wall
LDFLAGS=-std=c++11 -O3 -lm
SOURCES=src/main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fm
INCLUDES=
INPUT=testcase/input_0.dat
# INPUT=testcase/input_1.dat
# INPUT=testcase/input_2.dat
# INPUT=testcase/input_3.dat
# INPUT=testcase/input_4.dat
# INPUT=testcase/input_5.dat
# INPUT=testcase/input_0.3_20_30.dat
OUTPUT=$(subst input,output,$(INPUT))
SHELL := /bin/bash

all: bin/$(EXECUTABLE)

bin/$(EXECUTABLE): $(OBJECTS)
	mkdir -p $(@D)
	$(CC) $(LDFLAGS) $^ -o $@

%.o: %.cpp ${INCLUDES}
	$(CC) $(CFLAGS) $< -o $@

run: bin/$(EXECUTABLE)
	./$< '$(INPUT)' '$(OUTPUT)'

run_all: bin/$(EXECUTABLE)
	for i in 0 1 2 3 4 5; do \
		{ T[$$i]=$$( { TIMEFORMAT="%R"; time ./$< testcase/input_$$i.dat testcase/output_$$i.dat 1>&3 2>&4; } 2>&1 ); } 3>&1 4>&2; \
	done; \
	for i in 0 1 2 3 4 5; do \
		evaluator/evaluator.sh testcase/input_$$i.dat testcase/output_$$i.dat $${T[$$i]}; \
	done; \
	rm -f time_*.tmp

clean:
	rm -rf src/*.o bin