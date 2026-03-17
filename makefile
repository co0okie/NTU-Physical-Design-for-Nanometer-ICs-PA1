CC=g++
CFLAGS=-c -Wall
LDFLAGS=-std=c++11 -O3 -lm
SOURCES=src/main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fm
INCLUDES=
INPUT=testcase/input_0.5_5_9.dat
OUTPUT=$(subst input,output,$(INPUT))

all: bin/$(EXECUTABLE)

bin/$(EXECUTABLE): $(OBJECTS)
	mkdir -p $(@D)
	$(CC) $(LDFLAGS) $^ -o $@

%.o: %.cpp ${INCLUDES}
	$(CC) $(CFLAGS) $< -o $@

run: bin/$(EXECUTABLE)
	./$< '$(INPUT)' '$(OUTPUT)'

clean:
	rm -rf src/*.o bin