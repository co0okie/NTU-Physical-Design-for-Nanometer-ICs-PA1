CXX=g++
CXXFLAGS=-std=c++11 -O3 -Wall
LDFLAGS=-lm
SOURCES=src/main.cpp src/partition.cpp src/clusterning.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fm
INCLUDES=src/partition.h
# INPUT=testcase/input_0.dat
# INPUT=testcase/input_1.dat
# INPUT=testcase/input_2.dat
# INPUT=testcase/input_3.dat
# INPUT=testcase/input_4.dat
INPUT=testcase/input_5.dat
# INPUT=testcase/input_0.2_20_30.dat
# INPUT=testcase/input_0.5_5_9.dat
# INPUT=testcase/input_0.3_20_30.dat
OUTPUT=$(subst input,output,$(INPUT))
SHELL := /bin/bash

.PHONY: clean run run_all test evaluate tgz

all: bin/$(EXECUTABLE)

bin/$(EXECUTABLE): $(OBJECTS)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp ${INCLUDES}
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: bin/$(EXECUTABLE)
	./$< '$(INPUT)' '$(OUTPUT)'

run_all: bin/$(EXECUTABLE)
	for i in 0 1 2 3 4 5; do \
		{ T[$$i]=$$( { TIMEFORMAT="%R"; time ./$< testcase/input_$$i.dat testcase/output_$$i.dat 1>&3 2>&4; } 2>&1 ); } 3>&1 4>&2; \
	done; \
	for i in 0 1 2 3 4 5; do \
		echo; \
		echo "## input_$$i.dat"; \
		echo; \
		echo '```'; \
		evaluator/evaluator.sh testcase/input_$$i.dat testcase/output_$$i.dat $${T[$$i]} \
		| tee testcase/result_$$i.txt; \
		echo '```'; \
	done; \
	rm -f time_*.tmp

evaluate:
	for i in 0 1 2 3 4 5; do \
		echo; \
		echo "## input_$$i.dat"; \
		echo; \
		echo '```'; \
		evaluator/evaluator.sh testcase/input_$$i.dat testcase/output_$$i.dat 1 \
		| tee testcase/result_$$i.txt; \
		echo '```'; \
	done; \

bin/test: src/test.o $(filter-out src/main.o, $(OBJECTS))
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

test: bin/test
	./$<

tgz: b11107051_pa1.tgz

%.tgz: bin/fm src/*.cpp src/*.h makefile readme.txt report/report.pdf
	temp_dir=$$(mktemp -d); \
	mkdir -p "$$temp_dir"/$*/src "$$temp_dir"/$*/bin; \
	cp src/*.cpp src/*.h "$$temp_dir"/$*/src/; \
	cp bin/fm "$$temp_dir"/$*/bin/fm; \
	cp makefile readme.txt report/report.pdf "$$temp_dir"/$*; \
	tar zcvf $@ -C "$$temp_dir" $*; \
	rm -rf "$$temp_dir"

clean:
	rm -rf src/*.o testcase/output_*.dat testcase/result_*.txt bin log a.out b11107051_pa1.tgz report/report.aux report/report.fdb_latexmk report/report.fls report/report.log report/report.out report/report.synctex.gz report/report.xdv