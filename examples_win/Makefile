CXX = g++

TARGETS = \
test1_reservoir_with_ideal_gas \
test2_reservoir_valve_with_liquid \
test3_reservoir_pipe_with_ideal_gas\
test4_reservoir_valve_with_ideal_gas\
test5_2pipes_with_liquid \
test6_3pipes_with_liquid \
test7_2pipes_with_ideal_gas \
test8_frozen_mixture

INC = -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -I/Users/hoscsaba/program/CoolProp/include -I/Users/hoscsaba/program/CoolProp/externals/fmtlib
INC_CP = -I/Users/hoscsaba/program/CoolProp/include -I/Users/hoscsaba/program/CoolProp/externals/fmtlib
LINK = -L/Users/hoscsaba/program/CoolProp/build1 -lPSToolbox -lmy_tools -lpython2.7 -lCoolProp
LINK_CP = -L/Users/hoscsaba/program/CoolProp/build1 -lCoolProp
CFLAGS = -std=c++11 -pedantic -O3 -Wall -Wno-c++11-long-long

all:$(TARGETS)

test1_reservoir_with_ideal_gas: test1_reservoir_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test1_reservoir_with_ideal_gas.cpp
	./a.out

test2_reservoir_valve_with_liquid: test2_reservoir_valve_with_liquid.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test2_reservoir_valve_with_liquid.cpp
	./a.out

test3_reservoir_pipe_with_ideal_gas: test3_reservoir_pipe_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test3_reservoir_pipe_with_ideal_gas.cpp
	./a.out

test4_reservoir_valve_with_ideal_gas: test4_reservoir_valve_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test4_reservoir_valve_with_ideal_gas.cpp
	./a.out

test5_2pipes_with_liquid: test5_2pipes_with_liquid.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test5_2pipes_with_liquid.cpp
	./a.out

test6_3pipes_with_liquid: test6_3pipes_with_liquid.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test6_3pipes_with_liquid.cpp
	./a.out

test7_2pipes_with_ideal_gas: test7_2pipes_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test7_2pipes_with_ideal_gas.cpp
	./a.out

test8_frozen_mixture: test8_frozen_mixture.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test8_frozen_mixture.cpp
	./a.out

test9_CoolProp_base: test9_CoolProp_base.cpp
	clear
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) test9_CoolProp_base.cpp
	./a.out

test10_reservoir_valve_with_hydrogen: test10_reservoir_valve_with_hydrogen.cpp
	clear
	$(CXX) $(INC) $(INC_CP) $(LINK) $(LINK_CP) $(CFLAGS) test10_reservoir_valve_with_hydrogen.cpp
	./a.out

clean:
	$(RM) *.o *.a
