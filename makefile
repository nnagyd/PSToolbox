# Compiler
CXX := g++

# Compiler flags
CFLAGS = -g -Wall -std=c++11 -pedantic
INC = -I/usr/local/include/eigen3
LINK = -lmy_tools

INC_CP = -I/Users/hoscsaba/program/CoolProp/include -I/Users/hoscsaba/program/CoolProp/externals/fmtlib
LINK_CP = -L/Users/hoscsaba/program/CoolProp/build1 -lCoolProp

# Source files
SRCS := $(wildcard *.cpp)

# Object files
OBJS := $(SRCS:.cpp=.o)

# Executable
EXEC := output

# Build rule
$(EXEC): $(OBJS)
	libtool -static -o libmy_tools.a my_tools.o
	libtool -static -o libPSToolbox.a PSToolBoxRunner.o PSToolboxBaseEdge.o Gas.o IdealGas.o FrozenMixtureLiquidGas.o Units.o LWP.o SCP.o Reservoir.o Valve.o Connector.o CoolPropGas.o CoolPropHA.o Valve_with_Absorber.o
#$(CXX) $(CXXFLAGS) $^ -o $@

# Compile rule for each source file
my_tools.o: my_tools.cpp 
	$(CXX) $(INC) $(CFLAGS) my_tools.cpp -c -o my_tools.o

PSToolboxRunner.o: PSToolboxRunner.cpp 
	$(CXX) $(INC) $(CFLAGS) PSToolboxRunner.cpp -c -o PSToolboxRunner.o

PSToolboxBaseEdge.o: PSToolboxBaseEdge.cpp 
	$(CXX) $(INC) $(CFLAGS) PSToolboxBaseEdge.cpp -c -o PSToolboxBaseEdge.o

Gas.o: Gas.cpp 
	$(CXX) $(CFLAGS) Gas.cpp -c -o Gas.o

IdealGas.o: IdealGas.cpp 
	$(CXX) $(CFLAGS) IdealGas.cpp -c -o IdealGas.o

FrozenMixtureLiquidGas.o: FrozenMixtureLiquidGas.cpp 
	$(CXX) $(CFLAGS) FrozenMixtureLiquidGas.cpp -c -o FrozenMixtureLiquidGas.o

Units.o: Units.cpp
	$(CXX) $(CFLAGS) Units.cpp -c -o Units.o

LWP.o: LWP.cpp
	$(CXX) $(INC) $(CFLAGS) LWP.cpp -c -o LWP.o

SCP.o: SCP.cpp
	$(CXX) $(INC) $(CFLAGS) SCP.cpp -c -o SCP.o

Reservoir.o: Reservoir.cpp
	$(CXX) $(INC) $(CFLAGS) Reservoir.cpp -c -o Reservoir.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Reservoir.cpp -c -o Reservoir.o

Valve.o: Valve.cpp
	$(CXX) $(INC) $(CFLAGS) Valve.cpp -c -o Valve.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Valve.cpp -c -o Valve.o

Valve_with_Absorber.o: Valve_with_absorber.cpp
	$(CXX) $(INC) $(CFLAGS) Valve_with_Absorber.cpp -c -o Valve_with_Absorber.o

Connector.o: Connector.cpp
	$(CXX) $(INC) $(CFLAGS) Connector.cpp -c -o Connector.o

CoolPropGas.o: CoolPropGas.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) CoolPropGas.cpp -c -o CoolPropGas.o

CoolPropMixture.o: CoolPropMixture.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) CoolPropMixture.cpp -c -o CoolPropMixture.o

CoolPropHA.o: CoolPropHA.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) CoolPropHA.cpp -c -o CoolPropHA.o

# Phony target to clean the project
.PHONY: clean
clean:
	rm -f $(EXEC) $(OBJS)
