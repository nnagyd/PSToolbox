################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../SCP.cpp \
../my_tools.cpp 

O_SRCS += \
../my_tools.o 

OBJS += \
./SCP.o \
./my_tools.o 

CPP_DEPS += \
./SCP.d \
./my_tools.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O2 -g -pedantic -Wall -c -fmessage-length=0 -I/usr/local/include/eigen3 -v -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


