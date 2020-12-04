#g++ -I/usr/local/include/eigen3 -pedantic  -std=c++11 -c  SCP.cpp -o SCP.o
#ar rcs libSCP.a SCP.o

g++ -std=c++11 `pkg-config --cflags eigen3` -pedantic  -std=c++11 -c  SCP.cpp -o SCP.o
ar rcs libSCP.a SCP.o
