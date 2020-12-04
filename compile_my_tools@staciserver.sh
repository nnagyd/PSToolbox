#g++ -I/usr/local/include/eigen3 -pedantic -std=c++11 -c  my_tools.cpp -o my_tools.o
#ar rcs libmy_tools.a my_tools.o

g++ `pkg-config --cflags eigen3` -pedantic -std=c++11 -c  my_tools.cpp -o my_tools.o
ar rcs libmy_tools.a my_tools.o
mv my_tools.o staciserver/my_tools.o
