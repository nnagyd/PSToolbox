#include <iostream>
#include "../../CoolProp/include/CoolProp.h"

 //g++ -std=c++11 -Wall -O2 -o test9_CoolProp_base -I/usr/local/CoolProp/6.1.0/CoolProp/include -L/usr/local/CoolProp/6.1.0/CoolProp/build test9_CoolProp_base.cpp -lCoolProp

int main(){
   double T{293.15};
   double P{1e5};
   double res{CoolProp::PropsSI("D", "T", T, "P", P, "Water")};
   std::cout << "Density: " << res << std::endl;
}
