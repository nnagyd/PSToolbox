#include <iostream>
#include <string>
using namespace std;

// Base class
class Vehicle {
  public: 
    string brand = "Ford";
    virtual void honk()=0;
};

// Derived class
class Car: public Vehicle {
  public: 
    string model = "Mustang";
    void honk() {
      cout << "Tuut, tuut, says Mustang! \n" ;
      }
};

int main() {
  Car myCar;
  myCar.honk();
  cout << myCar.brand + " " + myCar.model;
  return 0;
}