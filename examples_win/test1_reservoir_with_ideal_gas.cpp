#define _USE_MATH_DEFINES

// MAC:
// g++ -IC:/ProgramData/chocolatey/lib/eigen/include -IC:/Users/hoscs/Documents/PSToolbox -LC:/Users/hoscs/Documents/PSToolbox -lPSToolbox.a -pedantic -std=c++11 test1_reservoir_with_ideal_gas.cpp
// Win
//  g++ -IC:/ProgramData/chocolatey/lib/eigen/include -IC:/Users/hoscs/Documents/PSToolbox -LC:/Users/hoscs/Documents/PSToolbox test1_reservoir_with_ideal_gas.cpp -lPSToolbox -pedantic -std=c++11

#include <stdio.h>
#include <iostream>
#include "my_tools.h"
#include "IdealGas.h"
#include "Reservoir.h"

using namespace std;

int main(int argc, char **argv) {

	IdealGas* gas = new IdealGas(1.4, 287.);
	Reservoir r1("tank", 1, gas, /*n_poly*/ 1.4, true);
	Reservoir r2("tank", 1, 1000., true);

	r1.Ini(1.e5, 283.);
	r2.Ini(1.e5);

	double t = 0, dt = 0.01;
	for (unsigned int i = 0; i < 10; i++) {
		printf("\n t=%5.3f s, " 
			"gas tank: p=%5.3f bar, rho=%5.3f kg/m3, a=%5.1f m/s, " 
			"liquid tank: p=%5.3f bar, a=%5.1f m/s",
		       t, 
		       r1.Get_dprop("p") / 1.e5, r1.gas->Get_rho(r1.Get_dprop("p"),r1.Get_dprop("T")),r1.Get_dprop("a"),
		       r2.Get_dprop("p") / 1.e5, r2.Get_dprop("a"));
		t += dt;

		r1.Update(dt, 1., 0.);
		r2.Update(dt, 1., 0.);

		//gas.Setp(r1.Get_dprop("p"));
		//gas.UpdatePolytropic(1.4);
	}
}
