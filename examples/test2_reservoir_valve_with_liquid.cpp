//#define _USE_MATH_DEFINES

// MAC:
// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lSCP -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long reservoir_valve_with_liquid.cpp

// SERVER:
// g++ `pkg-config --cflags eigen3` /home/cshos/git/PSToolbox/staciserver/SCP.o /home/cshos/git/PSToolbox/staciserver/my_tools.o -I/usr/include/python3.6m/ -lpython3.6m  -pedantic -O3 -Wall reservoir_valve_with_ideal_gas.cpp

#include <stdio.h>
#include "/Users/hoscsaba/program/PSToolbox/Reservoir.h"
#include "/Users/hoscsaba/program/PSToolbox/Valve.h"

using namespace std;

int main(int argc, char **argv) {

	double p_set = 3.e5, p0=1.e5, rho=1000., a=1400.;
	double D_bore=0.05, x_max=0.01;
	double A_nozzle=D_bore*D_bore*M_PI/4;
	double s=A_nozzle*0.1*p_set/x_max;

	Reservoir r1("tank", 100000., a,true);
	Valve v1("valve", 0.2, 0.05*2*sqrt(0.2*47300), s, D_bore, 0.7, rho, p_set, x_max, 0.8, 1.e5, true);

	r1.Ini(p0+p_set*1.1);
	v1.Ini(0.,0.);

	double t = 0, tmax=3., dt=1./v1.Get_dprop("freq")/50.;
	double mp,mpin=v1.Get_dprop("mp_nom");
	while (t<tmax) {
		mp=v1.Get_MassFlowRate_InCompressible(r1.Get_dprop("p"),p0,rho,v1.Get_dprop("x"));
		printf("\n t=%5.3f s, p=%5.3f barg, x=%5.3f mm, mp=%5.3f kg/s (imbalance: %5.3f%%)",
		       t, r1.Get_dprop("p")/1.e5-1., v1.Get_dprop("x_mm"),
		       mp, (mp-mpin)/mpin*100);
		t += dt;

		r1.Update(dt, mpin, mp);
		v1.Update(dt, r1.Get_dprop("p"), 1.e5,false);		
	}
}
