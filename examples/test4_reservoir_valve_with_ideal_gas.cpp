#define _USE_MATH_DEFINES

// MAC:
// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7  -lmy_tools -lPSToolbox -pedantic -O3 -Wall  -Wno-c++11-long-long reservoir_valve_with_ideal_gas.cpp

#include <stdio.h>

#include "/Users/hoscsaba/program/PSToolbox/my_tools.h"
#include "/Users/hoscsaba/program/PSToolbox/LWP.h"
#include "/Users/hoscsaba/program/PSToolbox/Reservoir.h"
#include "/Users/hoscsaba/program/PSToolbox/Valve.h"
#include "/Users/hoscsaba/program/PSToolbox/Connector.h"

#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"
void Plot(Reservoir* r, Valve* v);

using namespace std;

int main(int argc, char **argv) {

	double p_set = 3.e5, p0 = 1.e5, m = 0.2;
	double D_bore = 0.05, x_max = 0.01;
	double s = 47.3e3, k_crit = 2. * sqrt(s * m);

	IdealGas* gas = new IdealGas(1.4, 287.);
	Reservoir* r1 = new Reservoir("tank", 1., gas, /*n_poly*/ 1.0, true);
	Valve* v1 = new Valve("valve", m, 0.01 * k_crit, s, D_bore, 0.7, gas,
	                      p_set, x_max, /* r */ 0.8, 1.e5, true);

	r1->Ini(p0 + p_set * 0.95, 293.);
	v1->Ini(0., 0.);

	double t = 0, tmax = 10., dt = 1. / v1->Get_dprop("freq") / 10., tout = -1.e-10, dtout = tmax / 20.;
	double mp, mpin = v1->Get_dprop("mp_nom") * 1.0;
	while (t < tmax) {
		mp = v1->Get_MassFlowRate_Compressible_Choked(r1->Get_dprop("p"), r1->Get_dprop("T"), v1->Get_dprop("x"));
		if (t > tout) {
			printf("\n t=%5.3f s, p=%5.3f barg, x=%5.1f mm, mp=%5.3f kg/s (imbalance: %+5.1f%%)",
			       t, r1->Get_dprop("p") / 1.e5 - 1., v1->Get_dprop("x_mm"),
			       mp, (mp - mpin) / mpin * 100);
			tout += dtout;
		}

		t += dt;

		r1->Update(dt, mpin, mp);
		v1->Update(dt, r1->Get_dprop("p"), 1.e5, false);
	}
	Plot(r1, v1);
}

void Plot(Reservoir* r, Valve* v) {
	namespace plt = matplotlibcpp;
	vector<double> tr = r->Get_dvprop("t");
	vector<double> pr = r->Get_dvprop("p_bar");

	vector<double> tv = v->Get_dvprop("t");
	vector<double> xv = v->Get_dvprop("x_mm");
	vector<double> psetvec;
	for (unsigned int i = 0; i < tv.size(); i++)
		psetvec.push_back((v->Get_dprop("p_set_bar")) + 1.);

	plt::subplot(2, 1, 1);
	plt::plot(tr, pr, "b");
	plt::plot(tv, psetvec, "r--");
	plt::ylabel("p (abs), bar");

	plt::subplot(2, 1, 2);
	plt::plot(tv, xv, "m");
	plt::ylabel("lift, m");

	plt::show();
	//plt::save(fname.str());
}