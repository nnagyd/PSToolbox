#define _USE_MATH_DEFINES

// MAC:
// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7  -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long reservoir_pipe_with_ideal_gas.cpp

#include <stdio.h>
#include <iostream>
#include "/Users/hoscsaba/program/PSToolbox/my_tools.h"
#include "/Users/hoscsaba/program/PSToolbox/LWP.h"
#include "/Users/hoscsaba/program/PSToolbox/Reservoir.h"
#include "/Users/hoscsaba/program/PSToolbox/Valve.h"
#include "/Users/hoscsaba/program/PSToolbox/Connector.h"

#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"
void Plot(Reservoir* r, LWP* p);

using namespace std;

int main(int argc, char **argv) {

	double p0 = 1.e5, pt = p0 + 0.5e5, T0 = 293;

	IdealGas* gas = new IdealGas(1.4, 287.);
	Reservoir* r1 = new Reservoir("tank", 100000., gas, /*n_poly*/ 1.0, true);
	LWP* p1 = new LWP("pipe", "csp1", "csp2", 10., 0.05, 0.02, 0, 0, true, false, gas);
	Connector c;

	r1->Ini(pt, T0);
	p1->IniUniform(0., p0, T0, 20.);

	double t = 0, tmax = 1.0, dt = 0., pcon=1e5, Tcon=293;
	while (t < tmax) {
			//cout << r1->Info();
			cout << p1->Info(true);
			//cin.get();
		dt = p1->Get_dt() * 0.9;

		c.Connector_LWP_Reservoir_and_Pipe_Front(t+dt,r1,p1,pcon,Tcon);

		t += dt;

		r1->Update(dt, 0.0, p1->Get_dprop("mp_front"));
		//p1->Step("StaticPres_and_StaticTemp", pcon,Tcon,"Outlet", p0, T0, dt);
		p1->Step("StaticPres_and_StaticTemp", pcon,Tcon,"Wall", p0, T0, dt);

	}
	printf("\n\n t=%5.3f s, pt,p(0),p(L)=%5.3f, %5.3f, %5.3f bar"
	       ", v(0),v(L)= %+5.3f, %+5.3f m/s",
	       t, r1->Get_dprop("p") / 1.e5, p1->Get_dprop("p_front_bar"), p1->Get_dprop("p_back_bar"),
	       p1->Get_dprop("v_front"), p1->Get_dprop("v_back"));
	
	Plot(r1,p1);
}


void Plot(Reservoir* r, LWP* p){
namespace plt = matplotlibcpp;
	vector<double> tr = r->Get_dvprop("t");
	vector<double> pr = r->Get_dvprop("p_bar");

	vector<double> tp = p->Get_dvprop("t");
	vector<double> pf = p->Get_dvprop("p_front_bar");
	vector<double> pb = p->Get_dvprop("p_back_bar");
	vector<double> vf = p->Get_dvprop("v_front");
	vector<double> vb = p->Get_dvprop("v_back");

	for (int i=0; i<tp.size(); i++)
	printf("\n \t t=%6.4fs, %5.3f  %5.3f, %+5.2f %+5.2f",tp.at(i),pf.at(i),pb.at(i),vf.at(i),vb.at(i));
	
/*	stringstream fname(""), title("");
	title<<"air, L="<<(p->Get_dprop("L"))<<"m, (Izuchi:"<<Lcrit1<<", Hos: "<<Lcrit2<<")";
	fname<<"air_L_"<<(p->Get_dprop("L"))<<".png";
	cout<<endl<<"title: "<<title.str();
	cout<<endl<<"fname: "<<fname.str();
*/	
	plt::subplot(2,1,1);
	plt::plot(tr, pr,"b");
	//plt::plot(tr,pset,"r--");
	plt::plot(tp,pf,"k");
	plt::plot(tp,pb,"k--");
	plt::ylabel("p (abs), bar");
	//plt::title(title.str());

	plt::subplot(2,1,2);
	plt::plot(tp,vf,"k");
	plt::plot(tp,vb,"k--");
	plt::ylabel("v m/s");

	plt::show();
	//plt::save(fname.str());
}
