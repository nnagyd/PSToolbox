// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long 2pipes_with_ideal_gas.cpp

#include <stdio.h>
#include "/Users/hoscsaba/program/PSToolbox/LWP.h"
#include "/Users/hoscsaba/program/PSToolbox/IdealGas.h"
#include "/Users/hoscsaba/program/PSToolbox/Connector.h"

using namespace std;

int main(int argc, char **argv) {
	
	double p0=1.e5, dp=1.e5, T0=293;
	double L1=10., D1=0.1, lambda1=0.02;	
	double L2=10., D2=0.05, lambda2=0.02;

	IdealGas* gas = new IdealGas(1.4, 287.);	
	LWP *p1 = new LWP("p1", "csp1", "csp2", L1, D1, lambda1, 0, 0, false, false, gas);
	LWP *p2 = new LWP("p2", "csp2", "csp3", L2, D2, lambda2, 0, 0, false, false, gas);
	Connector c;

	p1->IniUniform(0., p0, T0, 20);
	p2->IniUniform(0., p0, T0, 20);

	double pL,pR,TL,TR,rhoL,rhoR,vL,vR;
	double t = 0., tmax=1., dt=1., dt2=dt;	
	
	while (t<tmax) {
		dt  = p1->Get_dt() * 0.5;
		dt2 = p2->Get_dt() * 0.5;
		if (dt>dt2)
			dt=dt2;

		t += dt;
		c.Connector_LWP_Pipes(t, p1, p2, pL, pR, TL, TR,rhoL,rhoR,vL,vR);
		p1->Step("TotalPres_and_TotalTemp",   p0+dp,T0,
		         "StaticPres_and_StaticTemp", pL,   TL, dt);
		p2->Step("StaticPres_and_StaticTemp", pR,TR,
		         "StaticPres_and_StaticTemp", p0, T0, dt);
		
		printf("\n t=%5.3f s, p=%5.3f barg, vL=%5.3f m/s, vR=%5.3f m/s",t, pL/1.e5,vL,vR); 
	}
}
