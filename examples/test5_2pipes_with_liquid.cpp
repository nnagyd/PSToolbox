// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long 2pipes_with_liquid.cpp

#include <stdio.h>
#include "/Users/hoscsaba/program/PSToolbox/SCP.h"
#include "/Users/hoscsaba/program/PSToolbox/Connector.h"

using namespace std;

int main(int argc, char **argv) {
	
	double p0=1.e5, dp=1.e5;
	double L1=10., D1=0.1, lambda1=0.02;	
	double L2=10., D2=0.05, lambda2=0.02;

	SCP *p1 = new SCP("p1", "csp1","csp2",1000.,1400,L1, D1, lambda1, 0., 0., true); 
	SCP *p2 = new SCP("p2", "csp2","csp3",1000.,1400,L2, D2, lambda2, 0., 0., true); 
	Connector c;

	p1->Ini(0., p0,1.);
	p2->Ini(0., p0,1.);

	double p,v1,v2;
	double t = 0., tmax=10.;
	double dt=p1->Get_dprop("dt");
	if (p2->Get_dprop("dt")<dt)
		dt=p2->Get_dprop("dt");
	
	while (t<tmax) {
		t += dt;
		c.Connector_SCP_Pipes(t, p1, false, p2, true, /*mpout*/ 0.,p,v1,v2); 
		p1->Step("Pressure",2e5,"Pressure",p);
		p2->Step("Pressure",p,"Pressure",1e5);
		
		printf("\n t=%5.3f s, p=%5.3f barg, v1=%5.3f m/s, v2=%5.3f m/s",t, p/1.e5,v1,v2); 
	}

	double A1=p1->Get_dprop("A");
	double A2=p2->Get_dprop("A");	
	double Qs = sqrt(dp/(lambda1*L1/D1/A1/A1+lambda2*L2/D2/A2/A2)/500.);
	printf("\n\n Steady-state solution: dp=%5.3f bar, v1s=%5.3f m/s, v2s=%5.3f m/s\n",dp/1.e5,Qs/A1,Qs/A2);
}
