// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test8_FrozenMixture.cpp

#include <stdio.h>
#include "/Users/hoscsaba/program/PSToolbox/IdealGas.h"
#include "/Users/hoscsaba/program/PSToolbox/FrozenMixtureLiquidGas.h"
//#include "/Users/hoscsaba/program/PSToolbox/Connector.h"

using namespace std;

int main(int argc, char **argv) {
	
	double xg=1.0;
	
	IdealGas* gas = new IdealGas(1.4,287.);
	FrozenMixtureLiquidGas* fm = new FrozenMixtureLiquidGas(xg,1400,1000,1e5,4000.,gas);
	
	printf("\n xg=%5.3f",xg);
	printf("\n a_mixture   = %5.2f m/s,   a_gas  =%5.2f m/s",fm->Get_SonicVel(293.,1e5),gas->Get_SonicVel(293,1.e5));
	printf("\n rho_mixture = %5.1f kg/m3, rho_gas =%5.1f kg/m3\n",fm->Get_rho(1e5,293),gas->Get_rho(1.e5,293));

	double T=293., pu=1.1e5, pd=1.0e5;
	printf("\n Forward flow, subsonic, pu/pd=%5.2f",pu/pd);
	printf("\n\t mass flux mixture   = %5.2f kg/s/m2",fm->Get_MassFlux(pu,T,pd,T));
	printf("\n\t mass flux ideal gas = %5.2f kg/s/m2\n",gas->Get_MassFlux(pu,T,pd,T));

	T=293., pu=3.e5, pd=1.0e5;
	printf("\n Forward flow, supersonic, pu/pd=%5.2f",pu/pd);
	printf("\n\t mass flux mixture   = %5.2f kg/s/m2",fm->Get_MassFlux(pu,T,pd,T));
	printf("\n\t mass flux ideal gas = %5.2f kg/s/m2\n",gas->Get_MassFlux(pu,T,pd,T));

	T=293., pu=1.0e5, pd=1.1e5;
	printf("\n Backward flow, subsonic, pu/pd=%5.2f",pu/pd);
	printf("\n\t mass flux mixture   = %5.2f kg/s/m2",fm->Get_MassFlux(pu,T,pd,T));
	printf("\n\t mass flux ideal gas = %5.2f kg/s/m2\n",gas->Get_MassFlux(pu,T,pd,T));

	T=293., pu=1.e5, pd=3.0e5;
	printf("\n Backward flow, supersonic, pu/pd=%5.2f",pu/pd);
	printf("\n\t mass flux mixture   = %5.2f kg/s/m2",fm->Get_MassFlux(pu,T,pd,T));
	printf("\n\t mass flux ideal gas = %5.2f kg/s/m2\n",gas->Get_MassFlux(pu,T,pd,T));


}
