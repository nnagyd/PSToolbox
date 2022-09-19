#include <math.h>
//#include <iostream>
#include "CoolProp.h"
#include "CoolPropGas.h"

using namespace std;

CoolPropGas::CoolPropGas(string _fluid_name) : Gas() {
		fluid_name = _fluid_name;
}

double CoolPropGas::Get_rho(double p, double T) {
		return CoolProp::PropsSI("D", "T", T, "P", p, fluid_name);
}

double CoolPropGas::Get_p(double rho, double T) {
		return CoolProp::PropsSI("P", "T", T, "D", rho, fluid_name);
}

double CoolPropGas::Get_T(double rho, double p) {
		return CoolProp::PropsSI("T", "D", rho, "P", p, fluid_name);
}

double CoolPropGas::Get_SonicVel(double T, double p) {
		return CoolProp::PropsSI("A", "T", T, "P", p, fluid_name);
}

double CoolPropGas::Get_T_from_ep(double e, double p) {
		return CoolProp::PropsSI("T", "U", e, "P", p, fluid_name);
}

double CoolPropGas::Get_T_from_erho(double e, double rho) {
		return CoolProp::PropsSI("T", "U", e, "D", rho, fluid_name);
}

double CoolPropGas::Get_e_from_Tp(double T, double p) {
		return CoolProp::PropsSI("U", "T", T, "P", p, fluid_name);
}

double CoolPropGas::Get_kappa_pv() {
		return 1.4;
}

double CoolPropGas::Get_kappa_Tv() {
		return 1.4;
}

double CoolPropGas::Get_kappa_Tp() {
		return 1.4;
}

double CoolPropGas::Get_cp(double T, double p) {
		return CoolProp::PropsSI("CPMass", "T", T, "P", p, fluid_name);
}

double CoolPropGas::Get_cV(double T, double p) {
		return CoolProp::PropsSI("CVMass", "T", T, "P", p, fluid_name);
}

double CoolPropGas::Get_cp() {
	cout<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!";
	cout<<endl<<"ERROR IN CoolPropGas::Get_cp() !!!!";
	cin.get();
		return 1.0;
}

double CoolPropGas::Get_cV() {
	cout<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!";
	cout<<endl<<"ERROR IN CoolPropGas::Get_cp() !!!!";
	cin.get();
		return 1.0;
}

double CoolPropGas::Get_MassFlux(double pu, double Tu, double pd, double Td) {

		//double rhou = Get_rho(pu,Tu);
		//double rhod = Get_rho(pd,Td);
		//double kappa = Get_kappa_Tp();
		double rhou = CoolProp::PropsSI("D", "T", Tu, "P", pu, fluid_name);
		double rhod = CoolProp::PropsSI("D", "T", Td, "P", pd, fluid_name);
		double kappa = 1.4;//Get_kappa_Tp();
		double massflux, dimless_massflux;

		// Change flow direction if necessary
		double dir_mul = 1.;
		if (pd > pu) {
				double p_tmp = pu;
				double rho_tmp = rhou;
				pu   = pd;
				rhou = rhod;
				pd   = p_tmp;	
				rhod = rho_tmp;
				dir_mul = -1.;
		}

		double pr = pd / pu; // pressure ration
		double eta_c = Get_eta_crit();

		if (pr < eta_c) {
				// Choked
				double exponent = (kappa + 1.) / (kappa - 1.);
				double tmp = kappa * pow(2. / (kappa + 1.), exponent);
				dimless_massflux =  sqrt(tmp);
		}
		else {
				// UnChoked
				double exp1 = 2. / kappa;
				double exp2 = (kappa + 1.) / kappa;
				double tmp = 2.*kappa / (kappa - 1);
				dimless_massflux = sqrt(tmp * (pow(pr, exp1) - pow(pr, exp2)));
		}
		return massflux =  dir_mul*sqrt(rhou * pu) * dimless_massflux;
}

double CoolPropGas::Get_eta_crit() {
	double kappa = 1.4;
		return pow(2. / (kappa + 1), kappa / (kappa - 1));
}


CoolPropGas::~CoolPropGas() {}
