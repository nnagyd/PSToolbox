#include <math.h>
//#include <iostream>
#include "CoolProp.h"
#include "HumidAirProp.h"
#include "CoolPropHA.h"

using namespace std;

// All CoolPropHA function wait for pressure in kPa!!!

CoolPropHA::CoolPropHA(double _R) : Gas() {
	R = _R;
}

double CoolPropHA::Get_rho(double p, double T) {
	//https://en.wikipedia.org/wiki/Density_of_air
	double TT=T-273;
	double psat=100*6.1078*pow(10,7.5*TT/(TT+237.3));
	double pv=R*psat;
	double pd=p-pv;
	double Md=0.0289652; // kg/mol
	double Mv=0.018016; // kg/mol
	double R=8.31446;
return (pd*Md+pv*Mv)/R/T;
}

double CoolPropHA::Get_p(double rho, double T) {
	return HumidAir::HAProps("P", "T", T, "D", rho, "R",R);
}

double CoolPropHA::Get_T(double rho, double p) {
	return HumidAir::HAProps("T", "D", rho, "P", p/1000., "R",R);
}

double CoolPropHA::Get_SonicVel(double T, double p) {
	double cp=HumidAir::HAProps("cp","T",T,"P",p/1000.,"R",R);
	double cV=HumidAir::HAProps("CV","T",T,"P",p/1000.,"R",R);
	double rho= Get_rho(p,T);
	//cout<<endl<<" R ="<<R;
	//cout<<endl<<" T ="<<T;
	//cout<<endl<<" p ="<<p;
	//cout<<endl<<" cp ="<<cp;
	//cout<<endl<<" cV ="<<cV;
	//cout<<endl<<" rho ="<<rho;
	return sqrt(cp/cV*p/rho);
}

double CoolPropHA::Get_T_from_ep(double e, double p) {
	return HumidAir::HAProps("T", "U", e, "P", p/1000., "R",R);
}

double CoolPropHA::Get_T_from_erho(double e, double rho) {
	return HumidAir::HAProps("T", "U", e, "D", rho, "R",R);
}

double CoolPropHA::Get_e_from_Tp(double T, double p) {
	return HumidAir::HAProps("U", "T", T, "P", p/1000., "R",R);
}

double CoolPropHA::Get_Prandtl_from_Tp(double T, double p) {
	return HumidAir::HAProps("Prandtl", "T", T, "P", p/1000., "R",R);
}

double CoolPropHA::Get_ThermalConductivity_from_Tp(double T, double p) {
	return HumidAir::HAProps("conductivity", "T", T, "P", p/1000., "R",R);
}

double CoolPropHA::Get_DynamicViscosity_from_Tp(double T, double p) {
	return HumidAir::HAProps("viscosity", "T", T, "P", p/1000., "R",R);
}

double CoolPropHA::Get_kappa_pv() {
	return 1.4;
}

double CoolPropHA::Get_kappa_Tv() {
	return 1.4;
}

double CoolPropHA::Get_kappa_Tp() {
	return 1.4;
}

double CoolPropHA::Get_cp(double T, double p) {
	return HumidAir::HAProps("cp", "T", T, "P", p/1000., "R",R);
}

double CoolPropHA::Get_cV(double T, double p) {
	return HumidAir::HAProps("CV", "T", T, "P", p/1000., "R",R);
}

double CoolPropHA::Get_cp() {
	cout<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!";
	cout<<endl<<"ERROR IN CoolPropHA::Get_cp() !!!!";
	cin.get();
	return 1.0;
}

double CoolPropHA::Get_cV() {
	cout<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!";
	cout<<endl<<"ERROR IN CoolPropHA::Get_cp() !!!!";
	cin.get();
	return 1.0;
}

double CoolPropHA::Get_MassFlux(double pu, double Tu, double pd, double Td) {

	//double rhou = Get_rho(pu,Tu);
	//double rhod = Get_rho(pd,Td);
	//double kappa = Get_kappa_Tp();
	double rhou = HumidAir::HAProps("D", "T", Tu, "P", pu/1000., "R",R);
	double rhod = HumidAir::HAProps("D", "T", Td, "P", pd/1000., "R",R);
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

double CoolPropHA::Get_eta_crit() {
	double kappa = 1.4;
	return pow(2. / (kappa + 1), kappa / (kappa - 1));
}


CoolPropHA::~CoolPropHA() {}
