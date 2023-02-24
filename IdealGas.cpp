#include <math.h>
#include <stdio.h>
#include <iostream>
#include "IdealGas.h"

using namespace std;

IdealGas::IdealGas() : Gas() {

	kappa = 1.4; //Adiabatic exponent
	R = 287.;
	cp = R * kappa / (kappa - 1.);
	cV = R / (kappa - 1.); //Heat capacity at constant volume
}

IdealGas::IdealGas(
		const double _kappa,
		const double _R) : Gas() {

	kappa = _kappa; //Adiabatic exponent
	R = _R;
	cp = R * kappa / (kappa - 1.);
	cV = R / (kappa - 1.); //Heat capacity at constant volume
}

double IdealGas::Get_rho(double p, double T) {
	return p / R / T;
}

double IdealGas::Get_p(double rho, double T) {
	return rho * R * T;
}

double IdealGas::Get_T(double rho, double p) {
	return p / rho / R;
}

double IdealGas::Get_SonicVel(double T, double p) {
	return sqrt(kappa * R * T);
}

double IdealGas::Get_T_from_ep(double e, double p) {
	return e / cV;
}

double IdealGas::Get_T_from_erho(double e, double rho) {
	return e / cV;
}


double IdealGas::Get_e_from_Tp(double T, double p) {
	return T * cV;
}

double IdealGas::Get_Prandtl_from_Tp(double T, double p){
	cout<<endl<<" WARNING! Get_Prandtl_from_Tp() is not implemented for class IdealGas!";
	cout<<endl<<"Returning 0.0!"<<endl;
	return 0.0;
}
double IdealGas::Get_ThermalConductivity_from_Tp(double T, double p){
	cout<<endl<<" WARNING! Get_Prandtl_from_Tp() is not implemented for class IdealGas!";
	cout<<endl<<"Returning 0.0!"<<endl;
	return 0.0;
}

double IdealGas::Get_DynamicViscosity_from_Tp(double T, double p){
	cout<<endl<<" WARNING! Get_Prandtl_from_Tp() is not implemented for class IdealGas!";
	cout<<endl<<"Returning 0.0!"<<endl;
	return 0.0;
}

double IdealGas::Get_kappa_pv() {
	return kappa;
}

double IdealGas::Get_kappa_Tv() {
	return kappa;
}

double IdealGas::Get_kappa_Tp() {
	return kappa;
}

double IdealGas::Get_cp(double p, double T) {
	return cp;
}

double IdealGas::Get_cV(double p, double T) {
	return cV;
}

double IdealGas::Get_MassFlux(double pu, double Tu, double pd, double Td) {

	double rhou = Get_rho(pu,Tu);
	double rhod = Get_rho(pd,Td);
	double kappa = Get_kappa_Tp();
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

double IdealGas::Get_eta_crit() {
	return pow(2. / (kappa + 1), kappa / (kappa - 1));
}


IdealGas::~IdealGas() {}
