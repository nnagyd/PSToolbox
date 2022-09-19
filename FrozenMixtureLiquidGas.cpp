#include <math.h>
#include <iostream>
#include "FrozenMixtureLiquidGas.h"

using namespace std;

FrozenMixtureLiquidGas::FrozenMixtureLiquidGas(
    const double _x_g,
    const double _aL,
    const double _rhoL_ref,
    const double _pL_ref,
    const double _cL,
    IdealGas* _gas)  : IdealGas() {

	x_g    = _x_g;
	rhoL_ref =  _rhoL_ref;
	pL_ref = _pL_ref;
	aL     = _aL;
	cL     = _cL;
	gas    = _gas;
	omega  = 1.;
}

double FrozenMixtureLiquidGas::Get_rho(double p, double T) {
	double rho_G = gas->Get_rho(p, T);
	double alpha_g = x_g / (x_g + (1 - x_g) * rho_G / rho_L(p));
	return alpha_g * rho_G + (1 - alpha_g) * rho_L(p);
}

double FrozenMixtureLiquidGas::Get_p(double rho, double T) {
	double pp = 1.;
	double R = gas->R;

	double b = -(pL_ref + R * T * x_g * rho + cL * cL * (rho * (1 - x_g) - rhoL_ref));
	double c = R * T * x_g * rho * (pL_ref - cL * cL * rhoL_ref);
	if (b * b - 4 * c < 0) {
		cout << endl << endl << " ERROR in FrozenMixtureLiquidGas::Get_p: negative square root argument!!!";
	} else
		pp = (-b + sqrt(b * b - 4 * c)) / 4.;
	return pp;
}

double FrozenMixtureLiquidGas::Get_T(double rho, double p) {
	cout << endl << endl << "!!!!!!! FrozenMixtureLiquidGas::Get_T(rho,T) not implemented!!!!!" << endl << endl;
	cin.get();
	return 293.;
}

double FrozenMixtureLiquidGas::Get_SonicVel(double T, double p) {
	double ppdp = 1.01 * p;
	double kappa = gas->Get_kappa_Tp();
	double TpdT = T * pow(ppdp / p, (kappa - 1) / kappa);
	double rho_1 = Get_rho(p, T);
	double rho_2 = Get_rho(ppdp, TpdT);
	return sqrt((ppdp - p) / (rho_2 - rho_1) );
}

double FrozenMixtureLiquidGas::Get_alpha_g(double p, double T) {
	double rho_G = gas->Get_rho(p, T);
	return  x_g / (x_g + (1 - x_g) * rho_G / rho_L(p));
}

double FrozenMixtureLiquidGas::rho_L(double p) {
	return (rhoL_ref + 1 / aL / aL * (p - pL_ref));
}

double FrozenMixtureLiquidGas::Get_T_from_e(double e) {
	double cV_mixture = x_g * gas->cV + (1 - x_g) * cL;
	return e / cV_mixture;
}

double FrozenMixtureLiquidGas::Get_e(double T) {
	double cV_mixture = x_g * gas->cV + (1 - x_g) * cL;
	return T * cV_mixture;
}

double FrozenMixtureLiquidGas::Get_kappa_pv() {
	return gas->kappa;
}

double FrozenMixtureLiquidGas::Get_kappa_Tv() {
	return gas->kappa;
}

double FrozenMixtureLiquidGas::Get_kappa_Tp() {
	return gas->kappa;
}

double FrozenMixtureLiquidGas::Get_cp(double p, double T) {
	return (x_g*gas->cV)+(1.-x_g)*cL;
}

double FrozenMixtureLiquidGas::Get_cV(double p, double T) {
	return (x_g*gas->cp)+(1.-x_g)*cL;
}

double FrozenMixtureLiquidGas::Get_MassFlux(double pu, double Tu, double pd, double Td) {

	double rhou = Get_rho(pu, Tu);
	double rhod = Get_rho(pd, Td);
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

	double pr = pd / pu; // pressure ratio
	omega = Get_alpha_g(pu, Tu) / kappa;
	double eta_c = Get_eta_crit();

	if (pr < eta_c) {
		// Choked
		dimless_massflux =  eta_c / sqrt(omega);
	}
	else {
		// UnChoked
		dimless_massflux = sqrt(-2.*(omega * log(pr) + (omega - 1) * (1 - pr))) / (omega * (1. / pr - 1.) + 1.);
	}
	return massflux =  dir_mul * sqrt(rhou * pu) * dimless_massflux;

	return 1.;
}

double FrozenMixtureLiquidGas::Get_eta_crit() {
	double x_g_old = x_g;
	bool x_g_too_small = false;
	if (x_g < 1.e-10) {
		x_g = 1.e-10;
		x_g_too_small = true;
	}

// Curve fitting for good initial guess, see Matlab code for details
	double eta_c = 0.5962 * pow(omega, 0.1671);
	double err = 1.e5, TOL = 1.e-3, f2, dx, dfdx;
	int iter = 0, ITER_MAX = 100;
	while (err > TOL) {
		err = eta_crit_fun(eta_c);
		//printf("\n iter: %2d -> eta_c=%5.3e, err=%+5.3e", iter, eta_c, err);
		dx = eta_c * 0.01;
		f2 = eta_crit_fun(eta_c + dx);
		dfdx = (f2 - err) / dx;
		eta_c -= err / dfdx ;
		iter++;
		if (iter > ITER_MAX) {
			cout << endl << "ERROR in FrozenMixtureLiquidGas::Get_eta_crit" << endl;
			cin.get();
		}
	}
	err = eta_crit_fun(eta_c);
	//printf("\n iter: %2d -> eta_c=%5.3e, err=%+5.3e", iter, eta_c, err);

	if (x_g_too_small)
		x_g = x_g_old;

	return eta_c;
}

double FrozenMixtureLiquidGas::eta_crit_fun(double x) {
	double val = x * x + (omega * omega - 2 * omega) * (1. - x) * (1. - x) + 2 * omega * omega * log(x) + 2 * omega * omega * (1 - x);
	return val;
}

FrozenMixtureLiquidGas::~FrozenMixtureLiquidGas() {}
