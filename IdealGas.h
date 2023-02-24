#pragma once
#include "Gas.h"

class IdealGas : public Gas
{
	public:
		IdealGas(
				const double _kappa, ///< [in] Adiabatic exponent*/
				const double _R /**< [in] temperature, T*/
			);
		IdealGas(IdealGas* g);
		IdealGas();
		~IdealGas();

		double Get_p(double rho, double T);
		double Get_rho(double p, double T);
		double Get_T(double p, double rho);
		double Get_e_from_Tp(double T, double p);
		double Get_SonicVel(double T, double p);
		double Get_T_from_ep(double e, double p);
		double Get_T_from_erho(double e, double rho);
		double Get_Prandtl_from_Tp(double T, double p);
		double Get_ThermalConductivity_from_Tp(double T, double p);
		double Get_DynamicViscosity_from_Tp(double T, double p);
		double Get_kappa_pv();
		double Get_kappa_Tv();
		double Get_kappa_Tp();
		double Get_pcrit();
		double Get_cp(double p, double T);
		double Get_cV(double p, double T);
		double Get_MassFlux(double pu, double Tu, double pd, double Td);
		double Get_eta_crit();

		double kappa;
		double R;
		double cp;
		double cV;
};
