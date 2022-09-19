#pragma once
#include "IdealGas.h"

class FrozenMixtureLiquidGas : public IdealGas
{
public:
  FrozenMixtureLiquidGas(
    const double x_g,
    const double _aL,
    const double _rhoL_ref,
    const double _pL_ref,
    const double _cV,
    IdealGas* gas);
  ~FrozenMixtureLiquidGas();

  double Get_p(double rho, double T);
  double Get_rho(double p, double T);
  double Get_T(double p, double rho);
  double Get_e(double T);
  double Get_SonicVel(double T, double p);
  double Get_T_from_e(double e);
  double Get_kappa_pv();
  double Get_kappa_Tv();
  double Get_kappa_Tp();
  double Get_pcrit(double pr, double omega);
  double Get_cp(double p, double T);
  double Get_cV(double p, double T);
  double Get_MassFlux(double pu, double Tu, double pd, double Td);

  double rho_L(double p);
  double Get_alpha_g(double p, double T);  

  double rhoL_ref, pL_ref, aL, cL;
  double x_g;
  IdealGas* gas;

private:
  double omega;
  double eta_crit_fun(double x);
  double Get_eta_crit();
};
