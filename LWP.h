#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "Units.h"
#include "IdealGas.h"

using namespace std;
using namespace Eigen;

class LWP: public Units
{
private:
  double t, L, D, A, lambda, he, hv, dt, lambda_p_2D, S0, g, dx;
  string name, node_from, node_to;
  double phi, alpha, gamma, mu;
  int Npts;
  VectorXd x, v, p, T, rho, q;
  bool ini_done;
  void UpdateInternalPoints();
  //VectorXd &p, VectorXd &v, VectorXd &T, VectorXd &rho);
  void BCLeft(string type, double val1, double val2);
  void BCRight(string type, double val1, double val2);

  double Source(int i);
  vector<vector<double> > data;
  vector<double> tmpvec;
  bool save_data, save_all_data;
  //bool do_plot_runtime;
  //bool ylim_isset;
  //double lim_pbottom, lim_pup, lim_vbottom, lim_vup;
  
  void Pack(bool is_half_step);
  void UnPackU(bool is_half_step);

  VectorXd phalf, vhalf, Thalf, rhohalf;
  VectorXd pnew, vnew, Tnew, rhonew;
  MatrixXd U, F, S, Uhalf, Fhalf, Shalf, Unew;
  void Ini(int Npts_mul);
  void UpdateTimeStep();

  void Add_data_row();
  double P_MIN, T_MIN, art_visc;

public:
  LWP(const string _nev,
    const string _cspe_nev,
    const string _cspv_nev,
      //const double _ro,
      //const double _a,
    const double _L,
    const double _D,
    const double _lambda,
    const double _he,
    const double _hv,
    const bool _save_data,
    const bool _save_all_data,
    Gas* gas);
  ~LWP();
  string Info(bool show_pts);
  void IniUniform(double _vini, double _pini, double _Tini, int Npts_mul);
  void IniUniform(double _vini, double _pini, double _Tini, double dt_target);
  void IniDistribution(const vector<double> _vini,
   const vector<double> _pini,
   const vector<double> _Tini);
  /*
  void Ini(double vini, double _pstart, double dt_target);
  void UpdateDimlessPars(double pref, double mp_nevl, double omega, double xref, double m);
  */
  double Get_dprop(string prop_string);
  void Set_dprop(string prop_string, double val);
  void Step(
    string BC_start_type, double BC_start_val1, double BC_start_val2,
    string BC_end_type, double BC_end_val1, double BC_end_val2,
    double dt_req);
  
 // double GetBetaAtFront(double t_target);
  bool GetC0AtFront(double t_target, double& p, double& v, double& T, double& rho);
  double GetBetaPrimitiveAtFront(double t_target);
  void GetAllPrimitiveAtFront(double t_target, double& pP, double& vP, double& TP, double& rhoP);

  //double GetAlphaAtEnd(double t_target);
  double GetC0AtEnd(double t_target);
  double GetAlphaPrimitiveAtEnd(double t_target);
  void GetAllPrimitiveAtEnd(double t_target, double& pP, double& vP, double& TP, double& rhoP);

  string GetName();
  double Get_dt();
  vector<double> Get_xgrid();
  vector<double> Get_p();
  vector<double> Get_v();
  vector<double> Get_dvprop(string prop_string);
  void Save_data();
  void UpdateDimlessPars(double pref, double mp_nevl, double omega, double xref, double m);

  string fname;
  Gas *gas;
};
