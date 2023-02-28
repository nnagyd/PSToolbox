#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "Units.h"


using namespace std;
using namespace Eigen;

//! Class to model a slightly compressible pipe
/*! A class that is used for modelling a Slightly Compressible Pipe. (The latter means that the
  liquid is treated as having constant density, and the speed of sound is much higher than the
  usual flow speeds in the pipe.)
  The solution of the system will be based on the method of characterisctics (MOC)
*/
class SCP: public Units
{
private:
  //Data
  double t, L, D, A, lambda, he, hv, ro, a, roa, dt, lambda_p_2D, S0, g;
  string name, node_from, node_to;
  double phi, alpha, gamma, mu;
  int Npts; //!< number of points the pipe is separated to during
  VectorXd x, v, p;
  bool ini_done; //!< whether the pipe has been initialised or not
  void BCLeft(string type, double val, double &pstart, double &vstart);
  void BCRight(string type, double val, double &pend, double &vend);
  double Source(int i);

  //Saving and plotting data
  vector< vector<double> > data;
  vector<double> tmpvec;
  vector<double> tstatus; //time data for the status
  vector< vector<double> > vstatus, pstatus;
  bool save_data; //!< whether to save data or not
  //bool do_plot_runtime;
  //bool ylim_isset;
  //double lim_pbottom, lim_pup, lim_vbottom, lim_vup;

public:
  SCP(const string _nev,
      const string _cspe_nev,
      const string _cspv_nev,
      const double _ro,
      const double _a,
      const double _L,
      const double _D,
      const double _lambda,
      const double _he,
      const double _hv,
      const bool save_data);
  ~SCP();
  string GetName();
  string Info(bool show_pts);
  void Ini(int Npts_mul);
  void Ini(double vini, double _pstart);
  void Ini(double vini, double _pstart, int Npts_mul);
  void Ini(double vini, double _pstart, double dt_target);
  void UpdateDimlessPars(double pref, double mp_nevl, double omega, double xref, double m);
  double Get_dprop(string prop_string);
  void Set_dprop(string prop_string, double val);
  void Step(
    string BC_start_type, double BC_start_val,
    string BC_end_type, double BC_end_val);
  double GetAlphaAtEnd(double t_target);
  double GetBetaAtFront(double t_target);
  double GetAlphaPrimitiveAtEnd(double t_target); // alias
  double GetBetaPrimitiveAtFront(double t_target); // alias
  void Save_data();
  void Save_status(bool newfile, bool atonce);
  float GetPenult(string what);
  vector<double> Get_dvprop(string prop_string);
  string fname;

};
