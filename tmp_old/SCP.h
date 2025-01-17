#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

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
      Ideal_Gas* gas);
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
  
  double GetBetaAtFront(double t_target);
  double GetBetaPrimitiveAtFront(double t_target);
  void GetAllPrimitiveAtFront(double t_target, double& pP, double& vP, double& TP, double& rhoP);

   double GetAlphaAtEnd(double t_target);
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
  Ideal_Gas *gas;
};


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
  void Save_data();
  void Save_status(bool newfile, bool atonce);
  float GetPenult(string what);
  vector<double> Get_dvprop(string prop_string);
  string fname;

};

//! Class for modeling ideal gases
/*!
  The class contains data from an ideal gas, and providing functions for the
*/
class Reservoir: public Units
{
private:
  bool save_data;
  string name;
  double vol, a, p, t;
  double beta;
  bool is_Ideal_Gas;
  VectorXd ReservoirODE(const double x, const VectorXd &y, const VectorXd &pars);

  void rk45(double x_0,
            double x_max,
            VectorXd y_0,
            vector<double> &x, vector<vector<double> > &y,
            const VectorXd & pars);

  vector<vector<double> > data;
  vector<double> tmpvec;
  bool Ini_done;
  double n_poly;

public:
  Reservoir(const string _nev,
            const double _vol,
            const double _a,
            const bool save_data);
  Reservoir(const string _name,
            const double _vol,
            Ideal_Gas* gas,
            double n_poly,
            const bool _save_data);
  ~Reservoir();
  string GetName();
  string Info();
  void Update(double delta_t, double mp_in, double mp_out);
  void UpdateDimlessPars(double pref, double mp_nevl, double omega);
  //void Ini();
  void Ini(double _pstart);
  void Ini(double _pstart, double _Tstart);
  double Get_dprop(string prop_string);
  void Set_dprop(string prop_string, double val);
  double GetP() {return p;}
  void List_data();
  void Save_data();
  vector<double> Get_dvprop(string prop_string);
  //void Plot();
  string fname;
  Ideal_Gas gas;
};

/*! \brief A class for handling both liquid and fluid valves
  A class for storing data regarding a valve and allowing for a solution of the differential equation
  describing the movement of the spring, mass damper system.
*/
class Valve: public Units
{
private:
  bool save_data; //!< whether to save the data of the valve or not
  string name; //!< name of the valve
  double m /*!< moving mass of the valve */, s /*!< spring rate in valve*/ , k  /*damping factor in valve*/ , A /*effective area of the valve */ , Cd /*!<discharge coefficient of the valve*/;
  double ro /*Desity of fluid in the valve */, Dbore /*Smallest cross section in valve*/, xe /*!< spring pre-compression */, p_set /* set pressure */;
  double mp_nevl /*Design pressure */, xmax /*maximum allowed displacement */, r /*!<restitution factor */;
  double t/*!< current time of the simulation*/, x /*!< displacement of the valve during the simulation */, v /*!< velocity of the valve head during the simulation*/, p /*inlet side pressure in the valve */;
  double pref /*!< refrence pressure*/, xref /*!< reference length*/, delta /*!< ratio of pre-compression and reference length */;
  bool ini_done; //!< whether the system ahs been initialized
  double x_simul_min /*! minimum x value reached during simulation*/, x_simul_max /*!<maximum x value reached during simulation */;
  double TINY_VELOCITY, TINY_OPENING; //!<small values to represent zero
  bool is_closed, is_fully_open; //!< describing the status of the valve for restitution

  VectorXd ValveODE(const double x, const VectorXd &y, const VectorXd &pars);

  void rk45(double x_0,
            double x_max,
            VectorXd y_0,
            vector<double> &x, vector<vector<double> > &y,
            const VectorXd & pars);

  double SumStatForce(const double x, const double pv, const double pb);
  double Aeff(const double x);
  double Aeffmax, a1, a2;


  vector<vector<double> > data;
  vector<double> tmpvec;
  Ideal_Gas gas;
  bool is_Ideal_Gas;

public:
  Valve(const string _nev,
        const double _m,
        const double _k,
        const double _s,
        const double _Dbore,
        const double _Cd,
        const double _ro,
        const double p_set,
        //const double _mp_nevl,
        const double _xmax,
        const double _r,
        const double _pref,
        const bool save_data);

  Valve(const string _nev,
        const double _m,
        const double _k,
        const double _s,
        const double _Dbore,
        const double _Cd,
        Ideal_Gas* _gas,
        const double p_set,
        //const double _mp_nevl,
        const double _xmax,
        const double _r,
        const double _pref,
        const bool save_data);
  ~Valve();
  string GetName();
  string Info();
  void Update(double delta_t, double p_valve, double p_back, bool monitor_min_max);
  void Ini();
  void Ini(double _xstart, double _vstart, double _pstart);
  void Ini(double _xstart, double _vstart);
  void UpdateDimlessPars();
  double Get_dprop(string prop_string);
  vector<double> Get_dvprop(string prop_string);
  void Set_dprop(string prop_string, double val);
  void List_data();
  void Save_data();
  //void Plot();
  string fname;
  void SetAeffCoeffs(const double _a1, const double _a2, const double Aeffmax);

  double Get_MassFlowRate_Compressible(double p_upstream, double T_upstream, double p_downstream, double T_downstream, double xx);
  double Get_MassFlowRate_Compressible_Choked(double p_upstream, double T_upstream, double x);
  double Get_MassFlowRate_Compressible_UnChoked(double p_upstream, double T_upstream, double p_downstream,double xx);

  double Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho, double x) ;

};
