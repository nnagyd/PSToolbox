#include <string>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


//! Class for modelling ideal gases
/*!
  The class contains data from an ideal gas, and providing functions for the
*/
class Ideal_Gas
{
public:
  //! A class constructor for containing Ideal gas data.
  /*! A class constructor used to calculate the various specific heats of an ideal gas and
    for retrieving parameters (density, pressure, sonic velocity) of the gas.
  */
  Ideal_Gas(
    const double _kappa, ///< [in] Adiabatic exponent*/
    const double _R /**< [in] Specific gas constant*/
  );
  Ideal_Gas(Ideal_Gas const & g);
  ~Ideal_Gas();

  //! Get the density of the gas
  /*! Calculate the density of an ideal gas at a given pressure and
    temperature. All units are to be in SI.
    \sa GetP()
    \param p Pressure in [Pa]
    \param T Temperature in [K]
    \return The density in [kg/m^3]
  */
  double GetRho(double p, double T);
  //! Get the pressure of the gas
  /*! Calculate the pressure of the ideal gas at a given
    density and temperature. All units in SI.
    \sa GetP()
    \param rho Density in [kg/m^3]
    \param T Absolute temperature [K]
    \return pressure in [Pa]
  */
  double GetP(double rho, double T);
  //! Return the sonic velocity
  /*! Returns the sonic velocity depentdent on the temperature
    \param T temperature in [K]
    \return sonic velocity in [m/s]
  */
  double GetSonicVel(double T);
  double kappa; //!< The ratio of specific heats
  double R; //!< The specific gas constant
  double cp; //!< Izobaric specific heat
  double cv; //!< Isochoric specific heat

};

//! Class to contain unit conversion values
/*! A class that contains the unit conversion values used in the various simulations
    Serves as a base class for the other classes s it easily provides conversion values that way.
  \sa LWP, SCP, Valve, Reservoir
*/
class Units
{
public:
  //! Constructor
  /*! Create an instance of the Units class. Values of the conversion
      are assigned during this phase. (and as such are not constants.)
  */
  Units();
  //! Destuctor
  ~Units();
  double ft_to_m; //!< Converts foot to meters
  double m_to_ft; //!< Convert meters to foot
  double inch_to_m; //!< Convert inches to meters
  double m_to_inch; //!< Convert metres to inches
  double lbm_to_kg; //!< Convert pounds to kilograms
  double kg_to_lbm; //!< Multiplier for conversion from kilogram to pound
  double lbf_to_N; //!< Multiplier for conversion of poundforce to Newtons
  double N_to_lbf; //!< Mulitplier for conversion of Newtons to poundforce
  double psi_to_Pa; //!<Mulitplier ofr conversion of pound per square inch to Pa
  double Pa_to_psi; //!< Multiplier for conversion of pascals to pound per square inch
  double psi_to_bar; //!< Mulitplier for conversion of psi to bar
  double bar_to_psi; //!< Multiplier for conversion of bar to psi
  double m3_to_inch3; //!< Multiplier for conversion of cubic meter to cumbic inch
  double inch3_to_m3; //!< Multiplier for conversion of cubic inches to cubic metres
};


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
  bool do_plot_runtime;
  bool ylim_isset;
  double lim_pbottom, lim_pup, lim_vbottom, lim_vup;
  Ideal_Gas *gas;
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
      Ideal_Gas *gas);
  ~LWP();
  string Info(bool show_pts);
  void IniUniform(double _vini, double _pini, double _Tini, int Npts_mul);
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
  double GetAlphaAtEnd(double t_target);
  double GetBetaAtFront(double t_target);
  string GetName();
  double Get_dt();
  vector<double> Get_xgrid();
  vector<double> Get_p();
  vector<double> Get_v();

  void Save_data();
  /*
  void Plot_p();
  void Plot_v();
  void Plot_T();
  void Plot_rho();
  void Plot_mp();

  void runtime_plot();
  
  void Set_plot_runtime(bool newval){do_plot_runtime = newval;};
  void Set_ylim(double pmin, double pmax, double vmin, double vmax);
  */
  string fname;
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
  bool do_plot_runtime;
  bool ylim_isset;
  double lim_pbottom, lim_pup, lim_vbottom, lim_vup;

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

class Reservoir: public Units
{
private:
  bool save_data;
  string name;
  double vol, a, p, t;
  double beta;
  VectorXd ReservoirODE(const double x, const VectorXd &y, const VectorXd &pars);

  void rk45(double x_0,
            double x_max,
            VectorXd y_0,
            vector<double> &x, vector<vector<double> > &y,
            const VectorXd & pars);

  vector<vector<double> > data;
  vector<double> tmpvec;

public:
  Reservoir(const string _nev,
            const double _vol,
            const double _a,
            const bool save_data);
  ~Reservoir();
  string GetName();
  string Info();
  void Update(double delta_t, double mp_in, double mp_out);
  void UpdateDimlessPars(double pref, double mp_nevl, double omega);
  void Ini();
  void Ini(double _pstart);
  double Get_dprop(string prop_string);
  void Set_dprop(string prop_string, double val);
  double GetP() {return p;}
  void List_data();
  void Save_data();
  vector<double> Get_dvprop(string prop_string);
  //void Plot();
  string fname;
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
  double Get_MassFlowRate_Compressible_Choked(double p_upstream, double rho_upstream, double kappa, double x);
  double Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho, double x) ;
  
};
