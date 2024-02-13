#define _USE_MATH_DEFINES
#include "SCP.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "my_tools.h"
#include "PSToolboxBaseEdge.h"

#include <cmath>

//! Constructor for the slightly compressible pipe class
SCP::SCP(const string _name, //!< [in] Name of the slightls compressible pipe
    const string _cspe_name, //!< [in] Name of the previous node
    const string _cspv_name, //!< [in] Name of the next node
    const double _ro, //!< [in] Density of the fluid in the system (Constant as it is incompressible)
    const double _a, //!< [in] Speed of sound in the system
    const double _L, //!< [in] Length of the pipe section
    const double _D, //!< [in] Diameter of the pipe
    const double _lambda, //!< [in] Friction loss factor
    const double _he, //!< [in] Height at the beginning [eleje] of the pipe
    const double _hv, //!< [in] Height at the end [v�ge] of the pipe
    const bool _save_data /**< [in] whether to save data*/ ) : PSToolboxBaseEdge(_name), Units() {
  save_data = _save_data;
  //name = _name;
  node_from = _cspe_name;
  node_to = _cspv_name;
  ro = _ro;
  a = _a; //Speed of sound;
  roa = ro * a;
  L = _L;
  D = _D;
  he = _he;
  hv = _hv;
  S0 = -(hv - he) / L; //slope
  A = D * D * M_PI / 4.; //Cross sectional area of the pipe
  lambda = _lambda;
  lambda_p_2D = lambda / 2 / D;
  ini_done = false;
  fname = name + ".dat"; //name of the file for saving data about the pipe
                         //do_plot_runtime = false;
                         //ylim_isset = false;
  g = 9.81;

  // Dummy initialization
  Npts = 1;
  gamma = 0.;
  alpha = 0.;
  phi = 0.;
  mu = 0.;
  t = 0.;
  dt = 0.;

}

//! Desctructor of the SCP class
SCP::~SCP() {}

/*! \brief Gives name of the pipe
  Returns the name of the SCP pipe
  \return name of the pipe
  */
string SCP::GetName() {
  return name;
}
/*! \brief Return the basic info of the pipe
  Returns the basic info about the pipe section. If so set, the data within the points
  is also displayed.
  \param show_pts Whether to show the data for each point or not
  \return a string containing all the pipe data, in a printable form.
  */
string SCP::Info() {
  bool show_pts=false;
  if (!ini_done) //If not initialized, it will do it
    Ini(1);

  std::ostringstream oss;
  oss << "\n\n Pipe name     : " << name;
  oss << "\n================================";
  oss << "\n      node_from : " << node_from;
  oss << "\n        node_to : " << node_to;
  oss << "\n              L : " << L << " m = " << L*m_to_inch << " in";
  oss << "\n              D : " << D << " m= " << D*m_to_inch << " in";
  oss << "\n             he : " << he << " m";
  oss << "\n             hv : " << hv << " m";
  oss << "\n             S0 : " << S0 << " -";
  oss << "\n             ro : " << ro << " kg/m^3";
  oss << "\n              a : " << a << " m/s";
  oss << "\n             dt : " << dt << " s";
  oss << "\n          f(2L) : " << a / (2. * L) << " Hz";
  oss << "\n         f(L/2) : " << a / (0.5 * L) << " Hz";
  oss << endl;
  oss << "\n            phi : " << phi << " -";
  oss << "\n          alpha : " << alpha << " -";
  oss << "\n          gamma : " << gamma << " -";
  oss << "\n             mu : " << mu << " -";
  //oss << "\n f_pipe/f_valve : " << 2.*L/a/(omega/(2.*M_PI)) << " -";
  oss << "\n f_pipe/f_valve : " << M_PI / gamma << " -";

  if (show_pts) { //In case of data needed from the data points, a
    oss << "\n\n # of grid pts.: " << Npts;

    oss << "\n       node #  : ";
    for (int i = 0; i < Npts; i++)
      oss << setw(6) << setfill(' ') << i << " ";

    oss << "\n       p (bar) : ";
    for (int i = 0; i < Npts; i++)
      oss << setw(6) << setprecision(3) << p(i) / 1.e5 << " ";

    oss << "\n       v (m/s) : ";
    for (int i = 0; i < Npts; i++)
      oss << setw(6) << setprecision(3) << v(i) << " ";
  }
  oss << endl;
  return oss.str();
}

void SCP::Ini(){
  Ini(1);
}

/*! \brief Initialize the SCP pipe
  Initializes the SCP pipe with a given number of points at 1 bar pressure and 0 velocity.
  \param Npts_mul Base number of points to create, actually 20 as many created
  */
void SCP::Ini(int Npts_mul) {
  t = 0.;

  Npts = 20 * Npts_mul; //20 times as many points created
  x = VectorXd::Zero(Npts);
  p = VectorXd::Zero(Npts);
  v = VectorXd::Zero(Npts);
  dt = L / (Npts - 1) / a; //timestep according to the CFL condition

  for (int i = 0; i < Npts; i++) {
    x(i) = i * L / (Npts - 1);
    p(i) = 1.e5;
    v(i) = 0.0;
  } //sets 1 bar pressure at 0 velocity in every location
  ini_done = true; //Save that initialization happened.

  tmpvec.push_back(t);
  tmpvec.push_back(p(0));
  tmpvec.push_back(p(Npts - 1));
  tmpvec.push_back(v(0));
  tmpvec.push_back(v(Npts - 1));
  tmpvec.push_back(v(0)*A * ro);
  tmpvec.push_back(v(Npts - 1)*A * ro);
  data.clear();
  data.reserve(100);
  data.push_back(tmpvec);
}

/*! \brief Initialize the SCP pipe
  Initializes the SCP pipe with 20 grid points and initial uniform flow velocity and pressure. (The generated pressure drops with friction.)
  \param vini Initial (uniform) velocity in the pipe
  \param pstart Pressure at the inlet, set up lambda to decrease with friction
  */
void SCP::Ini(double vini, double pstart) {
  t = 0.;

  Npts = 20;

  x = VectorXd::Zero(Npts);
  p = VectorXd::Zero(Npts);
  v = VectorXd::Zero(Npts);
  dt = L / (Npts - 1) / a; //Time step according to the CFL condition

  for (int i = 0; i < Npts; i++) {
    x(i) = i * L / (Npts - 1);
    p(i) = pstart - lambda * x(i) / D * ro / 2.* vini * abs(vini) + S0 * g * ro * x(i);
    v(i) = vini; //The velocity is the same as density is constant
  }
  ini_done = true;

  tmpvec.push_back(t);
  tmpvec.push_back(p(0));
  tmpvec.push_back(p(Npts - 1));
  tmpvec.push_back(v(0));
  tmpvec.push_back(v(Npts - 1));
  tmpvec.push_back(v(0)*A * ro);
  tmpvec.push_back(v(Npts - 1)*A * ro);
  data.clear();
  data.reserve(100);
  data.push_back(tmpvec);
}

/*! \brief Initialize the SCP pipe
  Initializes the SCP pipe with 20*Npts_mul grid points, initial uniform velocity and 
  pressure at the inlet. (The generated pressure drops with friction.)
  \param vini Initial (uniform) velocity in the pipe
  \param pstart Pressure at the inlet, set up lambda to decrease with friction
  \param Npts_mul Multiplier of 20 (e.g. for Npts_mul=2, Ngrid=40)
  */
void SCP::Ini(double vini, double pstart, int Npts_mul) {
  t = 0.;

  Npts = 20 * Npts_mul;

  x = VectorXd::Zero(Npts);
  p = VectorXd::Zero(Npts);
  v = VectorXd::Zero(Npts);
  dt = L / (Npts - 1) / a; //Time step according to the CFL condition

  for (int i = 0; i < Npts; i++) {
    x(i) = i * L / (Npts - 1);
    p(i) = pstart - lambda * x(i) / D * ro / 2.* vini * abs(vini) + S0 * g * ro * x(i);
    v(i) = vini; //The velocity is the same as density is constant
  }
  ini_done = true;

  tmpvec.push_back(t);
  tmpvec.push_back(p(0));
  tmpvec.push_back(p(Npts - 1));
  tmpvec.push_back(v(0));
  tmpvec.push_back(v(Npts - 1));
  tmpvec.push_back(v(0)*A * ro);
  tmpvec.push_back(v(Npts - 1)*A * ro);
  data.clear();
  data.reserve(100);
  data.push_back(tmpvec);
}

/*! \brief Initialize the SCP pipe
  Initializes the SCP pipe with number of grid points that satisfy the given timestep
  Initial uniform velocity distribution and a
  pressure at the inlet. (The generated pressure drops with friction.)
  \param vini Initial (uniform) velocity in the pipe
  \param pstart Pressure at the inlet, set up lambda to decrease with friction
  \param dt_target Target timestep, defines the grid
  */
void SCP::Ini(double vini, double pstart, double dt_target) {
  t = 0.;

  Npts = round(L / a / dt_target); // CFL condition reorganized
                                   //printf("\n\n L=%5.2f m, a=%5.1f m/s, dt_target=%5.3e s, Npts=%d ", L, a, dt_target, Npts);
  if (Npts < 20) {
    Npts = 20; // a minimum of 20 points is set
               //printf(" -> %d\n", Npts);
  }
  //cin.get();

  x = VectorXd::Zero(Npts);
  p = VectorXd::Zero(Npts);
  v = VectorXd::Zero(Npts);
  dt = L / (Npts - 1) / a;

  for (int i = 0; i < Npts; i++) {
    x(i) = i * L / (Npts - 1);
    p(i) = pstart - lambda * x(i) / D * ro / 2.* vini * abs(vini);
    v(i) = vini;
  }
  ini_done = true;

  tmpvec.push_back(t);
  tmpvec.push_back(p(0));
  tmpvec.push_back(p(Npts - 1));
  tmpvec.push_back(v(0));
  tmpvec.push_back(v(Npts - 1));
  tmpvec.push_back(v(0)*A * ro);
  tmpvec.push_back(v(Npts - 1)*A * ro);
  data.clear();
  data.reserve(100);
  data.push_back(tmpvec);
}

/*! \brief Calculates the dimensionless parameters
  Calculates the dimensionless parametwers of the pipe using the input parameters as a base.
  While these parametersa have a usual definitions, here values need to be passed, and as such any can be used
  for non-dimensionalization.
  \param pref [in] Reference pressure
  \param mp_nevl [in] Design mass flow rate
  \param omega [in] Referency frequency, ususally the natural frequency of the valve
  \param xref [in] Referecne length, usually equals the compession the reference (atmospheric) pressure would excert on the valve
  \param m [in] reference mass, ususally the moving mass of the spring-mass-damper system in the valve.
  */
void SCP::UpdateDimlessPars(double pref, double mp_nevl, double omega, double xref, double m) {
  phi   = lambda * xref / 2. / D;
  alpha = ro * A * a / m / omega;
  gamma = L * omega / a;
  mu    = ro * A * xref * omega / mp_nevl;
}

/*! \brief Return a property of the SCP.
  \param prop_string: A string describing the needed property ("L" | "L_feet" | "D" | "D_inch" | "A" | "dt" | "p_front" | "p_back" | "v_front" | "v_back" | "mp_front" | "mp_back" | "frek" | "tnext" | "lambda" | "phi" | "mu" | "rho" | "a")
  \return The value that was looked up.
  */
double SCP::Get_dprop(string prop_string) {
  double out=0.0;
  if (prop_string == "L")
    out = L;
  else if (prop_string == "L_feet")
    out = L * m_to_ft;
  else if (prop_string == "D")
    out = D;
  else if (prop_string == "D_inch")
    out = D * m_to_inch;
  else if (prop_string == "A")
    out = A;
  else if (prop_string == "dt")
    out = dt;
  else if (prop_string == "p_front")
    out = p(0);
  else if (prop_string == "p_back")
    out = p(Npts - 1);
  else if (prop_string == "v_front")
    out = v(0);
  else if (prop_string == "v_back")
    out = v(Npts - 1);
  else if (prop_string == "mp_front")
    out = v(0) * ro * A;
  else if (prop_string == "mp_back")
    out = v(Npts - 1) * ro * A;
  else if (prop_string == "frek")
    out = a / (2.*L);
  else if (prop_string == "tnext")
    out = t + dt;
  else if (prop_string == "lambda")
    out = lambda;
  else if (prop_string == "phi")
    out = phi;
  else if (prop_string == "alpha")
    out = alpha;
  else if (prop_string == "gamma")
    out = gamma;
  else if (prop_string == "mu")
    out = mu;
  else if ((prop_string == "ro") || (prop_string == "rho"))
    out = ro;
  else if (prop_string == "a")
    out = a;
  else {
    cout << endl
      << "ERROR! SCP::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl
      << endl;
    cout << endl << "Name of pipe: " << name << endl;
    cin.get();
  }
  return out;
}

/*! \brief Set a parameter
  Allows to reset the diameter or length the pipe. Updating the length will rescale the system.
  \param prop_string [in] The paramtere to change. Accepted values are "L" and "D"
  \param val [in] The new value of the set parameter, in SI units.
  \sa Get_dprop
  */
void SCP::Set_dprop(string prop_string, double val) {
  if (prop_string == "D") {
    D = val;
  } else if (prop_string == "L") {
    L = val;
    dt = L / (Npts - 1) / a;
  } else {
    cout << endl
      << "HIBA! Cso::Set_dprop(prop_string), ismeretlen bemenet: prop_string=" << prop_string << endl
      << endl;
  }
}

void SCP::Step(){

  double a, b;
  VectorXd pnew = VectorXd::Zero(Npts); //new pressure vector
  VectorXd vnew = VectorXd::Zero(Npts); //new velocity vector
  for (int i = 1; i < Npts - 1; i++) {
    a = (p(i - 1) + roa * v(i - 1)) + dt * roa * Source(i - 1); //
    b = (p(i + 1) - roa * v(i + 1)) - dt * roa * Source(i + 1);
    pnew(i) = (a + b) / 2.;
    vnew(i) = (a - b) / 2. / roa;
  }

  for (int i = 1; i < Npts-1; i++) {
    p(i) = pnew(i);
    v(i) = vnew(i);
  }

  t+=dt;

}

/*! \brief Take a timestep.
  Calculates a timestep in the system dependent on the boundary conditions. (Those are given by the same strings
  and values as the ones in BCLeft and BCRight) Accepted boundary condition types: "Pressure" & "Velocity"
  \param [in] BC_start_type the type of the boundary condition at the beginning of the pipe ("Pressure" | "Velocity")
  \param [in] BC_start_val The value associated with the boundary condition at the start of the pipe
  \param [in] BC_end_type the type of the boundary condition at the end of the pipe ("Pressure" | "Velocity")
  \param [in] BC_end_val The value associated with the boundary condition at the end of the pipe
  \sa BCLeft, BCRight
  */
void SCP::Step(
    string BC_start_type, double BC_start_val,
    string BC_end_type, double BC_end_val) {

  Step();

  BCLeft(BC_start_type, BC_start_val, p(0), v(0));

  BCRight(BC_end_type, BC_end_val, p(Npts - 1), v(Npts - 1));

  //	t += dt;

  if (save_data) {
    tmpvec.at(0) = t;
    tmpvec.at(1) = p(0);
    tmpvec.at(2) = p(Npts - 1);
    tmpvec.at(3) = v(0);
    tmpvec.at(4) = v(Npts - 1);
    tmpvec.at(5) = v(0) * A * ro;
    tmpvec.at(6) = v(Npts - 1) * A * ro;

    data.push_back(tmpvec);
  }
}
/*! \brief Left side bounbdary condition.
  The left side boundary condition, with multiple possible formations. Allows for a possibility of set
  pressure or velocity at the edge. If a wrong type is chosen an error message is displayed.
  \param type [in] Parameter to set at the boundary. Accepted values are "Pressure" and "Velocity".
  \param val [in] Value to set at the edge. Either in [Pa] or [m/s] dependent on the type.
  \param pp [out] The pressure at the left side. Modified in place.
  \param vv [out] The velocity at the left side. Modified in place.
  */
void SCP::BCLeft(string type, double val, double & pp, double & vv) {

  // pp - roa*vv = b
  double beta = (p(1) - roa * v(1)) - dt * roa * Source(1); //Constant value

  if (type == "Pressure") {
    pp = val;
    vv = (pp - beta) / roa;
  }
  else if (type == "Velocity") {
    vv = val;
    pp = beta + vv * roa;
  }
  else {
    cout << endl
      << "ERROR! SCP::BCLeft(), unknown BC type: " << type << endl
      << endl;
    cout << endl << "Name of pipe: " << name << endl;
    cin.get();
  }
}

void SCP::Set_BC_Left(string type, double val){
  BCLeft(type,val,p(0),v(0));
}

void SCP::Set_BC_Right(string type, double val){
  BCRight(type,val, p(Npts - 1), v(Npts - 1));
}

/*! \brief Right side boundary condition.
  The right side boundary condition, with multiple possible formations. Allows for a possibility of set
  pressure or velocity, or a valve at the edge. If a wrong type is chosen an error message is displayed.
  \param type [in] Parameter to set at the boundary. Accepted values are "Pressure" and "Velocity" and "Valve".
  \param val [in] Value to set at the edge. Either in [Pa] or [m/s] dependent on the type for pressure and velocity boundary conditions. For the valve it is Cd*Aflow(x)
  \param pp [out] The pressure at the right side. Modified in place.
  \param vv [out] The velocity at the right side. Modified in place.
  */
void SCP::BCRight(string type, double val, double & pp, double & vv) {

  // pp + roa*vv = alpha
  double alpha = (p(Npts - 2) + roa * v(Npts - 2)) + dt * roa * Source(Npts - 2);

  if (type == "Pressure") {
    pp = val;
    vv = (alpha - pp) / roa;
  }
  else if (type == "Velocity") {
    vv = val;
    // double a = (p(Npts - 2) + roa * v(Npts - 2)) + dt * roa * Source(Npts - 2);
    pp = alpha - vv * roa;
  }
  else if (type == "Valve") {
    double mul = val;
    //		double pend = p(Npts - 1);
    //		double vend = v(Npts - 1);
    // double alpha = (p(Npts - 2) + roa * v(Npts - 2)) + dt * roa * Source(Npts - 2);

    double a_ = ro * ro * A * A;
    double b_ = mul * mul * 2. * ro * roa;
    double c_ = -mul * mul * 2. * ro * alpha;
    if (b_ * b_ < 4.*a_ * c_) {
      cout << endl << "ERROR: SCP:BCRight(), Valve BC: Backflow!!!";
      //cin.get();
    }
    double vuj = 1. / 2. / a_ * (-b_ + sqrt(b_ * b_ - 4 * a_ * c_));
    double puj = alpha - roa * vuj;
    pp = puj;
    vv = vuj;
  }
  else {
    cout << endl
      << "ERROR! SCP::BCRight(), unknown BC type: " << type << endl
      << endl;
    cout << endl << "Name of pipe: " << name << endl;
    cin.get();
  }
}

/* \brief
   \param t_target the time we would like to step to. Needs to be set correctly,
   as if there is an incorrect setting value it will not work properly
   */
double SCP::GetAlphaAtEnd(double t_target) {
  double delta_t = t_target - t;
  double TOL = dt / 1000.;

  if (delta_t < 0) {
    if (fabs(delta_t) < TOL)
      delta_t = 0.;
    else {
      cout << endl
        << "ERROR! SCP::GetAlphaAtEnd(), delta_t = " << delta_t << " < 0 ! (TOL=" << TOL << ")" << endl;
      cout << endl << "Name of pipe: " << name << endl;
      cin.get();
    }
  }
  if (delta_t > dt) {
    if (delta_t - TOL < dt)
      delta_t = dt; //If it was within tolerance repair
    else { //if target time is too big, an error occurs
      cout << endl
        << "ERROR! SCP::GetAlphaAtEnd(), delta_t = " << delta_t << " > dt= " << dt << endl;
      cout << endl << "Name of pipe: " << name << endl;
      cin.get();
    }
  }

  double pp = p(Npts - 2) * delta_t / dt + p(Npts - 1) * (1. - delta_t / dt);
  double vv = v(Npts - 2) * delta_t / dt + v(Npts - 1) * (1. - delta_t / dt);
  double ss = Source(Npts - 2) * delta_t / dt + Source(Npts - 1) * (1. - delta_t / dt);
  double a = (pp + roa * vv) + dt * roa * ss;

  return a;

}

double SCP::GetAlphaPrimitiveAtEnd(double t_target){
  return GetAlphaAtEnd(t_target);
}

double SCP::GetBetaAtFront(double t_target) {
  double delta_t = t_target - t;
  double TOL = dt / 1000.;

  if (delta_t < 0) {
    if (fabs(delta_t) < TOL)
      delta_t = 0.;
    else {
      cout << endl << "Name of pipe: " << name << endl;
      cout << endl
        << "ERROR! SCP::GetBetaAtFront(), delta_t = " << delta_t << " < 0  ! (TOL=" << TOL << ")" << endl;
      cout << endl << " t_pipe = " << t << ", t_target=" << t_target << endl;
      cin.get();
    }
  }
  if (delta_t > dt) {
    if (delta_t - TOL < dt)
      delta_t = dt;
    else {
      cout << endl
        << "ERROR! SCP::GetBetaAtFront(), delta_t = " << delta_t << " > dt= " << dt << endl;
      cout << endl << "Name of pipe: " << name << endl;
      cin.get();
    }
  }

  double pp = p(1) * delta_t / dt + p(0) * (1. - delta_t / dt);
  double vv = v(1) * delta_t / dt + v(0) * (1. - delta_t / dt);
  double ss = Source(1) * delta_t / dt + Source(0) * (1. - delta_t / dt);
  double b = (pp - roa * vv) - dt * roa * ss;

  return b;

}
double SCP::GetBetaPrimitiveAtFront(double t_target) {
  return GetBetaAtFront(t_target);
}

/*! \brief Calculates local base pressure
  Calculates local base pressure from the slope and the
  \param i [in]
  \return Base pressure in the given location.
  */
double SCP::Source(int i) {
  return (g * S0 - lambda_p_2D * v(i) * abs(v(i)));
}

/*! \brief Exports savied data
  Exports data saved from previous settings (such as steps and initialization.) Only works if the
  save_data flag was set to true, otherwise an error message is displayed. Data is saved to the previously determined
  savefile name. No inputs/outputs.
  */
void SCP::Save_data() {
  //char fname [50];
  //sprintf (fname, "%s.dat", name.c_str());

  if (!save_data) {
    cout << endl << "WARNING! SCP: " << name;
    cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
  }
  else {
    cout << endl << "Saving to " << fname.c_str() << " ... ";

    FILE * pFile;
    pFile = fopen (fname.c_str(), "w");
    fprintf(pFile, "t (s); p(0) (bar); p(L) (bar); v(0) m/s; v(L) (m/s); mp(0) (kg/s), mp(L) (kg/s), L, D, lambda\n");
    for (int i = 0; i < data.size(); i++)
      fprintf(pFile, "%8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e\n",
          data.at(i).at(0),
          data.at(i).at(1) / 1.e5, data.at(i).at(2) / 1.e5,
          data.at(i).at(3), data.at(i).at(4),
          data.at(i).at(5), data.at(i).at(6),
          L, D, lambda);
    fclose (pFile);
    cout << " done. ";
  }
}

vector<double> SCP::Get_dvprop(string prop_string) {
  int Ntime = data.size();
  //	int Nvars = data.at(0).size();
  //cout<<endl<<"Ntime="<<Ntime<<", Nvars="<<Nvars<<endl;
  vector<double> out(Ntime);
  if (prop_string == "t")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(0);
  else if (prop_string == "p_front")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(1);
  else if (prop_string == "p_back")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(2);
  else if (prop_string == "p_front_bar")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(1) / 1.e5;
  else if (prop_string == "p_back_bar")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(2) / 1.e5;
  else if (prop_string == "v_front")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(3);
  else if (prop_string == "v_back")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(4);
  else if (prop_string == "mp_front")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(5);
  else if (prop_string == "mp_back")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(6);
  else {
    cout << endl
      << "ERROR! SCP::Get_dvprop(prop_string), unknown input: prop_string=" << prop_string << endl
      << endl;
    cout << endl << "Name of valve: " << name << endl;
    cin.get();
  }
  return out;
}

/*! \brief Export the penultimate values of a given return data
  Returns the penultimate value of the collected data. Only one value, the one specified
  by the input value is returned. Possible values: t,p_front, p_back, v_front, v_back, mp_front, mp_back

  \param The name of the value we would like tor return.
  */
float SCP::GetPenult(string what) {
  int penult = data.size() - 1 - 1; //Maximum allowed place is -1, the one before that
  int loc = 0;
  if (what == "t") { loc = 0; }
  else if (what == "p_front") { loc = 1; }
  else if (what == "p_back") { loc = 2; }
  else if (what == "v_front") { loc = 3; }
  else if (what == "v_back") { loc = 4; }
  else if (what == "mp_front") { loc = 5; }
  else if (what == "mp_back") { loc = 6; }
  else { loc = 0; }

  return data.at(penult).at(loc);
}



void SCP::Save_status(bool newfile, bool atonce) {
  if (!save_data) {
    cout << endl << "WARNING! LWP: " << name;
    cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
  }
  else if (atonce) {
    if (!newfile) { //No new file is created, stuff is added
      tstatus.push_back(t);
      vector<double> ptemp, vtemp, Ttemp;
      for (int i = 0; i < Npts; i++) {
        ptemp.push_back(p(i));
        vtemp.push_back(v(i));
      }
      pstatus.push_back(ptemp);
      vstatus.push_back(vtemp);
    }
    else {//save satus into a file
      FILE* pfile, * Tfile, * vfile;
      string pfname = name + "_pstatus.dat";
      string Tfname = name + "_Tstatus.dat";
      string vfname = name + "_vstatus.dat";
      pfile = fopen(pfname.c_str(), "w"); Tfile = fopen(Tfname.c_str(), "w"); vfile = fopen(vfname.c_str(), "w");

      fprintf(pfile, "# Npts = %d \n", Npts); fprintf(Tfile, "# Npts = %d \n", Npts); fprintf(vfile, "# Npts = %d \n", Npts); //Fist line
      for (int i = 0; i < tstatus.size(); i++) {
        fprintf(pfile, "%8.6e", tstatus.at(i)); fprintf(Tfile, "%8.6e", tstatus.at(i)); fprintf(vfile, "%8.6e", tstatus.at(i));
        for (int j = 0; j < Npts; j++) {
          fprintf(pfile, ";%8.6e", pstatus.at(i).at(j));
          fprintf(vfile, ";%8.6e", vstatus.at(i).at(j));
        }
        fprintf(pfile, "%-s \n", ""); fprintf(Tfile, "%-s \n", ""); fprintf(vfile, "%-s \n", ""); //close line
      }
      fclose(pfile); fclose(Tfile); fclose(vfile);
    }
  }
  else {
    std::cout << "Additional data saved" << std::endl;
    FILE* pfile, * Tfile=NULL, * vfile;
    string pfname = name + "_pstatus.dat";
    string vfname = name + "_vstatus.dat";
    if (newfile) {
      pfile = fopen(pfname.c_str(), "w"); vfile = fopen(vfname.c_str(), "w");
      fprintf(pfile, "# Npts = %d \n", Npts); fprintf(vfile, "# Npts = %d \n", Npts);
      fprintf(pfile, "%8.6e", t); fprintf(Tfile, "%8.6e", t); fprintf(vfile, "%8.6e", t);
      for (int i = 0; i < Npts; i++) { fprintf(pfile, ";%8.6e", p(i)); fprintf(vfile, ";%8.6e", v(i)); }
      fprintf(pfile, "%-s \n", "");  fprintf(vfile, "%-s \n", "");
    }
    else {
      pfile = fopen(pfname.c_str(), "a"); vfile = fopen(vfname.c_str(), "a");
      fprintf(pfile, "%8.6e", t); fprintf(vfile, "%8.6e", t);
      for (int i = 0; i < Npts; i++) { fprintf(pfile, ";%8.6e", p(i)); fprintf(vfile, ";%8.6e", v(i)); } //Avoids hanging semicolon at end
      fprintf(pfile, "%-s \n", ""); fprintf(vfile, "%-s \n", "");
    }
    fclose(pfile); fclose(vfile);
  }
}
void SCP::list_pv(){
  printf("\n\n Listing p (Pa) and v (m/s):");
  for (unsigned int i=0; i<Npts; i++){
    printf("\n\t #%2d: p=%5.3e, v=%5.3e",i,p(i),v(i));
  }
  printf("\n");
}

