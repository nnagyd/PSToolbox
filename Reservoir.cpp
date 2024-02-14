#define _USE_MATH_DEFINES
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

#include "Reservoir.h"
#include "my_tools.h"
//#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"

//! Class for modeling Reservoir 
/*!
  The class contains data for Reservoir
  */

Reservoir::Reservoir(const string _name,
    const double _vol,
    const double _a,
    const double _rho,
    const bool _save_data) : Units() {
  save_data = _save_data;
  name = _name;
  vol = _vol;
  a = _a;
  rho =_rho;
  fname = name + ".dat";
  beta = 0.;
  is_Gas = false;

  t = 0.;
  p = 0.;
  Ini_done = false;
  DEBUG=false;
}

Reservoir::Reservoir(const string _name,
    const double _vol,
    IdealGas* _gas,
    double _n_poly,
    const bool _save_data) : Units() {
  save_data = _save_data;
  name = _name;
  vol = _vol;
  fname = name + ".dat";
  beta = 0.;
  is_Gas = true;
  gas = _gas;

  n_poly = _n_poly;
  a = gas->Get_SonicVel(293.,1.e5);

  T=293.;
  t = 0;
  p = 0.;
  Ini_done = false;
  DEBUG=false;
}

Reservoir::Reservoir(const string _name,
    const double _vol,
    IdealGas* _gas,
    double _n_poly,
    const bool _save_data,
    const double Tt) : Units() {
  save_data = _save_data;
  name = _name;
  vol = _vol;
  fname = name + ".dat";
  beta = 0.;
  is_Gas = true;
  gas = _gas;

  n_poly = _n_poly;
  a = gas->Get_SonicVel(293.,1.e5);

  T=Tt;
  t = 0;
  p = 0.;
  Ini_done = false;
  DEBUG=false;
}



Reservoir::Reservoir(const string _name,
    const double _vol,
    Gas* _gas,
    double _n_poly,
    const bool _save_data) : Units() {
  save_data = _save_data;
  name = _name;
  vol = _vol;
  fname = name + ".dat";
  beta = 0.;
  is_Gas = true;
  gas = _gas;

  n_poly = _n_poly;
  a = gas->Get_SonicVel(293.,1.e5);

  t = 0;
  p = 0.;
  Ini_done = false;
  DEBUG=false;
}


Reservoir::~Reservoir() {}

void Reservoir::SetDebug(const bool _DEBUG){
  DEBUG=_DEBUG;
}

string Reservoir::GetName() {
  return name;
}

void Reservoir::Ini(double pstart) {
  if (is_Gas) {
    cout << endl << "ERROR! Reservoir::Ini(pstart) was called but the fluid is IdealGas" << endl;
    cout << endl << "    -> use Ini(pstart, Tstart) instead." << endl;
    cout << endl << "Name of pipe: " << name << endl;
    cin.get();
  }
  p = pstart;
  t = 0.;

  tmpvec.push_back(t);
  tmpvec.push_back(p);
  tmpvec.push_back(0.); //mp_in
  tmpvec.push_back(0.); //mp_out
  data.clear();
  data.reserve(100);
  data.push_back(tmpvec);
  Ini_done = true;
}

void Reservoir::Ini(double pstart, double Tstart) {
  if (!is_Gas) {
    cout << endl << "ERROR! Reservoir::Ini(pstart,Tstart) was called but the fluid is ** NOT ** IdealGas" << endl;
    cout << endl << "    -> use Ini(pstart) instead." << endl;
    cout << endl << "Name of pipe: " << name << endl;
    cin.get();
  }

  p = pstart;
  //	gas.Setp(pstart);
  //	gas.SetT(Tstart);
  a = gas->Get_SonicVel(Tstart,pstart);
  t = 0.;

  tmpvec.push_back(t);
  tmpvec.push_back(p);
  tmpvec.push_back(0.);
  tmpvec.push_back(0.);
  data.clear();
  data.reserve(100);
  data.push_back(tmpvec);

  Ini_done = true;
}


void Reservoir::UpdateDimlessPars(double pref, double mp_nevl, double omega) {
  if (is_Gas)
    a = gas->Get_SonicVel(293.,1.e5);
  beta = a * a * mp_nevl / vol / omega / pref;
}

void Reservoir::Update(double delta_t, double mp_in, double mp_out) {

  if (!Ini_done) {
    cout << endl << endl << "ERROR! Reservoir::Update() was called without initializing first";
    cout << endl << "    -> use Ini() before calling Update()." << endl;
    cout << endl << "Name of pipe: " << name << endl;
    cin.get();
  }

  vector<double> xvec;
  vector<vector<double> > yvec;
  VectorXd pars = VectorXd::Zero(2);
  pars(0) = mp_in;
  pars(1) = mp_out;
  VectorXd y_0 = VectorXd::Zero(1);
  y_0(0) = p;

  rk45(t, t + delta_t, y_0, xvec, yvec, pars);

  p = yvec.back().at(0);  
  t += delta_t;

  if (save_data) {
    // Save only endpoint
    tmpvec.at(0) = t;
    tmpvec.at(1) = p;
  data.push_back(tmpvec);
  }
}

VectorXd Reservoir::ReservoirODE(const double x, const VectorXd & y, const VectorXd & pars){ 
  VectorXd out(1);
  if (is_Gas)
    a = gas->Get_SonicVel(T,p);
  out(0) = a * a / vol * (pars(0) - pars(1));
  //printf("\n reservoir: p=%6.4f bar, a=%5.3e, vol=%5.3e, dq=%5.3e, dp/dt= %5.3f bar/s",y(0)/1.e5,a,vol,pars(0)-pars(1),out(0)/1.e5);
//cin.get();
  return out;
}


void Reservoir::rk45(
    double x_0, double x_max, VectorXd y_0,
    vector<double> &x, vector<vector<double> > &y,
    const VectorXd & pars ) {

  bool debug = false;
  bool is_rk = true;

  int Neq = y_0.size();

  x.clear();
  y.clear();
  x.push_back(x_0);
  std::vector<double> tmp(Neq);
  for (int i = 0; i < Neq; i++)
    tmp.at(i) = y_0[i];
  y.push_back(tmp);

  double c2 = 1. / 4., a21 = 1. / 4.;
  double c3 = 3. / 8., a31 = 3. / 32., a32 = 9. / 32.;
  double c4 = 12. / 13., a41 = 1932. / 2197., a42 = -7200. / 2197., a43 = 7296. / 2197.;
  double c5 = 1., a51 = 439. / 216., a52 = -8., a53 = 3680. / 513., a54 = -845. / 4104;
  double c6 = 1. / 2., a61 = -8. / 27., a62 = 2., a63 = -3544. / 2565., a64 = 1859. / 4104., a65 = -11. / 40.;
  double b1 = 16. / 135., b2 = 0., b3 = 6656. / 12825., b4 = 28561. / 56430., b5 = -9. / 50., b6 = 2. / 55.;
  double d1 = 25. / 216., d2 = 0., d3 = 1408. / 2565., d4 = 2197. / 4104., d5 = -1. / 5., d6 = 0.;

  //bool is_impact = false;
  bool last_step = false;
  double hiba_max = 1e-9, dx_min = 1.e-10;

  double hiba = 2.0 * hiba_max;
  double xx = x_0, dx = (x_max - x_0) / 100.;

  VectorXd yy = y_0;
  bool retake_step = true;

  if (debug) {
    printf("\n\nIntegrating...\n\t t=%7.5e", xx);
    for (int iii = 0; iii < Neq; iii++)
      printf(", y(%d)=%7.5e", iii, yy[iii]);
  }
  while (!last_step) {
    retake_step = true;

    VectorXd y1(Neq), y2(Neq), k1(Neq), k2(Neq), k3(Neq), k4(Neq), k5(Neq), k6(Neq);
    hiba = 2.0 * hiba_max;

    if (debug)
      cin.get();

    while (retake_step) {
      // Is this the last step?
      if ((dx > 0) && (xx + dx > x_max)) {
        dx = x_max - xx;
        retake_step = false;
        last_step = true;
      }
      if ((dx < 0) && (xx + dx < x_max)) {
        dx = (x_max - xx);
        retake_step = false;
        last_step = true;
      }

      if (is_rk) {
        // RK45
        k1 = ReservoirODE(xx, yy, pars);
        k2 = ReservoirODE(xx + c2 * dx, yy + a21 * dx * k1, pars);
        k3 = ReservoirODE(xx + c3 * dx, yy + a31 * dx * k1 + a32 * dx * k2, pars);
        k4 = ReservoirODE(xx + c4 * dx, yy + a41 * dx * k1 + a42 * dx * k2 + a43 * dx * k3, pars);
        k5 = ReservoirODE(xx + c5 * dx, yy + a51 * dx * k1 + a52 * dx * k2 + a53 * dx * k3 + a54 * dx * k4, pars);
        k6 = ReservoirODE(xx + c6 * dx, yy + a61 * dx * k1 + a62 * dx * k2 + a63 * dx * k3 + a64 * dx * k4 + a65 * dx * k5, pars);
        y1 = yy + dx * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6);
        y2 = yy + dx * (d1 * k1 + d2 * k2 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6);
      }
      else {
        // Heun
        k1 = ReservoirODE(xx, yy, pars);
        y1 = yy + dx * k1;
        y2 = yy + dx / 2 * (k1 + ReservoirODE(xx + dx, y1, pars));
      }

      // compare results
      hiba = (y1 - y2).norm() / y1.norm();
      if (debug) {
        printf("\n\nIntegrating...\n\t t=%7.5e, hiba=%7.5e", xx, hiba);
        for (int iii = 0; iii < Neq; iii++)
          printf(", y(%d)=%7.5e", iii, yy[iii]);
      }


      // Lepes beallitasa
      if (hiba > hiba_max) {
        dx /= 2.;
        retake_step = true;
        // printf("\n\t dx -> dx/2 (dx=%g)", dx);
        if (fabs(dx) < dx_min) {
          last_step = true;
          cout << endl << endl << "!!!! PANIC dxuj<dx_min !!!!" << endl << endl;
          exit(-1);
        }
      }
      else {
        retake_step = false;
        yy = y1;
        xx += dx;
        
        x.push_back(xx); //enter new time we stepped to
                         //Create a vector of arrays and save the y-s we had
        for (int iii = 0; iii < Neq; iii++)
          tmp.at(iii) = yy[iii];
        y.push_back(tmp);
        if (debug) //display succesful stepping
          printf(" ok.");
        if (hiba < hiba_max / 10.) { //if the error is really small, use a bigger timestep
          dx *= 2.;
          // printf("\n\t dx -> dx*2 (dxuj=%g)", dx);
        }
      }
    }
  }
  if (debug) { //report the y-s in the
    printf("\n\nIntegrating...\n\t t=%7.5e", xx);
    for (int iii = 0; iii < Neq; iii++)
      printf(", y(%d)=%7.5e", iii, yy[iii]);
    printf("\n Finished integration.\n");
  }
}

void Reservoir::List_data() {
  cout << endl << endl << "Reservoir " << name << endl;
  for (int i = 0; i < data.size(); i++)
    printf("t=%5.3f s, p=%5.3f bar\n", data.at(i).at(0), data.at(i).at(1) / 1.e5);
  cout << "Number of data points: " << data.size() << endl;
}

void Reservoir::Save_data() {
  //char fname [50];
  //sprintf (fname, "%s.dat", name.c_str());
  if (!save_data) {
    cout << endl << "WARNING! Reservoir: " << name;
    cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
  }
  else {
    cout << endl << "Saving to " << fname.c_str() << " ... ";

    FILE * pFile;
    pFile = fopen (fname.c_str(), "w");
    fprintf(pFile, "t (s); p (bar); mp_in (kg/s); mp_out (kg/s)\n");
    for (unsigned int i=0; i<data.size(); i++)
      fprintf(pFile, "%8.6e; %8.6e, %8.6e, %8.6e\n", data.at(i).at(0), data.at(i).at(1),data.at(i).at(2), data.at(i).at(3));

    fclose (pFile);

    cout << " done.";
  }
}

double Reservoir::Get_dprop(string prop_string) {
  double out = 0.;
  if ((prop_string == "Vres") || (prop_string == "V"))
    out = vol;
  else if (prop_string == "a") {
    if (is_Gas)
      out = gas->Get_SonicVel(293,p);
    else
      out = a;
  }
  else if (prop_string == "T") {
    if (is_Gas)
      out = 293.;
    else {
      if (DEBUG){
        cout << endl << "WARNING!";
        cout << endl << "Name of reservoir: " << name;
        cout << endl << "This reservoir is ** NOT ** defined with liquid content and not gas -> no temperature data is available, returning 293K." << endl;
      }
      out=293.;
    }
  }
  else if (prop_string == "p")
    out = p;
  else if (prop_string == "rho")
    out = gas->Get_rho(p,293.);
  else if (prop_string == "beta")
    out = beta;
  else {
    cout << endl
      << "ERROR! Reservoir::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl
      << endl;
    cout<<endl<<"Possible options: Vres,V | a | T | p | rho | beta";
    cout << endl << "Name of reservoir: " << name << endl;
    cin.get();
  }
  return out;
}

void Reservoir::Set_dprop(string prop_string, double val) {
  if (prop_string == "Vres")
    vol = val;
  else {
    cout << endl
      << "HIBA! Reservoir::Set_dprop(prop_string), ismeretlen bemenet: prop_string:" << prop_string << endl
      << endl;
  }
}

string Reservoir::Info() {

  std::ostringstream oss;
  oss << "\n\n Reservoir name : " << name;
  oss << "\n================================";
  oss << "\n          vol : " << vol << " m^3 = " << vol*m3_to_inch3 << " in^3";
  if (is_Gas)
    a = gas->Get_SonicVel(293,1.e5);

  oss << "\n            a : " << 	a << " m/s (T=293K, p=1bar)";
  oss << "\n            p : " << 	p/1.e5<<" bar";
  oss << "\n       is_Gas : " << 	is_Gas;
  //oss << "\n       beta : " << beta << " -";

  return oss.str();
}


vector<double> Reservoir::Get_dvprop(string prop_string) {
  int Ntime = data.size();
  //int Nvars = data.at(0).size();
  //cout<<endl<<"Ntime="<<Ntime<<", Nvars="<<Nvars<<endl;
  vector<double> out(Ntime);
  if (prop_string == "t")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(0);
  else if (prop_string == "p")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(1);
  else if (prop_string == "p_bar")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(1) / 1.e5;
  else if (prop_string == "mp_in")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(2);
  else if (prop_string == "mp_out")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(3);
  else {
    cout << endl;
    cout<<"Valid options: t | p | p_bar | mp_in | mp_out"<<endl;
    cout << endl << "Name of valve: " << name << endl;
    cin.get();
  }
  return out;
}

