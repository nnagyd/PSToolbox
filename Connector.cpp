#define _USE_MATH_DEFINES

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "my_tools.h"
#include "LWP.h"
#include "SCP.h"
#include "Valve.h"
#include "Reservoir.h"
#include "Connector.h"
#include <Eigen/Dense>

//#include <math.h>
#include <cmath>

// ! Default Constructor
Connector::Connector(bool _DEBUG) {
  DEBUG=_DEBUG;
}

// ! Default Destructor
Connector::~Connector() {}

//! Connect Reservoir to Valve
/*! Connect Reservoir to Valve
  \param t_target, s
  \param Reservoir* reservoir
  \param Valve* valve
  \param p_downstream
  \param &p pressure @ connection
  \param &T temperature @ connection
  */

void Connector::Connector_Reservoir_and_Valve_Inlet(double t_target, 
    Reservoir* r1, Valve *v1, bool INLET_PRESSURE_DROP, double p_downstream, double& p, double& mp){

  // Solve the following system for p,T,rho,v,a:
  // (1) mp = MassFlowRate(pt,Tt,pd,Tout,x)
  // (2) rho=rho(p,T)gt
  // (3) a  =a(p,T)
  // (4a) T/Tr=(pr/p)^(kappa-1)/kappa (gas)
  // (4b) pr=p+rho/2*v^2 (liquid)

  double pr    = r1->Get_dprop("p");
  double Tr    = r1->Get_dprop("T");
  double Abore = v1->Get_dprop("A");

  double rho, v, T, p_new, cp, cV, kappa;

  double err=1., TOL = 1e-3;
  int iter = 0, MAX_ITER = 10;
  p = r1->Get_dprop("p");
  if (r1->is_Gas)
    T = r1->Get_dprop("T");
  else 
    T=293.;

  if (DEBUG)
    cout<<endl<<"Connector::Connector_Reservoir_and_Pipe_Front -> entering debug mode"<<endl;

  while ((err > TOL) && (iter<=MAX_ITER)) {
    rho = r1->Get_dprop("rho");
    mp  = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
    v   = mp / rho / Abore;
    if (r1->is_Gas){
      cp    = r1->gas->Get_cp(p,T);
      cV    = r1->gas->Get_cV(p,T);
      kappa = cp/cV;
      T   = r1->gas->Get_T(p,rho);
      if (INLET_PRESSURE_DROP)
        p_new  = pr/pow(T/Tr,kappa/(kappa-1));
      else
        p_new=pr;
    }
    else{
      if (INLET_PRESSURE_DROP)
        p_new = pr-rho/2.*v*v;
      else p_new=pr;
    }
    if (DEBUG)
      printf("\n iter: %2d, p=%5.3e, rho=%5.3e, mp=%5.3e, v=%5.3e, p_new=%5.3e, err=%5.3e",
          iter,p,rho,mp,v,p_new,err);

    err = fabs(p-p_new) / 1.e5;
    p=p_new;
    if (iter==MAX_ITER){
      cout<<endl<<"Connector.Connector_Reservoir_and_Valve_Inlet() ERROR!!!"<<endl;
      cin.get();
    }
  }
}

//! Connect LWP end to Valve
/*! Connect LWP end to Valve
  \param t_target, s
  \param LWP* pipe
  \param Valve* valve
  \param p_downstream
  \param &p pressure @ connection
  \param &T temperature @ connection
  */

bool Connector::Connector_LWP_Pipe_Back_and_Valve(double t_target,
    LWP* p1, Valve* v1, double p_downstream) {
  // Solve the following system for p,T,rho,v:
  // (1) alpha=p+ro*a*v
  // (2) ro*Apipe*v=valve_mass_flow
  // (3) T=TL
  // (4) rho=rho(p,T)
  // (5) a  =a(T)

  //double kappa=p1->gas->kappa;
  //double cp   =p1->gas->cp;

  double update_OK=true;
  double rho, v, a, alpha, Apipe;

  Apipe = p1->Get_dprop("A");
  //p1->GetAllPrimitiveAtEnd(t_target, pL, v, T, rho);
  alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

  double err1, err2, err = 1.e5, mp, TOL = 1e-5;
  int iter = 0, MAX_ITER = 10;
  double p = p1->Get_dprop("p_back");
  double T = p1->Get_dprop("T_back");
  while ((err > TOL) && (iter<=MAX_ITER)) {
    a   = p1->gas->Get_SonicVel(T,p);
    rho = p1->gas->Get_rho(p, T);
    mp  = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
    v   = mp / rho / Apipe;
    p   = alpha - rho * a * v;
    //T   = p1->gas->Get_T(p,rho);
    err1 = (alpha - (p + rho * a * v)) / 1.e5;
    err2 = rho * v * Apipe - mp;
    err = sqrt(err1 * err1 + err2 * err2);

    iter++;
    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_Pipe_Back_and_Valve() -> MAX_ITER reached, stopping iteration.";
      printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, a=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e", iter, alpha, p, T, rho, a, v, err1, err2);
      //cin.get();
      update_OK=false;
    }
  }

  p1->BCRight("StaticPres_Outlet",p,T,true);

  return update_OK;
}

bool Connector::Connector_LWP_Pipe_Back_and_Valve_with_Absorber(double t_target,
    LWP* p1, Valve_with_Absorber* v1, double p_downstream, double& p, double& T) {
  // Solve the following system for p,T,rho,v:
  // (1) alpha=p+ro*a*v
  // (2) ro*Apipe*v=valve_mass_flow
  // (3) T=TL
  // (4) rho=rho(p,T)
  // (5) a  =a(T)

  //double kappa=p1->gas->kappa;
  //double cp   =p1->gas->cp;

  double update_OK=true;
  double rho, v, a, alpha, Apipe;

  Apipe = p1->Get_dprop("A");
  //p1->GetAllPrimitiveAtEnd(t_target, pL, v, T, rho);
  alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

  double err1, err2, err = 1.e5, mp, TOL = 1e-5;
  int iter = 0, MAX_ITER = 10;
  p = p1->Get_dprop("p_back");
  while ((err > TOL) && (iter<=MAX_ITER)) {
    a   = p1->gas->Get_SonicVel(T,p);
    rho = p1->gas->Get_rho(p, T);
    mp  = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
    v   = mp / rho / Apipe;
    p   = alpha - rho * a * v;
    //T   = p1->gas->Get_T(p,rho);
    err1 = (alpha - (p + rho * a * v)) / 1.e5;
    err2 = rho * v * Apipe - mp;
    err = sqrt(err1 * err1 + err2 * err2);

    iter++;
    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_Pipe_Back_and_Valve() -> MAX_ITER reached, stopping iteration.";
      printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, a=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e", iter, alpha, p, T, rho, a, v, err1, err2);
      //cin.get();
      update_OK=false;
    }
  }

  //cout<<endl<<"\n\n exiting Connector...";
  //printf("\n\t p=%5.3e, T=%5.3e, rho=%5.3e, a=%5.3e, v=%5.3e, mp=%+5.3e",p, T, rho, a, v, mp);
  return update_OK;
}


void Connector::Connector_LWP_Reservoir_and_Pipe_Front(double t_target,
    Reservoir* r1, LWP* p1, bool inlet_pressure_drop){
  double pt = r1->Get_dprop("p");
  double Tt = r1->Get_dprop("T");

  double pin,Tin;
  bool is_inflow=false, is_outflow=false;
  if (!inlet_pressure_drop){
    is_inflow=p1->BCLeft("StaticPres_and_StaticTemp_Inlet",pt,Tt,false);
    pin=pt; Tin=Tt;
  }
  else{
    is_inflow=Connector_LWP_Reservoir_and_Pipe_Front_inlet(t_target,r1,p1,inlet_pressure_drop,pt,Tt);
    pin=pt;
    Tin=Tt;
    // Old pt, Tt values were overwritten load them again
    pt = r1->Get_dprop("p");
    Tt = r1->Get_dprop("T");
  }

  if (is_inflow)
    is_inflow=p1->BCLeft("StaticPres_and_StaticTemp_Inlet",pin,Tin,true);
  else{
    is_outflow=p1->BCLeft("StaticPres_Outlet",pt,0,false);

    if (is_outflow)
      is_outflow=p1->BCLeft("StaticPres_Outlet",pt,0,true);
    else{
      p1->BCLeft("Wall",0,0,true);
      cout<<endl<<"WARNING: Connector_LWP_Reservoir_and_Pipe_Front() -> Neither inflow, nor outflow! Applying wall BC.";
      if (DEBUG)
        cin.get();
    }

  }
}

bool Connector::Connector_LWP_Reservoir_and_Pipe_Front_inlet(double t_target,
    Reservoir* r1, LWP* p1, bool inlet_pressure_drop,
    double& p, double& T) {
  // Solve the following system for p,T,rho,v:
  // (1)  Tt=T+v^2/2/cp            Isentropic flow from the reservoir to the pipe
  // (2)  p-rho*a*v = beta_front   Charactersistic equation
  // (3)  rho = rho(p,T)
  // (4a) p=pt;
  // (4b) pt/rhot^k = p/rho^k

  if (DEBUG)
    cout<<endl<<endl<< "Connector::Connector_LWP_Reservoir_and_Pipe_Front: entering DEBUG mode"<<endl;

  double cp=1000., kappa_Tv=1.4;
  double pt = r1->Get_dprop("p");
  double Tt = r1->Get_dprop("T");
  //double rhot = p1->gas->Get_rho(pt, Tt);

  double rho, v, a, beta;
  //p1->GetAllPrimitiveAtFront(t_target, p, v, T, rho);
  beta = p1->GetBetaPrimitiveAtFront(t_target);
  //bool is_ok = false;

  // Assume inflow
  v = p1->Get_dprop("v_front");
  if (v<0)
    v=0.;

  double err = 1.e5, TOL = 0.001 /*m/s*/, dv, f, f1, df, vnew;
  int iter = 0, MAX_ITER = 50;
  while (fabs(err) > TOL) {		

    f   = Connector_LWP_Reservoir_and_Pipe_Front_fun(v, beta, r1, p1, inlet_pressure_drop);

    if (fabs(v)>0.001)
      dv = v * 0.01;
    else
      dv=0.001;
    f1 = Connector_LWP_Reservoir_and_Pipe_Front_fun(v + dv, beta, r1, p1,inlet_pressure_drop);

    df = (f1 - f) / dv;
    vnew = v - f / df;

    err = f;

    if (DEBUG)  
      printf("\n (front, inflow) inner iter #%2d: pt=%5.3e, beta=%5.3e, v=%5.3e, err=%5.3e",iter, pt, beta, v, err);

    iter++; 

    double RELAX=0.5;
    v = (1.-RELAX)*v + RELAX*vnew;

    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (1)";
      exit(-1);
    }
  }

  double Tnew;
  iter=0; err=1.e5; T=293;
  while (fabs(err)>0.1){
    cp = p1->gas->Get_cp(p,T);
    kappa_Tv = p1->gas->Get_kappa_Tv();
    Tnew     = Tt-v*v/2/cp;
    if (inlet_pressure_drop)
      p   = pt * pow(Tnew / Tt, kappa_Tv/(kappa_Tv - 1.));
    else 
      p=pt;
    rho = p1->gas->Get_rho(p, Tnew);
    a   = p1->gas->Get_SonicVel(Tnew,p);

    err=T-Tnew;
    T=Tnew;

    if (DEBUG)  
      printf("\n (front, inflow) outer iter #%2d: cp=%5.3e, kappa_Tv=%5.3e, p=%5.3e, T=%5.3e, v=%5.3e, err=%5.3e",iter, cp, kappa_Tv, p, T, v, err);

    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (2)";
      exit(-1);
    }

    iter++;
  }

  bool is_inflow=false;
  if (v>0)
    is_inflow=true;
  return is_inflow;
}

void Connector::Connector_LWP_Reservoir_and_Pipe_Front(double t_target,
    Reservoir* r1, LWP* p1, bool inlet_pressure_drop,
    double& p, double& T) {
  // Solve the following system for p,T,rho,v:
  // (1)  Tt=T+v^2/2/cp            Isentropic flow from the reservoir to the pipe
  // (2)  p-rho*a*v = beta_front   Charactersistic equation
  // (3)  rho = rho(p,T)
  // (4a) p=pt;
  // (4b) pt/rhot^k = p/rho^k

  if (DEBUG)
    cout<<endl<<endl<< "Connector::Connector_LWP_Reservoir_and_Pipe_Front: entering DEBUG mode"<<endl;

  double cp=1000., kappa_Tv=1.4;
  double pt = r1->Get_dprop("p");
  double Tt = r1->Get_dprop("T");
  //double rhot = p1->gas->Get_rho(pt, Tt);

  double rho, v, a, beta;
  //p1->GetAllPrimitiveAtFront(t_target, p, v, T, rho);
  beta = p1->GetBetaPrimitiveAtFront(t_target);
  bool is_ok = false;

  // Assume inflow
  v = 1.;

  double err = 1.e5, TOL = 0.001 /*m/s*/, dv, f, f1, df, vnew;
  int iter = 0, MAX_ITER = 50;
  while (fabs(err) > TOL) {		

    f   = Connector_LWP_Reservoir_and_Pipe_Front_fun(v, beta, r1, p1, inlet_pressure_drop);

    if (fabs(v)>0.001)
      dv = v * 0.01;
    else
      dv=0.001;
    f1 = Connector_LWP_Reservoir_and_Pipe_Front_fun(v + dv, beta, r1, p1,inlet_pressure_drop);

    df = (f1 - f) / dv;
    vnew = v - f / df;

    err = f;

    if (DEBUG)  
      printf("\n (front, inflow) inner iter #%2d: pt=%5.3e, beta=%5.3e, v=%5.3e, err=%5.3e",iter, pt, beta, v, err);

    iter++; 

    double RELAX=0.5;
    v = (1.-RELAX)*v + RELAX*vnew;

    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_LWP_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (1)";
      exit(-1);
    }
  }

  double Tnew;
  iter=0; err=1.e5; T=293;
  while (fabs(err)>0.1){
    cp = p1->gas->Get_cp(p,T);
    kappa_Tv = p1->gas->Get_kappa_Tv();
    Tnew     = Tt-v*v/2/cp;
    if (inlet_pressure_drop)
      p   = pt * pow(Tnew / Tt, kappa_Tv/(kappa_Tv - 1.));
    else 
      p=pt;
    rho = p1->gas->Get_rho(p, Tnew);
    a   = p1->gas->Get_SonicVel(Tnew,p);

    err=T-Tnew;
    T=Tnew;

    if (DEBUG)  
      printf("\n (front, inflow) outer iter #%2d: cp=%5.3e, kappa_Tv=%5.3e, p=%5.3e, T=%5.3e, v=%5.3e, err=%5.3e",iter, cp, kappa_Tv, p, T, v, err);

    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (2)";
      exit(-1);
    }

    iter++;
  }

  if (v < 0) {
    is_ok=false;
    if (DEBUG){
      cout << endl<<endl << " Connector_Reservoir_and_Pipe_Front() -> Assumed inflow, but computed v<0";
      //printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err=%5.3e",
      //	iter, pt, rhot, beta, p, T, rho, v, err);
      cout<<endl<<"Trying outflow";		
    }
  }
  else
    is_ok = true;

  // Need to try again with outflow
  if (!is_ok) {

    //double pM, vM;
    double TM=p1->Get_dprop("T_front");
    double pM=p1->Get_dprop("p_front");
    double rhoM=p1->gas->Get_rho(pM,TM);
    bool is_C0_ok = p1->GetC0AtFront(t_target);
    if (is_C0_ok) {
      rho = rhoM;
      T = TM * pow(rho / rhoM, kappa_Tv - 1.);
      p   = p1->gas->Get_p(rho, T);
      a   = p1->gas->Get_SonicVel(T,p);
      v = (p - beta) / rho / a;
      if (DEBUG){
        printf("\n (front) OUTFLOW beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e",beta, p, T, rho, v);
        if (v>0){
          cout<<endl<<"ERROR!!! Assumed outflow but computed inflow."<<endl;
          //	exit(-1);
        }
      }
    }
    else {
      // Still inflow, set velocity to 0.
      T = TM;
      rho = rhoM;
      p   = p1->gas->Get_p(rho, T);
      a   = p1->gas->Get_SonicVel(T,p);
      v   = (p - beta) / rho / a;
      if (DEBUG){
        printf("\n (front) ERROR !! UNDECIDED beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e",beta, p, T, rho, v);
        //exit(-1);
      }
    }
  }
}

double Connector::Connector_LWP_Reservoir_and_Pipe_Front_fun(double v, double beta, Reservoir* r1, LWP* p1,bool inlet_pressure_drop) {

  double pt = r1->Get_dprop("p");
  double Tt = r1->Get_dprop("T");
  //double rhot = p1->gas->Get_rho(pt, Tt);
  double cp   = p1->gas->Get_cp(pt,Tt);
  double kappa_Tv = p1->gas->Get_kappa_Tv();

  // The absolute value is very important!
  double T   = Tt-v*fabs(v)/2./cp;
  double p=pt;
  if (inlet_pressure_drop)
    p   = pt * pow(T / Tt, kappa_Tv/(kappa_Tv - 1.));
  double rho = p1->gas->Get_rho(p, T);
  double a   = p1->gas->Get_SonicVel(T,p);
  return (p - rho * a * v) - beta;
}

double Connector::signed_sqrt(double x) {
  if (x > 0)
    return sqrt(x);
  else
    return (-sqrt(-x));
}

bool Connector::Connector_SCP_Pipe_Back_and_Valve(double t_target,
    SCP* p1, Valve* v1, double p_downstream, double rho, double a, double& p) {
  // Solve the following system for p,T,rho,v:
  // (1) alpha=p+ro*a*v
  // (2) ro*Apipe*v=valve_mass_flow

  if (DEBUG)
    cout<<endl<<endl<< "Connector::Connector_SCP_Pipe_Back_and_Valve: entering DEBUG mode"<<endl;

  bool update_OK=true;

  double f, f1, df, dp, v, Apipe = p1->Get_dprop("A");
  double x = v1->Get_dprop("x");
  double alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

  double err1, err2, err = 1.e5, mp, TOL = 10., pnew;
  int iter = 0, MAX_ITER = 100;
  double RELAX=0.8;
  while ((fabs(err) > TOL) && (iter<=MAX_ITER)) {

    // 1. sol
    //mp  = v1->Get_MassFlowRate_InCompressible(p, p_downstream, rho, x);
    //v   = mp / rho / Apipe;
    //p   = alpha - rho * a * v;
    //err1 = (alpha - (p + rho * a * v)) / 1.e5;
    //err2 = rho * v * Apipe - (v1->Get_MassFlowRate_InCompressible(p, p_downstream, rho, x));

    // Newton (this is faaar quicker)
    mp  = v1->Get_MassFlowRate(p, 293., p_downstream, 293., x);
    v   = mp / rho / Apipe;
    f = p - (alpha - rho * a * v);

    dp = 0.01 * p;
    if (fabs(dp)<10.)
      dp = 10.;

    mp  = v1->Get_MassFlowRate(p + dp, 293., p_downstream, 293., x);
    v   = mp / rho / Apipe;
    f1  = (p + dp) - (alpha - rho * a * v);
    df  = (f1 - f) / dp;

    pnew = p - f / df;
    p=RELAX*pnew+(1-RELAX)*p;

    err1 = f;
    err2 = 0;
    if (DEBUG)
      printf("\n\t mp=%5.3e, x=%5.3e, v=%5.3e, f=%5.3e, f1=%5.3e, puj=%5.3e",mp,x,v,f,f1,p/1e5);

    err = sqrt(err1 * err1 + err2 * err2);

    //printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
    //				iter, alpha, p, v, err1, err2);

    iter++;
    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_SCP_Pipe_Back_and_Valve() -> MAX_ITER reached, exiting";
      printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
          iter, alpha, p, v, err1, err2);
      update_OK=false;
    }
  }

  return update_OK;
}

bool Connector::Connector_SCP_Pipe_Front_and_Valve_Outlet(double t_target,
    SCP* p1, Valve* v1, double p_upstream, double rho, double a, double& p) {
  // Solve the following system for p,T,rho,v:
  // (1) beta=p-ro*a*v
  // (2) ro*Apipe*v=valve_mass_flow

  if (DEBUG)
    cout<<endl<<endl<< "Connector::Connector_SCP_Pipe_Back_and_Valve: entering DEBUG mode"<<endl;

  bool update_OK=true;

  double f, f1, df, dp, v, Apipe = p1->Get_dprop("A");
  double x = v1->Get_dprop("x");
  double beta = p1->GetBetaPrimitiveAtFront(t_target);

  double err = 1.e5, mp, TOL = 10., pnew;
  int iter = 0, MAX_ITER = 100;
  double RELAX=0.8;
  while ((fabs(err) > TOL) && (iter<=MAX_ITER)) {

    // 1. sol
    //mp  = v1->Get_MassFlowRate_InCompressible(p_upstream, p, rho, x);
    //v   = mp / rho / Apipe;
    //p   = beta + rho * a * v;
    //err1 = (beta - (p - rho * a * v)) / 1.e5;
    //err2 = rho * v * Apipe - (v1->Get_MassFlowRate_InCompressible(p_upstream, p, rho, x));

    // Newton (this is faaar quicker)
    mp  = v1->Get_MassFlowRate(p_upstream, 293., p, 293., x);
    v   = mp / rho / Apipe;
    f = p - (beta + rho * a * v);

    dp = 0.01 * p;
    if (fabs(dp)<10.)
      dp = 10.;

    mp  = v1->Get_MassFlowRate(p_upstream, 293., p+dp, 293., x);
    v   = mp / rho / Apipe;
    f1  = (p + dp) - (beta + rho * a * v);
    df  = (f1 - f) / dp;

    pnew = p - f / df;
    p=RELAX*pnew+(1-RELAX)*p;

    err = fabs(f);
    if (DEBUG)
      printf("\n\t pu=%5.3e, pd=%5.3e, mp=%5.3e, x=%5.3e, v=%5.3e, err=%5.3e",p_upstream,p,mp,x,v,err);

    iter++;
    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_SCP_Pipe_Front_and_Valve() -> MAX_ITER reached, exiting";
      printf("\n (back) iter #%2d: beta=%5.3e, p=%5.3e, v=%5.3e, err=%+5.3e",
          iter, beta, p, v, err);
      update_OK=false;
    }
  }

  return update_OK;
}


/*!
  Connects SCP front to reservoir
  \param t_target [in] target time value
  \param *r1 [in] pointer to reservoir
  \param *p1 [in] pointer to pipe
  \param rho [in] density, kg/m3
  \param a [in] dsonic velocity, m/s
  \param inlet pressure drop [in] bool, if true, inlet pressure drop issumed, i.e. pr=p+rho/2*v^2
  \param p [out] pressure at the first node of the pipe, result of computations
  */

void Connector::Connector_SCP_Reservoir_and_Pipe_Front(double t_target,
    Reservoir* r1, SCP* p1, double rho, double a, bool inlet_pressure_drop, double& p) {
  // Solve the following system for p,T,rho,v:
  // (1) p-rho*a*v = beta_front   Charactersistic equation
  // (2a) p=pt-rho/2*v^2       if inlet_pressure_drop=true
  // (2b) p=pt                 if inlet_pressure_drop=false

  double IPD_mul=0;
  if (inlet_pressure_drop)
    IPD_mul=1.;
  double pt = r1->Get_dprop("p");
  double v = 0;
  double beta = p1->GetBetaPrimitiveAtFront(t_target);

  p = pt;
  double err1, err2, err = 1.e5, TOL = 1e-3;
  int iter = 0, MAX_ITER = 50;
  while (fabs(err) > TOL) {
    v   = (p - beta) / rho / a; // "Corrector"
    p   = pt - IPD_mul*rho / 2 * v * v;
    err1 = pt - (p + IPD_mul*rho * v * v / 2.);
    err2 = (p - rho * a * v - beta) / 1.e5;
    err = sqrt(err1 * err1 + err2 * err2);

    //printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);

    iter++;
    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached.";
      printf("\n (front) iter #%2d: pt=%5.3e, p=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e",
          iter, pt, p, v, err1, err2);
      cin.get();
    }
  }
  //cin.get();
}


void Connector::Connector_SCP_Pipes(double t_target,
    SCP *p1, bool is_front1,
    SCP *p2, bool is_front2,
    double mpout,
    double& p, double& v1, double& v2) {

  // Solves:
  // p + d1*rho*a*v1 = alpha1
  // p + d2*rho*a*v2 = alpha2
  // d1*A1*v1 + d2*A2*v2 = mpout/rho
  //
  // Solution: p=(sum Ai*alphai-mpout*a)/(sum Ai)
  //
  // is_front = false -> d1=+1, alpha=pL+rho*a*vL
  // is_front = true  -> d1=-1, alpha=pR-rho*a*vR

  double alpha1, alpha2;
  double A1  = p1->Get_dprop("A");
  double A2  = p2->Get_dprop("A");
  double rho = p1->Get_dprop("rho");
  double a   = p1->Get_dprop("a");
  int d1 = 1, d2 = 1;
  if (is_front1) {
    alpha1 = p1->GetBetaAtFront(t_target);
    d1 = -1;
  }
  else
    alpha1 = p1->GetAlphaAtEnd(t_target);

  if (is_front2) {
    alpha2 = p2->GetBetaAtFront(t_target);
    d2 = -1;
  }
  else
    alpha2 = p2->GetAlphaAtEnd(t_target);

  p  = (A1 * alpha1 + A2 * alpha2 - mpout * a) / (A1 + A2);
  v1 = (alpha1 - p) / (d1 * rho * a);
  v2 = (alpha2 - p) / (d2 * rho * a);

}

void Connector::Connector_SCP_Pipes(double t_target,
    SCP *p1, bool is_front1,
    SCP *p2, bool is_front2,
    SCP *p3, bool is_front3,
    double mpout,
    double& p, double& v1, double& v2, double& v3) {

  // Solves:
  // p + d1*rho*a*v1 = alpha1
  // p + d2*rho*a*v2 = alpha2
  // p + d3*rho*a*v3 = alpha3
  // d1*A1*v1 + d2*A2*v2 + d3*A3*v3 = mpout/rho
  //
  // Solution: p=(sum Ai*alphai-mpout*a)/(sum Ai)
  //
  // is_front = false -> d1=+1, alpha=pL+rho*a*vL
  // is_front = true  -> d1=-1, alpha=pR-rho*a*vR

  double alpha1, alpha2, alpha3;
  double A1  = p1->Get_dprop("A");
  double A2  = p2->Get_dprop("A");
  double A3  = p3->Get_dprop("A");
  double rho = p1->Get_dprop("rho");
  double a   = p1->Get_dprop("a");
  int d1 = 1, d2 = 1, d3 = 1;
  if (is_front1) {
    alpha1 = p1->GetBetaAtFront(t_target);
    d1 = -1;
  }
  else
    alpha1 = p1->GetAlphaAtEnd(t_target);

  if (is_front2) {
    alpha2 = p2->GetBetaAtFront(t_target);
    d2 = -1;
  }
  else
    alpha2 = p2->GetAlphaAtEnd(t_target);

  if (is_front3) {
    alpha3 = p3->GetBetaAtFront(t_target);
    d3 = -1;
  }
  else
    alpha3 = p3->GetAlphaAtEnd(t_target);

  p  = (A1 * alpha1 + A2 * alpha2 + A3 * alpha3 - mpout * a) / (A1 + A2 + A3);
  v1 = (alpha1 - p) / (d1 * rho * a);
  v2 = (alpha2 - p) / (d2 * rho * a);
  v3 = (alpha3 - p) / (d3 * rho * a);
}

void Connector::Connector_LWP_Pipes(double t_target,
    LWP *p1, bool is_front1,
    LWP *p2, bool is_front2,
    LWP *p3, bool is_front3,
    const double mpout, const double T){

  // Solves:
  // p + d1*rho*a*v1 = alpha1
  // p + d2*rho*a*v2 = alpha2
  // p + d3*rho*a*v3 = alpha3
  // d1*A1*v1 + d2*A2*v2 + d3*A3*v3 = mpout/rho
  //
  // Solution: p=(sum Ai*alphai-mpout*a)/(sum Ai)
  //
  // is_front = false -> d1=+1, alpha=pL+rho*a*vL
  // is_front = true  -> d1=-1, alpha=pR-rho*a*vR

  double alpha1, alpha2, alpha3;
  double A1  = p1->Get_dprop("A");
  double A2  = p2->Get_dprop("A");
  double A3  = p3->Get_dprop("A");
  double p=0., p_new, rho, a, v1, v2, v3;

  int d1 = 1, d2 = 1, d3 = 1;
  if (is_front1) {
    alpha1 = p1->GetBetaPrimitiveAtFront(t_target);
    d1 = -1;
  }
  else
    alpha1 = p1->GetAlphaPrimitiveAtEnd(t_target);

  if (is_front2) {
    alpha2 = p2->GetBetaPrimitiveAtFront(t_target);
    d2 = -1;
  }
  else
    alpha2 = p2->GetAlphaPrimitiveAtEnd(t_target);

  if (is_front3) {
    alpha3 = p3->GetBetaPrimitiveAtFront(t_target);
    d3 = -1;
  }
  else
    alpha3 = p3->GetAlphaPrimitiveAtEnd(t_target);

  try{
    double err_p=1.e5, ERR_P_MAX=10.;
    int iter=0, ITER_MAX=20;
    while (err_p>ERR_P_MAX){
      p_new  = (A1 * alpha1 + A2 * alpha2 + A3 * alpha3 - mpout * a) / (A1 + A2 + A3);	
      a = p1->gas->Get_SonicVel(T,p);
      err_p=sqrt((p-p_new)*(p-p_new));
      p=p_new;
      if (DEBUG)
        printf("\n\t iter #%d/%d, p=%5.3e, a=%5.3e, err_p=%5.3e",iter,ITER_MAX,p_new,a,err_p);

      iter++;
      if (iter==ITER_MAX)
        throw(iter);
    }
  }
  catch (int iter){
    cout<<endl<<"!!! ERROR !!!";
    cout<<endl<<"Connector::LWP_Pipes() -> too many iterations!!!";
    cout<<endl<<"Exiting..."<<endl<<endl;
    exit(-1);
  }

  p=p_new;
  rho = p1->gas->Get_rho(p,T);
  v1 = (alpha1 - p) / (d1 * rho * a);
  v2 = (alpha2 - p) / (d2 * rho * a);
  v3 = (alpha3 - p) / (d3 * rho * a);

  Set_LWP_BC(p1,is_front1,p,T,v1);
  Set_LWP_BC(p2,is_front2,p,T,v2);
  Set_LWP_BC(p3,is_front3,p,T,v3);
}

void Connector::Set_LWP_BC(LWP* p1, const bool is_front, const double p, const double T, const double v){

  double is_inflow,is_outflow;
  if (is_front){
    if (v>0)
      is_inflow=p1->BCLeft("StaticPres_and_StaticTemp_Inlet",p,T,true);
    else
      is_outflow=p1->BCLeft("StaticPres_Outlet",p,0,true);
  }
  else{
    if (v<0)
      is_inflow=p1->BCRight("StaticPres_and_StaticTemp_Inlet",p,T,true);
    else
      is_outflow=p1->BCRight("StaticPres_Outlet",p,0,true);
  }
}

void Connector::Connector_LWP_Pipes(double t_target,
    LWP *p1, LWP *p2){
  /*
     void Connector::Connector_LWP_Pipes(double t_target,
     LWP *p1, LWP *p2,
     double& pL, double& pR,
     double& TL, double& TR,
     double& rhoL, double& rhoR,
     double& vL, double& vR) {
     */
  //bool DEBUG=true;

  // Solves:
  // eq. (1) pL + rhoL*aL*vL = alpha_L
  // eq. (2) pR - rhoR*aR*vR = beta_R
  // eq. (3) rhoL*vL*AL   = rhoR*vR*AR (continuity)
  // eq. (4) rhoL*vL^2*AL = rhoR*vR^2*AR+(pR-pL)*max(AL,AR) (impulse)
  // eq. (5) (rhoL*vL*eL + pL*vL)*AL = (rhoR*vR*eR + pR*vR)*AR (energy)
  // eq. (6-7) a = sonic_vel(T,p)
  // eq. (8-9) e = internal_energy(T,p)
  // eq. (10-11) gas law
  // eq. (12) rhoL = rho0 (C0 characteristic)
  //
  // 12 unkowns: p, T, rho, v, a, e @ R(ight) and L(eft)

  bool is_inflow, is_outflow;

  double AR = p1->Get_dprop("A");
  double AL = p2->Get_dprop("A");
  double A = AR;
  if (AL > AR)
    A = AL;

  double pL   = p1->Get_dprop("p_back");
  double TL   = 293.;//p1->Get_dprop("T_back");
  double vL   = p1->Get_dprop("v_back");
  double	rhoL = p1->gas->Get_rho(pL, TL);
  double aL   = p1->gas->Get_SonicVel(TL,pL);
  double eL   = p1->gas->Get_e_from_Tp(TL,pL) + vL * vL / 2;

  double pR   = p2->Get_dprop("p_front");
  double TR   = 293.;//p2->Get_dprop("T_front");
  double vR   = p2->Get_dprop("v_front");
  double rhoR = p2->gas->Get_rho(pR, TR);
  double aR   = p2->gas->Get_SonicVel(TR,pR);
  double eR   = p1->gas->Get_e_from_Tp(TR,pR) + vR * vR / 2;

  double pR_new,pL_new;//,rhoR_new,rhoL_new,vL_new,vR_new;

  double alphaL = p1->GetAlphaPrimitiveAtEnd(t_target);
  double betaR  = p2->GetBetaPrimitiveAtFront(t_target);
  if (DEBUG){
    printf("\n\n Connector::Connector_LWP_Pipes()");
    printf("\n pL=%5.3e, TL=%5.3e, vL=%5.3e, rhoL=%5.3e, aL=%5.3e, eL=%5.3e, alphaL=%5.3e",
        pL,TL,vL,rhoL,aL,eL,alphaL);
    printf("\n pR=%5.3e, TR=%5.3e, vR=%5.3e, rhoR=%5.3e, aR=%5.3e, eR=%5.3e, betaR =%5.3e",
        pR,TR,vR,rhoR,aR,eR,betaR);
    cin.get();
  }

  // Preliminary esrimate assuming constant temperature
  double aa=1./rhoR/AR-1./rhoL/AL;
  double bb=aR+aL;
  double cc=betaR*AR-alphaL*AL;
  double DD=bb*bb-4*aa*cc;
  double mp=0.0;
  if (DD>0){
    mp=(-bb+sqrt(DD))/2./aa;
  }
  else{
    cout<<endl<<endl<<"ERRORRRRRR!"<<endl;
    cin.get();
  }

  pL = alphaL-aL/AL*mp;
  pR = betaR +aR/AR*mp;
  rhoL = p1->gas->Get_rho(pL,TL);
  rhoR = p1->GetC0AtEnd(t_target);
  vL=mp/rhoL/AL;
  vR=mp/rhoR/AR;

  if (DEBUG){
    printf("\n\nEstimate:");
    printf("\n mp=%5.3e",mp);
    printf("\n pL=%5.3e, TL=%5.3e, vL=%5.3e, rhoL=%5.3e",pL,TL,vL,rhoL);
    printf("\n pR=%5.3e, TR=%5.3e, vR=%5.3e, rhoR=%5.3e",pR,TR,vR,rhoR);
    cin.get();
  }

  double err_p = 1e5, TOL_p = 10.;
  double err_T = 1e5, TOL_T = .1;
  int iter = 0, MAX_ITER = 1000;

  while ((fabs(err_p) > TOL_p) || (fabs(err_T) > TOL_T)) {
    iter++;
    if (iter == MAX_ITER) {
      cout << endl << "!!!ERROR!!! Connector::Connector_LWP_Pipes() -> MAX_ITER reached. Exiting" << endl;
      //printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);
      exit(-1);
    }
    // eq. (1) & (2)
    double tmpL=rhoL*vL*vL*AL;
    double tmpR=rhoR*vR*vR*AR;
    double Atmp = max(AL,AR);
    //if (iter<3)
    //pL_new = alphaL - rhoL * aL * vL;
    pL_new = tmpR/Atmp + pR - tmpL/Atmp; 

    pR_new = betaR + rhoR * aR * vR;


    double v_TRESHOLD=0.01;
    vL = (alphaL-pL_new)/rhoL/aL; 
    // eq. (12)
    if (fabs(vL)>v_TRESHOLD){
      rhoL = p1->GetC0AtEnd(t_target);
      mp=rhoL*vL*AL;
      if (rhoL<0){
        cout<<endl<<"!!! ERROR !!!! rhoL="<<rhoL<<" < 0 !!!"<<endl;
        cin.get();
      }
      vR = mp / AR / rhoR;
      rhoR = (pL+rhoL*eL-pR)/eR; 
    }
    else{
      pL_new = alphaL - rhoL * aL * vL;
      vL=0.;
      vR=0.;
      rhoR = p2->GetC0AtFront(t_target);
      rhoL = rhoR; 
    }

    // eq. (6-7)
    aL = p1->gas->Get_SonicVel(TL,pL_new);
    aR = p2->gas->Get_SonicVel(TR,pR_new);
    // eq. (8-9)
    eL = p1->gas->Get_e_from_Tp(TL,pL_new) + vL * vL / 2.;
    eR = p2->gas->Get_e_from_Tp(TR,pR_new) + vR * vR / 2.;
    // eq. (10-11)
    double TL_new = p1->gas->Get_T(rhoL, pL_new);
    double TR_new = p2->gas->Get_T(rhoR, pR_new);

    err_T = sqrt(pow(TL - TL_new, 2.) + pow(TR - TR_new, 2.));
    err_p = sqrt(pow(pL - pL_new, 2.) + pow(pR - pR_new, 2.));
    pL = (pL+pL_new)/2.;
    pR = (pR+pR_new)/2.;
    TL = (TL+TL_new)/2.;
    TR = (TR+TR_new)/2.;
    mp=vL*rhoL*AL;
    double mpR=vR*rhoR*AR;
    if (DEBUG){
      printf("\n\t iter=%2d, (L,R) p=%5.3f,%5.3f, T=%5.3f,%5.3f, rho=%5.3f,%5.3f, v=%5.3f,%5.3f, mp=%5.3f, %5.3f, err=%5.3e,%5.3e",
          iter,pL_new/1.e5,pR_new/1.e5,TL,TR,rhoL,rhoR,vL,vR,mp,mpR,err_T,err_p);
    }
  }
  if (vR>0){
    is_outflow=p1->BCRight("StaticPres_Outlet",pL,0,true);
    is_inflow=p2->BCLeft("StaticPres_and_StaticTemp_Inlet",pR,TR,true);
  }
  else{
    is_outflow=p2->BCRight("StaticPres_Outlet",pR,0,true);
    is_inflow=p1->BCLeft("StaticPres_and_StaticTemp_Inlet",pL,TL,true);
  }
}
