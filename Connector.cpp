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

void Connector_Reservoir_and_Valve(double t_target, Reservoir* r1, Valve *v1, double& p, double& mp){

	// Solve the following system for p,T,rho,v,a:
	// (1) mp = MassFlowRate(pt,Tt,pout,Tout,x)
	// (2) rho=rho(p,T)gt
	// (4) a  =a(p,T)
	// (5) T/Tt=(pt/p)^(kappa-1)/kappa

	//double kappa=p1->gas->kappa;
	//double cp   =p1->gas->cp;

	//double update_OK=true;
	//double rho, v, a, alpha;
	//Apipe = p1->Get_dprop("A");
	//p1->GetAllPrimitiveAtEnd(t_target, pL, v, T, rho);
	//alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

	//double err1, err2, err = 1.e5, mp, TOL = 1e-3;
	//int iter = 0, MAX_ITER = 10;
	//p = pL;
	//while ((err > TOL) && (iter<=MAX_ITER)) {
	//	a   = p1->gas->Get_SonicVel(T,p);
	//	rho = p1->gas->Get_rho(p, T);
	//	mp  = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
	//	v   = mp / rho / Apipe;
	//	p   = alpha - rho * a * v;
	//T   = gas->Get_T(p,rho);
	//	err1 = (alpha - (p + rho * a * v)) / 1.e5;
	//	err2 = rho * v * Apipe - mp;
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

	p1->BCRight("StaticPres_and_StaticTemp",p,T,true);

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
	double rho, v, a, alpha, Apipe, pL;

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
		LWP *p1, LWP *p2,
		double& pL, double& pR,
		double& TL, double& TR,
		double& rhoL, double& rhoR,
		double& vL, double& vR) {

	bool DEBUG=true;

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

	double AR = p1->Get_dprop("A");
	double AL = p2->Get_dprop("A");
	double A = AR;
	if (AL > AR)
		A = AL;

	pL   = p1->Get_dprop("p_back");
	TL   = p1->Get_dprop("T_back");
	vL   = p1->Get_dprop("v_back");
	rhoL = p1->gas->Get_rho(pL, TL);
	double aL   = p1->gas->Get_SonicVel(TL,pL);
	double eL   = p1->gas->Get_e_from_Tp(TL,pL) + vL * vL / 2;

	pR   = p2->Get_dprop("p_front");
	TR   = p2->Get_dprop("T_front");
	vR   = p2->Get_dprop("v_front");
	rhoR = p2->gas->Get_rho(pR, TR);
	double aR   = p2->gas->Get_SonicVel(TR,pR);
	double eR   = p1->gas->Get_e_from_Tp(TR,pR) + vR * vR / 2;

	double alphaL = p1->GetAlphaPrimitiveAtEnd(t_target);
	double betaR  = p2->GetBetaPrimitiveAtFront(t_target);
	if (DEBUG){
		printf("\n\n Connector::Connector_LWP_Pipes()");
		printf("\n pL=%5.3e, TL=%5.3e, vL=%5.3e, rhoL=%5.3e, aL=%5.3e, eL=%5.3e, alphaL=%5.3e",
				pL,TL,vL,rhoL,aL,eL,alphaL);
		printf("\n pR=%5.3e, TR=%5.3e, vR=%5.3e, rhoR=%5.3e, aR=%5.3e, eR=%5.3e, betaR=%5.3e",
				pR,TR,vR,rhoR,aR,eR,betaR);
		cin.get();
	}
	double err_p = 1e5, TOL_p = 10.;
	double err_T = 1e5, TOL_T = .1;
	int iter = 0, MAX_ITER = 20;
	while ((fabs(err_p) > TOL_p) || (fabs(err_T) > TOL_T)) {
		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector::Connector_LWP_Pipes() -> MAX_ITER reached. Exiting" << endl;
			//printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);
			exit(-1);
		}
		// eq. (1) & (2)
		double pL_new = alphaL - rhoL * aL * vL;
		double pR_new = betaR + rhoR * aR * vR;
		err_p = sqrt(pow(pL - pL_new, 2.) + pow(pR - pR_new, 2.));
		pL = pL_new;
		pR = pR_new;

		double v_TRESHOLD=0.1;
		// eq. (12)
		if (vL>v_TRESHOLD){
			rhoL = p1->GetC0AtEnd(t_target);			
			// eq. (3)
			vL = rhoR * AR * vR / AL / rhoL;
			// eq. (4)
			double tmp = rhoL * vL * vL * AL - (pR - pL) * A;
			tmp = fabs(tmp);
			if (tmp < 0) {
				cout << endl << "!!!ERROR!!! Connector::Connector_LWP_Pipes() -> negative value at sqrt() of vR" << endl;
				//printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);
				exit(-1);
			}
			else
				vR = sqrt(tmp / rhoR / AR);
			// eq. (5)
			if (fabs(vR>0.001))
				rhoR = ((rhoL * vL * eL + pL * vL) * AL / AR - pR*vR) / eR/vR;
			else
				rhoR=rhoL;
		}
		else{
			rhoR = p2->GetC0AtFront(t_target);			
			// eq. (3)
			vR = rhoL * AL * vL / AR / rhoR;
			// eq. (4)
			// rhoL*vL^2*AL = rhoR*vR^2*AR+(pR-pL)*max(AL,AR) (impulse)
			double tmp = rhoR * vR * vR * AR + (pR - pL) * A;
			tmp = fabs(tmp);
			if (tmp < 0) {
				cout << endl << "!!!ERROR!!! Connector::Connector_LWP_Pipes() -> negative value at sqrt() of vL" << endl;
				//printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);
				exit(-1);
			}
			else
				vL = sqrt(tmp / rhoL / AL);
			// eq. (5)
			// (rhoL*vL*eL + pL*vL)*AL = (rhoR*vR*eR + pR*vR)*AR (energy)
			if (fabs(vL>0.001))
				rhoL = ((rhoR * vR * eR + pR * vR) * AR / AL - pL*vL) / eL/vL;
			else
				rhoL=rhoR;
		}

		// eq. (6-7)
		aL = p1->gas->Get_SonicVel(TL,pL);
		aR = p2->gas->Get_SonicVel(TR,pR);
		// eq. (8-9)
		eL = p1->gas->Get_e_from_Tp(TL,pL) + vL * vL / 2;
		eR = p2->gas->Get_e_from_Tp(TR,pR) + vR * vR / 2;
		// eq. (10-11)
		double 	TL_new = p1->gas->Get_T(rhoL, pL);
		double TR_new = p2->gas->Get_T(rhoR, pR);
		err_T = sqrt(pow(TL - TL_new, 2.) + pow(TR - TR_new, 2.));
		TL = TL_new;
		TR = TR_new;
		if (DEBUG){
			printf("\n\t iter=%2d, pL_new=%5.3e, pR_new=%5.3e, rhoL=%5.3e, vL=%5.3e,vR=%5.3e,rhoR=%5.3e",
					iter,pL_new,pR_new,rhoL,vL,vR,rhoR);
		}
	}
}
