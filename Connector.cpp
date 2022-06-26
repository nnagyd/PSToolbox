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
Connector::Connector() {}

// ! Default Destructor
Connector::~Connector() {}

void Connector_Reservoir_and_Valve(double t_target, Reservoir* r, Valve *v, double& p, double& mp){

	// Solve the following system for p,T,rho,v,a:
	// (1) mp = MassFlowRate(pt,Tt,pout,Tout,x)
	// (2) rho=rho(p,T)gt
	// (4) a  =a(p,T)
	// (5) T/Tt=(pt/p)^(kappa-1)/kappa

	//double kappa=p1->gas->kappa;
	//double cp   =p1->gas->cp;

	double update_OK=true;
	double rho, v, a, alpha;
	Apipe = p1->Get_dprop("A");
	p1->GetAllPrimitiveAtEnd(t_target, pL, v, T, rho);
	alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

	double err1, err2, err = 1.e5, mp, TOL = 1e-3;
	int iter = 0, MAX_ITER = 10;
	p = pL;
	while ((err > TOL) && (iter<=MAX_ITER)) {
		a   = p1->gas->Get_SonicVel(T,p);
		rho = p1->gas->Get_rho(p, T);
		mp  = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
		v   = mp / rho / Apipe;
		p   = alpha - rho * a * v;
		//T   = gas->Get_T(p,rho);
		err1 = (alpha - (p + rho * a * v)) / 1.e5;
		err2 = rho * v * Apipe - mp;
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
		LWP* p1, Valve* v1, double p_downstream, double& p, double& T) {
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
	p1->GetAllPrimitiveAtEnd(t_target, pL, v, T, rho);
	alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

	double err1, err2, err = 1.e5, mp, TOL = 1e-3;
	int iter = 0, MAX_ITER = 10;
	p = pL;
	while ((err > TOL) && (iter<=MAX_ITER)) {
		a   = p1->gas->Get_SonicVel(T,p);
		rho = p1->gas->Get_rho(p, T);
		mp  = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
		v   = mp / rho / Apipe;
		p   = alpha - rho * a * v;
		//T   = gas->Get_T(p,rho);
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
		Reservoir* r1, LWP* p1, double& p, double& T) {
	// Solve the following system for p,T,rho,v:
	// (1)  Tt=T+v^2/2/cp            Isentropic flow from the reservoir to the pipe
	// (2)  p-rho*a*v = beta_front   Charactersistic equation
	// (3)  rho = rho(p,T)
	// (4a) p=pt;
	// (4b) pt/rhot^k = p/rho^k

	double cp   = p1->gas->Get_cp();
	double kappa_Tv = p1->gas->Get_kappa_Tv();
	double pt = r1->Get_dprop("p");
	double Tt = r1->Get_dprop("T");
	double rhot = p1->gas->Get_rho(pt, Tt);

	double rho, v, a, beta;
	p1->GetAllPrimitiveAtFront(t_target, p, v, T, rho);
	beta = p1->GetBetaPrimitiveAtFront(t_target);

	bool is_ok = false;

	// Assume inflow
	rho = rhot;

	double err = 1.e5, TOL = 1e-3, drho, f, f1, df, rhonew;
	int iter = 0, MAX_ITER = 50;
	while (err > TOL) {
		f   = Connector_LWP_Reservoir_and_Pipe_Front_fun(rho, beta, r1, p1);

		if (fabs(rho) < 0.01)
			drho = 0.01;
		else
			drho = rho * 0.01;
		f1 = Connector_LWP_Reservoir_and_Pipe_Front_fun(rho + drho, beta, r1, p1);

		df = (f1 - f) / drho;
		rhonew = rho - f / df;

		err = f;

		/*T = Tt * pow(rho / rhot,kappa_Tv - 1.);
		  a   = p1->gas->Get_SonicVel(T);
		  p   = p1->gas->Get_p(rho,T);
		  v = signed_sqrt(2 * cp * (Tt - T));
		  printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e",
		  iter, pt, rhot, beta, p, T, rho, v, err, err1, err2, err3);
		  cin.get();*/

		iter++; rho = rhonew;
		if (iter == MAX_ITER) {
			T = Tt * pow(rho / rhot, kappa_Tv - 1.);
			p   = p1->gas->Get_p(rho, T);
			a   = p1->gas->Get_SonicVel(T,p);
			v = signed_sqrt(2 * cp * (Tt - T));
			cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached.";
			/*printf("\n (front) INFLOW iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err=%5.3e",
			  iter, pt, rhot, beta, p, T, rho, v, err);
			  cin.get();*/
		}
	}
	T = Tt * pow(rho / rhot, kappa_Tv - 1.);
	p   = p1->gas->Get_p(rho, T);
	a   = p1->gas->Get_SonicVel(T,p);
	v = signed_sqrt(2 * cp * (Tt - T));

	if (v < 0) {
		is_ok=false;
		/*cout << endl << " Connector_Reservoir_and_Pipe_Front() -> Assumed inflow, but computed v<0";
		  printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err=%5.3e",
		  iter, pt, rhot, beta, p, T, rho, v, err);*/
		//cin.get();
	}
	else
		is_ok = true;

	// Need to try again with outflow
	if (!is_ok) {

		double pM, vM, TM, rhoM;
		bool is_C0_ok = p1->GetC0AtFront(t_target, pM, vM, TM, rhoM);
		if (is_C0_ok) {
			rho = rhoM;
			T = TM * pow(rho / rhoM, kappa_Tv - 1.);
			p   = p1->gas->Get_p(rho, T);
			a   = p1->gas->Get_SonicVel(T,p);
			v = (p - beta) / rho / a;
			/*printf("\n (front) OUTFLOW pM=%5.3e, rhoM=%5.3e, TM=%5.3e, vM=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e",
			  pM, rhoM, TM, vM, beta, p, T, rho, v);
			  cin.get();*/
		}
		else {
			// Still inflow, set velocity to 0.
			T = TM;
			rho = rhoM;
			p   = p1->gas->Get_p(rho, T);
			a   = p1->gas->Get_SonicVel(T,p);
			v   = (p - beta) / rho / a;
			//printf("\n (front) (!) UNDECIDED beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e",beta, p, T, rho, v);
			/*cin.get();*/
		}
	}
}

double Connector::Connector_LWP_Reservoir_and_Pipe_Front_fun(double rho, double beta, Reservoir* r1, LWP* p1) {

	double cp   = p1->gas->Get_cp();
	double kappa_Tv = p1->gas->Get_kappa_Tv();
	double pt = r1->Get_dprop("p");
	double Tt = r1->Get_dprop("T");
	double rhot = p1->gas->Get_rho(pt, Tt);
	//double beta = p1->GetBetaPrimitiveAtFront(t_target);

	double T = Tt * pow(rho / rhot, kappa_Tv - 1.);
	double p   = p1->gas->Get_p(rho, T);
	double a   = p1->gas->Get_SonicVel(T,p);
	double v = signed_sqrt(2 * cp * (Tt - T));
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

	bool update_OK=true;

	double f, f1, df, dp, v, Apipe = p1->Get_dprop("A"), x = v1->Get_dprop("x");
	double alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

	double err1, err2, err = 1.e5, mp, TOL = 1e-3;
	int iter = 0, MAX_ITER = 100;
	while ((err > TOL) && (iter<=MAX_ITER)) {

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
		mp  = v1->Get_MassFlowRate(p + dp, 293., p_downstream, 293., x);
		v   = mp / rho / Apipe;
		f1  = (p + dp) - (alpha - rho * a * v);
		df  = (f1 - f) / dp;

		p = p - f / df;

		err1 = f;
		err2 = 0;
		//printf("\n\t\t puj=%5.3e",p/1e5);
		//cin.get();

		err = sqrt(err1 * err1 + err2 * err2);

		//printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
		//				iter, alpha, p, v, err1, err2);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_SCP_Pipe_Back_and_Valve() -> MAX_ITER reached, exiting";
			printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
					iter, alpha, p, v, err1, err2);
			//cin.get();
			update_OK=false;
		}
	}
	return update_OK;
}


void Connector::Connector_SCP_Reservoir_and_Pipe_Front(double t_target,
		Reservoir* r1, SCP* p1, double rho, double a, double& p) {
	// Solve the following system for p,T,rho,v:
	// (1) p-rho*a*v = beta_front   Charactersistic equation
	// (1) p=pt-rho/2*v^2;

	double pt = r1->Get_dprop("p");
	double v = 0;
	double beta = p1->GetBetaPrimitiveAtFront(t_target);

	p = pt;
	double err1, err2, err = 1.e5, TOL = 1e-3;
	int iter = 0, MAX_ITER = 50;
	while (err > TOL) {
		v   = (p - beta) / rho / a; // "Corrector"
		p   = pt - rho / 2 * v * v;
		err1 = pt - (p + rho * v * v / 2.);
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

	// Solves:
	// eq. (1) pL + rhoL*aL*vL = alpha_L
	// eq. (2) pR - rhoR*aR*vR = beta_R
	// eq. (3) rhoL*vL*AL   = rhoR*vR*AR (continuity)
	// eq. (4) rhoL*vL^2*AL = rhoR*vR^2*AR+(pR-pL)*max(AL,AR) (impulse)
	// eq. (5) (rhoL*vL*eL + pL*vL)*AL = (rhoR*vR*eR + pR*vR)*AR (energy)
	// eq. (6-7) a = sonic_vel(T)
	// eq. (8-9) e = internal_energy(T)
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
	double eL   = p1->gas->Get_e(TL) + vL * vL / 2;

	pR   = p2->Get_dprop("p_front");
	TR   = p2->Get_dprop("T_front");
	vR   = p2->Get_dprop("v_front");
	rhoR = p2->gas->Get_rho(pR, TR);
	double aR   = p2->gas->Get_SonicVel(TR,pR);
	double eR   = p1->gas->Get_e(TR) + vR * vR / 2;

	double alphaL = p1->GetAlphaPrimitiveAtEnd(t_target);
	double betaR  = p2->GetBetaPrimitiveAtFront(t_target);

	double err_p = 1e5, TOL_p = 10.;
	double err_T = 1e5, TOL_T = .1;
	int iter = 0, MAX_ITER = 20;
	while ((err_p > TOL_p) || (err_T > TOL_T)) {
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

		// eq. (12)
		rhoL = p1->GetC0AtEnd(t_target);
		// eq. (3)
		vL = rhoR * AR * vR / AL / rhoL;
		// eq. (4)
		double tmp = rhoL * vL * vL * AL - (pR - pL) * A;
		if (tmp < 0) {
			cout << endl << "!!!ERROR!!! Connector::Connector_LWP_Pipes() -> nexative value at sqrt() of vR" << endl;
			//printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);
			exit(-1);
		}
		else
			vR = sqrt(tmp / rhoR / AR);
		// eq. (5)
		double vL_p_vR = rhoL * AL / rhoR / AR;
		rhoR = ((rhoL * vL_p_vR * eL + pL * vL_p_vR) * AL / AR - pR) / eR;
		// eq. (6-7)
		aL = p1->gas->Get_SonicVel(TL,pL);
		aR = p2->gas->Get_SonicVel(TR,pR);
		// eq. (8-9)
		eL = p1->gas->Get_e(TL) + vL * vL / 2;
		eR = p2->gas->Get_e(TR) + vR * vR / 2;
		// eq. (10-11)
		double 	TL_new = p1->gas->Get_T(rhoL, pL);
		double TR_new = p2->gas->Get_T(rhoR, pR);
		err_T = sqrt(pow(TL - TL_new, 2.) + pow(TR - TR_new, 2.));
		TL = TL_new;
		TR = TR_new;

	}
}
