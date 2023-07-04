#define _USE_MATH_DEFINES
#include "LWP.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

LWP::LWP(const string _name,
		const string _cspe_name,
		const string _cspv_name,
		//const double _ro,
		//const double _a,
		const double _L,
		const double _D,
		const double _lambda,
		const double _he,
		const double _hv,
		const bool _save_data,
		const bool _save_all_data,
		Gas *_gas) : Units() {

	save_data = _save_data;
	save_all_data = _save_all_data;
	name = _name;
	node_from = _cspe_name;
	node_to = _cspv_name;
	L = _L;
	D = _D;
	he = _he;
	hv = _hv;
	S0 = -(hv - he) / L;
	A = D * D * M_PI / 4.;
	lambda = _lambda;
	lambda_p_2D = lambda / 2 / D;
	ini_done = false;
	fname = name + ".dat";
	//do_plot_runtime = false;
	//ylim_isset = false;
	g = 9.81;
	gas = _gas;

	P_MIN = 0.01e5;
	T_MIN = 273.15 - 100.;
	art_visc = 0.6;

	// Dummy initialization
	Npts = 1;
	gamma = 0.;
	alpha = 0.;
	phi = 0.;
	mu = 0.;
	dx = 1.;
	t = 0.;
	dt = 0.;
	Tw=293;
	HEAT_TRANSFER_ON=false;
	HEAT_TRANSFER_MODEL="adiabatic";

	left_boundary_points_already_updated=false;
	right_boundary_points_already_updated=false;
	is_BCRight_Wall=false;

	if (save_all_data) {
		string pfname;
		FILE * pFile;
		pfname = name + "_p.dat";
		pFile = fopen(pfname.c_str(), "w");
		fclose(pFile);

		pfname = name + "_v.dat";
		pFile = fopen(pfname.c_str(), "w");
		fclose(pFile);

		pfname = name + "_T.dat";
		pFile = fopen(pfname.c_str(), "w");
		fclose(pFile);

		pfname = name + "_rho.dat";
		pFile = fopen(pfname.c_str(), "w");
		fclose(pFile);

		pfname = name + "_mdot.dat";
		pFile = fopen(pfname.c_str(), "w");
		fclose(pFile);
	}
}

LWP::~LWP() {}

string LWP::Info(bool show_pts) {

	if (!ini_done)
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
	oss << "\n             dt : " << dt << " s";
	//oss << "\n          f(2L) : " << a / (2. * L) << " Hz";
	//oss << "\n         f(L/2) : " << a / (0.5 * L) << " Hz";
	oss << endl;
	oss << "\n        mean(p) : " << p.mean() / 1.e5 << " bar";
	oss << "\n        mean(v) : " << v.mean() << " m/s";
	oss << "\n        mean(T) : " << T.mean() - 273.15 << " C";
	oss << "\n      mean(rho) : " << rho.mean() << " kg/m3";
	oss << "\n      actual dt : " << dt << " s";
	oss << endl;
	/*oss << "\n            phi : " << phi << " -";
	  oss << "\n          alpha : " << alpha << " -";
	  oss << "\n          gamma : " << gamma << " -";
	  oss << "\n             mu : " << mu << " -";*/
	//oss << "\n f_pipe/f_valve : " << 2.*L/a/(omega/(2.*M_PI)) << " -";
	//oss << "\n f_pipe/f_valve : " << M_PI / gamma << " -";

	if (show_pts) {
		oss << "\n\n t = " << t << " s, dt = " << dt << " s";
		oss << "\n # of grid pts.: " << Npts << std::fixed;

		oss << "\n       node #     : ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setfill(' ') << i << " ";

		oss << "\n       p (bar)    : ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(3) << p(i) / 1.e5 << " ";

		oss << "\n       v (m/s)    : ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(1) << v(i) << " ";

		oss << "\n       T (C)      : ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(1) << T(i) - 273.15 << " ";

		oss << "\n       rho (kg/m3): ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(3) << rho(i) << " ";

		oss << "\n            Ma (-): ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(3) << fabs(v(i)) / gas->Get_SonicVel(T(i),p(i)) << " ";

		oss << "\n        mp (kg/s) : ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(3) << rho(i)*v(i)*A << " ";

		/*	double mp_mean = rho.mean() * v.mean() * A;
			oss << "\n       mp_imp (%) : ";
			for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(1) << fabs((rho(i)*v(i)*A - mp_mean) / mp_mean) * 100 << " ";*/
	}
	oss << endl;
	return oss.str();
}

/*! \brief Initialize the LWP pipe
  Initializes the LWP pipe with as many points to statisfy a given timestep, an initial uniform (due to incompressibility) speed and a
  pressure at the inlet. (The generated pressure drops with friction.)
  \param vini Initial (uniform) velocity in the pipe
  \param pstart Pressure at the inlet, set up to decrease with friction
  \param dt_target Target timestep, defines the grid
  */
void LWP::IniUniform(double _vini, double _pini, double _Tini, double dt_target) {
	t = 0.;

	double a = gas->Get_SonicVel(_Tini,_pini);
	Npts = round(L / a / dt_target); // CFL condition reorganized
					 //printf("\n\n L=%5.2f m, a=%5.1f m/s, dt_target=%5.3e s, Npts=%d ", L, a, dt_target, Npts);
	if (Npts < 20) {
		Npts = 20; // a minimum of 20 points is set
			   //printf(" -> %d\n", Npts);
	}
	Ini(Npts);

	for (int i = 0; i < Npts; i++) {
		x(i) = i * L / (Npts - 1);
		v(i) = _vini;
		p(i) = _pini;
		T(i) = _Tini;
		rho(i) = gas->Get_rho(p(i), T(i));
	}

	UpdateTimeStep();

	tmpvec.push_back(0);
	tmpvec.push_back(p(0));
	tmpvec.push_back(p(Npts - 1));
	tmpvec.push_back(v(0));
	tmpvec.push_back(v(Npts - 1));
	tmpvec.push_back(T(0));
	tmpvec.push_back(T(Npts - 1));
	tmpvec.push_back(rho(0));
	tmpvec.push_back(rho(Npts - 1));
	tmpvec.push_back(v(0)*A * rho(0));
	tmpvec.push_back(v(Npts - 1)*A * rho(Npts - 1));
	data.clear();
	data.reserve(100);
	data.push_back(tmpvec);
}

// Direct, user-defined initialization of the primitive variables
// with uniform distribution
void LWP::IniUniform(double _vini, double _pini, double _Tini, int _Npts_mul) {
	Npts = _Npts_mul*20;
	Ini(Npts);
	for (int i = 0; i < Npts; i++) {
		x(i) = i * L / (Npts - 1);
		v(i) = _vini;
		p(i) = _pini;
		T(i) = _Tini;
		rho(i) = gas->Get_rho(p(i), T(i));
	}
	UpdateTimeStep();

	tmpvec.push_back(0);
	tmpvec.push_back(p(0));
	tmpvec.push_back(p(Npts - 1));
	tmpvec.push_back(v(0));
	tmpvec.push_back(v(Npts - 1));
	tmpvec.push_back(T(0));
	tmpvec.push_back(T(Npts - 1));
	tmpvec.push_back(rho(0));
	tmpvec.push_back(rho(Npts - 1));
	tmpvec.push_back(v(0)*A * rho(0));
	tmpvec.push_back(v(Npts - 1)*A * rho(Npts - 1));
	data.clear();
	data.reserve(100);
	data.push_back(tmpvec);

}


// Direct, user-defined initialization of the primitive variables
// with point-wise quantities
void LWP::IniDistribution(const vector<double> _vini,
		const vector<double> _pini,
		const vector<double> _Tini) {
	Npts = _vini.size();

	if (_pini.size() != Npts) {
		cout << endl << endl << "ERROR!!! LWP::IniDistribution";
		cout << endl << "\t vini.size() = " << _vini.size() << " != ";
		cout << endl << "\t pini.size() = " << _pini.size();
	}

	if (_Tini.size() != Npts) {
		cout << endl << endl << "ERROR!!! LWP::IniDistribution";
		cout << endl << "\t vini.size() = " << _vini.size() << " != ";
		cout << endl << "\t Tini.size() = " << _Tini.size();
	}

	Ini(Npts);
	for (int i = 0; i < Npts; i++) {
		x(i) = i * L / (Npts - 1);
		v(i) = _vini.at(i);
		p(i) = _pini.at(i);
		T(i) = _Tini.at(i);
		rho(i) = gas->Get_rho(p(i), T(i));
	}

	UpdateTimeStep();

	tmpvec.push_back(0);
	tmpvec.push_back(p(0));
	tmpvec.push_back(p(Npts - 1));
	tmpvec.push_back(v(0));
	tmpvec.push_back(v(Npts - 1));
	tmpvec.push_back(T(0));
	tmpvec.push_back(T(Npts - 1));
	tmpvec.push_back(rho(0));
	tmpvec.push_back(rho(Npts - 1));
	tmpvec.push_back(v(0)*A * rho(0));
	tmpvec.push_back(v(Npts - 1)*A * rho(Npts - 1));
	data.clear();
	data.reserve(100);
	data.push_back(tmpvec);
}

// This is a private function that needs to be called after public Initialization
// This builds the internal variables:
// 	phalf,...Thalf; pnew...Tnew;
//  Uhalf, Fhalf, Shalf; Unew
void LWP::Ini(int _Npts) {
	t = 0.;

	Npts = _Npts;
	dx = L / (Npts - 1.);

	x = VectorXd::Zero(Npts);
	p = VectorXd::Zero(Npts);
	v = VectorXd::Zero(Npts);
	rho = VectorXd::Zero(Npts);
	T = VectorXd::Zero(Npts);
	q = VectorXd::Zero(Npts);

	phalf = VectorXd::Zero(Npts - 1);
	vhalf = VectorXd::Zero(Npts - 1);
	rhohalf = VectorXd::Zero(Npts - 1);
	Thalf = VectorXd::Zero(Npts - 1);

	pnew = VectorXd::Zero(Npts);
	vnew = VectorXd::Zero(Npts);
	Tnew = VectorXd::Zero(Npts);
	rhonew = VectorXd::Zero(Npts);


	U = MatrixXd::Zero(Npts, 3);
	F = MatrixXd::Zero(Npts, 3);
	S = MatrixXd::Zero(Npts, 3);

	Uhalf = MatrixXd::Zero(Npts - 1, 3);
	Fhalf = MatrixXd::Zero(Npts - 1, 3);
	Shalf = MatrixXd::Zero(Npts - 1, 3);

	Unew = MatrixXd::Zero(Npts, 3);
	ini_done = true;
}

void LWP::UpdateTimeStep() {
	dt = 1.e5;
	for (int i = 0; i < Npts; i++) {
		double a = gas->Get_SonicVel(T(i),p(i));
		double dt_new = dx / (fabs(v(i)) + a);
		if (dt_new < dt)
			dt = dt_new;
	}
}

double LWP::Get_dt() {
	return dt;
}

vector<double> LWP::Get_xgrid() {
	std::vector<double> out(Npts);
	for (unsigned int i = 0; i < Npts; i++)
		out.at(i) = i * L / (Npts - 1);
	return out;
}

vector<double> LWP::Get_p() {
	std::vector<double> out(Npts);
	for (unsigned int i = 0; i < Npts; i++)
		out.at(i) = p(i);
	return out;
}

vector<double> LWP::Get_v() {
	std::vector<double> out(Npts);
	for (unsigned int i = 0; i < Npts; i++)
		out.at(i) = v(1);
	return out;
}

void LWP::Step(double dt_req) {

	// Check if required timestep is too big.
	//UpdateTimeStep();

	if (dt_req > dt) {
		cout << endl << endl << " ERROR!! LWP pipe: " << name;
		cout << endl << endl << " LWP::Step(dt_req): dt_req=" << dt_req;
		cout << " > dt=" << dt << " !!!!";
		cout << endl << "Exiting..." << endl;
		exit(-1);
	}
	else
		dt = dt_req; // It is OK to take a smaller timestep

	// Go ahead with the update
	UpdateInternalPoints();

if (is_BCRight_Wall)
	BCRight("Wall",0,0,true);

	// We assume that the boundary points were already updated by the Connector class
	if (!left_boundary_points_already_updated || !right_boundary_points_already_updated){
		cout<<endl<<endl<<"ERROR!!!";
		cout<<endl<<endl<<"LWP::Step() boundary points were not yet updated";
		cout<<endl<<"This typically happens if the Connector class was NOT used earlier.";
		cout<<endl<<"In such cases, use the LWP::Step(string,double,double,string,double,double,double) function instead,";
		cout<<endl<<"which updates the internal points AND the boundaries.";
		cout<<endl<<"Alternatively, call the appropriate Connector class member function before this.";
		cout<<endl<<endl;
		exit(-1);
	}

	// Close timestep
	for (int i = 0; i < Npts; i++) {
		if (pnew(i) < P_MIN)
			pnew(i) = P_MIN;
		if (Tnew(i) < T_MIN)
			Tnew(i) = T_MIN;
		p(i) = pnew(i);
		v(i) = vnew(i);
		T(i) = Tnew(i);
		rho(i) = rhonew(i);
	}
	t += dt;
	UpdateTimeStep();

	left_boundary_points_already_updated=false;
	right_boundary_points_already_updated=false;

	if (save_data) {
		tmpvec.at(0) = t;
		tmpvec.at(1) = p(0);
		tmpvec.at(2) = p(Npts - 1);
		tmpvec.at(3) = v(0);
		tmpvec.at(4) = v(Npts - 1);
		tmpvec.at(5) = T(0);
		tmpvec.at(6) = T(Npts - 1);
		tmpvec.at(7) = rho(0);
		tmpvec.at(8) = rho(Npts - 1);
		tmpvec.at(9) = v(0) * A * rho(0);
		tmpvec.at(10) = v(Npts - 1) * A * rho(Npts - 1);

		data.push_back(tmpvec);
	}

	if (save_all_data) {
		Add_data_row();
	}
}


void LWP::Step(
		string BC_start_type, double BC_start_val1, double BC_start_val2,
		string BC_end_type, double BC_end_val1, double BC_end_val2,
		double dt_req) {

	// Check if required timestep is too big.
	//UpdateTimeStep();

	if (dt_req > dt) {
		cout << endl << endl << " ERROR!! LWP pipe: " << name;
		cout << endl << endl << " LWP::Step(dt_req): dt_req=" << dt_req;
		cout << " > dt=" << dt << " !!!!";
		cout << endl << "Exiting..." << endl;
		exit(-1);
	}
	else
		dt = dt_req; // It is OK to take a smaller timestep

	// Go ahead with the update
	UpdateInternalPoints();

	if(left_boundary_points_already_updated || right_boundary_points_already_updated){
		cout<<endl<<endl<<"ERROR!!!";
		cout<<endl<<endl<<"LWP::Step() boundary points were already updated";
		cout<<endl<<"This typically happens if the Connector class was used earlier.";
		cout<<endl<<"In such cases, use the LWP::Step(dt_req) function insteady,";
		cout<<endl<<"which updates the internal points only."<<endl<<endl;
		exit(-1);
	}
	else{
		BCLeft(BC_start_type, BC_start_val1, BC_start_val2,true);
		left_boundary_points_already_updated=true;
		BCRight(BC_end_type, BC_end_val1, BC_end_val2,true);
		right_boundary_points_already_updated=true;
	}

	// Close timestep
	for (int i = 0; i < Npts; i++) {
		if (pnew(i) < P_MIN)
			pnew(i) = P_MIN;
		if (Tnew(i) < T_MIN)
			Tnew(i) = T_MIN;
		p(i) = pnew(i);
		v(i) = vnew(i);
		T(i) = Tnew(i);
		rho(i) = rhonew(i);
	}
	t += dt;
	UpdateTimeStep();

	left_boundary_points_already_updated=false;
	right_boundary_points_already_updated=false;

	if (save_data) {
		tmpvec.at(0) = t;
		tmpvec.at(1) = p(0);
		tmpvec.at(2) = p(Npts - 1);
		tmpvec.at(3) = v(0);
		tmpvec.at(4) = v(Npts - 1);
		tmpvec.at(5) = T(0);
		tmpvec.at(6) = T(Npts - 1);
		tmpvec.at(7) = rho(0);
		tmpvec.at(8) = rho(Npts - 1);
		tmpvec.at(9) = v(0) * A * rho(0);
		tmpvec.at(10) = v(Npts - 1) * A * rho(Npts - 1);

		data.push_back(tmpvec);
	}

	if (save_all_data) {
		Add_data_row();
	}
}

void LWP::Add_data_row() {
	string pfname;

	FILE * pFile;
	pfname = name + "_p.dat";
	pFile = fopen(pfname.c_str(), "a");
	fprintf(pFile, "%+8.6e", p(0));
	for (int i = 1; i < p.size(); i++)
		fprintf(pFile, "; %+8.6e", p(i));
	fprintf(pFile, "\n");
	fclose (pFile);

	pfname = name + "_v.dat";
	pFile = fopen(pfname.c_str(), "a");
	fprintf(pFile, "%+8.6e", v(0));
	for (int i = 1; i < v.size(); i++)
		fprintf(pFile, ";%+8.6e ", v(i));
	fprintf(pFile, "\n");
	fclose (pFile);

	pfname = name + "_T.dat";
	pFile = fopen(pfname.c_str(), "a");
	fprintf(pFile, "%+8.6e", T(0));
	for (int i = 1; i < T.size(); i++)
		fprintf(pFile, ";%+8.6e ", T(i));
	fprintf(pFile, "\n");
	fclose (pFile);

	pfname = name + "_rho.dat";
	pFile = fopen(pfname.c_str(), "a");
	fprintf(pFile, "%+8.6e", rho(0));
	for (int i = 1; i < rho.size(); i++)
		fprintf(pFile, ";%+8.6e; ", rho(i));
	fprintf(pFile, "\n");
	fclose (pFile);

	pfname = name + "_mdot.dat";
	pFile = fopen(pfname.c_str(), "a");
	fprintf(pFile, "%+8.6e", rho(0)*v(0)*A);
	for (int i = 1; i < T.size(); i++)
		fprintf(pFile, ";%+8.6e ", rho(i)*v(i)*A);
	fprintf(pFile, "\n");
	fclose (pFile);
}

void LWP::UpdateInternalPoints() {

	Pack(/*is_half_step*/false);

	// Predictor step
	double Umean, diff_F, Smean;

	for (unsigned int j = 0; j < 3; j++) {
		for (unsigned int i = 0; i < Npts - 1; i++) {
			Umean  = (U(i, j) + U(i + 1, j)) / 2.;
			diff_F = (F(i + 1, j) - F(i, j)) / dx;
			Smean  = (S(i, j) + S(i + 1, j)) / 2.;
			Uhalf(i, j) = Umean + (Smean - diff_F) * dt / 2.;
		}
	}

	// Unpack primitive variables from Uhalf
	UnPackU(/*is_half_step*/true);
	Pack(/*is_half_step*/true);

	// Full step
	for (unsigned int j = 0; j < 3; j++) {
		for (unsigned int i = 1; i < Npts - 1; i++) {
			diff_F = (Fhalf(i, j) - Fhalf(i - 1, j)) / dx;
			Smean  = (Shalf(i, j) + Shalf(i - 1, j)) / 2.;
			Unew(i, j) = U(i, j) + (Smean - diff_F) * dt;
		}
	}

	// Unpack Unew
	UnPackU(/*is_half_step*/false);
}

void LWP::Pack(bool is_half_step) {
	double e;
	if (is_half_step) {
		for (unsigned int i = 0; i < Npts - 1; i++) {
			e = gas->Get_e_from_Tp(Thalf(i),phalf(i));
			Fhalf(i, 0) = rhohalf(i) * vhalf(i) * A;
			Fhalf(i, 1) = ( rhohalf(i) * vhalf(i) * vhalf(i) + phalf(i) ) * A;
			Fhalf(i, 2) = ( rhohalf(i) * vhalf(i) * e + phalf(i) * vhalf(i) ) * A;

			Shalf(i, 0) = 0.;
			Shalf(i, 1) = -A * rhohalf(i) / 2. * lambda / D * vhalf(i) * fabs(vhalf(i));
			Shalf(i, 2) = HeatFluxThroughWall(Thalf(i),vhalf(i),rhohalf(i),phalf(i));
		}
	}
	else {
		for (unsigned int i = 0; i < Npts; i++) {
			e = gas->Get_e_from_Tp(T(i),p(i));
			U(i, 0) = rho(i) * A;
			U(i, 1) = rho(i) * v(i) * A;
			U(i, 2) = rho(i) * e * A;

			F(i, 0) = rho(i) * v(i) * A;
			F(i, 1) = ( rho(i) * v(i) * v(i) + p(i) ) * A;
			F(i, 2) = ( rho(i) * v(i) * e + p(i) * v(i) ) * A;

			S(i, 0) = 0.;
			S(i, 1) = -A * rho(i) / 2. * lambda / D * v(i) * fabs(v(i));
			S(i, 2) = HeatFluxThroughWall(T(i),v(i),rho(i),p(i));
		}
	}
}

void LWP::UnPackU(bool is_half_step) {

	if (is_half_step) {
		for (unsigned int i = 0; i < Npts - 1; i++) {
			rhohalf(i) = Uhalf(i, 0) / A;
			vhalf(i) = Uhalf(i, 1) / Uhalf(i, 0);
			Thalf(i) = gas->Get_T_from_erho(Uhalf(i, 2) / Uhalf(i, 0),rhohalf(i));
			phalf(i) = gas->Get_p(rhohalf(i), Thalf(i));
		}
	}
	else {
		for (unsigned int i = 1; i < Npts-1; i++) {
			rhonew(i) = Unew(i, 0) / A;
			vnew(i) = Unew(i, 1) / Unew(i, 0);
			Tnew(i) = gas->Get_T_from_erho(Unew(i, 2) / Unew(i, 0),rhonew(i));
			pnew(i) = gas->Get_p(rhonew(i), Tnew(i));
		}
		// artificial viscosity
		double dv;

		//double pp = 0.;
		for (unsigned int i = 1; i < Npts - 1; i++) {
			dv = vnew(i - 1) - vnew(i);
			if (dv > 0)
				q(i) = -art_visc * art_visc * rhonew(i) * dv * dv;
			else
				q(i) = 0.;
			pnew(i) += q(i);
		}
	}
}

bool LWP::BCLeft(string type, double val1, double val2,bool write_to_node) {

	double beta_primitive = GetBetaPrimitiveAtFront(t + dt);
	double vBC, pBC, TBC,rhoBC;
	bool BC_result_is_consistent=false;
	bool ok = false;

	if (type == "Wall") {
		vBC = 0;
		pBC = beta_primitive;
		double p_old = p(0);
		//double r_old = rho(0);
		double T_old = T(0);
		// Korrigalni kell majd meg a sebesseggel
		double kappa = gas->Get_kappa_Tp();
		TBC = T_old * pow(pBC / p_old, (kappa - 1) / kappa);
		rhoBC = gas->Get_rho(pBC, TBC);
		ok = true;
		BC_result_is_consistent=true;
	}

	/*if (type == "MassFlowIn_and_T") {
	  double mp = val1;
	  Tnew(0) = val2;
	  double  a = gas->Get_SonicVel( Tnew(0));
	  vnew(0) = (a - beta) * 2. / (kappa - 1);
	  if (vnew(0) < 0) {
	  cout << endl << endl << " ERROR!! LWP pipe: " << name;
	  cout << endl << endl << " LWP::BCLeft(): MassFlowIn_and_T -> negative velocity";
	  cout << endl << "  prescribed mass flow rate: " << mp << " kg/s";
	  cout << endl << "  prescribed temperature   : " << Tnew(0) - 273.15 << " C";
	  cout << endl << "  computed flow velocity   : " << v << " m/s";
	  cout << endl << " Something is wrong. Try increasing the node number and/or decreasing the timestep.";
	  cout << endl << "Exiting..." << endl;
	  exit(-1);
	  }
	  rhonew(0) = mp / A / vnew(0);
	  pnew(0) = gas->GetP(rhonew(0), Tnew(0));
	  ok = true;
	  }*/

	if (type == "StaticPres_and_StaticTemp_Inlet") {
		pBC  = val1;
		TBC  = val2;
		rhoBC = gas->Get_rho(pBC, TBC);
		vBC  = (pBC - beta_primitive) / rhoBC / gas->Get_SonicVel(TBC,pBC);
		BC_result_is_consistent=true;
		if (vBC<0){
			BC_result_is_consistent=false;
			vBC=0;
if (write_to_node){
			cout << endl << endl << " WARNING!! LWP pipe: " << name;
			cout << endl << "\tLWP::BCLeft(): StaticPres_and_StaticTemp_Inlet -> negative velocity";
			cout << endl << "\tUsing Wall BC instead. If the situation persists, switch to other boundary condition";
			cout<<endl;
}
			// Just copied the Wall BC here
			pBC = beta_primitive;
			double p_old = p(0);
			double T_old = T(0);
			double kappa = gas->Get_kappa_Tp();
			TBC = T_old * pow(pBC / p_old, (kappa - 1) / kappa);
			rhoBC = gas->Get_rho(pBC, TBC);

			//printf("\tvBC=%5.3e, pBC=%5.3e, TBC=%5.3e, rhoBC=%5.3e\n",vBC,pBC,TBC,rhoBC);
			//cin.get();
		}
		ok = true;
	}

	if (type == "StaticPres_Outlet") {
		pBC = val1;
		rhoBC=GetC0AtFront(t+dt);
		TBC = gas->Get_T(rhoBC,pBC);
		double aBC=gas->Get_SonicVel(TBC,pBC);
		vBC = (pBC-beta_primitive)/aBC/rhoBC;
//		cout<<endl<<"pBC   = "<<pBC;
//		cout<<endl<<"rhoBC = "<<rhoBC;
//		cout<<endl<<"TBC   = "<<TBC;
//		cout<<endl<<"aBC      = "<<aBC;
//		cout<<endl<<"vBC   = "<<vBC;
//		cin.get();
		BC_result_is_consistent=true;
		if (vBC > 0) {
			BC_result_is_consistent=false;
			cout<<endl<<"LWP -> BCLeft -> StaticPres_Outlet";
			cout << endl << "WARNING: vBC = " << vBC << ", should be negative";
		}
		ok = true;
	}

	if (!ok) {
		cout << endl << endl
			<< "ERROR! LWP::BCLeft(), unknown BC type: " << type << endl;
		cout << "Possible choices:" << endl;
		cout << "\t Wall" << endl;
		//cout << "\t MassFlowIn_and_T" << endl;
		cout << "\t StaticPres_and_StaticTemp_Inlet" << endl;
		cout << "\t StaticPres_Outlet" << endl;
		//cout << "\t TotalPres_and_TotalTemp" << endl;
		//cout << "\t TotalPres_and_TotalTemp_Isentropic" << endl;
		//cout << "\t Opening" << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}

	if (write_to_node){
		vnew(0)=vBC;
		pnew(0)=pBC;
		rhonew(0)=rhoBC;
		Tnew(0)=TBC;
		left_boundary_points_already_updated=true;
	}

	return BC_result_is_consistent;
}

bool LWP::BCRight(string type, double val1, double val2,bool write_to_node) {


	double alpha_primitive = GetAlphaPrimitiveAtEnd(t + dt);
	double vBC,pBC,rhoBC,TBC;
	bool ok = false;
	bool BC_result_is_consistent=true;

	if (type == "Wall") {
		vBC = 0.;
		pBC = alpha_primitive;
		double T_old = T(Npts - 1);
		double p_old = p(Npts - 1);
		double kappa_Tp = gas->Get_kappa_Tp();
		TBC = T_old * pow(pBC / p_old, (kappa_Tp - 1) / (kappa_Tp));
		rhoBC = gas->Get_rho(pBC, TBC);
		ok = true;
	}

	if (type == "StaticPres_and_StaticTemp") {
		pBC   = val1;
		TBC   = val2;
		double a = gas->Get_SonicVel(TBC,pBC);
		alpha = GetAlphaPrimitiveAtEnd(t + dt);
		rhoBC = gas->Get_rho(pBC, TBC);
		vBC   = (alpha - pBC) / rhoBC / a;
		ok = true;
	}

	if (!ok) {
		cout << endl << endl
			<< "ERROR! LWP::BCRight(), unknown BC type: " << type << endl;
		cout << "Possible choices:" << endl;
		cout << "\t Wall" << endl;
		//cout << "\t MassFlowIn_and_T" << endl;
		cout << "\t StaticPres_and_StaticTemp" << endl;
		//cout << "\t TotalPres_and_TotalTemp" << endl;
		//cout << "\t TotalPres_and_TotalTemp_Isentropic" << endl;
		//cout << "\t Opening" << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}

	if (write_to_node){
		vnew(Npts-1)=vBC;
		pnew(Npts-1)=pBC;
		rhonew(Npts-1)=rhoBC;
		Tnew(Npts-1)=TBC;
		right_boundary_points_already_updated=true;
	}

	return BC_result_is_consistent;

}

/*double LWP::GetBetaAtFront(double t_target) {

  double delta_t = t_target - t;
  double TOL = dt / 1000.;

  if (delta_t < 0) {
  if (fabs(delta_t) < TOL)
  delta_t = 0.;
  else {
  cout << endl
  << "ERROR! LWP::GetBetaAtFront(), delta_t = " << delta_t << " < 0 ! (TOL=" << TOL << ")" << endl;
  cout << endl << "Name of pipe: " << name << endl;
  cin.get();
  }
  }
  if (delta_t > dt) {
  if (delta_t - TOL < dt)
  delta_t = dt;
  else {
  cout << endl
  << "ERROR! LWP::GetBetaAtFront(), delta_t = " << delta_t << " > dt= " << dt << endl;
  cout << endl << "Name of pipe: " << name << endl;
  cin.get();
  }
  }

  double kR = (gas->kappa) * (gas->R);
  double a0 = sqrt(kR * T(0));
  double v0 = v(0);
  double a1 = sqrt(kR * T(1));
  double v1 = v(1);
  double dxL = dx * (a0 - v0) / (a0 - a1 - v0 + v1 + dx / delta_t);

  if (dxL < 0) {
  cout << endl
  << "ERROR! LWP::GetAlphaAtEnd(), dxL = " << dxL << " < 0 !!!" << endl;
  cout << endl << "Name of pipe: " << name << endl;
  cin.get();
  }
  if (dxL > dx) {
  cout << endl
  << "ERROR! LWP::GetAlphaAtEnd(), dxL/dx = " << dxL / dx << " > 1 !!!" << endl;
  cout << endl << "Name of pipe: " << name << endl;
  cin.get();
  }

  double vv = v0 * (dx - dxL) / dx + v1 * dxL / dx;
  double aa = a0 * (dx - dxL) / dx + a1 * dxL / dx;
  double beta = aa - (gas->kappa - 1.) / 2.*vv;

  return beta;

  }*/

void LWP::GetAllPrimitiveAtFront(double t_target, double& pR, double& vR, double& TR, double& rhoR) {

	double delta_t = t_target - t;
	double TOL = dt / 1000.;

	if (delta_t < 0) {
		if (fabs(delta_t) < TOL)
			delta_t = 0.;
		else {
			cout << endl
				<< "ERROR! LWP::GetBetaPrimitiveAtFront(), delta_t = " << delta_t << " < 0 ! (TOL=" << TOL << ")" << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
	}
	if (delta_t > dt) {
		if (delta_t - TOL < dt)
			delta_t = dt;
		else {
			cout << endl
				<< "ERROR! LWP::GetBetaPrimitiveAtFront(), delta_t = " << delta_t << " > dt= " << dt << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
	}

	//double kR = (gas->kappa) * (gas->R);
	double a0 = gas->Get_SonicVel(T(0),p(0));
	double v0 = v(0);
	double a1 = gas->Get_SonicVel(T(1),p(1));
	double v1 = v(1);
	double dxL = dx * (a0 - v0) / (a0 - a1 - v0 + v1 + dx / delta_t);

	if (dxL < 0) {
		cout << endl
			<< "ERROR! LWP::GetBetaPrimitiveAtFront(), dxL = " << dxL << " < 0 !!!" << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}
	if (dxL > dx) {
		cout << endl
			<< "ERROR! LWP::GetBetaPrimitiveAtFront(), dxL/dx = " << dxL / dx << " > 1 !!!" << endl;
		cout << endl << "Name of pipe: " << name << endl;
		//cin.get();
	}

	vR   = v(0)   * (dx - dxL) / dx + v(1) * dxL / dx;
	TR   = T(0)   * (dx - dxL) / dx + T(1) * dxL / dx;
	pR   = p(0)   * (dx - dxL) / dx + p(1) * dxL / dx;
	rhoR = rho(0) * (dx - dxL) / dx + rho(1) * dxL / dx;

}

double LWP::GetBetaPrimitiveAtFront(double t_target) {
	double pR, vR, TR, rhoR;
	GetAllPrimitiveAtFront(t_target, pR, vR, TR, rhoR);
	return pR - rhoR * gas->Get_SonicVel(TR,pR) * vR;
}

//bool LWP::GetC0AtFront(double t_target, double& pM, double& vM, double& TM, double& rhoM) {
double LWP::GetC0AtFront(double t_target) {
	bool DEBUG=false;
	double delta_t = t_target - t;
	double TOL = dt / 1000.;
	bool is_ok = true;

	if (delta_t < 0) {
		if (fabs(delta_t) < TOL)
			delta_t = 0.;
		else {
			cout << endl
				<< "ERROR! LWP::GetC0AtFront(), delta_t = " << delta_t << " < 0 ! (TOL=" << TOL << ")" << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
	}
	if (delta_t > dt) {
		if (delta_t - TOL < dt)
			delta_t = dt;
		else {
			cout << endl
				<< "ERROR! LWP::GetC0AtFront(), delta_t = " << delta_t << " > dt= " << dt << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
	}

	//double a0 = gas->Get_SonicVel(T(0));
	double v0 = v(0);
	//double a1 = gas->Get_SonicVel(T(1));
	double v1 = v(1);
	double dxM = dx * v0 / (v0 - v1 - dx / delta_t);

	if (dxM < 0.) {
		is_ok = false;
		if (DEBUG){
			cout << endl << "WARNING! LWP::GetC0AtFront(), dxM/dx = " << dxM / dx << " < 0 !!!";
			if (v(0) < 0) {
				cout << endl << "dx=" << dx;
				cout << endl << "dt=" << delta_t;
				cout << endl << "v(0)=" << v(0) << " < 0 ";
				cout << endl << "v(1)=" << v(1);
				cout << endl << "we should have v(0)< v1+dx/dt = " << v(1) << " + " << dx / delta_t << " = " << v(1) + dx / delta_t;
			}
			if (v(0) > 0) {
				cout << endl << "dx=" << dx;
				cout << endl << "dt=" << delta_t;
				cout << endl << "v(0)=" << v(0) << " > 0 ";
				cout << endl << "v(1)=" << v(1);
				cout << endl << "we should have v(0)> v1+dx/dt = " << v(1) << " + " << dx / delta_t << " = " << v(1) + dx / delta_t;
			}
			cout << endl << "Name of pipe: " << name << " --> is_ok set to false";
			//cin.get();
		}
	}
	if (dxM > dx) {
		is_ok = false;
		if (DEBUG){
			cout << endl
				<< "ERROR! LWP::GetC0AtFront(), dxM/dx = " << dxM / dx << " > 1 !!!" << endl;
			cout << endl << "Name of pipe: " << name << endl;
			//cin.get();
		}
	}


	double vM   = v(0)   * (1. - dxM / dx) + v(1) * dxM / dx;
	double TM   = T(0)   * (1. - dxM / dx) + T(1) * dxM / dx;
	double pM   = p(0)   * (1. - dxM / dx) + p(1) * dxM / dx;
	double rhoM = rho(0) * (1. - dxM / dx) + rho(1) * dxM / dx;

	if (DEBUG){
		if (!is_ok) {
			printf("\n\n--------------------\n dxM/dx=%5.3f", dxM / dx);
			printf("\n         0           M           1");
			printf("\n p  =   %5.3e    %5.3e    %5.3e", p(0) / 1.e5, pM / 1.e5, p(1) / 1.e5);
			printf("\n rho=   %5.3e    %5.3e    %5.3e", rho(0), rhoM, rho(1));
			printf("\n T  =   %5.3e    %5.3e    %5.3e", T(0), TM, T(1));
			printf("\n v  =  %+5.3e   %+5.3e   %+5.3e\n", v(0), vM, v(1));
		}
	}
	return rhoM;

}

/*double LWP::GetAlphaAtEnd(double t_target) {
  double delta_t = t_target - t;
  double TOL = dt / 1000.;
  */
/*cout<<endl<<"entering GetAlphaAtEnd():";
  cout<<endl<<"t_target   : "<<t_target;
  cout<<endl<<"delta_t/dt : "<<delta_t/dt;
  cin.get();*/
/*
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

   double kR = (gas->kappa) * (gas->R);
   double aNm1 = sqrt(kR * T(Npts - 2));
   double vNm1 = v(Npts - 2);
   double aN = sqrt(kR * T(Npts - 1));
   double vN = v(Npts - 1);
   double dxR = dx * (aN + vN) / (aN - aNm1 + vN - vNm1 + dx / delta_t);

*/
/*cout<<endl<<" aNm1   = :"<<aNm1;
  cout<<endl<<" aN     = :"<<aN;
  cout<<endl<<" vNm1   = :"<<vNm1;
  cout<<endl<<" vN     = :"<<vN;
  cout<<endl<<" dxR/dx = :"<<dxR/dx<<endl;
  cin.get();*/
/*
   if (dxR < 0) {
   cout << endl
   << "ERROR! LWP::GetAlphaAtEnd(), dxR = " << dxR << " < 0 !!!" << endl;
   cout << endl << "Name of pipe: " << name << endl;
   cin.get();
   }
   if (dxR > dx) {
   cout << endl
   << "ERROR! LWP::GetAlphaAtEnd(), dxR/dx = " << dxR / dx << " > 1 !!!" << endl;
   cout << endl << "Name of pipe: " << name << endl;
   cin.get();
   }

   double vv = vN * (dx - dxR) / dx + vNm1 * dxR / dx;
   double aa = aN * (dx - dxR) / dx + aNm1 * dxR / dx;
   double alpha = aa + (gas->kappa - 1.) / 2.*vv;

   return alpha;

   }*/

double LWP::GetC0AtEnd(double t_target) {
	double delta_t = t_target - t;
	double TOL = dt / 1000.;

	/*cout<<endl<<"entering GetAlphaAtEnd():";
	  cout<<endl<<"t_target   : "<<t_target;
	  cout<<endl<<"delta_t/dt : "<<delta_t/dt;
	  cin.get();*/

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
				<< "ERROR! SCP::GetC0AtEnd(), delta_t = " << delta_t << " > dt= " << dt << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
	}

	//double kR = (gas->kappa) * (gas->R);
	double aNm1 = gas->Get_SonicVel(T(Npts - 2),p(Npts - 1));
	double vNm1 = v(Npts - 2);
	double rhoNm1 = rho(Npts - 2);
	double aN = gas->Get_SonicVel(T(Npts - 1),p(Npts - 1));
	double vN = v(Npts - 1);
	double rhoN = rho(Npts - 1);
	double dxR = dx * (0 * aN + vN) / (0 * aN - 0 * aNm1 + vN - vNm1 + dx / delta_t);


	/*cout<<endl<<" aNm1   = :"<<aNm1;
	  cout<<endl<<" aN     = :"<<aN;
	  cout<<endl<<" vNm1   = :"<<vNm1;
	  cout<<endl<<" vN     = :"<<vN;
	  cout<<endl<<" dxR/dx = :"<<dxR/dx<<endl;
	  cin.get();*/

	if (dxR < 0) {
		cout << endl
			<< "ERROR! LWP::GetC0AtEnd, dxR = " << dxR << " < 0 !!!" << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}
	if (dxR > dx) {
		cout << endl
			<< "ERROR! LWP::GetC0AtEnd, dxR/dx = " << dxR / dx << " > 1 !!!" << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}

	double rho0 = rhoN * (dx - dxR) / dx + rhoNm1 * dxR / dx;
	//double aa = aN * (dx - dxR) / dx + aNm1 * dxR / dx;
	//double alpha = aa + (gas->kappa - 1.) / 2.*vv;

	return rho0;

}

void LWP::GetAllPrimitiveAtEnd(double t_target, double& pP, double& vP, double& TP, double& rhoP) {
	double delta_t = t_target - t;
	double TOL = dt / 1000.;

	if (delta_t < 0) {
		if (fabs(delta_t) < TOL)
			delta_t = 0.;
		else {
			cout << endl << "Name of pipe: " << name << endl;
			cout << endl
				<< "ERROR! LWP::GetAllPrimitiveAtEnd(), delta_t = " << delta_t << " < 0  ! (TOL=" << TOL << ")" << endl;
			cout << endl << " t_pipe = " << t << ", t_target=" << t_target << endl;
			cin.get();
		}
	}
	if (delta_t > dt) {
		if (delta_t - TOL < dt)
			delta_t = dt;
		else {
			cout << endl
				<< "ERROR! LWP::GetAllPrimitiveAtEnd(), delta_t = " << delta_t << " > dt= " << dt << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
	}

	//	double kR = (gas->kappa) * (gas->R);
	double aNm1 = gas->Get_SonicVel(T(Npts - 2),T(Npts-2));
	double vNm1 = v(Npts - 2);
	double aN = gas->Get_SonicVel(T(Npts - 1),p(Npts-1));
	double vN = v(Npts - 1);
	double dxR = dx * (aN + vN) / (aN - aNm1 + vN - vNm1 + dx / delta_t);


	/*cout<<endl<<" aNm1   = :"<<aNm1;
	  cout<<endl<<" aN     = :"<<aN;
	  cout<<endl<<" vNm1   = :"<<vNm1;
	  cout<<endl<<" vN     = :"<<vN;
	  cout<<endl<<" dxR/dx = :"<<dxR/dx<<endl;
	  cin.get();*/

	double TOL_dx_rel = 0.1 / 100.;
	if (dxR < 0.) {
		if (fabs(dxR / dx) > TOL_dx_rel) {
			cout << endl
				<< "ERROR! LWP::GetAllPrimitiveAtEnd(), dxR = " << dxR << " < 0 !!! (dx=" << dx << ")" << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
		else
			dxR = 0.;
	}
	if (dxR > dx) {
		if (fabs((dxR - dx) / dx) > TOL_dx_rel) {
			cout << endl
				<< "ERROR! LWP::GetAllPrimitiveAtEnd(), dxR/dx = " << dxR / dx << " > 1 !!!" << endl;
			cout << endl << "Name of pipe: " << name << endl;
			cin.get();
		}
		else
			dxR = dx;
	}

	vP = v(Npts - 1) * (dx - dxR) / dx + v(Npts - 2) * dxR / dx;
	TP = T(Npts - 1) * (dx - dxR) / dx + T(Npts - 2) * dxR / dx;
	pP = p(Npts - 1) * (dx - dxR) / dx + p(Npts - 2) * dxR / dx;
	rhoP = rho(Npts - 1) * (dx - dxR) / dx + rho(Npts - 2) * dxR / dx;
	//double rhoL = pL+roL*aL*vL;

	//return rhoP;

}


double LWP::GetAlphaPrimitiveAtEnd(double t_target) {
	double pL, vL, TL, rhoL;
	GetAllPrimitiveAtEnd(t_target, pL, vL, TL, rhoL);

	return pL + rhoL * gas->Get_SonicVel(TL,pL) * vL;

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
void LWP::UpdateDimlessPars(double pref, double mp_nevl, double omega, double xref, double m) {
	if (!ini_done) {
		IniUniform(0, 1e5, 293, 20);
		cout << endl;
		cout << "WARNING! Trying to call LWP::UpdateDimlessPars without initializing!";
		cout << endl << "Name of pipe: " << name << endl;
		cout << "Initializing with v=0, p=1bar, T=283K, Npts=20" << endl;
	}

	phi   = lambda * xref / 2. / D;
	double ro = gas->Get_rho(p(0), T(0));
	double a  = gas->Get_SonicVel(T(0),p(0));
	alpha = ro * A * a / m / omega;
	gamma = L * omega / a;
	mu    = ro * A * xref * omega / mp_nevl;
}

void LWP::Set_dprop(string prop_string, double val) {
	if (prop_string == "art_visc") {
		art_visc = val;
	} else {
		cout << endl
			<< "HIBA! LWP::Set_dprop(prop_string,val), unknown property: prop_string=" << prop_string << endl
			<< endl;
	}
}

vector<double> LWP::Get_dvprop(string prop_string) {
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
			out.at(i) = data.at(i).at(9);
	else if (prop_string == "mp_back")
		for (unsigned int i = 0; i < Ntime; i++)
			out.at(i) = data.at(i).at(10);
	else {
		cout << endl
			<< "ERROR! LWP::Get_dvprop(prop_string), unknown input: prop_string=" << prop_string << endl;
		cout<<"options: t | p_front | p_front_bar | p_back | p_back_bar | v_front | v_back | mp_front | mp_back";
		cout << endl << "Name of the LWP: " << name << endl;
		cin.get();
	}
	return out;
}


double LWP::Get_dprop(string prop_string) {
	double out = 0.0;
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
	else if (prop_string == "t")
		out = t;
	else if (prop_string == "p_front")
		out = p(0);
	else if (prop_string == "p_back")
		out = p(Npts - 1);
	else if (prop_string == "p_front_bar")
		out = p(0) / 1.e5;
	else if (prop_string == "p_back_bar")
		out = p(Npts - 1) / 1.e5;
	else if (prop_string == "v_front")
		out = v(0);
	else if (prop_string == "v_back")
		out = v(Npts - 1);
	else if (prop_string == "rho_front")
		out = rho(0);
	else if (prop_string == "rho_back")
		out = rho(Npts - 1);
	else if (prop_string == "T_front")
		out = T(0);
	else if (prop_string == "T_back")
		out = T(Npts - 1);
	else if (prop_string == "a_front")
		out = gas->Get_SonicVel(T(0),p(0));
	else if (prop_string == "a_back")
		out = gas->Get_SonicVel(T(Npts - 1),p(Npts - 1));
	else if (prop_string == "M_front")
		out = fabs(v(0)) / gas->Get_SonicVel(T(0),p(0));
	else if (prop_string == "M_back")
		out = fabs(v(Npts - 1)) / gas->Get_SonicVel(T(Npts - 1),p(Npts - 1));
	else if (prop_string == "mp_front")
		out = v(0) * gas->Get_rho(p(0), T(0)) * A;
	else if (prop_string == "mp_back")
		out = v(Npts - 1) * gas->Get_rho(p(Npts - 1), T(Npts - 1)) * A;
	else if (prop_string == "frek")
		out = gas->Get_SonicVel(T.mean(),p.mean()) / (2.*L);
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
	else if (prop_string == "art_visc")
		out = art_visc;
	else if (prop_string == "rho_mean")
		out = gas->Get_rho(p.mean(), T.mean());
	else if (prop_string == "a_mean")
		out = gas->Get_SonicVel(T.mean(),p.mean());
	else {
		cout << endl
			<< "ERROR! LWP::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl
			<< endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}
	return out;
}

/*! Save time history at the two ends of the pipe.
*/

void LWP::Save_data() {
	//char fname [50];
	//sprintf (fname, "%s.dat", name.c_str());
	if (!save_data) {
		cout << endl << "WARNING! LWP: " << name;
		cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
	}
	else {
		cout << endl << "Saving to " << fname.c_str() << " ... ";

		FILE * pFile;
		pFile = fopen (fname.c_str(), "w");
		fprintf(pFile, "t (s);") ;					// 0
		fprintf(pFile, "p(0) (bar); p(L) (bar);"); 	// 1-2
		fprintf(pFile, "v(0) m/s; v(L) (m/s);"); 	// 3-4
		fprintf(pFile, "T(0) (C); T(L) (C);"); 		// 5-6
		fprintf(pFile, "rho(0) (kg/m3); rho(L) (kg/m3);");	// 7-8
		fprintf(pFile, "mp(0) (kg/s), mp(L) (kg/s)\n");		// 9-10

		for (int i = 0; i < data.size(); i++)
			fprintf(pFile, "%8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e\n",
					data.at(i).at(0),
					data.at(i).at(1) / 1.e5, data.at(i).at(2) / 1.e5,
					data.at(i).at(3), data.at(i).at(4),
					data.at(i).at(5), data.at(i).at(6),
					data.at(i).at(7), data.at(i).at(8),
					data.at(i).at(9), data.at(i).at(10));
		fclose (pFile);
		cout << " done. ";
	}
}

/*! Save current distribution of all variables to file.
*/

void LWP::Save_distribution() {
	cout << endl << "Saving distribution to " << name+"_distr.dat" << " ... ";

	FILE * pFile;
	string tmp=name+"_distr.dat";
	pFile = fopen (tmp.c_str(), "w");
	fprintf(pFile, " t=%5.3e s",t) ;	
	fprintf(pFile, "\n x (m), p (Pa), v (m/s), T (K), rho (kg/m3)");
	if (HEAT_TRANSFER_ON)
		fprintf(pFile,", alpha (W/m^2/K), Q (W)");
	for (int i = 0; i < p.size(); i++){
		fprintf(pFile, "\n%8.6e; %8.6e; %8.6e; %8.6e; %8.6e",x(i),p(i),v(i),T(i),rho(i));
		if (HEAT_TRANSFER_ON){
			double alpha_heat;
			double Q = HeatFluxThroughWall(T(i),v(i),rho(i),p(i),alpha_heat);
			fprintf(pFile,"; %8.6e; %8.6e",alpha_heat,Q);
		}
	}
	fclose (pFile);
	cout << " done. ";
}
/*! \brief Sets up heat transfer through wall 
  \param Tw [in] Wall (reference) temperature
  \param model [in] "adiabatic", "constant_alpha" or "Dittus-Boelter"
  \param tmp [in] if model="constant_alpha", this sets the heat transfer coeff alpha
  If heat transfer is switched on, Save_distribution() will add alpha and Q values for each node.
  */

void LWP::SetHeatTransfer(double _Tw,string model, double tmp){
	HEAT_TRANSFER_ON=true;
	Tw=_Tw;
	HEAT_TRANSFER_MODEL=model;
	heat_transfer_par=tmp;
}

double LWP::HeatFluxThroughWall(double Tgas, double v, double rho, double p){
	double alpha_heat;
	return LWP::HeatFluxThroughWall(Tgas, v, rho, p, alpha_heat);
}

double LWP::HeatFluxThroughWall(double Tgas, double v, double rho, double p, double &alpha_heat){
	double q=0.;
	if (HEAT_TRANSFER_MODEL=="adiabatic"){
		q=0.; alpha_heat=0.0;
	}
	else if (HEAT_TRANSFER_MODEL=="constant_alpha"){
		// convective heat transfer coeff.
		alpha_heat=heat_transfer_par;
		// heat flux, W/s on length element dx
		q=-alpha_heat*D*M_PI*dx*(Tgas-Tw);
	}
	else if (HEAT_TRANSFER_MODEL=="Dittus-Boelter"){
		// Heat transfer through pipe wall for turbulent flow.
		// Dittus-Boelter, for L/D>60
		// Q=alpha*(Tgas-Tw)*A (W/s)
		// Nu=alpha*L/lambda
		// Nu=0.023*Re^0.8*Pr^0.4
		// Re=v*D/nu
		// Pr=<CoolProp par>

		// Input data
		//cout<<endl<<"------------------------";
		//cout<<endl<<"Tgas="<<Tgas;
		//cout<<endl<<"p   ="<<p;
		//cout<<endl<<"rho ="<<rho;
		//cout<<endl<<"v   ="<<v;

		// Thermal conductivity
		double lambda_gas = gas->Get_ThermalConductivity_from_Tp(Tgas,p);
		//cout<<endl<<"lambda_gas="<<lambda_gas;

		// Prandtl number
		double Pr=gas->Get_Prandtl_from_Tp(Tgas,p);
		//cout<<endl<<"Pr="<<Pr;

		// Kinematic viscosity
		double nu = gas->Get_DynamicViscosity_from_Tp(Tgas,p)/rho;
		//cout<<endl<<"nu="<<nu;

		// Reynolds number
		double Re = fabs(v)*D/nu;
		//cout<<endl<<"Re="<<Re;

		// Nu correlation
		double Nu=0.023*pow(Re,0.8)*pow(Pr,0.4);
		//cout<<endl<<"Nu="<<Nu;

		// convective heat transfer coeff.
		alpha_heat=Nu*lambda_gas/dx;
		//cout<<endl<<"alpha="<<alpha_heat;
		// heat flux, W/s on length element dx
		q=-alpha_heat*D*M_PI*dx*(Tgas-Tw);
	}
	else{
		cout << endl
			<< "ERROR! LWP::HeatFluxThroughWall(), unknown heat transfer model: prop_string=" << HEAT_TRANSFER_MODEL << endl;
		cout<<"Valid options: adiabatic | constant_alpha | Dittus-Boelter"<<endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}
	return q;
}

