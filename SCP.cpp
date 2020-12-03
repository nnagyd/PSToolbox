#define _USE_MATH_DEFINES
#include "SCP.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "my_tools.h"

#include <math.h>

Ideal_Gas::Ideal_Gas(
    const double _kappa,
    const double _R ) {

	kappa = _kappa; //Adiabatic exponent
	R = _R;
	cp = R * kappa / (kappa - 1.);
	cv = R / (kappa - 1.); //Heat capacity at constant volume
}
//! Copy contructor
/*!  Copies data from an already existing class*/
Ideal_Gas::Ideal_Gas(Ideal_Gas const & g) {
	kappa = g.kappa;
	R = g.R;
}


//! Get the density of the gas
/*! Calculate the density of an ideal gas at a given pressure and
	temperature. All units are to be in SI.
	\sa GetP()
	\param p Pressure in [Pa]
	\param T Temperature in [K]
	\return The density in [kg/m^3]
*/
double Ideal_Gas::GetRho(double p, double T) {
	double rho = p / R / T;
	return rho;
}

/*! \brief Calculate the density.
  Calculate the density using the ideal gas equation.
  \param p [in] Density
  \param T [in] Temperature
  \return The pressure
*/
double Ideal_Gas::GetP(double rho, double T) {
	double p = rho * R * T;
	return p;
}


/*! \brief Get the speed of sound
	Calculate the speed of sound based on the temperature of the gas.
	\param [in] T The temperature
	\return The speed of sound
*/
double Ideal_Gas::GetSonicVel(double T) {
	return sqrt(kappa * R * T);
}
//! Destructor of the Ideal gas function
Ideal_Gas::~Ideal_Gas() {}

//==================================================

//! Unit conversion constants
/*! Class filled up with conversion coefficients between various Imperial
	and SI measurements. */
Units::Units() {
	ft_to_m = 304.8 / 1000.;
	m_to_ft = 1 / ft_to_m;
	inch_to_m = 2.54 / 100.;
	m_to_inch = 1. / inch_to_m;
	lbm_to_kg = 0.45359237;
	kg_to_lbm = 1. / lbm_to_kg;
	lbf_to_N = 4.44822;
	N_to_lbf = 1. / lbf_to_N;
	psi_to_Pa = 6894.75;
	Pa_to_psi = 1. / psi_to_Pa;
	psi_to_bar = psi_to_Pa / 1e5;
	bar_to_psi = 1. / psi_to_bar;
	inch3_to_m3 = inch_to_m * inch_to_m * inch_to_m;
	m3_to_inch3 = 1. / inch3_to_m3;
}

Units::~Units() {}

//==================================================

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
         Ideal_Gas *_gas) : Units() {

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
	do_plot_runtime = false;
	//ylim_isset = false;
	g = 9.81;
	gas = _gas;

	P_MIN = 0.01e5;
	T_MIN = 273.15 - 100.;
	art_visc = 0.6;


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
			oss << setw(6) << setprecision(3) << v(i) << " ";

		oss << "\n       T (C)      : ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(3) << T(i) - 273.15 << " ";

		oss << "\n       rho (kg/m3): ";
		for (int i = 0; i < Npts; i++)
			oss << setw(6) << setprecision(3) << rho(i) << " ";
	}
	oss << endl;
	return oss.str();
}

// Direct, user-defined initialization of the primitive variables
// with uniform distribution
void LWP::IniUniform(double _vini, double _pini, double _Tini, int _Npts) {
	Npts = _Npts;
	Ini(Npts);
	for (int i = 0; i < Npts; i++) {
		x(i) = i * L / (Npts - 1);
		v(i) = _vini;
		p(i) = _pini;
		T(i) = _Tini;
		rho(i) = gas->GetRho(p(i), T(i));
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
		rho(i) = gas->GetRho(p(i), T(i));
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
		double a = gas->GetSonicVel(T(i));
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
	BCLeft(BC_start_type, BC_start_val1, BC_start_val2);
	BCRight(BC_end_type, BC_end_val1, BC_end_val2);

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
			e = (gas->cv) * Thalf(i);
			Fhalf(i, 0) = rhohalf(i) * vhalf(i) * A;
			Fhalf(i, 1) = ( rhohalf(i) * vhalf(i) * vhalf(i) + phalf(i) ) * A;
			Fhalf(i, 2) = ( rhohalf(i) * vhalf(i) * e + phalf(i) * vhalf(i) ) * A;

			Shalf(i, 0) = 0.;
			Shalf(i, 1) = -A * rhohalf(i) / 2. * lambda / D * vhalf(i) * fabs(vhalf(i));
			Shalf(i, 2) = 0;
		}
	}
	else {
		for (unsigned int i = 0; i < Npts; i++) {
			e = (gas->cv) * T(i);
			U(i, 0) = rho(i) * A;
			U(i, 1) = rho(i) * v(i) * A;
			U(i, 2) = rho(i) * e * A;

			F(i, 0) = rho(i) * v(i) * A;
			F(i, 1) = ( rho(i) * v(i) * v(i) + p(i) ) * A;
			F(i, 2) = ( rho(i) * v(i) * e + p(i) * v(i) ) * A;

			S(i, 0) = 0.;
			S(i, 1) = -A * rho(i) / 2. * lambda / D * v(i) * fabs(v(i));
			S(i, 2) = 0.;
		}
	}
}

void LWP::UnPackU(bool is_half_step) {

	if (is_half_step) {
		for (unsigned int i = 0; i < Npts - 1; i++) {
			rhohalf(i) = Uhalf(i, 0) / A;
			vhalf(i) = Uhalf(i, 1) / Uhalf(i, 0);
			Thalf(i) = Uhalf(i, 2) / Uhalf(i, 0) / (gas->cv);
			phalf(i) = gas->GetP(rhohalf(i), Thalf(i));
		}
	}
	else {
		for (unsigned int i = 0; i < Npts; i++) {
			rhonew(i) = Unew(i, 0) / A;
			vnew(i) = Unew(i, 1) / Unew(i, 0);
			Tnew(i) = Unew(i, 2) / Unew(i, 0) / (gas->cv);
			pnew(i) = gas->GetP(rhonew(i), Tnew(i));
		}
		// artificial viscosity
		double dv;

		double pp = 0.;
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

/*
void LWP::runtime_plot() {
	//printf("\n\n Plotting...");
	int ll = p.size();
	vector<double> xdv(ll), vdv(ll), pdv(ll), qdv(ll), Tdv(ll);
	for (int i = 0; i < ll; i++) {
		xdv.at(i) = x(i);
		vdv.at(i) = v(i);
		pdv.at(i) = p(i) / 1.e5;
		qdv.at(i) = q(i) / 1.e5;
		Tdv.at(i) = T(i)-273.15;
		//printf("\n x=%5.3f m, p=%5.3f bar, v=%5.3f m/s",xdv.at(i),pdv.at(i),vdv.at(i));
	}

	plt::clf();
	plt::subplot(3, 1, 1);
	plt::plot(xdv, pdv,"k",xdv, qdv,"r");
	//plt::xlabel("x (m)");
	plt::ylabel("p (bar)");
	stringstream s;
	s.str("");
	s << "Pipe name: " << name << ", t = " << t << " s";
	plt::title(s.str());
	//if (ylim_isset)
	//	plt::ylim(lim_pbottom, lim_pup);
	plt::grid(true);

	plt::subplot(3, 1, 2);
	plt::plot(xdv, vdv);
	//if (ylim_isset)
	//	plt::ylim(lim_vbottom, lim_vup);
	plt::xlabel(" x (m)");
	plt::ylabel("v (m/s)");
	//plt::title(name);
	plt::grid(true);

	plt::subplot(3, 1, 3);
	plt::plot(xdv, Tdv);
	//if (ylim_isset)
	//	plt::ylim(lim_vbottom, lim_vup);
	plt::xlabel(" x (m)");
	plt::ylabel(" T (C)");
	//plt::title(name);
	plt::grid(true);


	//plt::show();
	plt::pause(0.01);
}
*/

void LWP::BCLeft(string type, double val1, double val2) {

	double beta = GetBetaAtFront(t + dt);
	double kappa = gas->kappa;
	double R = gas->R;
	double cp = gas->cp;

	bool ok = false;

	if (type == "Wall") {
		vnew(0) = 0;
		double a = beta;
		Tnew(0) = a * a / kappa / R;
		double p_old = p(0);
		double r_old = rho(0);
		// Korrigalni kell majd meg a sebesseggel
		pnew(0) = R * Tnew(0) * pow(R * Tnew(0) * pow(r_old, kappa) / p_old, 1. / (kappa - 1));
		rhonew(0) = gas->GetRho(pnew(0), Tnew(0));
		ok = true;
	}

	if (type == "MassFlowIn_and_T") {
		double mp = val1;
		Tnew(0) = val2;
		double  a = sqrt(kappa * R * Tnew(0));
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
	}

	if (type == "StaticPres_and_StaticTemp") {
		pnew(0)  = val1;
		Tnew(0)  = val2;
		double a = sqrt(kappa * R * Tnew(0));
		vnew(0)  = (a - beta) * 2. / (kappa - 1.);
		rhonew(0) = gas->GetRho(pnew(0), Tnew(0));
		ok = true;
	}

	if (type == "Opening") {
		double pt = val1;
		double Tt = val2;
		double at = sqrt(kappa * R * Tt);
		double a0 = sqrt(kappa * R * T(0));

		double aa = pow((kappa - 1.) / 2., 2.) + kappa * R / 2. / cp;
		double bb = beta * (kappa - 1.);
		double cc = beta * beta - kappa * R * Tt;
		double D  = bb * bb - 4.*aa * cc;
		double vv = (-bb + sqrt(D)) / 2. / aa;

		cout << endl << "beta=" << beta << " ?< at=" << at;
		cout << endl << " vv = " << vv << endl;
		if (beta < at) {
			double aa = pow((kappa - 1.) / 2., 2.) + kappa * R / 2. / cp;
			double bb = beta * (kappa - 1.);
			double cc = beta * beta - kappa * R * Tt;
			double D  = bb * bb - 4.*aa * cc;
			if (beta / a0 < (3 - kappa) / 2) { //We have a supersonic flow entering the system
				//Does not use anything from the inside, both temperature and pressure of tank used
				//vnew(0) = at; //sonic flow, we have no Laval nozzle
				//Tnew(0) = Tt - (vnew(0) * vnew(0)) / (2 * gas->cp);
				Tnew(0) = Tt / (1 + kappa * R / 2 / gas->cp);
				vnew(0) = sqrt(kappa * R * Tnew(0));
				pnew(0) = pt * pow(Tnew(0) / Tt, kappa / (kappa - 1.));
			}
			else { //Subsonic flow, all is well
				vnew(0) = (-bb + sqrt(D)) / 2. / aa;
				Tnew(0) = pow(beta + (kappa - 1) / 2.0 * vnew(0), 2.0) / kappa / R;
				pnew(0) = pt * pow(Tnew(0) / Tt, kappa / (kappa - 1.));
			}

			cout << endl << "LWP::BCLeft Opening -> SUBSONIC INLET" << endl;
			printf("\n pt = %5.3f bar, Tt = %5.3f K", pt / 1e5, Tt);
			printf("\n p  = %5.3f bar, T  = %5.3f K, v = %5.3f m/s", pnew(0) / 1e5, Tnew(0), vnew(0));
			printf("\n Tt = %5.3f K ?= %5.3f K", Tt, Tnew(0) + vnew(0)*vnew(0) / 2. / cp);
			printf("\n T/Tt = %5.3f ?= %5.3f = (p/pt)^((k-1)/k)", Tnew(0) / Tt, pow(pnew(0) / pt, (kappa - 1.) / kappa));
			cout << endl;
		}
		else {
			vnew(0) = (-bb + sqrt(D)) / 2. / aa;
			double aa = beta + (kappa - 1.) / 2.*vnew(0);
			Tnew(0) = aa * aa / kappa / R;
			pnew(0) = pt;

			/*double dxB = -v(0) * dt / (1 + (v(1) - v(0)) * dt / dx);
			if (dxB<0){
				cout<<endl<<endl<<"v(0)="<<v(0)<<", v(1)="<<v(1);
				cout<<endl<<"ERROR! LWP::BCLeft:: Opening -> dxB<0, dxB/dx= "<<dxB/dx<<" !!!"<<endl;
				exit(-1);
			}
			if (dxB>dx){
				cout<<endl<<"ERROR! LWP::BCLeft:: Opening -> dxB>dx (dxB/dx)="<<dxB/dx<<" !!!"<<endl;
				exit(-1);
			}
			double pB = p(0) * (1 - dxB / dx) + p(1) * dxB / dx;
			double TB = T(0) * (1 - dxB / dx) + T(1) * dxB / dx;
			double tmp = pow(pB / pt, (kappa - 1) / kappa);
			Tnew(0) = TB / tmp;
			double aa = sqrt(kappa * R * Tnew(0));
			vnew(0) = (beta-aa) * 2. / (kappa - 1.);*/

			cout << endl << "LWP::BCLeft Opening -> SUBSONIC OUTLET" << endl;
			printf("\n p  = %5.3f bar, T  = %5.3f K, v = %5.3f m/s", pnew(0) / 1e5, Tnew(0), vnew(0));
			printf("\n pt = %5.3f bar, Tt = %5.3f K", pt / 1e5, Tt);
			cout << endl;
		}
		rhonew(0) = gas->GetRho(pnew(0), Tnew(0));

//cin.get();

		ok = true;
	}

	if (!ok) {
		cout << endl
		     << "ERROR! LWP::BCLeft(), unknown BC type: " << type << endl
		     << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}
}

void LWP::BCRight(string type, double val1, double val2) {

	double alpha = GetAlphaAtEnd(t + dt);
	double kappa = gas->kappa;
	double R = gas->R;
	double cp = gas->cp;

	/*cout<<endl<<"GetAlphaAtEnd() finished."<<endl;
	cin.get();*/
	bool ok = false;

	if (type == "Wall") {
		vnew(Npts - 1) = 0.;
		double a = alpha;
		Tnew(Npts - 1) = a * a / kappa / R;
		double p_old = p(Npts - 1);
		double r_old = rho(Npts - 1);
		// Korrigalni kell majd meg a sebesseggel
		pnew(Npts - 1) = R * Tnew(Npts - 1) * pow(R * Tnew(Npts - 1) * pow(r_old, kappa) / p_old, 1. / (kappa - 1));
		rhonew(Npts - 1) = gas->GetRho(pnew(Npts - 1), Tnew(Npts - 1));
		ok = true;
	}

	if (type == "StaticPres_and_StaticTemp") {
		pnew(Npts - 1)   = val1;
		Tnew(Npts - 1)   = val2;
		double a = sqrt(kappa * R * Tnew(Npts - 1));
		vnew(Npts - 1)   = (alpha - a) * 2. / (kappa - 1);
		rhonew(Npts - 1) = gas->GetRho(pnew(Npts - 1), Tnew(Npts - 1));
		ok = true;
	}

	if (type == "Opening") {
		double pt = val1;
		double Tt = val2;
		double at = sqrt(kappa * R * Tt);
		double D  = at * at / alpha / alpha * (kappa + 1.) / 2. - 1.;

		if (D < 0) {
			// Limit outlet sonic velocity
			vnew(Npts - 1) = at * sqrt(2. / (kappa + 1));
			// TODO!!!
			Tnew(Npts - 1) = vnew(Npts - 1) * vnew(Npts - 1) / R / kappa;
			pnew(Npts - 1)  = pnew(Npts - 2);

			cout << endl << "LWP::BCRight Opening -> SONIC OUTLET" << endl;
			cout << endl << " WARNING! THIS BC IS NOT READY YET." << endl;
			printf("\n pt = %5.3f bar, Tt = %5.3f K, at = %5.3f m/s", pt / 1e5, Tt, at);
			printf("\n p  = %5.3f bar, T  = %5.3f K, v  = %5.3f m/s", pnew(Npts - 1) / 1e5, Tnew(Npts - 1), vnew(Npts - 1));
			cout << endl;
		}
		else {
			// Subsonic inlet/outlet
			vnew(Npts - 1) = 2. / (kappa + 1.) * (alpha - sqrt(D));

			if (vnew(Npts - 1) < 0) {
				// subsonic inlet
				Tnew(Npts - 1) = Tt - v(Npts - 1) * v(Npts - 1) / 2. / cp;
				pnew(Npts - 1)  = pt * pow(T(Npts - 1) / Tt, kappa / (kappa - 1.));

				cout << endl << "LWP::BCRight Opening -> SUBSONIC INLET" << endl;
				printf("\n pt = %5.3f bar, Tt = %5.3f K", pt / 1e5, Tt);
				printf("\n p  = %5.3f bar, T  = %5.3f K, v = %5.3f m/s", pnew(Npts - 1) / 1e5, Tnew(Npts - 1), vnew(Npts - 1));
				printf("\n Tt = %5.3f K ?= %5.3f K", Tt, Tnew(Npts - 1) + vnew(Npts - 1)*vnew(Npts - 1) / 2. / cp);
				printf("\n T/Tt = %5.3f ?= %5.3f = (p/pt)^((k-1)/k)", Tnew(Npts - 1) / Tt, pow(pnew(Npts - 1) / pt, (kappa - 1.) / kappa));
				cout << endl;
			}
			else {
				// subsonic outlet
				pnew(Npts - 1) = pt;
				Tnew(Npts - 1) = T(Npts - 2);
				double aa = sqrt(kappa * R * Tnew(Npts - 1));
				vnew(Npts - 1) = (alpha - aa) * 2. / (kappa - 1.);

				cout << endl << "LWP::BCRight Opening  -> SUBSONIC OUTLET" << endl;
				printf("\n p  = %5.3f bar, T  = %5.3f K, v = %5.3f m/s", pnew(0) / 1e5, Tnew(0), vnew(0));
				printf("\n pt = %5.3f bar, Tt = %5.3f K", pt / 1e5, Tt);
				cout << endl;
			}
		}
		rhonew(Npts - 1) = gas->GetRho(pnew(Npts - 1), Tnew(Npts - 1));

//cin.get();

		ok = true;
	}

	if (type == "Valve") {
		double Afx = val1; //flow through area dependent on the displacement. This should not be considered outside info, as it is geometry, not status.
		double pback = val2; //the back pressure of the system
		double alpha = GetAlphaAtEnd(t + dt);
		double kappa = gas->kappa;
		double R = gas->R;
		double psic = 1.;//gas->psi_c;
		double aNm1 = sqrt(kappa * R * T(Npts - 2));

		if (Afx == 0) { //It is closed, then we need no other BC
			vnew(Npts - 1) = 0.;
			double a = alpha;
			Tnew(Npts - 1) = a * a / kappa / R;
			// Korrigalni kell majd meg a sebesseggel
			/*pnew(Npts - 1) = R * Tnew(Npts - 1) * pow(R * Tnew(Npts - 1) * pow(r_old, kappa) / p_old, 1. / (kappa - 1));
			rhonew(Npts - 1) = gas->GetRho(pnew(Npts - 1), Tnew(Npts - 1));*/ //Csaba's version
			rhonew(Npts - 1) = pow(Tnew(Npts - 1) / T(Npts - 1), 1.0 / (kappa - 1)) * rho(Npts - 1);
			pnew(Npts - 1) = gas->GetP(rhonew(Npts - 1), Tnew(Npts - 1));
		}
		//Now the open considerations!
		else if (alpha / aNm1 >= (kappa + 1) / 2) { // if it is faster than the speed of sound
			/*double zeta = 0;
			if (v(Npts - 1) - v(Npts - 2) != 0.0) {
				zeta = (pow(M_E, (v(Npts - 2) - v(Npts - 1)) / dx * dt) - 1.0) * v(Npts - 1) * dx / (v(Npts - 2) - v(Npts - 1)); //Location of the fluid leaving
			}
			zeta = abs(zeta);
			double Tzeta = T(Npts - 1) + (T(Npts - 2) - T(Npts - 1)) / dx * zeta;
			double rho_zeta = rho(Npts - 1) + (rho(Npts - 2) - rho(Npts - 1)) / dx * zeta; //velocity at leaving point

			vnew(Npts - 1) = (Afx * psic / sqrt(kappa) * alpha) / (Afx * psic / sqrt(kappa) * (kappa - 1) / 2.0 + A);
			double anew = alpha - (kappa - 1.0) / 2.0 * vnew(Npts-1);
			Tnew(Npts - 1) = anew * anew / kappa / R;
			rhonew(Npts - 1) = rho_zeta * pow(Tnew(Npts - 1) / Tzeta, 1 / (kappa - 1));
			pnew(Npts - 1) = gas->GetP(rhonew(Npts - 1), Tnew(Npts - 1));*/ //Old remaining code
			//In case of supersonic flow, no outside effects are taken into account
			//Same setup as in the case of the normal outflow
			double Unp1[3] = { 0.0, 0.0, 0.0 };
			for (int j = 0; j < 3; j++) { Unp1[j] = U(Npts - 1, j) + dt * ((S(Npts - 1, j) + S(Npts - 2, j)) / 2.0 - (F(Npts - 1, j) - F(Npts - 2, j)) / dx); } //Do the calculation
			//Unpacking
			rhonew(Npts - 1) = Unp1[0] / A;
			vnew(Npts - 1) = Unp1[1] / Unp1[0];
			Tnew(Npts - 1) = Unp1[2] / Unp1[0] / gas->cv;
			pnew(Npts - 1) = gas->GetP(rhonew(Npts - 1), Tnew(Npts - 1));
		}
		else { //subsonic outlet, the valve can not function as an inlet
			double zeta = 0; //Location of exiting fluid droplet
			if (v(Npts - 1) - v(Npts - 2) != 0.0) {
				zeta = (pow(M_E, (v(Npts - 2) - v(Npts - 1)) / dx * dt) - 1.0) * v(Npts - 1) * dx / (v(Npts - 2) - v(Npts - 1)); //Location of the fluid leaving
			}
			zeta = abs(zeta);
			double a0 = sqrt(kappa * R * T(Npts - 1)); double a1 = sqrt(kappa * R * T(Npts - 2)); //Speed of sounds necessary
			double T_zeta = T(Npts - 1) + (T(Npts - 2) - T(Npts - 1)) / dx * zeta; //density at leaving point
			double p_zeta = p(Npts - 1) + (p(Npts - 2) - p(Npts - 1)) / dx * zeta; //pressure at leaving point

			//Starting to solve a three element non-linear system. --With old equations
			/*double epsilon = 1.0;
			int steps = 0;
			double valdisp[3]; double fdisp[3];
			Eigen::VectorXd xnew(3), xold(3), fx(3), deltax(3);
			Eigen::MatrixXd Jac(3, 3);
			xold(0) = v(Npts - 1); xold(1) = p(Npts - 1); xold(2) = rho(Npts - 1); //fill up with last values
			for (int i = 0; i < 3; i++) { valdisp[i] = xold(i); }
			while (/*epsilon > 0.00001 && steps < 255) { //Numerical solution
				fx(0) = xold(0) - Afx / A * sqrt(2 * (xold(1) - pback) / xold(2));
				fx(1) = alpha - sqrt(xold(1) / xold(2) * kappa) - (kappa - 1) / 2 * xold(0);
				fx(2) = xold(1) / p_zeta - pow(xold(2) / rho_zeta,kappa);
				for (int i = 0; i < 3; i++) { fdisp[i] = fx(i); }
				//Now come elements of the jacobian
				Jac(0, 0) = 1.0; //df0/dx = 1.0
				Jac(0, 1) = Afx / A * 2 / xold(2) * 1 / (2 * sqrt(2*(xold(1) - pback) / xold(2)));
				Jac(0, 2) = Afx / A * 1 / (2 * sqrt(2 * (xold(1) - pback) / xold(2))) * -2*(xold(1) - pback) / xold(2) / xold(2);
				Jac(1, 0) = -1.0 * (kappa - 1) / 2.0;
				Jac(1, 1) = sqrt(1 / xold(2) * kappa) * 1 / (2 * sqrt(xold(1)));
				Jac(1, 2) = sqrt(xold(1) * kappa) * -1 / (2 * pow(xold(2), 3 / 2));
				Jac(2, 0) = 0.0;
				Jac(2, 1) = 1 / p_zeta;
				Jac(2, 2) = -kappa * pow(xold(2), kappa - 1) * pow(1 / rho_zeta, kappa);
				//Do the step
				xnew = xold - Jac.inverse() * fx;
				for (int i = 0; i < 3; i++) { deltax(i) = (xnew(i) - xold(i)) / (xold(i) + 1e-5); valdisp[i] = xnew(i); } //create relative errors
				epsilon = deltax.norm();
				steps++;
				xold = xnew;
			}
			//Now unload results
			vnew(Npts - 1) = xold(0);
			pnew(Npts - 1) = xold(1);
			rhonew(Npts - 1) = xold(2);
			Tnew(Npts - 1) = pnew(Npts - 1) / rhonew(Npts - 1) / gas->R;*/
			double anew = (2 * alpha / (kappa - 1)) * 1 / (Afx / A * psic * sqrt(1 / kappa) + 2 / (kappa - 1));
			vnew(Npts - 1) = 2 / (kappa - 1) * (alpha - anew);
			Tnew(Npts - 1) = anew * anew / kappa / R;
			pnew(Npts - 1) = p_zeta * pow(Tnew(Npts - 1) / T_zeta, kappa / (kappa - 1));
			rhonew(Npts - 1) = gas->GetRho(pnew(Npts - 1), Tnew(Npts - 1));
		}
		ok = true;
	}

	if (!ok) {
		cout << endl
		     << "ERROR! LWP::BCRight(), unknown BC type: " << type << endl
		     << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}
}

double LWP::GetBetaAtFront(double t_target) {

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

}

double LWP::GetAlphaAtEnd(double t_target) {
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


	/*cout<<endl<<" aNm1   = :"<<aNm1;
	cout<<endl<<" aN     = :"<<aN;
	cout<<endl<<" vNm1   = :"<<vNm1;
	cout<<endl<<" vN     = :"<<vN;
	cout<<endl<<" dxR/dx = :"<<dxR/dx<<endl;
	cin.get();*/

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

double LWP::Get_dprop(string prop_string) {
	double out;
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
	else if (prop_string == "M_front")
		out = fabs(v(0)) / sqrt(gas->kappa * gas->R * T(0));
	else if (prop_string == "M_back")
		out = fabs(v(Npts - 1)) / sqrt(gas->kappa * gas->R * T(Npts - 1));
	else if (prop_string == "mp_front")
		out = v(0) * gas->GetRho(p(0), T(0)) * A;
	else if (prop_string == "mp_back")
		out = v(Npts - 1) * gas->GetRho(p(Npts - 1), T(Npts - 1)) * A;
	else if (prop_string == "frek")
		out = gas->GetSonicVel(T.mean()) / (2.*L);
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
	else if (prop_string == "ro")
		out = gas->GetRho(p.mean(), T.mean());
	else if (prop_string == "a")
		out = gas->GetSonicVel(T.mean());
	else {
		cout << endl
		     << "ERROR! LWP::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl
		     << endl;
		cout << endl << "Name of pipe: " << name << endl;
		cin.get();
	}
	return out;
}

/*! Save data from the simulation.
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
/*
void LWP::Plot_p() {
	printf("\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), pe(ll), pv(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		pe.at(i) = data.at(i).at(1) / 1.e5;
		pv.at(i) = data.at(i).at(2) / 1.e5;
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, pe, "k", t, pv, "b");
	plt::xlabel("t (s)");
	plt::ylabel("p (bar)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}

void LWP::Plot_v() {
	printf("\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), xe(ll), xv(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		xe.at(i) = data.at(i).at(3);
		xv.at(i) = data.at(i).at(4);
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, xe, "k", t, xv, "b");
	plt::xlabel("t (s)");
	plt::ylabel("v (m/s)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}

void LWP::Plot_T() {
	printf("\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), xe(ll), xv(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		xe.at(i) = data.at(i).at(5) -273.15;
		xv.at(i) = data.at(i).at(6) -273.15;
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, xe, "k", t, xv, "b");
	plt::xlabel("t (s)");
	plt::ylabel("T (C)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}

void LWP::Plot_rho() {
	printf("\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), xe(ll), xv(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		xe.at(i) = data.at(i).at(7) / 1.e5;
		xv.at(i) = data.at(i).at(8) / 1.e5;
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, xe, "k", t, xv, "b");
	plt::xlabel("t (s)");
	plt::ylabel("rho (kg/m3)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}

void LWP::Plot_mp() {
	printf("\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), xe(ll), xv(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		xe.at(i) = data.at(i).at(9) / 1.e5;
		xv.at(i) = data.at(i).at(10) / 1.e5;
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, xe, "k", t, xv, "b");
	plt::xlabel("t (s)");
	plt::ylabel("mp (kg/s)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}

*/

//==================================================


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
         const bool _save_data /**< [in] whether to save data*/ ) : Units() {
	save_data = _save_data;
	name = _name;
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
	do_plot_runtime = false;
	ylim_isset = false;
	g = 9.81;

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
string SCP::Info(bool show_pts) {

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
	Initializes the SCP pipe with a given number of points, an initial uniform (due to incompressibility) speed and a
	pressure at the inlet. (The generated pressure drops with friction.)
	\param vini Initial (uniform) velocity in the pipe
	\param pstart Pressure at the inlet, set up to decrease with friction
	\param Npts_mul Base number of points to create, actually 20 as many created
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
		v(i) = vini; //The valocity is the same as density is constant
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
	Initializes the SCP pipe with as many points to statisfy a given timestep, an initial uniform (due to incompressibility) speed and a
	pressure at the inlet. (The generated pressure drops with friction.)
	\param vini Initial (uniform) velocity in the pipe
	\param pstart Pressure at the inlet, set up to decrease with friction
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
	\param prop_string: A string describing the needed property
	\return The value that was looked up.
*/
double SCP::Get_dprop(string prop_string) {
	double out;
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
	else if (prop_string == "ro")
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

/*! \brief Do a timestep.
	Calculates a timestep in the system dependent on the boundary conditions. (Those are given by the same strings
	and values as the ones in BCLeft and BCRight) Accepted boundary condition types: "Pressure" & "Velocity"
	\param [in] BC_start_type the type of the boundary condition at the beginning of the pipe
	\param [in] BC_start_val The value associated with the boundary condition at the start of the pipe
	\param [in] BC_end_type the type of the boundary condition at the end of the pipe
	\param [in] BC_end_val The value associated with the boundary condition at the end of the pipe
	\sa BCLeft, BCRight
*/
void SCP::Step(
    string BC_start_type, double BC_start_val,
    string BC_end_type, double BC_end_val) {

	double a, b;
	VectorXd pnew = VectorXd::Zero(Npts); //new pressure vector
	VectorXd vnew = VectorXd::Zero(Npts); //new velocity vector
	for (int i = 1; i < Npts - 1; i++) {
		a = (p(i - 1) + roa * v(i - 1)) + dt * roa * Source(i - 1); //
		b = (p(i + 1) - roa * v(i + 1)) - dt * roa * Source(i + 1);
		pnew(i) = (a + b) / 2.;
		vnew(i) = (a - b) / 2. / roa;
	}

	BCLeft(BC_start_type, BC_start_val, pnew(0), vnew(0));

	BCRight(BC_end_type, BC_end_val, pnew(Npts - 1), vnew(Npts - 1));

	for (int i = 0; i < Npts; i++) {
		p(i) = pnew(i);
		v(i) = vnew(i);
	}
	t += dt;

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


/*! \brief Right side bounbdary condition.
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
		double pend = p(Npts - 1);
		double vend = v(Npts - 1);
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
	int Nvars = data.at(0).size();
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
		FILE* pfile, * Tfile, * vfile;
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

/*
void SCP::Plot_p() {
	printf("\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), pe(ll), pv(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		pe.at(i) = data.at(i).at(1) / 1.e5;
		pv.at(i) = data.at(i).at(2) / 1.e5;
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, pe, "k", t, pv, "b");
	plt::xlabel("t (s)");
	plt::ylabel("p (bar)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}

void SCP::Set_ylim(double pmin, double pmax, double vmin, double vmax) {
	lim_pbottom = pmin;
	lim_pup = pmax;
	lim_vbottom = vmin;
	lim_vup = vmax;
	ylim_isset = true;
}

void SCP::runtime_plot(const double t) {
	//printf("\n\n Plotting...");
	int ll = p.size();
	vector<double> xdv(ll), vdv(ll), pdv(ll);
	for (int i = 0; i < ll; i++) {
		xdv.at(i) = x(i);
		vdv.at(i) = v(i);
		pdv.at(i) = p(i) / 1.e5;
		//printf("\n x=%5.3f m, p=%5.3f bar, v=%5.3f m/s",xdv.at(i),pdv.at(i),vdv.at(i));
	}

	plt::clf();
	plt::subplot(2, 1, 1);
	plt::plot(xdv, pdv);
	//plt::xlabel("x (m)");
	plt::ylabel("p (bar)");
	stringstream s;
	s.str("");
	s << "Pipe name: " << name << ", t = " << t << " s";
	plt::title(s.str());
	if (ylim_isset)
		plt::ylim(lim_pbottom, lim_pup);
	plt::grid(true);

	plt::subplot(2, 1, 2);
	plt::plot(xdv, vdv);
	if (ylim_isset)
		plt::ylim(lim_vbottom, lim_vup);
	plt::xlabel(" x (m)");
	plt::ylabel("v (m/s)");
	//plt::title(name);
	plt::grid(true);

	//plt::show();
	plt::pause(0.01);
}

void SCP::Plot_v() {
	printf("\n\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), ve(ll), vv(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		ve.at(i) = data.at(i).at(3);
		vv.at(i) = data.at(i).at(4);
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, ve, "k", t, vv, "b");
	plt::xlabel("t (s)");
	plt::ylabel("v (m/s)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}
void SCP::Plot_distr() {
	printf("\n\n Plotting...");
	int ll = x.size();
	vector<double> xx(ll), pp(ll), vv(ll);
	for (int i = 0; i < ll; i++) {
		xx.at(i) = x(i);
		pp.at(i) = p(i) / 1.e5;
		vv.at(i) = v(i);
		printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::subplot(2,1,1)
	plt::plot(xx, pp);
	//plt::xlabel("t (s)");
	//plt::ylabel("v (m/s)");
	//plt::title(name);
	plt::grid(true);
	plt::show();
}
*/
// ===========================================================

Reservoir::Reservoir(const string _name,
                     const double _vol,
                     const double _a,
                     const bool _save_data) : Units() {
	save_data = _save_data;
	name = _name;
	vol = _vol;
	a = _a;
	fname = name + ".dat";
}

Reservoir::~Reservoir() {}

string Reservoir::GetName() {
	return name;
}


void Reservoir::Ini(double pstart) {
	p = pstart;
	t = 0.;

	tmpvec.push_back(t);
	tmpvec.push_back(p);
	data.clear();
	data.reserve(100);
	data.push_back(tmpvec);
}

void Reservoir::UpdateDimlessPars(double pref, double mp_nevl, double omega) {
	beta = a * a * mp_nevl / vol / omega / pref;
}

void Reservoir::Update(double delta_t, double mp_in, double mp_out) {
	vector<double> xvec;
	vector<vector<double> > yvec;
	VectorXd pars = VectorXd::Zero(2);
	pars(0) = mp_in;
	pars(1) = mp_out;
	VectorXd y_0 = VectorXd::Zero(1);
	y_0(0) = p;

	rk45(t, t + delta_t, y_0, xvec, yvec, pars);

	p = yvec.back().back();
	t += delta_t;

	if (save_data) {
		tmpvec.at(0) = t;
		tmpvec.at(1) = p;

		data.push_back(tmpvec);
	}
}

VectorXd Reservoir::ReservoirODE(const double x, const VectorXd & y, const VectorXd & pars) {
	VectorXd out(1);
	out(0) = a * a / vol * (pars(0) - pars(1));
	// cout << "\n dfdt=" << out(0);
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

	bool is_impact = false;
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

				x.push_back(xx);
				for (int iii = 0; iii < Neq; iii++)
					tmp.at(iii) = yy[iii];

				y.push_back(tmp);
				if (debug)
					printf(" ok.");
				if (hiba < hiba_max / 10.) {
					dx *= 2.;
					// printf("\n\t dx -> dx*2 (dxuj=%g)", dx);
				}
			}
		}
	}
	if (debug) {
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
		fprintf(pFile, "t (s); p (bar)\n");
		for (int i = 0; i < data.size(); i++)
			fprintf(pFile, "%8.6e; %8.6e\n", data.at(i).at(0), data.at(i).at(1) / 1.e5);
		fclose (pFile);

		cout << " done.";
	}
}

/*void Reservoir::Plot() {
	printf("\n\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), p(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		p.at(i) = data.at(i).at(1) / 1.e5;
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, p);
	plt::xlabel("t (s)");
	plt::ylabel("p (bar)");
	plt::title(name);
	plt::grid(true);
	plt::show();
}*/


double Reservoir::Get_dprop(string prop_string) {
	double out;
	if (prop_string == "Vres")
		out = vol;
	else if (prop_string == "a")
		out = a;
	else if (prop_string == "p")
		out = p;
	else if (prop_string == "beta")
		out = beta;
	else {
		cout << endl
		     << "ERROR! Reservoir::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl
		     << endl;
		cout << endl << "Name of pipe: " << name << endl;
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
	oss << "\n        vol : " << vol << " m^3 = " << vol*m3_to_inch3 << " in^3";
	oss << "\n          a : " << a << " m/s";
	oss << "\n       beta : " << beta << " -";

	return oss.str();
}

vector<double> Reservoir::Get_dvprop(string prop_string) {
	int Ntime = data.size();
	int Nvars = data.at(0).size();
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
	else {
		cout << endl
		     << "ERROR! Valve::Get_dvprop(prop_string), unknown input: prop_string=" << prop_string << endl
		     << endl;
		cout << endl << "Name of valve: " << name << endl;
		cin.get();
	}
	return out;
}

// ===========================================================

Valve::Valve(const string _name,
             const double _m,
             const double _k,
             const double _s,
             const double _Dbore,
             const double _Cd,
             const double _ro,
             const double _p_set,
             //const double _mp_nevl,
             const double _xmax,
             const double _r,
             const double _pref,
             const bool _save_data) : Units() {
	save_data = _save_data;
	name = _name;
	m = _m;
	k = _k;
	s = _s;
	Dbore = _Dbore;
	Cd = _Cd;
	ro = _ro;
	A = Dbore * Dbore * M_PI / 4.;
	p_set = _p_set;
	xe = p_set * A / s;
	xmax = _xmax;

	mp_nevl = Get_MassFlowRate_InCompressible(1.1 * p_set, 0, ro, xmax); //Cd * xmax * M_PI * Dbore * sqrt(2.*ro * 1.1 * p_set);

	r = _r;
	pref = _pref;

	ini_done = false;
	fname = name + ".dat";
	x_simul_min = xmax;
	x_simul_max = 0.;

	TINY_VELOCITY = 1e-3;
	TINY_OPENING = 1e-6;

	// AeffCoeffs
	a1 = 0.;
	a2 = 0.;
	Aeffmax = 1.;
}

void Valve::UpdateDimlessPars() {
	xref = A * pref / s;
	delta = xe / xref;
}


Valve::~Valve() {}

string Valve::GetName() {
	return name;
}


void Valve::Ini() {
	x = 1. / 1000.;
	v = 0.;
	t = 0.;

	tmpvec.push_back(t);
	tmpvec.push_back(x);
	tmpvec.push_back(v);
	tmpvec.push_back(0);
	data.clear();
	data.reserve(100);
	data.push_back(tmpvec);
}

void Valve::Ini(double _xstart, double _vstart, double _pstart) {
	x = _xstart;
	v = _vstart;
	p = _pstart;
	t = 0.;

	if (x < TINY_OPENING)
		is_closed = true;
	else
		is_closed = false;

	if (x > (xmax - TINY_OPENING))
		is_fully_open = true;
	else
		is_fully_open = false;

	tmpvec.push_back(t);
	tmpvec.push_back(x);
	tmpvec.push_back(v);
	tmpvec.push_back(0);
	data.clear();
	data.reserve(100);
	data.push_back(tmpvec);
}

void Valve::Ini(double _xstart, double _vstart) {
	x = _xstart;
	v = _vstart;
	p = 0.;
	t = 0.;

	if (x < TINY_OPENING)
		is_closed = true;
	else
		is_closed = false;

	if (x > (xmax - TINY_OPENING))
		is_fully_open = true;
	else
		is_fully_open = false;

	tmpvec.push_back(t);
	tmpvec.push_back(x);
	tmpvec.push_back(v);
	tmpvec.push_back(0);
	data.clear();
	data.reserve(100);
	data.push_back(tmpvec);
}


void Valve::Update(double delta_t, double pv, double pb, bool monitor_min_max) {
	vector<double> xvec;
	vector<vector<double> > yvec;
	VectorXd pars = VectorXd::Zero(2);
	pars(0) = pv;
	pars(1) = pb;
	p = pv;
	VectorXd y_0 = VectorXd::Zero(2);
	y_0(0) = x;
	y_0(1) = v;

	bool debug = false;

	if (debug) {
		printf("\n t = %g s, is_closed = %d, is_fully_open = %d, x=%g mm, v=%g, (pv-pb)/pset=%g %%",
		       t, is_closed, is_fully_open, x * 1000., v, (pv - pb) / p_set * 100);
		cin.get();
	}

	if (is_closed) {
		double Fsum = SumStatForce(0., pars(0), pars(1));

		//if (debug) printf("\n\t Fsum = %g", Fsum);

		x = 0.;
		v = 0.;
		if (Fsum > 0.) {
			// Open only in the next timstep
			is_closed = false;
			if (debug) {
				printf("\n\t t=%5.3e, x=0, Fsum = %g -> Valve opens!", t, Fsum);
				cin.get();
			}
		}
		else
			// Valve stays closed
			is_closed = true;

	}
	else if (is_fully_open) {
		double Fsum = SumStatForce(xmax, pars(0), pars(1));

		x = xmax;
		v = 0.;
		if (Fsum < 0.) {
			// Start closing only in the next timstep
			is_fully_open = false;
			if (debug) {
				printf("\n\t t=%5.3e, x=xmax, Fsum = %g -> Valve starts closing!", t, Fsum);
				cin.get();
			}
		}
		else
			// Valve stays fully open
			is_fully_open = true;

	}
	else {
		// Valve is open, take a step
		rk45(t, t + delta_t, y_0, xvec, yvec, pars);

		x = yvec.back().at(0);
		v = yvec.back().at(1);

		// Handle hitting the seat
		if (x < 0.) {
			x = 0.;
			if (debug) printf("\n x = %g < 0", x);
			if (fabs(v) < TINY_VELOCITY) {

				is_closed = true;
				double Fsum = SumStatForce(0, pars(0), pars(1));

				if (debug) {
					printf("\n\t v=%g -> valve shuts ! (Fsum = %g)",
					       v, Fsum);
					cin.get();
				}
				v = 0.;
			}
			else {
				double vnew = r * fabs(v);
				if (debug) {
					printf("\n\t -> valve bounces! ");
					printf("\n\t v+ = %g", vnew);
					//	cin.get();
				}
				v = vnew;
			}
		}

		// Handle hitting the upper stopper
		if (x > xmax) {
			x = xmax;
			if (debug) printf("\n x/xmax = %g > 1", x / xmax);
			if (fabs(v) < TINY_VELOCITY) {

				is_fully_open = true;
				double Fsum = SumStatForce(xmax, pars(0), pars(1));

				if (debug) {
					printf("\n\t v=%g -> valve fully opens ! (Fsum = %g)",
					       v, Fsum);
					cin.get();
				}
				v = 0.;
			}
			else {
				double vnew = -r * fabs(v);
				if (debug) {
					printf("\n\t -> valve bounces! ");
					printf("\n\t v+ = %g", vnew);
					//	cin.get();
				}
				v = vnew;
			}
		}
	}


	t += delta_t;
	if (save_data) {
		// Save only endpoint
		tmpvec.at(0) = t;
		tmpvec.at(1) = x;
		tmpvec.at(2) = v;
		tmpvec.at(3) = Get_MassFlowRate_InCompressible(pv, pb, ro, x); // Cd * x * M_PI * Dbore * sqrt(2 * ro * (pv - pb));
		data.push_back(tmpvec);
		// Save all data
		//for (int i=0; i<yvec.size(); i++){
		//		tmpvec.at(0) = xvec.at(i);
		//		tmpvec.at(1) = yvec.at(i).at(0);
		//		tmpvec.at(2) = yvec.at(i).at(1);
		//
		//	}
		//data.push_back(tmpvec);
	}

	if (monitor_min_max) {
		if (x < x_simul_min)
			x_simul_min = x;
		if (x > x_simul_max)
			x_simul_max = x;
	}
}
/* \brief Summ the static forces of the system
	Calcualte the static forces acting on the valve (pressure and compression dependent)
	The equation is $$ A_{eff} (x) \left(  p_v - p_b \right) + s (x+x_e)
	where s is the sriong compressuion, p a re the pressure, x the displacement, and xe a pre-compression

	\param x [in] the displacement
	\param pv [in] the pressure in the upstream sid of tghe valve
	\param pb [in] the back pressure of the valve
	\return the static forces (displacement dependent one) of the valve
*/
double Valve::SumStatForce(const double x, const double pv, const double pb) {

	return (Aeff(x) * (pv - pb)  - s * (x + xe));
}

/*! \brief set teh coefficients forthe function describing the effective area
	Set two elements of a second order function (and a maximum value) for the relative multiplicator
	of the effective area.

	\param _a1 [in] first order parameter
	\param _a2 [in] second order parameter
	\param _Aeffmax [in] maximum area relative to valve cross section
*/
void Valve::SetAeffCoeffs(const double _a1, const double _a2, const double _AeffMax) {
	a1 = _a1;
	a2 = _a2;
	Aeffmax = _AeffMax;
}

/*! \brief Calculate the effective area
	Calcualte the effective are based on the effective area coefficients set previously. Uses a second order correlation
	\param x [in] the valve disk displacement
	\return The effective area
	\sa SetACoeffs
*/
double Valve::Aeff(const double x) {
	//double a1=1.0;
	//double a2=0.0;
	//double des_Aeffmax = s * (xmax + xe) / (1.1 * p_set + p0 - 1e5) / A; // 2J3: 1.55
	double des_Aeffmax = s * (xmax + xe) / (1.1 * p_set) / A; // 2J3: 1.55
	//cout<<endl<<"\n Aeffmax required = "<<des_Aeffmax<<endl;
	//cin.get();
	double Aeff = (1. + a1 * (x / xmax) + a2 * pow(x / xmax, 2.));
	if (Aeff > Aeffmax)
		Aeff = Aeffmax;
	return Aeff * A;
}
/* \brief The differential equation describing the
	An equation that returns the zero solved equations of the valve ODE
	Equations:
	$$o(0) = v$$
	$$o(1) = \left( A_{eff}(p_v - p_b) + s(xe+x) - k v \right) \cdot \frac{1}{m} $$
	The first equation is simply returning the speed of the valve stem, and the second is the expression
	of the second derivative of x

	\param x [in] the time variable in the equation
	\param y [in] the vector of the different variables y=[x,dx/dt]
	\param pars [in] the parameters for the differential equation
	\return the first and second derivative of the displacement of the valve
*/
VectorXd Valve::ValveODE(const double x, const VectorXd & y, const VectorXd & pars) {
	// Warning: x is time!
	VectorXd out(2);
	out(0) = y(1);
	//out(1) = (A * pars(0) - k * y(1) - s * (y(0) + xe)) / m;
	double x_v = y(0);
	double v_v = y(1);
	double pv = pars(0);
	double pb = pars(1);
	out(1) = (SumStatForce(x_v, pv, pb) - k * v_v) / m;

	/*printf("\n\t s     = %g",s);
	printf("\n\t xe    = %g",xe);
	printf("\n\t Apres = %g",A);
	printf("\n\t x     = %g",y(0));
	printf("\n\t v     = %g",y(1));
	printf("\n\t p     = %g",pars(0));
	 cout << "\n dfdt=" << out.transpose();
	 cin.get();*/
	return out;
}

/* \brief RK45 method for solving the valve differential equation
	Uses a Runge-Kutta-Fehlberg (RK-45) adaptable timestep method to move the differential equation of the
	valve further in time.
	!!! x variables denote time !!!

	\param x_0 [in] the current time of the system
	\param x_max [in] the time amount we actually want to go (maximum we can move forward in time)
	\param y_0 [in] initial displacement and its derivatives
	\param x [out] vector of the times the method stepped to while working
	\param y [out]
	\param pars [in] the parameters of the valve
*/

void Valve::rk45(
    double x_0, double x_max, VectorXd y_0,
    vector<double> &x, vector<vector<double> > &y,
    const VectorXd & pars ) {

	bool debug = false; //allows for turning debug mode on
	bool is_rk = true; //If true a RK-45 method is used, otherwise the Heun method.

	int Neq = y_0.size(); //the number of equations is the size of the equation vector
	//clear the x and y vectors
	x.clear();
	y.clear();
	x.push_back(x_0); //insert x0 to the first place
	std::vector<double> tmp(Neq);
	for (int i = 0; i < Neq; i++)
		tmp.at(i) = y_0[i];
	y.push_back(tmp);

	//Set up the Runge-Kutta-Fehlberg constantns
	double c2 = 1. / 4., a21 = 1. / 4.;
	double c3 = 3. / 8., a31 = 3. / 32., a32 = 9. / 32.;
	double c4 = 12. / 13., a41 = 1932. / 2197., a42 = -7200. / 2197., a43 = 7296. / 2197.;
	double c5 = 1., a51 = 439. / 216., a52 = -8., a53 = 3680. / 513., a54 = -845. / 4104;
	double c6 = 1. / 2., a61 = -8. / 27., a62 = 2., a63 = -3544. / 2565., a64 = 1859. / 4104., a65 = -11. / 40.;
	double b1 = 16. / 135., b2 = 0., b3 = 6656. / 12825., b4 = 28561. / 56430., b5 = -9. / 50., b6 = 2. / 55.; //Fifth order approximation constants
	double d1 = 25. / 216., d2 = 0., d3 = 1408. / 2565., d4 = 2197. / 4104., d5 = -1. / 5., d6 = 0.; //Fourth order approximation constants

	bool is_impact = false;
	bool last_step = false;
	double hiba_max = 1e-9, dx_min = 1.e-10;

	double hiba = 2.0 * hiba_max;
	double xx = x_0; //the time during the RK simulation
	double dx = (x_max - x_0) / 100.; //size of the timestep in the RK method

	VectorXd yy = y_0;
	bool retake_step = true;

	if (debug) {
		printf("\n\nIntegrating...\n\t t=%7.5e", xx);
		for (int iii = 0; iii < Neq; iii++)
			printf(", y(%d)=%7.5e", iii, yy[iii]); //print every value during debug mode
	}
	while (!last_step) {
		retake_step = true;

		VectorXd y1(Neq), y2(Neq), k1(Neq), k2(Neq), k3(Neq), k4(Neq), k5(Neq), k6(Neq); //create vectors for the components
		hiba = 2.0 * hiba_max;

		if (debug) //wait for an enter if dubuggung -- alows readability
			cin.get();

		while (retake_step) {
			// Is this the last step? (To reach the x_max outer timestep I want to get to)
			if ((dx > 0) && (xx + dx > x_max)) { //we will reach the end in one step
				dx = x_max - xx;
				retake_step = false;
				last_step = true;
			}
			if ((dx < 0) && (xx + dx < x_max)) { //we have passed the maximum time
				dx = (x_max - xx);
				retake_step = false;
				last_step = true;
			}

			if (is_rk) { //This swithc is affected by changing the source code
				// RK45
				k1 = ValveODE(xx, yy, pars);
				k2 = ValveODE(xx + c2 * dx, yy + a21 * dx * k1, pars);
				k3 = ValveODE(xx + c3 * dx, yy + a31 * dx * k1 + a32 * dx * k2, pars);
				k4 = ValveODE(xx + c4 * dx, yy + a41 * dx * k1 + a42 * dx * k2 + a43 * dx * k3, pars);
				k5 = ValveODE(xx + c5 * dx, yy + a51 * dx * k1 + a52 * dx * k2 + a53 * dx * k3 + a54 * dx * k4, pars);
				k6 = ValveODE(xx + c6 * dx, yy + a61 * dx * k1 + a62 * dx * k2 + a63 * dx * k3 + a64 * dx * k4 + a65 * dx * k5, pars);
				y1 = yy + dx * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6); //fifth order approximation
				y2 = yy + dx * (d1 * k1 + d2 * k2 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6); //fourth order approximation
			}
			else {
				// Heun
				k1 = ValveODE(xx, yy, pars);
				y1 = yy + dx * k1;
				y2 = yy + dx / 2 * (k1 + ValveODE(xx + dx, y1, pars));
			}

			// compare results
			hiba = (y1 - y2).norm() / y1.norm();
			if (debug) { //Display the errors in debug mode
				printf("\n\nIntegrating...\n\t t=%7.5e, hiba=%7.5e", xx, hiba);
				for (int iii = 0; iii < Neq; iii++)
					printf(", y(%d)=%7.5e", iii, yy[iii]);
			}


			// Lepes beallitasa
			if (hiba > hiba_max) { //If the accuracy is not sufficient
				dx /= 2.; // a smaller timestep
				retake_step = true; // is used for retaking
				// printf("\n\t dx -> dx/2 (dx=%g)", dx);
				if (fabs(dx) < dx_min) { //Dispaly error message in case the necessary timestep is smaller than the minimum allowed one
					last_step = true;
					cout << endl << endl << "!!!! PANIC dxuj<dx_min !!!!" << endl << endl;
					exit(-1);
				}
			}
			else { //acceptable
				retake_step = false;
				yy = y1; //use the fifth order approximation
				xx += dx; //timestep

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

void Valve::List_data() {
	cout << endl << endl << "Valve " << name << endl;
	for (int i = 0; i < data.size(); i++)
		printf("t=%5.3f s, x=%5.3f mm, v=%5.3f m/s, mp=%5.3f kg/s\n",
		       data.at(i).at(0),
		       data.at(i).at(1) * 1000.,
		       data.at(i).at(2),
		       data.at(i).at(3));
	cout << "Number of data points: " << data.size() << endl;
}


void Valve::Save_data() {
//	char fname [50];
//	sprintf (fname, "%s.dat", name.c_str());
	if (!save_data) {
		cout << endl << "WARNING! Valve: " << name;
		cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
	}
	else {
		cout << endl << "Saving to " << fname.c_str() << " ... ";

		FILE * pFile;
		pFile = fopen (fname.c_str(), "w");
		fprintf(pFile, "t (s); x (mm); v (m/s); mdot (kg/s)\n");
		for (int i = 0; i < data.size(); i++)
			fprintf(pFile, "%8.6e; %8.6e; %+8.6e; %+8.6e\n",
			        data.at(i).at(0),
			        data.at(i).at(1) * 1000.,
			        data.at(i).at(2),
			        data.at(i).at(3));
		fclose (pFile);
		cout << " done.";
	}
}

/*void Valve::Plot() {
	printf("\n\n Plotting...");
	int ll = data.size();
	vector<double> t(ll), x(ll);
	for (int i = 0; i < ll; i++) {
		t.at(i) = data.at(i).at(0);
		x.at(i) = data.at(i).at(1) * 1000.;
		//printf("\n #%3d: t=%5.3f s, x=%5.3f mm",i,t.at(i),x.at(i));
	}
	plt::plot(t, x);
	plt::xlabel("t (s)");
	plt::ylabel("x (mm)");
	plt::grid(true);
	plt::show();
}*/

double Valve::Get_MassFlowRate_Compressible_Choked(double p_upstream, double rho_upstream, double kappa, double xx) {
	double A_flowthrough = Dbore * M_PI * xx;
	double exponent = (kappa + 1.) / (kappa - 1.);
	double tmp = kappa * pow(2. / (kappa + 1.), exponent);
	return Cd * A_flowthrough * sqrt(rho_upstream * p_upstream * tmp);
}

double Valve::Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho, double xx) {
	double A_flowthrough = Dbore * M_PI * xx;
	double dp = p_upstream - p_downstream;
	double DP_MIN = 0.1e5;
	double mp;
	if (fabs(dp) < DP_MIN) {
		double mp_lin = Cd * A_flowthrough * sqrt(2. * ro * DP_MIN);
		mp = dp / DP_MIN * mp_lin;
	}
	else
		mp = Cd * A_flowthrough * sqrt(2. * ro * fabs(dp)) * signum(dp);

	return mp;
}

vector<double> Valve::Get_dvprop(string prop_string) {
	int Ntime = data.size();
	int Nvars = data.at(0).size();
	//cout<<endl<<"Ntime="<<Ntime<<", Nvars="<<Nvars<<endl;
	vector<double> out(Ntime);
	if (prop_string == "t")
		for (unsigned int i = 0; i < Ntime; i++)
			out.at(i) = data.at(i).at(0);
	else if (prop_string == "x")
		for (unsigned int i = 0; i < Ntime; i++)
			out.at(i) = data.at(i).at(1);
	else if (prop_string == "x_mm")
		for (unsigned int i = 0; i < Ntime; i++)
			out.at(i) = data.at(i).at(1) * 1000.;
	else if (prop_string == "x_rel_percent")
		for (unsigned int i = 0; i < Ntime; i++)
			out.at(i) = data.at(i).at(1) / xmax * 100.;
	else if (prop_string == "v")
		for (unsigned int i = 0; i < Ntime; i++)
			out.at(i) = data.at(i).at(2);
	else {
		cout << endl
		     << "ERROR! Valve::Get_dvprop(prop_string), unknown input: prop_string=" << prop_string << endl
		     << endl;
		cout << endl << "Name of valve: " << name << endl;
		cin.get();
	}
	return out;
}

double Valve::Get_dprop(string prop_string) {
	double out;
	if (prop_string == "Aft")
		out = x * M_PI * Dbore;
	else if (prop_string == "Cd")
		out = Cd;
	else if (prop_string == "Dbore")
		out = Dbore;
	else if (prop_string == "Dbore_inch")
		out = Dbore * m_to_inch;
	else if (prop_string == "x")
		out = x;
	else if (prop_string == "x_mm")
		out = x * 1000.;
	else if (prop_string == "x_rel_percent")
		out = x / xmax * 100.;
	else if (prop_string == "xmax")
		out = xmax;
	else if (prop_string == "v")
		out = v;
	else if (prop_string == "c1")
		out = Dbore * M_PI;
	else if (prop_string == "A_flowthrough")
		out = Dbore * M_PI * x;
	else if (prop_string == "Frugo")
		out = s * (x + xe);
	else if (prop_string == "A")
		out = A;
	else if (prop_string == "p_set")
		out = p_set;
	else if (prop_string == "p_set_barg")
		out = p_set / 1.e5;
	else if (prop_string == "p_set_psig")
		out = p_set / 1.e5 * bar_to_psi;
	else if (prop_string == "pref")
		out = pref;
	else if (prop_string == "pref_bar")
		out = pref / 1.e5;
	else if (prop_string == "frek")
		out = sqrt(s / m) / 2. / M_PI;
	else if (prop_string == "omega")
		out = sqrt(s / m);
	else if (prop_string == "xe")
		out = xe;
	else if (prop_string == "s")
		out = s;
	else if (prop_string == "s_kN_p_m")
		out = s / 1000.;
	else if (prop_string == "s_lbf_p_inch")
		out = s * N_to_lbf / m_to_inch;
	else if ((prop_string == "mp_nevl") || (prop_string == "mp_nom") || (prop_string == "mp_cap"))
		out = mp_nevl;
	else if (prop_string == "mp_nevl_lbm_p_s")
		out = mp_nevl * kg_to_lbm;
	else if (prop_string == "m")
		out = m;
	else if (prop_string == "m_lbm")
		out = m * kg_to_lbm;
	else if (prop_string == "mp")
		out = Cd * Dbore * M_PI * x * sqrt(2 * ro * p);
	else if (prop_string == "k")
		out = k;
	else if (prop_string == "kcrit")
		out = 2.*sqrt(s * m);
	else if (prop_string == "k_p_kcrit_percent")
		out = k / (2.*sqrt(s * m)) * 100.;
	else if (prop_string == "xref")
		out = xref;
	else if (prop_string == "delta")
		out = delta;
	else if (prop_string == "x_simul_min")
		out = x_simul_min;
	else if (prop_string == "x_simul_max")
		out = x_simul_max;
	else {
		cout << endl
		     << "ERROR! Valve::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl
		     << endl;
		cout << endl << "Name of valve: " << name << endl;
		cin.get();
	}
	return out;
}

void Valve::Set_dprop(string prop_string, double val) {
	if (prop_string == "x")
		x = val;
	else if (prop_string == "v")
		v = val;
	else {
		cout << endl
		     << "HIBA! Valve::Set_dprop(prop_string), ismeretlen bemenet: prop_string:" << prop_string << endl
		     << endl;
	}
}

string Valve::Info() {

	if (!ini_done)
		Ini();

	std::ostringstream oss;
	oss << "\n\n Valve name : " << name;
	oss << "\n================================";
	oss << "\n       mass : " << m << " kg = " << m*kg_to_lbm << "lbm";
	oss << "\n          k : " << k << " N/(m/s)";
	oss << "\n          s : " << s << " N/m = " << s*N_to_lbf / m_to_inch << "lbf/in";
	oss << "\n      omega : " << sqrt(s / m) << " rad/s";
	oss << "\n          f : " << sqrt(s / m) / 2. / M_PI << " Hz";
	oss << "\n      Dbore : " << Dbore * 1000. << " mm = " << Dbore*m_to_inch << " in";
	oss << "\n         Cd : " << Cd;
	oss << "\n       xmax : " << xmax * 1000. << " mm = " << xmax*m_to_inch << " in";
	oss << "\n         ro : " << ro << " kg/m3";
	oss << "\n      p_set : " << p_set / 1.e5 << " barg = " << p_set / 1.e5*bar_to_psi << " psig";
	oss << "\n         xe : " << xe * 1000. << " mm = " << xe*m_to_inch << " in";
	oss << "\n     F(x=0) : " << s*xe << " N ( = s*xe )";
	oss << "\n   F(p_set) : " << p_set*A << " N ( = p_set * A )";
	oss << endl;
	/*oss << "\n          x : " << x * 1000. << " mm";
	oss << "\n          p : " << p / 1.e5 << " bar";
	oss << "\n       F(x) : " << s*(x + xe) << " N";
	oss << "\n       F(p) : " << p*A << " N";
	oss << endl;*/
	oss << "\n   capacity : " << mp_nevl << " kg/s = " << mp_nevl*kg_to_lbm << "lbm/s";
	oss << "\n       pref : " << pref / 1.e5 << " bar = " << pref*Pa_to_psi << " barg";
	oss << "\n       xref : " << xref * 1000. << " mm";
	oss << "\n  xref/xmax : " << xref / xmax << " -";
	oss << "\n      delta : " << delta << " -";

	return oss.str();
}