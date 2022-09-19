#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

#include "Valve.h"
#include "my_tools.h"
//#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"

//! Class for modeling valves
/*!
  The class contains data for direct spring operated pressure reliev valve
  */

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
	xmax_unrestricted = _xmax;
	RT = 0.;
	xmax = (1.-RT/100.)*xmax_unrestricted;
	
	pref = _pref;
	is_Gas = false;
	mp_nevl = Get_MassFlowRate(1.1 * p_set+pref, 293., pref, 293., xmax_unrestricted); //Cd * xmax * M_PI * Dbore * sqrt(2.*ro * 1.1 * p_set);
	r = _r;
	//pref = _pref;

	ini_done = false;
	fname = name + ".dat";
	x_simul_min = xmax;
	x_simul_max = 0.;

	TINY_VELOCITY = 1e-3;
	TINY_OPENING = 1e-6;

	//is_Gas = false;

	// AeffCoeffs
	a1 = 0.;
	a2 = 0.;
	Aeffmax = 1.;
}

Valve::Valve(const string _name,
		const double _m,
		const double _k,
		const double _s,
		const double _Dbore,
		const double _Cd,
		Gas* _gas,
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
	//ro = _ro;
	A = Dbore * Dbore * M_PI / 4.;
	p_set = _p_set;
	xe = p_set * A / s;
	xmax_unrestricted = _xmax;
	RT=0.;
	xmax = (1.-RT/100.)*xmax_unrestricted;

	is_Gas = true;

	gas = _gas;//Ideal_Gas(_gas);

	pref = _pref;
	//double Temp_at_NominalMassFlowRate = 293.;
	//mp_nevl = Get_MassFlowRate_Compressible_Choked(1.1 * p_set + 1.e5, Temp_at_NominalMassFlowRate, xmax);
	//mp_nevl = Get_MassFlowRate(1.1 * p_set + 1.e5, Temp_at_NominalMassFlowRate, xmax);
	mp_nevl = Get_MassFlowRate(1.1 * p_set + pref,293., pref, 293.,xmax_unrestricted);

	r = _r;
	//pref = _pref;

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
		tmpvec.at(3) = Get_MassFlowRate(pv, 293., pb, 293., x); // Cd * x * M_PI * Dbore * sqrt(2 * ro * (pv - pb));
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
	//printf("\n \t Valve: Aeff*dp=%5.3f, s*(x+xe)=%5.3f N, +++++ %5.3e, %5.3e, %5.3e",Aeff(x) * (pv - pb),s * (x + xe),pv,pb,x);
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
	//	double des_Aeffmax = s * (xmax + xe) / (1.1 * p_set) / A; // 2J3: 1.55
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

	//bool is_impact = false;
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


double Valve::Get_MassFlowRate(double p_upstream, double T_upstream, double p_downstream, double T_downstream, double xx){
	if (is_Gas)
		return Get_MassFlowRate_Compressible(p_upstream, T_upstream, p_downstream, T_downstream, xx);
	else
		return Get_MassFlowRate_InCompressible(p_upstream, p_downstream, ro, xx);
}

double Valve::Get_MassFlowRate_Compressible(double p_upstream, double T_upstream, double p_downstream, double T_downstream, double xx) {
	/*double dir_mul = 1.;
	  if (p_downstream > p_upstream) {
	  double p_tmp = p_upstream;
	  double T_tmp = T_upstream;
	  p_upstream = p_downstream;
	  T_upstream = T_downstream;
	  p_downstream = p_tmp;
	  T_downstream = T_tmp;
	  dir_mul = -1.;
	  }

	  double mp = 0.;
	  if (p_downstream / p_upstream < gas->Get_pcrit())
	  mp = Get_MassFlowRate_Compressible_Choked(p_upstream, T_upstream, xx) ;
	  else
	  mp = Get_MassFlowRate_Compressible_UnChoked(p_upstream, T_upstream, p_downstream, xx);

	  mp *= dir_mul;*/
	double A_flowthrough = Dbore * M_PI * xx;
	double massflux = gas->Get_MassFlux(p_upstream, T_upstream, p_downstream, T_downstream);
	return Cd*A_flowthrough*massflux;
}


double Valve::Get_MassFlowRate_Compressible_Choked(double p_upstream, double T_upstream, double xx) {
	double rho_upstream = gas->Get_rho(p_upstream, T_upstream);
	double kappa = gas->Get_kappa_Tp();

	double A_flowthrough = Dbore * M_PI * xx;
	double exponent = (kappa + 1.) / (kappa - 1.);
	double tmp = kappa * pow(2. / (kappa + 1.), exponent);
	return Cd * A_flowthrough * sqrt(rho_upstream * p_upstream * tmp);
}

double Valve::Get_MassFlowRate_Compressible_UnChoked(double p_upstream, double T_upstream, double p_downstream, double xx) {
	double rho_upstream = gas->Get_rho(p_upstream, T_upstream);
	double kappa = gas->Get_kappa_Tp();

	double A_flowthrough = Dbore * M_PI * xx;
	double exp1 = 1. / kappa;
	double exp2 = (kappa + 1.) / kappa;
	double tmp = 2.*rho_upstream * p_upstream * kappa / (kappa - 1);
	double pr = p_downstream / p_upstream;
	return Cd * A_flowthrough * sqrt(tmp * (pow(pr, exp1) - pow(pr, exp2)));
}

double Valve::Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho, double xx) {
	double A_flowthrough = Dbore * M_PI * xx;
	double dp = p_upstream - p_downstream;
	double DP_MIN = 0.01e5;
	double G0, mp;
	if (fabs(dp) < DP_MIN) {
		G0=sqrt(2. * ro * DP_MIN)*dp / DP_MIN;
		mp = Cd * A_flowthrough * G0;		
	}
	else{
		G0=sqrt(2. * ro * fabs(dp)) * signum(dp);
		mp = Cd * A_flowthrough * G0;
	}

	return mp;
}


vector<double> Valve::Get_dvprop(string prop_string) {
	int Ntime = data.size();
	//	int Nvars = data.at(0).size();
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
	double out = 0.0;
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
	else if (prop_string == "RT")
		out=RT;
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
	else if (prop_string == "p_set_bar")
		out = p_set / 1.e5;
	else if (prop_string == "p_set_psi_abs")
		out = p_set / 1.e5 * bar_to_psi;
	else if (prop_string == "pref")
		out = pref;
	else if (prop_string == "pref_bar")
		out = pref / 1.e5;
	else if ((prop_string == "frek") || (prop_string == "freq"))
		out = sqrt(s / m) / 2. / M_PI;
	else if (prop_string =="dt")
		out = 2*M_PI/sqrt(s/m)/100.;
	else if (prop_string == "omega")
		out = sqrt(s / m);
	else if ((prop_string == "xe") || (prop_string == "x0"))
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
	/*	else if (prop_string == "mp"){
		if (is_Gas){
		}
		else{
		out = Get_MassFlowRate_InCompressible(p,);
	//out = Cd * Dbore * M_PI * x * sqrt(2 * ro * p);
	}*/
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
	else if (prop_string == "RT"){
		RT=val;
		xmax = xmax_unrestricted*(1.-RT/100.);
	}
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
	oss << "\n        mass : " << m << " kg = " << m*kg_to_lbm << "lbm";
	oss << "\n           k : " << k << " N/(m/s)";
	oss << "\n           s : " << s << " N/m = " << s*N_to_lbf / m_to_inch << "lbf/in";
	oss << "\n       omega : " << sqrt(s / m) << " rad/s";
	oss << "\n           f : " << sqrt(s / m) / 2. / M_PI << " Hz";
	oss << "\n       Dbore : " << Dbore * 1000. << " mm = " << Dbore*m_to_inch << " in (nozzle diameter)";
	oss << "\n xmax unrest.: " << xmax_unrestricted * 1000. << " mm = " << xmax_unrestricted*m_to_inch << " in";
	oss << "\n          RT : " << RT << "% -> xmax = "<<xmax*1000<<" mm";
	oss << "\n A orif. area: " << A*1.e6<<" mm^2 = Dbore^2*pi/4";
	oss << "\n      A@xmax : " << Dbore*M_PI*xmax*1.e6<<" mm^2 = Dbore*pi*xmax";
	oss << "\n          Cd : " << Cd;
	if (is_Gas)	
		oss << "\n          ro : " << gas->Get_rho(1.e5, 273+16) << " kg/m3 @ 1 bar, 16C";
	else
		oss << "\n          ro : " << ro << " kg/m3";
	oss << "\n       p_set : " << p_set / 1.e5 << " barg = " << p_set / 1.e5*bar_to_psi << " psig";
	oss << "\n          xe : " << xe * 1000. << " mm = " << xe*m_to_inch << " in";
	oss << "\n      F(x=0) : " << s*xe << " N ( = s*xe )";
	oss << "\n    F(p_set) : " << p_set*A << " N ( = p_set * A )";
	oss << "\n   F(x=xmax) : " << s*(xe+xmax_unrestricted) << " N ( = s*(xe+xmax_unrestricted )";
	oss << "\nF(1.1*p_set) : " << 1.1*p_set*A << " N ( = 1.1*p_set * A )";
	//oss << "\n "
	oss << endl;

	UpdateDimlessPars();

	oss << "\n   capacity : " << mp_nevl << " kg/s = " << mp_nevl*kg_to_lbm << "lbm/s";
	double mp_RT=Get_MassFlowRate(1.e5+p_set*1.1,293,1e5,293,xmax);
	oss << "\n  cap. @ RT : " << mp_RT << " kg/s = " << mp_RT*kg_to_lbm << "lbm/s, @ 110% of pset and xmax_restricted (pback=1.e5, Tin=293K)";
	oss << "\n sum(F) @ RT: " << SumStatForce(xmax,1.e5+1.1*p_set, 1.e5) <<" N, restoring force @ xmax_restricted, 1.1*pset"; 
	oss << "\n       pref : " << pref / 1.e5 << " bar = " << pref*Pa_to_psi << " barg";
	oss << "\n       xref : " << xref * 1000. << " mm";
	oss << "\n  xref/xmax : " << xref / xmax << " -";
	oss << "\n      delta : " << delta << " -";

	return oss.str();
}


/*void Valve::Plot(bool show, bool save) {
  namespace plt = matplotlibcpp;
  int n = data.size();
  std::vector<double> t(n), x(n), v(n), mp(n);
  for (int i = 0; i < n; ++i) {
  t.at(i) = data.at(i).at(0);
  x.at(i) = data.at(i).at(1) * 1000.;
  v.at(i) = data.at(i).at(2);
  mp.at(i) = data.at(i).at(3);
  }

// Plot line from given x and y data. Color is selected automatically.
plt::subplot(3, 1, 1);
plt::plot(t, x);
plt::ylabel("x, mm");

string title = "Valve: " + name;
plt::title(fname);

plt::subplot(3, 1, 2);
plt::plot(t, v);
plt::ylabel("v, m/s");

plt::subplot(3, 1, 3);
plt::plot(t, mp);
plt::ylabel("mp, kg/s");
plt::xlabel("t, s");

if (show)
plt::show();

// Enable legend.
//plt::legend();
// save figure
if (save) {
string fname = name + ".png";
plt::save(fname);
}
}*/
