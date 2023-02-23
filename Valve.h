#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "Units.h"
#include "IdealGas.h"
//#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"

using namespace std;
using namespace Eigen;

/*! \brief A class for handling both liquid and fluid valves
  A class for storing data regarding a valve and allowing for a solution of the differential equation
  describing the movement of the spring, mass damper system.
  */
class Valve: public Units
{
	private:
		bool save_data; //!< whether to save the data of the valve or not
		string name; //!< name of the valve
		double m; /*!< moving mass of the valve */
		double s; /*!< spring rate in valve*/ 
		double k;  /*damping factor in valve*/
		double A; /*effective area of the valve */
		double Cd; /*!<discharge coefficient of the valve*/
		double ro /*Desity of fluid in the valve */, Dbore /*Smallest cross section in valve*/, xe /*!< spring pre-compression */, p_set /* set pressure */;
		double mp_nevl; /*Design pressure */
		double xmax_unrestricted; /*maximum allowed displacement */
		double RT; /* restriction, percent, 0=unrestricted*/
		double xmax; /* restricted (actual) max. lift*/
		double r /*!<restitution factor */;
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

		//double SumStatForce(const double x, const double pv, const double pb);
		// double Aeff(const double x);
		double Aeffmax, a1, a2, a3;
		vector<vector<double> > data;
		vector<double> tmpvec;
		Gas* gas;
		bool is_Gas;


		double Get_MassFlowRate_Compressible(double p_upstream, double T_upstream, double p_downstream, double T_downstream, double xx);
		double Get_MassFlowRate_Compressible_Choked(double p_upstream, double T_upstream, double x);
		double Get_MassFlowRate_Compressible_UnChoked(double p_upstream, double T_upstream, double p_downstream, double xx);
		double Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho, double x) ;
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
				Gas* _gas,
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

		double Aeff(const double x);
		void SetAeffCoeffs(const double _a1, const double _a2, const double _a3, const double Aeffmax);

		double Get_MassFlowRate(double p_upstream, double T_upstream, double p_downstream, double T_downstream, double xx);

//		double Get_MassFlowRate_Compressible(double p_upstream, double T_upstream, double p_downstream, double T_downstream, double xx);
//		double Get_MassFlowRate_Compressible_Choked(double p_upstream, double T_upstream, double x);
//		double Get_MassFlowRate_Compressible_UnChoked(double p_upstream, double T_upstream, double p_downstream, double xx);
//		double Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho, double x) ;
		//void Plot(bool show, bool save);
		double SumStatForce(const double x, const double pv, const double pb);
};
