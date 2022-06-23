#define _USE_MATH_DEFINES
#include <stdio.h>
#include <vector>
#include <iostream>
#include <utility>  // pair, make_pair
#include <random>
#include <numeric>

#include <Eigen/Dense>

#include "my_tools.h"
/*#include </usr/include/eigen3/Eigen/Dense>*/
//#include </usr/local/Cellar/eigen/3.3.4/include/eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

//Compare two evaluated multidimensional functions
/*! Compare two evaluated multidimensional functions where the data is stored in a pair, 
	of which the first one is the evaluation value (double) and the second one is the evaluation location
	\param i First pair of value and vector
	\param j Second pair of value and vector
	\return The boolean comparing the evaluated values (".first" in the pair)
*/
bool compare(const pair<double, VectorXd>&i, const pair<double, VectorXd>&j)
{
	return i.first < j.first;
}

//double signum(double x);

double mean(vector<double> v){
	/*int N=x.size();
	double mean = 0;
	for (unsigned int i=0; i<N; i++)
		mean+=x.at(i);
	return (mean/(double)N);*/
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	return sum / v.size();

/*double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
double stdev = std::sqrt(sq_sum / v.size() - mean * mean);*/
}

double stdev(vector<double> v){
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
double mean = sum / v.size();

double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
double stdev = std::sqrt(fabs(sq_sum / v.size() - mean * mean));
return stdev;
}

//=============================================================================

//! An implementation of the multidimensional Newton-Raphson method.
/*! This function finds the root of a function from a given starting point using the multi-dimensional Newton-Raphson method. 
	The implementation also allows for the use of a VectorXd of parameters.
	The maximum available error and number of iterations also need to be declared.
	\return Whether the search forthe solution was succesful */
bool NR(VectorXd (*fun)(VectorXd & x, const VectorXd & pars), /**< [in] pointer to the function that the solution is searched for*/
        VectorXd & x, /**< [in,out] The starting location and the solution */
        const VectorXd & pars, /**< [in] Input parameters for the function */
        const bool show_iter, /**<  [in] Display the number of iterations*/
        const double err_max, /**< [in] Maximum allowed error */
        const int iter_max, /**< [in] Maximum number of iterations*/
        const bool stop_on_failure /**< [in] Exit the program if iteration was succesful */) {

	int Neq = x.size(); //gather number of equations

	VectorXd f(Neq), xnew;
	double err;//, err_max = 1.e-9;
	int iter = 0;//, iter_max = 500;
	bool success = true; //assume success

	f = (*fun)(x, pars); //Evaluate function at start point

	err = f.norm();
	if (show_iter) {
		printf("\n\t iter=%2d, |f|=%5.3e", iter, err);
		// cin.get();
	}

	double RELAX = 0.1;
	while ((err > err_max) && (iter<iter_max)) {
		iter++;
		// cout << endl << "x=" << x << endl;
		// cout << endl << "jac=" << num_jac((*fun), x, pars) << endl;
		// cout << endl << "jac.inverse()=" << num_jac((*fun), x, pars).inverse() << endl;
		// cin.get();

		xnew = x - num_jac((*fun), x, pars).inverse() * f;
		x = (1.-RELAX)*xnew + RELAX*x;
		f = (*fun)(x, pars); 

		err = f.norm();
		if (show_iter)
			printf("\n\t iter=%2d, |f|=%5.3e", iter, err);
		if (iter == iter_max) {
			success = false;
			if (stop_on_failure){
				printf("\n\t iter=%2d, |f|=%5.3e", iter, err);
				cout << endl << "TOO MANY ITERATIONS!" << endl << "NOW EXITING..." << endl;
				exit(-1);
			}			
		}
	}
	if (show_iter) {
		cout << "\n\t Solution vector: x=(";
		printf("%+5.3f", x(0));
		for (int i = 1; i < Neq; i++)
			printf(", %+5.3f", x(i));
		printf(")\n");
	}
	return success;
}

//! Calculation of the Jacobian of a multi dimensional parametric function
/*! Calculates the Jacobian of a  given parametric function while being
	provided the location and the parameters.
	\return The Jacobian matrix*/
MatrixXd num_jac(VectorXd (*fun)(VectorXd & x, const VectorXd & pars), /**< The function to derive  */
                 VectorXd & x, /**< [in] The starting location and the solution */
                 const VectorXd & pars /**< [in] Input parameters for the function */ ) {

	int Neq = x.size();
	MatrixXd out(Neq, Neq); //Establish the Jacobian to be filled
	VectorXd f0(Neq);
	VectorXd f1(Neq);
	VectorXd x0, dx;
	f0 = (*fun)(x, pars);

	double TINY = 0.00001;

	x0 = x;
	for (unsigned int i = 0; i < Neq; i++) {
		dx = VectorXd::Zero(Neq);
		if (fabs(x(i)) > TINY)
			dx(i) = x(i) * 0.0001;
		else
			dx(i) = TINY; //Established step size
		x = x0 + dx;
		f1 = (*fun)(x , pars);
		f1 = (f1 - f0) / dx(i); //numerical difference
		for (unsigned int j = 0; j < Neq; j++)
			out(j, i) = f1(j); //filling up the columns
	}
	x = x0;

	return out;
}

//! Calculation of the gradient of a function
/*! Calcultates the gradient of a parametric multidimensional scalar function at a
	given location and set of parameters.
	\return A vector containing the derivatives
	\sa num_jac

\return vector of the derivative in the base coordinates*/
VectorXd num_grad(double (*fun)(VectorXd & x, const VectorXd & pars), /**< [in] parametric function to derive */
                  VectorXd & x, /**< [in] derivation location*/
                  const VectorXd & pars /**< [in] paramters of the parametric function*/ ) {

	int Neq = x.size();
	VectorXd out(Neq); //results vector
	double f0, f1;
	VectorXd x0, dx;

	double TINY = 0.0001;
	x0 = x;
	f0 = (*fun)(x, pars);
	for (unsigned int i = 0; i < Neq; i++) {
		dx = VectorXd::Zero(Neq);
		if (fabs(x(i)) > TINY)
			dx(i) = x(i) * TINY;
		else
			dx(i) = TINY;
		x = x0 + dx;
		f1 = (*fun)(x , pars);
		f1 = (f1 - f0) / dx(i);
		out(i) = f1;
	}
	x = x0;

	return out;
}

//! Create a (lin)space
/*! Creates a distribution between a minimum and maximum value (both included). The
	distribution has a fixed number of points, and can be either linear or exponential.
	\param x_min The lower bound
	\param x_max The upper bound
	\param Nsteps Number of steps (elements in the return matrix)
	\param is_linear If true the distribution is linear, else it is exponential
	\return A vector with the linspace
*/
VectorXd linspace(double x_min, double x_max, int Nsteps, bool is_linear) {
	double d_ap;
	VectorXd out = VectorXd::Zero(Nsteps);
	out(0) = x_min;
	if (is_linear) {
		d_ap = (x_max - x_min) / ((double) Nsteps - 1.);
		for (unsigned int i = 1; i < Nsteps; i++)
			out(i) = out(i - 1) + d_ap;
	}
	else {
		d_ap = pow(x_max / x_min, 1. / Nsteps);
		for (unsigned int i = 1; i < Nsteps; i++)
			out(i) = out(i - 1) * d_ap;
	}

	return out;
}

vector<double> linspace_vec(double x_min, double x_max, int Nsteps, bool is_linear) {
	double d_ap;
	vector<double> out(Nsteps);
	out.at(0) = x_min;
	if (is_linear) {
		d_ap = (x_max - x_min) / ((double) Nsteps - 1.);
		for (unsigned int i = 1; i < Nsteps; i++)
			out.at(i) = out.at(i - 1) + d_ap;
	}
	else {
		d_ap = pow(x_max / x_min, 1. / Nsteps);
		for (unsigned int i = 1; i < Nsteps; i++)
			out.at(i) = out.at(i - 1) * d_ap;
	}

	return out;
}


//! Determine the sign of the number
/*! Determine the sign of a given double and return 1.0 if positive and -1.0 if negative
	\param x the value to evaluate
	\return -1.0 if negative and 1.0 if positive 
*/
double signum(double x) {
	// if (x > 0.) return 1.;
	// if (x < 0.) return -1.;
	// return 0.;
	if (x > 0.)
		return 1.;
	else
		return -1.;
}

//! Use the Gradient Descent to find the minimum of a function
/*! Determines the minimum location of a multi-dimensional parametric function, using 
	given parameters and a starting method. Uses the Gradient descent (also known as
	the steepest descent.)
*/
void SteepestDescent(
    double (*fun)(VectorXd & x, const VectorXd & pars), /**< [in] The function to evaluate */
    VectorXd & x, /**< [in, out] The start and end point*/
    const VectorXd & pars, /**< [in] The parameters of the evaluated function*/
    const bool show_iter /**< Whether to print out the current iterations or not */) {

	int Nvar = x.size();

	VectorXd grad(Nvar), gradp(Nvar), dx(Nvar), dxp(Nvar);
	double gamma, f, err = 1.e5, err_max = 1.e-6;
	int iter = 0, iter_max = 500;
	
	gamma = 0.0001; //step multiplier


	while (err > err_max) {
		f = (*fun)(x, pars); //Evaluate the function in place

		grad = num_grad(fun, x, pars); //Gather the gradient as well

		if (show_iter)
			cout << "\n step " << iter << ", x=(" << x.transpose() << "), err=" << err << ", gamma = " << gamma << ",  | f |= " << f << ", grad = (" << grad.transpose() << ")";

		dx = gamma * grad; //stepsize
		err = dx.norm(); //The "error" here is the size of the step.
		x -= dx; //do the stepping
		iter++;
		if (iter == iter_max) {
			cout << "\n\n ERROR: my_tools.SteepestDescent() -> Too many steps!";
			exit(-1); //Program exits if too many steps have been done
		}

		// Update stepsize
		if (iter > 1) {
			gamma = -dxp.transpose() * (grad - gradp);
			double tmp = (grad - gradp).norm();
			gamma = gamma / tmp / tmp; //Calculates gamma with the Barzilai-Borwein method to ensure that the solution converges (doi:10.1093/imanum/8.1.141)
			// if (show_iter)
			// cout << "\n\t gamma=" << gamma;
		}
		gradp = grad;
		dxp = dx; //set up having a previous step
		// cin.get();
	}
}

//! A 2D Nelder-Mead optimum searcher (for illustratory purposes) DO NOT USE!
/*! A 2D Nelder-Mead optimum searcher. Accepts an any dimensional (but will be treated as 2) 
	parametric function, along with a starting point, radius and of course the parameters.
	There is no real output, only the convergence is displayed on the console, with the solution
	always running for a 100 iterations. 
*/
void NelderMead2D(
    double (*fun)(VectorXd & x, const VectorXd & pars), /**< [in] The function to optimize */
    VectorXd & x, /**< [in] The starting (centre) point*/
    double radius, /*< [in] The radius od the starting simplex around the centre point*/
    const VectorXd & pars, /**< [in] The parameters necessary for evaluating the function */
    const bool show_iter) /**< [in] Show iterations? (not used, iterations always shown)*/ {

//	int Nvar = x.size();
	int Npts = 3;

	double alpha = 1.; //reflection coefficient
	double alpha_min = 0.001;
	double gamma = 2.; // gamma > 1, expansion coefficient
	double ro = 0.5; // 0<ro<0.5, contraction coefficient
	double sigma = 0.5; //shrink coefficient

	// Initialize
	double f1, f2, f3, fr, fe, fc, ftmp;
	VectorXd x1(2), x2(2), x3(2), xo(2), xr(2), xe(2), xc(2), xtmp(2);
	// double radius = 0.1;
	x1(0) = x(0) + radius * sin(0.);
	x2(0) = x(0) + radius * sin(2. * M_PI / 3.);
	x3(0) = x(0) + radius * sin(2. * 2.* M_PI / 3.);
	x1(1) = x(1) + radius * cos(0.);
	x2(1) = x(1) + radius * cos(2. * M_PI / 3.);
	x3(1) = x(1) + radius * cos(2. *2 * M_PI / 3.); 
	//Starting simplex established as a regular triangle around the given starting point with the predetermined radius

	// Eval
	f1 = (*fun)(x1, pars);
	f2 = (*fun)(x2, pars);
	f3 = (*fun)(x3, pars);

	vector< pair<double, VectorXd> >pts;
	pts.push_back(make_pair(f1, x1));
	pts.push_back(make_pair(f2, x2));
	pts.push_back(make_pair(f3, x3));
	//Collects the data points and their evaluation into pairs, which then are inserted into a std::vector
	//First value is the evaluation value and the second one is the location VectorXd


	/* for (int i = 0; i < Npts; i++)
	 	printf("\n x%d: (%+5.3f,%5.3f), f=%5.3e",
	 	       i, pts.at(i).second(0), pts.at(i).second(1), pts.at(i).first); */

	//Begins the iterative process
	for (int step = 0; step < 100; step++) {

		printf("\n\n STEP #%d\n==============", step);

		// Sort the three points according to evaluation value such as f(u)<f(v)<f(w)
		// where f(w) is the "worst" point, while the others are "good"
		sort(pts.begin(), pts.end(), compare);

		printf("\n Sorted points:");
		for (int i = 0; i < Npts; i++)
			printf("\n\t x%d: (%+5.3f,%5.3f), f=%5.3e",
			       i, pts.at(i).second(0), pts.at(i).second(1), pts.at(i).first);

		// Centroid calculated between the (two) "good"
		xo = (pts.at(0).second + pts.at(1).second) / 2.;
		// printf("\n\t    xo = (%+5.3f,%5.3f)", xo(0), xo(1));

		// Reflection of the "worst" point through this centroid
		xr = xo + alpha * (xo - pts.at(2).second);
		fr = (*fun)(xr, pars);
		printf("\n Reflected point:");
		printf("\n\t xr: (%+5.3f,%5.3f), f=%5.3e", xr(0), xr(1), fr);

		if (fr < pts.at(2).first) { //If the reflected value is better than the "worst"
			printf("\n fr<f3, ");
			if (fr < pts.at(1).first) { //If it is also better than the second
				printf(" fr<f2, ");
				if (fr > pts.at(0).first) { //but worse than the best
					printf(" fr>f1 -> Replacing x3 by xr");
					pts.at(2) = make_pair(fr, xr);  //Worst point simply needs to be replaced
				}
				else { //If the new point is the best, Expansion is performed
					printf(" fr<f1, the best so far -> performing Expansion");
					xe = xo + gamma * (xr - xo); //expanded point computed
					fe = (*fun)(xe, pars); //evaluate the expanded point
					printf("\n After expansion:");
					printf("\n\t xe: (%+5.3f,%5.3f), f=%5.3e", xe(0), xe(1), fe);
					if (fe < fr) {
						printf("\n\t fe<fr, -> Replacing x3 by xe");
						pts.at(2) = make_pair(fe, xe); //If the expanded is better, the worst point is replaced with it
					}
					else {
						printf("\n\t fe>fr, -> Replacing x3 by xr");
						pts.at(2) = make_pair(fr, xr); //If the expanded is worse then we replace the worst with the reflected
					}
				}
			}
			else { // fr > f2
				printf(" fr>f2 -> performing contraction ");
				xc = xo + ro * (pts.at(2).second - xo);
				fc = (*fun)(xc, pars); //contracted point calculated and evaluated
				printf("\n After contraction:");
				printf("\n\t xc: (%+5.3f,%5.3f), f=%5.3e", xc(0), xc(1), fc);
				if (fc < pts.at(2).first) {
					printf("\n\t fc<f3, -> Replacing x3 by xc");
					pts.at(2) = make_pair(fc, xc); //If contracted is better than the worst, replace it with it
				}
				else { //All points are shrunk around the best one (except itself of course)
					printf("\n\t fc>f3, -> Shrink");
					xtmp = pts.at(0).second + sigma * (pts.at(1).second - pts.at(0).second);
					ftmp = (*fun)(xtmp, pars);
					pts.at(1) = make_pair(ftmp, xtmp);

					xtmp = pts.at(0).second + sigma * (pts.at(2).second - pts.at(0).second);
					ftmp = (*fun)(xtmp, pars);
					pts.at(2) = make_pair(ftmp, xtmp);
				}
			}
		}
		else { //If the reflection did not give us a better point
			alpha /= 2.;
			printf("\n fr>f3 ??? Retrying with smaller alpha = %g", alpha);
			if (alpha < alpha_min) {
				printf("\n\n ERROR!!! alpha<alpha_min, exiting...\n\n");
				break;
			}
		}
		cin.get();
	}
}



void NelderMead(
    double (*fun)(VectorXd & x, const VectorXd & pars), /**< [in] The function to evaluate */
    VectorXd & x, /**< [in,out] The starting vector*/
    double radius, /*[in] Radius of the first simplex*/
    const VectorXd & pars, /**< The parameters necessary for the evaluation of  */
    const double x_TOL,
    const double fun_TOL,
    const bool show_iter) {

//	int Nvar = x.size();
	int Npts = 3;

//	double alpha = 1.; //Reflectiomn coefficient
//	double alpha_min = 0.001; //Minimum allowed reflection coefficient
//	double gamma = 2.; // gamma > 1, expansion coefficient
//	double ro = 0.5; // 0<ro<0.5 contraction coefficient
//	double sigma = 0.5; //Shrink coefficient

	// Initialize
	double f1, f2, f3, fR, fC1, fC2, fC, fS, fM, fE;
	VectorXd x1(2), x2(2), x3(2), xR(2), xC1(2), xC2(2), xC(2), xS(2), xM(2), xE(2);
	//Calculate the points in a triangle arount the centre point 
	x1(0) = x(0) + radius * sin(0);
	x2(0) = x(0) + radius * sin(2. * M_PI / 3.);
	x3(0) = x(0) + radius * sin(2. * 2.* M_PI / 3.);
	x1(1) = x(1) + radius * cos(0.);
	x2(1) = x(1) + radius * cos(2. * M_PI / 3.);
	x3(1) = x(1) + radius * cos(2. * 2. * M_PI / 3.);

	// Evaluate all three starting points
	f1 = (*fun)(x1, pars);
	f2 = (*fun)(x2, pars);
	f3 = (*fun)(x3, pars);

	//A list of them is created, where the list elements are pairs. Of these, the first value is
	//a double containing the evaluated value while the second is the location vector.
	vector< pair<double, VectorXd> >pts;
	pts.push_back(make_pair(f1, x1));
	pts.push_back(make_pair(f2, x2));
	pts.push_back(make_pair(f3, x3)); 

	double diff_x = 1e5;
	//double diff_x_max = 1.e-5;
	int iter = 0, iter_max = 1000;

	//Main itartions starting if convergence still happening and
	// function values are not low enough
	while ((diff_x > x_TOL) && (pts.at(0).first > fun_TOL)) {

		// Sort the points by the evaluation value
		sort(pts.begin(), pts.end(), compare); //First two are good, while the largest one is the "bad" point

		if (show_iter) {
			printf("\n\n STEP #%d\n==============", iter);

			printf("\n Sorted points:");
			for (int i = 0; i < Npts; i++)
				printf("\n\t x%d: (%+5.3e,%+5.3e), f=%5.3e",
				       i, pts.at(i).second(0), pts.at(i).second(1), pts.at(i).first);
			// cin.get();
		}

		// Centroid is calculated from the two good points, and evaluated
		xM = (pts.at(0).second + pts.at(1).second) / 2.;
		fM = (*fun)(xM, pars);

		// Reflection is done trough the centre point
		xR = 2.*xM - pts.at(2).second;
		fR = (*fun)(xR, pars);

		if (fR < pts.at(1).first) { 
			// printf("\n\n fR < fG, Case i.");
			if (pts.at(0).first < fR) { 
				//if the reflection point between the the good ones, it replaces the worst
				// printf("\n\n fB< fR");
				pts.at(2) = make_pair(fR, xR); 
				
			}
			else {
				//if the reflection point is the best, an expanded point is calculated
				//The better of the refleced/expanded point is used to replece the worst point
				xE = 2 * xR - xM;
				fE = (*fun)(xE, pars);
				if (fE < pts.at(0).first)
					pts.at(2) = make_pair(fE, xE);
				else
					pts.at(2) = make_pair(fR, xR);
			}
		}
		else {
			// printf("\n\n fR > fG, Case ii.");
			if (fR < pts.at(2).first)
				pts.at(2) = make_pair(fR, xR); //if the expansion is between the ok and worst point, it replaces the worst
			xC1 = (pts.at(2).second + xM) / 2.;
			fC1 = (*fun)(xC1, pars);
			xC2 = (xR + xM) / 2.;
			fC2 = (*fun)(xC2, pars);
			if (fC1 < fC2) {
				xC = xC1;
				fC = fC1;
			}
			else {
				xC = xC2;
				fC = fC2;
			}
			if (fC < pts.at(2).first)
				pts.at(2) = make_pair(fC, xC);
			else {
				xS = (pts.at(2).second + pts.at(0).second) / 2;
				fS = (*fun)(xS, pars);
				pts.at(2) = make_pair(fS, xS);
				pts.at(1) = make_pair(fM, xM);
			}
		}
		diff_x = (pts.at(0).second - pts.at(1).second).norm();
		iter++;
		if (iter == iter_max) { //no convergence within the maximum allowed number of iterations
			printf("\n\n !!! ERROR !!! NelderMead, iter=%d, EXITING!", iter);
			printf("\nLast point: x=(%+6.4e, %6.4f), f=%6.4e\n\n",
			       pts.at(0).second(0), pts.at(0).second(1), pts.at(0).first);
			cin.get();
		}
	}
	x = pts.at(0).second;
	if (show_iter)
		printf("\n xmin = (%+6.4e, %+6.4e), f=%6.4e\n\n",
		       pts.at(0).second(0), pts.at(0).second(1), pts.at(0).first);
}

//Generate a random number in a given range
/*! Generates a random double precision in a range (inclusive of the edges) as
	an uniform distribution
	\param MIN_RAND Minimum allowed value
	\param MAX_RAND Maximum allowed value
	\return the random number
*/
double gen_random_double(double MIN_RAND, double MAX_RAND)
{
	double range = MAX_RAND - MIN_RAND;
	double random = range * ((((double) rand()) / (double) RAND_MAX)) + MIN_RAND ;
	return random;
}

//Swaps two vectors
/*! Swaps the values of two vectors
	\param x The first vector
	\param y The second vector
*/
void swap(VectorXd x, VectorXd y) {
	VectorXd tmp;
	tmp = x;
	x = y;
	y = tmp;
}

//The data reading sections begin.

void Input::add_datapair(Datapair _dp)
{
	elements++;
	dataarray[elements - 1] = _dp;
}

double Input::extract_value(std::string name)
{
	int i = 0; //counter to avoud calling non-existent elements
	double res = 0.0; //value to return
	bool over = false; //measures if value has been found

	while (i < elements && !over) {
		if (dataarray[i].name == name) { res = dataarray[i].value; over = true; }
		i++;

	}
	if (!over) { std::cout << name << " was not found in the system. 0 was passed as the value. Please check input" << std::endl; }
	return res;
}

void Input::parse_file(std::string filename)
{
	std::string line; //The line to read
	std::string name, val;
	std::size_t splitplace1, splitplace2; //The placement of the
	Datapair temp_pair;//{ "a", 0.0 };

	std::ifstream infile;
	infile.open(filename, std::ifstream::in);

	while (std::getline(infile, line)) {
		/*for (int i = 0; i < line.size(); i++) {
			if (line[i] == ';') { splitplace = i; }
		} //collect location of the ;*/
		splitplace1 = line.find(";");
		name = line.substr(0, splitplace1); //first should always be found
		splitplace2 = line.find(";", splitplace1 + 1);
		if (splitplace2 == std::string::npos) { val = line.substr(splitplace1 + 1, line.length() - splitplace1); }
		else { val = val = line.substr(splitplace1 + 1, splitplace2 - splitplace1 - 1); }
		/*std::cout << "Name: " << name << std::endl;
		std::cout << "Value: " << val << std::endl;*/

		elements++;
		dataarray[elements - 1].name = name;
		dataarray[elements - 1].value = atof(val.c_str());
	}
}
