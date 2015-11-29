/*  Legendre polynomials

	An in-class exercise using equations for the Polynomials and solvers as
	shown in the lecture/manuscript for Numerical Methods for CSE
	by Prof. R. Hiptmair, ETH Zürich

	Include the Eigen3 library as shown in documentation for Eigen3.

	use piping to store the .m file. Example call:
	legendre >legendre.m

 	This program calculates
 	- Legendre Polynomials P0 to P8 and their derivatives in interval [-1,1]
 	  -> plots them using MatlabPlotter
 	- Gauss points / zero points for P1 to P8
 	  - using secant method as solver
 	  - using secant falsi method as solver
 	  -> plots them using MatlabPlotter
 	- uses these Gauss points to calculate the weights for GL quadrature
 	- applies this GL quadrature to a function f(x) = e^(x^2) over an interval [a,b] = [3,6]
 	- plots the relative error (comparison of P1 to P8 vs reference result by Wolfram|Alpha)

	v1.0.3 2015-11-22 / 2015-11-29 Pirmin Schmid
*/

//#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "matlab_plotter.h"

using namespace std;
using namespace Eigen;

//------------------------------------------------------------------------------

// 3-term recursion for sequences (Pn)n  and (Pn')n for n in 0 to (N-1)
// evaluation of multiple x_i values in parallel
// input:        vector x in R^n with x_0 to x_{n-1}
// input/output: Lx and DLx in R^Nxn
//               thus, number of given rows N eq. number of elements desired
//               of the sequences (Pn)n  and (Pn')n
void legvals(const VectorXd& x, MatrixXd& Lx, MatrixXd& DLx) {
	// check input
	long n = x.size();
	long N = Lx.rows();
	if(Lx.cols() != n || DLx.rows() != N || DLx.cols() != n) {
		cout << "Error in legvals(): Dimensions mismatch." << endl;
		exit(1);
	}

	double denominator_inv = 0.0;
	double numerator1 = 0.0;
	double numerator2 = 0.0;

	// we need a row vector for the calculations below (otherwise Eigen assertion fails)
	// this seemed the easiest way for me to get one.
	RowVectorXd xr = x;

	DLx.row(0) = MatrixXd::Zero(1, n);
	Lx.row(0) = MatrixXd::Ones(1, n);

	DLx.row(1) = MatrixXd::Ones(1, n);
	Lx.row(1) = x;

	for(long i = 2; i < N; i++) {
		// these values are calculated fresh for each iteration to avoid any accumulating rounding error
		// when calculated iteratively from prior values
		denominator_inv = 1.0 / (double)i; // one division here -> many much faster multiplications than divisions later
		numerator1 = (double)(2 * i - 1);
		numerator2 = (double)(i - 1);

		DLx.row(i) = denominator_inv * (numerator1 * (Lx.row(i-1) + DLx.row(i-1).cwiseProduct(xr)) - numerator2 * DLx.row(i-2));
		Lx.row(i)  = denominator_inv * (numerator1 * Lx.row(i-1).cwiseProduct(xr) - numerator2 * Lx.row(i-2));
	}
}

//------------------------------------------------------------------------------

using EvalFunction = function<double (const double, const int)>;

// computes Pk(x) for scalar x
double Pkx(const double x, const int k) {
	if(k < 2) {
		if(k == 0) {
			return 1.0;
		}

		if(k == 1) {
			return x;
		}

		cout << "Error in Pkx(): Negative k values are not valid." << endl;
		exit(1);
	}

	double denominator = 0.0;
	double numerator1 = 0.0;
	double numerator2 = 0.0;

	double result_minus2 = 1.0;
	double result_minus1 = x;
	double result = 0.0;

	for(int i = 2; i <= k; i++) {
		denominator = (double)i;
		numerator1 = (double)(2 * i - 1);
		numerator2 = (double)(i - 1);

		result = (numerator1 * result_minus1 * x - numerator2 * result_minus2) / denominator;
		result_minus2 = result_minus1;
		result_minus1 = result;
	}
	return result;
}

//------------------------------------------------------------------------------

using Solver = function<double (double, double, EvalFunction, int, const double, const double, const int)>;

// translation of the Matlab function secant 2.3.25 in the manuscript
// modified to include the additional parameter k
double secant(double x0, double x1, EvalFunction f, int k, const double rtol, const double atol, const int maxIterations) {
	double f0 = f(x0, k);
	double fn = 0.0;
	double s = 0.0;
	for(int i=0; i < maxIterations; i++) {
		fn = f(x1, k);
		s = fn * (x1-x0) / (fn-f0); // correction
		x0 = x1;
		x1 = x1 - s;
		if( abs(s) < max(atol, rtol * min(abs(x0), abs(x1))) ) {
			return x1;
		}
		f0 = fn;
	}

	// default, best guess after maxIterations
	return x1;
}

// translation of the Matlab function secant_falsi on the exercise sheet
// modified to include the additional parameter k
double secant_falsi(double x0, double x1, EvalFunction f, int k, const double rtol, const double atol, const int maxIterations) {
	double f0 = f(x0, k);
	double fn = 0.0;
	double s = 0.0;
	for(int i=0; i < maxIterations; i++) {
		fn = f(x1, k);
		s = fn * (x1-x0) / (fn-f0); // correction
		if(f(x1 - s, k) * fn < 0.0) {
			x0 = x1;
			f0 = fn;
		}
		x1 = x1 - s;
		if( abs(s) < max(atol, rtol * min(abs(x0), abs(x1))) ) {
			return x1;
		}
	}

	// default, best guess after maxIterations
	return x1;
}

//------------------------------------------------------------------------------

#define MAX_ITERATIONS 100

// calculate zeros of Pk, k in 1 to n using the secant rule for end points {-1, 1} of the interval [-1, 1] and the zeros
// of the previous Legendre polynomial as initial guesses. Correction based termination criterion
// input:  n size
//         rtol and atol relative and absolute tolerance
// return: nxn upper triangular matrix, to actually get such an upper triangular
//         I assume that row j indicates the j-th zero for j in 1 to k
//                       column k indicates the solutions for Pk (thus, we will have column vectors of solutions)
//         note: for C++, index 0 will refer to j=1 and k=1 respectively, and so on.
MatrixXd gaussPts(const int n, Solver z, const double rtol = 1e-10, const double atol = 1e-12) {
	MatrixXd zeros = MatrixXd::Zero(n, n);

	// find the first for P1 -> will be in [0,0]
	zeros(0,0) = z(-1.0, 1.0, Pkx, 1, rtol, atol, MAX_ITERATIONS);

	// get the zeros for P2 to Pn (will be in columns 1 to (n-1)
	for(int i = 1; i < n; i++) {
		// get first zero
		zeros(0, i) = z(-1.0, zeros(0, i-1), Pkx, i+1, rtol, atol, MAX_ITERATIONS);

		// get last zero
		zeros(i, i) = z(zeros(i-1, i-1), 1.0, Pkx, i+1, rtol, atol, MAX_ITERATIONS);

		// get the zeros in-between
		for(int j = 1; j < i; j++) {
			zeros(j, i) = z(zeros(j-1, i-1), zeros(j, i-1), Pkx, i+1, rtol, atol, MAX_ITERATIONS);
		}
	}

	return zeros;
}

//------------------------------------------------------------------------------

#define A 3
#define B 6

using Function = function<double (const double)>;

// just a simple function that does not have a primitive (Stammfunktion) that can be expressed
// in R space. thus: suitable for numerical integration / quadrature
double test_function_for_quadrature(const double x) {
	return exp(x * x);
}

// Wolfram|Alpha calculated the quadrature of this function in [3,6] to be
// 3.644831077835569048422984645481051411815484722480248338949090926023254915628803401963716304967305392 10^14
// change this reference value if you change the function or A, B

#define REFERENCE_RESULT 3.6448310778355690e14

// applies the Gauss-Legendre quadrature for the given function over a defined interval [a,b]
// input:  f     a function that fits the type definition of Function
//         a, b  define interval [a,b]
//         w, x  weights and Gauss points for the given Legendre Polynomial in standard interval [-1,1]
//               size of both arrays must match, of course
// return: quadrature approximation for this function in interval [a,b]
double GLquadrature(const Function f, const double a, const double b, const ArrayXd& w, const ArrayXd& x) {
	int n = w.size();
	if(n != x.size()) {
		cout << "vectors of weights and Gauss points must have the same size" << endl;
		exit(1);
	}

	double half_delta = 0.5 * (b-a);
	double avg = 0.5 * (a+b);
	ArrayXd weights = half_delta * w;
	ArrayXd xs = half_delta * x;
	xs += avg;

	ArrayXd ys = xs.unaryExpr(f);
	ys *= weights;
	return ys.sum();
}

//------------------------------------------------------------------------------

#define MIN -1
#define MAX 1
#define MAX_K 8
#define N_X 600

int main() {
	// initialization
	MatlabPlotter p;
	p.comment("Legendre polynomials");
	p.comment("Code generated by legendre.cpp");

	vector<string> colors = {"k-", "b-", "g-", "r-", "c-", "m-", "y-"};
	vector<string> colors2 = {"ko", "bo", "go", "ro", "co", "mo", "yo"};
	int n_colors = colors.size();

	vector<string> description = {"P_{0}", "P_{1}", "P_{2}", "P_{3}", "P_{4}", "P_{5}", "P_{6}", "P_{7}", "P_{8}"};


	// (1) get a visual impression -> plot the Legendre polynomials up to k = MAX_K
	VectorXd xx = VectorXd::LinSpaced(N_X, MIN, MAX);
	MatrixXd Lx(MAX_K+1, N_X);
	MatrixXd DLx(MAX_K+1, N_X);
	legvals(xx, Lx, DLx);

	p.figure("Legendre polynomials 0 to 8");
	VectorXd yy = Lx.row(0);
	p.plot(xx, yy, colors[0]);
	p.hold();
	p.title("Legendre polynomials 0 to 8");
	for(int i=1; i <= MAX_K; i++) {
		yy = Lx.row(i);
		p.plot(xx, yy, colors[i % n_colors]);
	}
	p.xylabels("x", "P_{i}(x) for i={0, 1, ..., 8}");
	p.legend("P_{0}", "P_{1}", "P_{2}", "P_{3}", "P_{4}", "P_{5}", "P_{6}", "P_{7}", "P_{8}");
	p.hold(false);

	// (2) get a visual impression -> plot the derivatives of the Legendre polynomials up to k = MAX_K
	p.figure("Derivatives of Legendre polynomials 0 to 8");
	yy = DLx.row(0);
	p.plot(xx, yy, colors[0]);
	p.hold();
	p.title("Derivatives of Legendre polynomials 0 to 8");
	for(int i=1; i <= MAX_K; i++) {
		yy = DLx.row(i);
		p.plot(xx, yy, colors[i % n_colors]);
	}
	p.xylabels("x", "dP_{i}(x)/dx for i={0, 1, ..., 8}");
	p.legend("dP_{0}/dx", "dP_{1}/dx", "dP_{2}/dx", "dP_{3}/dx", "dP_{4}/dx", "dP_{5}/dx", "dP_{6}/dx", "dP_{7}/dx", "dP_{8}/dx");
	p.hold(false);

	// (3) Find zeros with secant and secant falsi methods
	//     and since we are iterating thru the Legendre polynomials and their Gauss points
	//     -> use the gained insight to apply GL quadrature to a given function and interval
	//        and measure the error of approximation
	MatrixXd zeros_x = gaussPts(MAX_K, secant);
	MatrixXd zeros_x_falsi = gaussPts(MAX_K, secant_falsi);
	VectorXd zeros_y = VectorXd::Zero(MAX_K);

	vector<double> x = {MIN, MAX};
	vector<double> y = {0.0, 0.0};

	ArrayXd i_values(MAX_K+1);
	i_values[0] = 0.0;
	ArrayXd areas(MAX_K+1);
	areas[0] = 0.0; // P0 is not tested

	for(int i=1; i <= MAX_K; i++) {
		p.figure("Legendre polynomial " + description[i] + " Zeros / Gauss points by secant and secant falsi methods.");
		yy = Lx.row(i);
		p.plot(xx, yy, colors[i % n_colors]);
		p.hold();
		p.title("Legendre polynomial " + description[i] + " Zeros / Gauss points by secant and secant falsi methods.");
		// secant
		VectorXd gyy = zeros_y.head(i);
		VectorXd gxx = zeros_x.col(i-1);
		gxx = gxx.head(i);
		p.plot(gxx, gyy, "ko");
		// true position of these "zeros"
		MatrixXd Lgx(i+1, i);
		MatrixXd DLgx(i+1, i);
		legvals(gxx, Lgx, DLgx);
		VectorXd true_gyy = Lgx.row(i);
		p.plot(gxx, true_gyy, "b*"); // note: This may not be visible if secant falsi gets the same result
		                             // but P8 shows a clear difference
		p.comment("zeros for " + description[i] + " by secant method");
		for(int j=0; j < i; j++) {
			p.comment("x = " + to_string(gxx[j]) + " y = " + to_string(true_gyy[j]) + " expected 0.0");
		}

		// falsi
		gxx = zeros_x_falsi.col(i-1);
		gxx = gxx.head(i);
		p.plot(gxx, gyy, "r*");
		// true position of these "zeros"
		legvals(gxx, Lgx, DLgx);
		true_gyy = Lgx.row(i);
		p.comment("zeros for " + description[i] + " by secant falsi method");
		for(int j=0; j < i; j++) {
			p.comment("x = " + to_string(gxx[j]) + " y = " + to_string(true_gyy[j]) + " expected 0.0");
		}
		// calculation of weights for GL-quadrature on standard interval [-1,1]
		p.comment("");
		p.comment(description[i] + ": Weights w_i and gauss points x_i needed for Gauss-Legendre quadrature");
		p.comment("(integration approximation) on interval [-1, 1]. Use appropriate scaling for other intervals.");
		ArrayXd DLgx_squared = DLgx.row(i).array();
		DLgx_squared *= DLgx_squared;
		ArrayXd gww = gxx.cwiseProduct(gxx).array();
		gww = -gww; // intermediary step since - seems not to be defined in combination with scalars
		gww = gww + 1.0;
		gww = gww.cwiseProduct(DLgx_squared);
		gww = 2.0 / gww;
		for(int j=0; j < i; j++) {
			string index = to_string(j);
			p.comment("w_" + index + " = " + to_string(gww[j]) +" x_" + index + " = " + to_string(gxx[j]) );
		}

		// let's test this on an effective quadrature
		double area = GLquadrature(test_function_for_quadrature, A, B, gww, gxx.array());
		p.comment("-> quadrature of f(x)=e^(x^2) in [" + to_string(A) + "," + to_string(B) + "] approx. = " + to_string(area));
		areas[i] = area;
		i_values[i] = (double)i;

		// null line
		p.plot(x, y, "k:");
		// info
		p.xylabels("x", description[i] + "(x)");
		p.legend(description[i], "zeros by secant", "true y value of these zeros", "zeros by secant falsi");
		p.hold(false);
	}

	// show relative errors (P1 to P8 vs reference result)
	double negReference = -REFERENCE_RESULT;
	ArrayXd error = areas + negReference; // workaround since - is not accepted with scalars (while + is)
	error = error.cwiseAbs() / REFERENCE_RESULT;

	// do not show P0
	i_values = i_values.tail(MAX_K);
	error = error.tail(MAX_K);

	p.figure("Relative errors of P1 to P8 vs reference result from Wolfram|Alpha", MatlabPlotter::LINEAR);
	p.plot(i_values, error, "r*-");
	p.xylabels("P_{i}", "Relative error vs reference result (lin scale)");
	p.title("Relative errors of P1 to P8 vs reference result from Wolfram|Alpha");

	p.figure("Relative errors of P1 to P8 vs reference result from Wolfram|Alpha", MatlabPlotter::SEMILOGY);
	p.plot(i_values, error, "r*-");
	p.xylabels("P_{i}", "Relative error vs reference result (log scale)");
	p.title("Relative errors of P1 to P8 vs reference result from Wolfram|Alpha");

	return 0;
}