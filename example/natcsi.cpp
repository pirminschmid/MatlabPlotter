/*  example for MatlabPlotter: Natural cubic splines (core problem)

	An in-class exercise using equations for the cubic splines as
	shown in the lecture/manuscript for Numerical Methods for CSE
	by Prof. R. Hiptmair, ETH Zürich

	Include the Eigen3 library as shown in documentation for Eigen3.

	use piping to store the .m file. Example call:
	natcsi >example.m

	given nodes (t_i,y_i) for i = 0,...,n

	and using eq 3.5.5, the calculation of y = s(x) boils down to
	- finding the correct interval [t_j-1,t_j] for x
	- apply s|[t_j-1,t_j](x) = s_j(x) =
		y_j-1 * (1 - 3 * tau^2 + 2 * tau^3) +
		y_j   * (3 * tau^2 - 2 * tau^3) +
		m_j-1 * (tau - 2 * tau^2 + tau^3) +
		m_j   * (-tau^2 + tau^3)

		with
		m1_j = m_j-1 = h_j * c_j-1,
		m2_j = m_j   = h_j * c_j
		
		c_j   = s'(t_j)
		
		h_j   = t_j - t_j-1

		tau   = (t - t_j-1) / h_j

	current runtime:
	- initialization: O(n), Thomas algorithm for solving the tridiagonal system matrix
	- evaluation:
		* worst case:        O(log n)
		* best/typical case: O(1)
	
	space requirements: O(n) all the time             

	output: is given again as a Matlab script that
	allows plotting all data in Matlab.

	v1.0 2015-11-05 / 2015-11-07 Pirmin Schmid 
*/

#include <cassert>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "matlab_plotter.h"

using namespace std;
using namespace Eigen;

// cubic spline interpolant with natural boundaries
class NatCSI {
public:
	
	// build cubic spline interpolant with natural boundaries
	// input:  t  nodes of the grid for pairs (t_i, y_i), sorted!
	//         y  values y_i at t_i for pairs (t_i, y_i), sorted!
	// input in std::vector<double> type
	NatCSI(const vector<double> &t, const vector<double> &y) {
		// map std::vector<double> to a VectorXd
		VectorXd tmp_t = VectorXd::Map(t.data(), t.size());
		VectorXd tmp_y = VectorXd::Map(y.data(), y.size());
		
		// store a copy of the data
		this->t = tmp_t;
		this->y = tmp_y;
		
		compute_coefficients();
	}
	
	// build cubic spline interpolant with natural boundaries
	// input:  t  nodes of the grid for pairs (t_i, y_i), sorted!
	//         y  values y_i at t_i for pairs (t_i, y_i), sorted!
	// input in Eigen::VectorXd type
	NatCSI(const VectorXd &t, const VectorXd &y) {
		// store a copy of the input vectors
		this->t = t;
		this->y = y;
		
		compute_coefficients();
	}
	
	// interpolant evaluation at x
	// input:   x value where to evaluate the spline
	// return:  y = s(x) if x in defined range of s(x)
	//          or 0.0 otherwise
	double operator() (double x) {
		// out of defined interval
		if(x < min || max < x) {
			return 0.0;
		}
		
		// get proper j, using a heuristic
		// 1) check whether x is in last range, O(1)
		// 2) check whether x is in next range, O(1)
		// 3) use binary search tree, O(log n)
		if(!in_range(x, last_j)) {
			if(last_j < n && in_range(x, last_j+1)) {
				// note: the second test is only applied if last_j actually allows it
				last_j++;
			}
			else {
				find_index(x);
			}
		}
		
		// apply function in interval j
		int j = last_j; // just for convenience
		
		double tau  = (x - t[j-1]) / h[j];
		double tau2 = tau * tau;
		double tau3 = tau2 * tau;
		
		double result = y[j-1] * (1.0 - 3.0 * tau2 + 2.0 * tau3);
		result += y[j] * (3.0 * tau2 - 2.0 * tau3);
		result += m1[j] * (tau - 2.0 * tau2 + tau3);
		result += m2[j] * (-tau2 + tau3);
		return result;
	}
	
private:
	// we need to store the vectors t and y
	// and a vector for m and h as described above
	// m is calculated after having solved a LSE
	// A * c = d as described below (and in the manuscript)
	
	VectorXd t, y, h, m1, m2;
	
	// number of segments of the spline function
	// note: n is a valid index for the vectors. Their size is (n+1).
	int n = 0;
	
	// valid range for x
	double min = 0.0;
	double max = 0.0;
	
	// since most requests for x calculations are in vicinity
	// of the prior one (or even sorted), remembering
	// the last requested section will help to bring
	// evaluation time to O(1) compared to O(log n) if only
	// a binary search tree is used.
	// the binary search tree is currently used as a fallback
	// in case the current partition/segment or the next one
	// has not been successful.
	int last_j = 1;
	
	// returns true if x is in interval j [t_j-1, t_j]
	bool in_range(const double x, const int j) const {
		assert(1 <= j && j <= n);
		return t[j-1] <= x && x <= t[j];
	}
	
	// sets last_j to  index of matching interval j [t_j-1, t_j] if x in range
	// note: x had to be tested before whether it is actually within the range of the spline function
	// note: optimizations (first test whether last interval was already ok or test the next
	// interval have already been done when one calls this function)
	// finds next j in O(log n) time
	void find_index(const double x) {
		assert(min <= x && x <= max);
		
		int j = last_j;
		int a = 1;
		int b = n;
		bool found = false;
		while(!found) {
			if(x < t[j-1]) {
				b = j;
				j = (a+b)/2;
			}
			else if(t[j] < x) {              
				if(a == b-1) {
					// special case needed due to integer division rounding down
					j = b;
					found = true;
				}
				else {
					a = j;
					j = (a+b)/2;
				}
			}
			else {
				// match
				found = true;
			}
		}
		
		last_j = j;
	}
	
	void compute_coefficients() {
		if(t.size() != y.size()) {
			// only minimalistic error management
			cout << "Error: size of t and y must be identical." << endl;
			exit(1);
		}
		
		// n
		n = (int)t.size() - 1;
		
		// min/max
		min = t[0];
		max = t[n];
		
		// h_j = t_j - t_j-i
		h.resize(n + 1);
		h.segment(1, n) = t.tail(n) - t.head(n);
		
		// build a virtual (n+1)x(n+1) system matrix A; however, only the diagonal vectors will actually be built
		// A = diag(a) + diag(b,1) + diag(b,-1), with size(a) = n+1 and size(b) = n
		
		// b_i = 1/h_i+1 for i=0 to n-1
		VectorXd b = h.tail(n).cwiseInverse();
		
		// a_i = 2/h_i + 2/h_i+1 for i=1 to n-1
		// a_0 = 2/h_1 (see natural cubic spline)
		// a_n = 2/h_n (see natural cubic spline)
		VectorXd a(n + 1);
		a.segment(1, n-1) = 2.0 * (h.segment(1, n-1).cwiseInverse() + h.tail(n-1).cwiseInverse());
		a[0] = 2.0 / h[1];
		a[n] = 2.0 / h[n];
		
		// right hand side d for A*c = d
		// for i=1 to n-1
		// d_i = 3.0 * ( (y_i - y_i-1) / h2_i  +  (y_i+1 - y_i) / h2_i+1 )
		// d_0 = 3.0 * (y_1 - y_0) / h2_1
		// d_n = 3.0 * (y_n - y_n-1) / h2_n
		
		VectorXd h2 = h.cwiseProduct(h); // h^2
		
		// note: y_i - y_i-1 is stored in delta_y[i-1]
		VectorXd delta_y = y.tail(n) - y.head(n);
		
		// quotient is either left or right side inside the brackets
		VectorXd quotient = delta_y.cwiseQuotient(h2.tail(n));
		
		// build d
		VectorXd d(n + 1);
		d.segment(1, n-1) = quotient.head(n-1) + quotient.tail(n-1);
		d[0] = delta_y[0] / h2[1];
		d[n] = delta_y[n-1] / h2[n];
		d = 3.0 * d;
		
		// Solution of this LSE with n+1 equations and n+1 unknowns c0 to cn
		// using the Thomas algorithm for tridiagonal matrices, needs O(n)
		// typical description with vectors b1..bn, a2..an, c1..cn-1
		// 'b' is a here
		// 'a' and 'c' are b here
		// 'x' is our c
		// note: different indices, too; also we have size n+1 instead of n
		
		// forward preparation
		VectorXd bprime(n);
		VectorXd dprime(n+1);
		
		// cprime -> bprime
		bprime[0] = b[0] / a[0];
		for(int i = 1; i < (n-1); ++i) {
			bprime[i] = b[i] / (a[i] - bprime[i-1] * b[i]);
		}
		
		// dprime
		dprime[0] = d[0] / a[0];
		for(int i = 1; i < n; ++i) {
			dprime[i] = (d[i] - dprime[i-1] * b[i]) / (a[i] - bprime[i-1] * b[i]);
		}
		
		// back substitution to get x, or c in our case
		VectorXd c(n + 1);
		c[n] = dprime[n];
		for(int i = n-1; i >= 0; --i) {
			c[i] = dprime[i] - bprime[i] * c[i+1];
		}
		
		// we do not actually need c but
		// m1_j = m_j-1 = h_j * c_j-1,
		// m2_j = m_j   = h_j * c_j
		m1.resize(n+1);
		m2.resize(n+1);
		m1.segment(1, n) = h.tail(n).cwiseProduct( c.head(n) );
		m2.segment(1, n) = h.tail(n).cwiseProduct( c.tail(n) );
	}
	
};

// just a quick implementation of the function
double f(const double x) {
	return exp(sin(2.0 * M_PI * x));
}

#define MIN -3
#define MAX 3
#define N_T 61
#define N_X 600

int main() {
	// initialization
	VectorXd t = VectorXd::LinSpaced(N_T, MIN, MAX);
	VectorXd y(N_T);
	for(int i = 0; i < N_T; ++i) {
		y[i] = f(t[i]);
	}
	
	NatCSI ncsi(t, y);
	
	// test 1
	VectorXd xx = VectorXd::LinSpaced(N_X, MIN, MAX);
	VectorXd yy(N_X);
	for(int i = 0; i < N_X; ++i) {
		yy[i] = ncsi(xx[i]);
	}
	
	// test 2. test the binary search in find_index()
	VectorXd y2(N_T);
	for(int i = N_T-1; i >= 0; --i) {
		y2[i] = ncsi(t[i]);
	}
	
	// test 3. test border case (last node)
	// make sure that the 1st AND the last node are highlighted
	// in the plot (last node: see special case handling due to
	// integer division in find_index() )
	VectorXd x3(2);
	VectorXd y3(2);
	x3 << t[0], t[N_T-1];
	for(int i = 0; i < x3.size(); ++i) {
		y3[i] = ncsi(x3[i]);
	}
	
	// test 4. test std::vector<double> (no legend item)
	vector<double> x4 = {MIN, MAX};
	vector<double> y4 = {1.0, 1.0};

	// plot results
	MatlabPlotter p;
	p.comment("Example: natural cubic splines");
	p.comment("Code generated by natcsi.cpp");
	p.figure("Natural CSI");
	p.plot_fx("exp(sin(2.*pi.*x))", MIN, MAX, N_X, "b-");
	p.hold();
	p.plot(t, y, "b*");
	p.plot(xx, yy, "r--");
	p.plot(t, y2, "ro");
	p.plot(x3, y3, "g*");
	p.plot(x4, y4, "k:");
	p.legend(5, "exp(sin(2*pi*x))", "t0 to tn", "Natural CSI", "check findindex()", "check corner cases");
	p.hold(false);
}
