/*  matlab_plotter.h

	A simple interface to plot C++ results by Matlab. Basically, a .m file is composed.
	This is an include-only library without dependency except the C++ standard library.
	If you are using the Eigen3 library, specific plot methods are offered for
	Eigen vectors, too.

	The 'program' is printed to stdout. Use piping to store it in a .m file.
	Do not write other output to stdout when you intend to use this method.

	Yes, it's trivial, but it fit the bill when needed.

	Feedback welcome: mailbox@pirmin-schmid.ch

	v0.1 2015-11-05, first draft; interface may change

	(c) 2015 Pirmin Schmid, MIT License

	license text: https://github.com/pirminschmid/MatlabPlotter/tree/master/LICENSE
*/


#ifndef _MATLAB_PLOTTER_H
#define _MATLAB_PLOTTER_H


#include <cstdarg>
#include <iostream>


class MatlabPlotter {
public:
	MatlabPlotter() {
		// currently nothing
	}

	~MatlabPlotter() {
		// currently nothing
	}

	void comment(const char *comment) const {
		std::cout << "% " << comment << std::endl;
	}

	// starts a new figure with name title
	void figure(const char *title) const {
		std::cout << "figure('Name','" << title << "');" << std::endl;
	}

	// writes the Matlab commands hold on or hold off
	void hold(const bool on = true) const {
		if(on) {
			std::cout << "hold on;" << std::endl;
		}
		else {
			std::cout << "hold off;" << std::endl;
		}
	}


	// prints the legend. This is a variadic function
	// input:  count number of strings afterwards
	//         strings (type char *) 
	void legend(const int count, ...) const {
		std::cout << "legend(";
		va_list args;
		va_start(args, count);
		for(int i = 0; i < count; ++i) {
			if(i > 0) {
				std::cout << ",";
			}
			char *item = va_arg(args, char *);
			std::cout << "'" << item << "'";
		}
		va_end(args);
		std::cout << ");" << std::endl;
	}

	// plots values stored in standard library vectors x and y using Matlab style description style
	void plot(const std::vector<double> &x, const std::vector<double> &y, const char *style) const {
		if(x.size() != y.size()) {
			std::cout << "ERROR - plot(): vectors x and y must have equal size." << std::endl;
			return;
		}

		print_row_vector("x", x);
		print_row_vector("y", y);
		std::cout << "plot(x,y,'" << style << "');" << std::endl;
	}


#ifdef EIGEN_CORE_H
	// conditional compilation in case the Eigen library is used

	// plots values stored in Eigen3 vectors x and y using Matlab style description style
	// This function is only available if this include detects Eigen to be loaded before
	void plot(const Eigen::VectorXd &x, const Eigen::VectorXd &y, const char *style) const {
		if(x.size() != y.size()) {
			std::cout << "ERROR - plot(): vectors x and y must have equal size." << std::endl;
			return;
		}

		print_row_vector("x", x);
		print_row_vector("y", y);
		std::cout << "plot(x,y,'" << style << "');" << std::endl;
	}

#endif // EIGEN_CORE_H	


	// plots the function f using linspace for x with min, max and n=number of points
	// note: f must be described in Matlab syntax that allows vector-wise operations
	//       therefore, use .* ./ .^ instead of * / and ^, respectively
	//       you must use x as the variable name
	void plot_fx(const char *f, const double min, const double max, const int n, const char *style) const {
		std::cout << "f = @(x) " << f << ";" << std::endl;
		std::cout << "x = linspace(" << min << "," << max << "," << n << ");" << std::endl;
		std::cout << "plot(x,f(x),'" << style << "');" << std::endl;
	}

	// prints any string text as is
	void raw(const char *text) const {
		std::cout << text << std::endl;
	}

private:
	// prints a given std::vector as Matlab row vector
	void print_row_vector(const char *name, const std::vector<double> &values) const {
		bool first = true;
		std::cout << name << " = [";
		for(double v : values) {
			if(first) {
				first = false;
				std::cout << v;
			}
			else {
				std::cout << ", " << v;
			}
		}
		std::cout << "];" << std::endl;		
	}


#ifdef EIGEN_CORE_H
	// conditional compilation in case the Eigen library is used

	// prints a given std::vector as Matlab row vector
	void print_row_vector(const char *name, const Eigen::VectorXd &values) const {
		bool first = true;
		std::cout << name << " = [";
		for(int i=0; i < values.size(); ++i) {
			if(first) {
				first = false;
				std::cout << values[i];
			}
			else {
				std::cout << ", " << values[i];
			}
		}
		std::cout << "];" << std::endl;
	}

#endif // EIGEN_CORE_H

};

#endif // _MATLAB_PLOTTER_H
