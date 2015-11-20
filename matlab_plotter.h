/*  matlab_plotter.h

	A simple interface to plot C++ results by Matlab. Basically, a .m file is composed.
	This is an include-only library without dependency except the C++ standard library
	(C++11). If you are using the Eigen3 library, specific plot methods are offered for
	Eigen vectors, too.

	The 'program' is printed to stdout. Use piping to store it in a .m file.
	Do not write other output to stdout when you intend to use this method.

	Yes, it's trivial, but it fit the bill when needed.

	Feedback welcome: mailbox@pirmin-schmid.ch

	v0.3 2015-11-05 / 2015-11-20; early development, interface may change

	get latest version from: https://github.com/pirminschmid/MatlabPlotter

	(c) 2015 Pirmin Schmid, MIT License

	license text: https://github.com/pirminschmid/MatlabPlotter/tree/master/LICENSE
*/


#ifndef _MATLAB_PLOTTER_H
#define _MATLAB_PLOTTER_H

#include <iostream>
#include <string>
#include <type_traits>

class MatlabPlotter {
public:
	MatlabPlotter() {
		// currently nothing
	}

	~MatlabPlotter() {
		// currently nothing
	}


	//--- public types

	enum PlotType {
		LINEAR,
		SEMILOGY,
		SEMILOGX,
		LOGLOG
	};


	//--- public methods

	// writes a comment
	void comment(const std::string &comment) const {
		std::cout << "% " << comment << std::endl;
	}


	// starts a new figure with name title and currently set plotType
	void figure(const std::string &title) const {
		std::cout << "figure('Name','" << title << "');" << std::endl;
	}


	// starts a new figure with name title and sets plotType for all future plots
	void figure(const std::string &title, const PlotType newType) {
		plotType = newType;
		figure(title);
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


	// prints the legend. This is a variadic function using C++ variadic template functionality
	// input: one or more variables, typically strings (type std::string)
	// NOTICE: This method does not have the first parameter 'int count' anymore (see v0.1).
	//         Please adjust your code if needed. Sorry for the inconvenience.
	template <typename T>
	void legend(T s) {
		if(lps == INACTIVE) {
			// only one argument
			std::cout << "legend('" << s << "');" << std::endl;
		}
		else {
			// last of a series of arguments
			std::cout << "'" << s << "');" << std::endl;
			lps = INACTIVE;
		}
	}

	template<typename T, typename... Args>
	void legend(T first, Args... more) {
		if(lps == INACTIVE) {
			// prefix
			std::cout << "legend(";
			lps = ACTIVE;
		}

		// first / or next argument
		std::cout << "'" << first << "',";

		// handle more arguments
		legend(more...);
	}


	// plots values stored in std::vector<double> or Eigen::VectorXd vectors x and y as Matlab row vectors
	template<typename Vector>
	void plot(const Vector &x, const Vector &y, const std::string &style) const {
#ifdef EIGEN_CORE_H
		static_assert(std::is_base_of<std::vector<double>, Vector>::value ||
					  std::is_base_of<Eigen::VectorXd, Vector>::value, "Vector type mismatch! Use std::vector<double> or Eigen::VectorXd.");
#else
		static_assert(std::is_base_of<std::vector<double>, Vector>::value, "Vector type mismatch! Use std::vector<double>.");
#endif
		if(x.size() != y.size()) {
			std::cout << "ERROR - plot(): vectors x and y must have equal size." << std::endl;
			return;
		}

		print_row_vector("x", x);
		print_row_vector("y", y);
		std::cout << plotCommand() << "(x,y,'" << style << "');" << std::endl;
	}


	// prints a given std::vector<double> or Eigen::VectorXd as Matlab row vector
	template<typename Vector>
	void print_row_vector(const std::string &name, const Vector &values) const {
#ifdef EIGEN_CORE_H
		static_assert(std::is_base_of<std::vector<double>, Vector>::value ||
					  std::is_base_of<Eigen::VectorXd, Vector>::value, "Vector type mismatch! Use std::vector<double> or Eigen::VectorXd.");
#else
		static_assert(std::is_base_of<std::vector<double>, Vector>::value, "Vector type mismatch! Use std::vector<double>.");
#endif
		bool first = true;
		std::cout << name << " = [";
		int n = values.size();
		for(int i = 0; i < n; ++i) {
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


	// plots the function f using linspace for x with min, max and n=number of points
	// note: f must be described in Matlab syntax that allows vector-wise operations
	//       therefore, use .* ./ .^ instead of * / and ^, respectively
	//       you must use x as the variable name
	void plot_fx(const std::string &f, const double min, const double max, const int n, const std::string &style) const {
		std::cout << "f = @(x) " << f << ";" << std::endl;
		std::cout << "x = linspace(" << min << "," << max << "," << n << ");" << std::endl;
		std::cout << plotCommand() << "(x,f(x),'" << style << "');" << std::endl;
	}


	// prints any string text as is
	void raw(const std::string &text) const {
		std::cout << text << std::endl;
	}


	// subplot(m,n,p) command as defined by Matlab
	void subplot(const int m, const int n, const int p) const {
		std::cout << "subplot(" << m << "," << n << "," << p << ");" << std::endl;
	}


	// adds x and y axis labels
	void xylabels(const std::string &xlabel, const std::string &ylabel) const {
		std::cout << "xlabel('" << xlabel << "');" << std::endl;
		std::cout << "ylabel('" << ylabel << "');" << std::endl;
	}

private:

	// use a tiny finite state machine to help controlling the variadic template methods for legend()
	enum LegendPrintState {
		INACTIVE,
		ACTIVE
	};

	LegendPrintState lps = INACTIVE;


	// control the plotting type/scale for the figure
	PlotType plotType = LINEAR;

	std::string plotCommand() const {
		switch(plotType) {
			case LINEAR:
				return "plot";
			case SEMILOGY:
				return "semilogy";
			case SEMILOGX:
				return "semilogx";
			case LOGLOG:
				return "loglog";
			default:
				return "error_wrong_plot_type";
		}
	}

};

#endif // _MATLAB_PLOTTER_H
