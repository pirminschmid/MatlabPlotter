
CHANGELOG
=========

v0.4.1 (2015-11-27) to v0.4.3 (2015-11-28)
-----------------------------------------
- minor changes in example2
- improved README
- added automatic testing by travis-ci

v0.4 (2015-11-23)
-----------------
- added title()
- removed restrictive type assertion at the beginning of plot() and print_row_vector().
  If these functions are used with wrong input type, the error messages will be more cryptic
  but this gives more flexibility for additional types that can use these templates
- new example2: Legendre polynomials and Gauss-Legendre quadrature  

v0.3 (2015-11-20)
-----------------
- plot() rewritten as templated method (no API change)
- print_row_vector() is now public
- added subplot()

v0.2 (2015-11-16)
-----------------
- added semilog and loglog plots
- added axis labels
- std::strings
- legend() changed to variadic template style (notice: API change)

v0.1 (2015-11-05)
-----------------
- Hello world!
- first draft -- warning: API may/will change
- figures
- plotting for C++ and Eigen vectors
- plot function in given interval
- classic style variadic function signature used for legend()
- thus: C style strings (char *)
- first example: Natural cubic spline interplolation