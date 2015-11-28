MatlabPlotter
=============
[![Build Status](https://travis-ci.org/pirminschmid/MatlabPlotter.svg?branch=master)](https://travis-ci.org/pirminschmid/MatlabPlotter)

A simple C++ class that writes Matlab .m programs as output of calculations e.g. using the Eigen C++ library. 

This may be useful if you do not want to add additional dependencies to your program or write the plotting code by yourself, but you know that all involved parties have Matlab installed on their computers.

Note: Due to the use of templates, any container type (e.g. vector, array,...) of any library can be used as input vector for plot() and print\_row\_vector(), as long as the container supports the standard methods size() and operator\[\]\(\), and the embedded scalar type can be printed by cout << and is suitable for plotting in Matlab. 

Current version v0.4 (2015-11-23). See [changelog][changelog] for details. [Feedback][feedback] welcome.


Usage
-----

Include matlab\_plotter.h into your C++ program. It needs the standard library of C++11. Use methods as shown in the examples [natcsi.cpp][example] ([folder][folder]) and [legendre.cpp][example2] ([folder][folder2]).

Note: the examples need [Eigen][eigen] to be installed on your system. You may use cmake for building.
However, the plotter can be used independently of it.


License
-------

Copyright (c) 2015 Pirmin Schmid, [MIT license][license].


[changelog]:https://github.com/pirminschmid/MatlabPlotter/tree/master/CHANGELOG.md
[example]:https://github.com/pirminschmid/MatlabPlotter/tree/master/example/natcsi.cpp
[folder]:https://github.com/pirminschmid/MatlabPlotter/tree/master/example
[example2]:https://github.com/pirminschmid/MatlabPlotter/tree/master/example2/legendre.cpp
[folder2]:https://github.com/pirminschmid/MatlabPlotter/tree/master/example2
[eigen]:http://eigen.tuxfamily.org
[license]:https://github.com/pirminschmid/MatlabPlotter/tree/master/LICENSE
[feedback]:mailto:mailbox@pirmin-schmid.ch?subject=MatlabPlotter
