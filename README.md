MatlabPlotter
=============

A simple C++ class that writes Matlab .m programs as output of calculations e.g. using the Eigen C++ library.

This may be useful if you do not want to add additional dependencies to your program or write the plotting code by yourself, but you know that all involved parties have Matlab installed on their computers.

Current version v0.2 (2015-11-16). See [changelog][changelog] for details.


Usage
-----

Include matlab_plotter.h into your C++ program. Use methods as shown in the [example natcsi.cpp][example] ([folder][folder]).

Note: the example needs [Eigen][eigen] to be installed on your system. You may use cmake for building.
However, the plotter can be used independently of it.


License
-------

Copyright (c) 2015 Pirmin Schmid, [MIT license][license].


[changelog]:https://github.com/pirminschmid/MatlabPlotter/tree/master/CHANGELOG.md
[example]:https://github.com/pirminschmid/MatlabPlotter/tree/master/example/natcsi.cpp
[folder]:https://github.com/pirminschmid/MatlabPlotter/tree/master/example
[eigen]:http://eigen.tuxfamily.org
[license]:https://github.com/pirminschmid/MatlabPlotter/tree/master/LICENSE