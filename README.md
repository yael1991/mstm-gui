This software is designed to be used with CMake, so we recommend using it to build the application.

The MSTM application requires a Fortran compiler and OpenMP.

In addition, the GUI requires:

Python -- http://www.python.org/

Qt4 (user interface library) -- http://qt-project.org/

PyQt4 (Python library for Qt4) -- http://wiki.python.org/moin/PyQt4

matlibplot (Plotting library for Python) -- http://matplotlib.org/

These should be available easily through most Linux package managers and can be downloaded individually for Windows systems.  If they are installed through Linux, CMake should be able to find them without any problems.  Windows users may have to point CMake to the installed locations.