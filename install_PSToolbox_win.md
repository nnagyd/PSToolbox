# Installing and using PSToolbox on Windows

## Install
First, you need a couple of utility programs:

* Use Chocolatex for instaaling Unix/Mac-like programs: https://chocolatey.org/

* Next, install 'make' command by typing: ``choco install make``

* You'll also need MinGW (Minimalist Gnu for Windows): ``choco install mingw``

* We will also need the eigen package: ``choco install eigen``


Next you need to clone the repository from GitHub. Once you're done, you need to tailor the Windows makefile   ``makefile_win.mk``:

* Specify the path to Eigen, e.g. ``INC = -IC:/ProgramData/chocolatey/lib/eigen/include``

* On Windows systems, you use the ``ar rcs`` utility to build static libraries, e.g. ``ar rcs libmy_tools.a my_tools.o`` (instead of ``libtool -static -o libmy_tools.a my_tools.o``). This is already set in the sample ``makefile_win.mk`` file.

At this point, you should be able to compile the toolbox apart from the classes depending on the ``Coolprop`` library. 

## Using the library

Navigate to the ````

