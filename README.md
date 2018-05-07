# PHY905_Projects
This is the repository for projects of PHY905: Computational Physics.

In each project, there are bascially three repositorys:

In "Report", my report gives a comprehensive summary of the project.

In "src", source codes and the makefile are given. Type "make" in command line to compile codes.

In "benchmark", some output files are stored for reference.

2018.3.5 Update for Project2

Executbales in "src" include armdiallo diagonalization; classical & cyclic Jacobi; bisection method; Jacobi and armadillo for one-electron system; and armadillo for two-electron system. Last executable require several inputs as size n, rmax and omega. Also, two unit tests folder are included, they are tests for armadillo and classical Jacobi method.

Partial resluts are in "benchmark", I don't want to include too much to cause chaos.

2018.4.1 Update for Project3

Executbales in "src" include the Earth-Sun system, Earth-Sun-Jupiter system and precission of Mercury. The Euler forward and velocity-Verlet method are implemented in the Earth-Sun system. The rests use the VV method. The first one can also be used to calculate the escape velocity. Execitables start form t =0 and require an input of t_end.
Total number of steps need to be adjusted in codes.

Test results for the Euler and VV-method are in "benmark".

Some plots of orbits and data are included in "plot".

See final summarized report in "report".

2018.5.7 Update for Project4

Executbales in "src" include "vmc.x" for the first trail wave function and "vmc2.x" for the second trail wave function. In "vmc.x", the omega and step is required for input. Only omega is required for "vmc2.x", so step need to be adjusted by hand.

Test results for both (omega=1) are in "benmark".

See final summarized report in "report".

