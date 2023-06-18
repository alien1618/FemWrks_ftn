----------------------------------------------------------------------------------
FemWrks
----------------------------------------------------------------------------------
FemWrks is a 2D/3D finite element solver written in modern fortran
and is used to solve the general transport equation, incompressible
navier-stokes equations and can be used for stress/deformation analysis.
It can use first and second order elements and can solve the finite
element equations using element-by-element technique without the 
need to form global matrices.
---------------------------------------------------------------------------------
FemWrks has been tested to run on Debian and FreeBSD systems. To run FemWrks perform
the following steps:

0. If you haven't compiled the code in your system, compile it by:

chmod u+x mk.sh run.sh genvid.sh

./mk.sh

or if using FreeBSD use:
bash mk.sh

This step can be skipped once code is compiled.

1. Inside the FemWrks directory create a folder called "sim"

2. Copy a case study from the "examples" folder and rename folder to "in" and place it inside the "sim" directory

3. edit the input parameters for the simulation as needed for all the files inside "in" folder and run the simulation by:

chmod u+x run.sh 

./run.sh

or if using FreeBSD use:
bash run.sh

4. To use the 2D python plotter:

python3 pltmsh.py <scale val> <show mesh flag 0/1> <variable 1> variable 2> 

for example:
python3 pltmsh.py 5 0 u phi

3D data is best viewed by opening the generated vtk files in Paraview or VisIt or Salome-platform.

5. To export the generated simulation pics to a video run script:

chmod u+x genvid.sh

./genvid.sh <variable name>

for example:
./genvid.sh u
