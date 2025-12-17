# TEM

This is TEM, an algorithm for modeling time-domain seismic electromagnetic radiation.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

AUTHOR: Xiaodong Yang, Qinhua Yu

EMAIL: xdyang@ustc.edu.cn, yuqinhua1127@mail.ustc.edu.cn

ADDRESS: School of Earth and Space Sciences, University of Science and Technology of China, Hefei 230026, China

AIM: An algorithm that carries out time-domain seismic electromagnetic radiation modeling. The algorithm is based on the 3-D time-domain finite element method with unstructured edge-based meshes, and an unconditional stable time-stepping method is applied to acquire stable solutions. Several models are provided along with the program. The program is compatible with PC, workstation, and cluster usage with OpenMP and MPI parallelization.

# Version history

1.0     2025/12/16      initial release

# 1 Environment

Linux operating system

Install Intel OneAPI basekit and hpckit packages, refer to following website:

https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit?operatingsystem=linux?operatingsystem=linux

Set environmental variables BEFORE running example codes as:

source /opt/intel/oneapi/setvars.sh --force

Set stacksize as:

ulimit -s unlimited

# 2 Run example models

## example: Single-fault model

./TEM -f sph1
./TEM -f sph2
./TEM -f sph3
./TEM -f sph4
./TEM -f sph5
./TEM -f sph6
./TEM -f sph7
./TEM -f sph8
./TEM -f sph9
./TEM -f sph10
./TEM -f sph11
./TEM -f sph12
./TEM -f sph13
./TEM -f sph14
./TEM -f sph15
./TEM -f sph16
./TEM -f sph17
./TEM -f sph18
./TEM -f sph19

# 3 Build your own model

## 3.1 Gmsh script

Refer to https://gmsh.info/ for Gmsh scripting info.

Here is a sample Gmsh file sph17.geo that could be used by FETD program, you can find it within the package:

// set the min and max mesh size

SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 600;

// set the physical volumes with physical ID 11, 12, 13, ..., the physical domains must start with number 11, and no more than 50

Sphere(1)={0,0,0,3470};
Sphere(2)={0,0,0,5000};
Sphere(3)={0,0,0,6371};
Sphere(4)={0,0,0,6391};
Sphere(5)={0,0,0,6411};
Sphere(6)={0,0,0,6431};
Sphere(7)={0,0,0,6451};
Sphere(8)={0,0,0,6471};
f() = Abs(Boundary{ Volume{2}; });
Coherence;
Physical Volume ( "core",11 ) = {1};
Physical Volume ( "mantle1",12 ) = {2};
Physical Volume ( "mantle2",13 ) = {3};
Physical Volume ( "air1",14 ) = {4};
Physical Volume ( "air2",15 ) = {5};
Physical Volume ( "air3",16 ) = {6};
Physical Volume ( "air4",17 ) = {7};
Physical Volume ( "air5",18 ) = {8};

// set line sources, each source has its own physical ID 51, 52, 53, ..., the physical number of lines must start with number 51

pt1=newp;Point(pt1)={-34.9947,-39.9705,6337.6066};
pt2=newp;Point(pt2)={-34.1912,-40.928,6335.4415};
l1=newreg;Line(l1)={pt1,pt2};
Physical Line ("lsource1",51)={l1};
...

// this section is optional, you can set near which receiver points the mesh is going to have a smaller size than the background

pr1=newp;Point(pr1)={-129.2655,-42.2436,6369.5384};
...

// field, needs to be changed according to the model

Field[1] = Distance;
Field[1].Sampling = 1;
Field[1].CurvesList = {l1};
Field[1].PointsList = {pt1,pt2};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 1;
Field[2].SizeMax = 600;
Field[2].DistMin = 0;
Field[2].DistMax = 3000;
Field[2].StopAtDistMax = 1;
Field[3] = Distance;
Field[3].PointsList = {pr1};
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 0.1;
Field[4].SizeMax = 600;
Field[4].DistMin = 0;
Field[4].DistMax = 3000;
Field[4].StopAtDistMax = 1;
Field[5] = Distance;
Field[5].SurfacesList =  {f()};
Field[6] = Threshold;
Field[6].InField = 5;
Field[6].LcMin = 400;
Field[6].LcMax = 600;
Field[6].DistMin = 0;
Field[6].DistMax = 3000;
Field[6].StopAtDistMax = 1;
Field[7] = Min;
Field[7].FieldsList = {2,4,6};
Background Field = 7;

## 3.2 Meshing

The mesh should be saved as 2 seperate files:

1. The gmsh default .msh format file;
3. The gmsh INRIA MEDIT .mesh format file with "physical entities" selected.

## 3.3 .cntl file

This file should have the same prefix as the mesh files. Here is a sample sph17.cntl file that could be used by FETD program:

! set the receiver locations with x-, y- and z- coordinates

receivers: 35
-129.2655 -42.2436 6369.5384
...

! set the media chatacteristics, each record should contain the physical ID as the same as .geo file,

! followed by the resistivity, the relative permittivity and the relative magnetic conductivity.

isotropic: 8
11 1e-2 1.0 1.0
12 1e-2 1.0 1.0
13 1e-2 1.0 1.0
14 1.3e-12 1.0 1.0
15 1.4e-11 1.0 1.0
16 1.6e-10 1.0 1.0
17 1.7e-08 1.0 1.0
18 1.7e-05 1.0 1.0

! set source condiguration, the first line defines the number of source lines

! the following lines define the pulse type (2:Delta pulse), the pulse width of the source and the maximum current of each source

source: 176
51(physical ID) 2(should not be changed) 0(should not be changed) 0.01(t_0 as mentioned in the paper) 1(A as mentioned in the paper) -34.9947 -39.9705 6337.6066(set the sorce locations with x-, y- and z- coordinates) -34.1912 -40.928 6335.4415(set the sorce locations with x-, y- and z- coordinates)
...

! number to divide pulse width, the initial time step is thus calculated by t_0/sp_division

sp_division: 50

! the maximum time of the computation

time_maximum: 100

! number of steps to double the time step size

sp_double: 50

! n_multi as mentioned in the paper

sp_dblesize: 2

! Select the electromagnetic field components for output(1-Yes, 0-No)
! Ex, Ey, Ez, dBx, dBy, dBz, Bx, By, Bz

output:
0 0 0 0 0 0 1 1 1

## 3.4 Run and plot

./TEM -f sph1

After running the FETD program, the required field componets at receiver locations will be ouput to .field file.