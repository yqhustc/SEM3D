# SA

This is SA, an algorithm used for inverting finite fault electromagnetic sources.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

AUTHOR: Qinhua Yu, Xiaodong Yang

EMAIL: yuqinhua1127@mail.ustc.edu.cn, xdyang@ustc.edu.cn

ADDRESS: School of Earth and Space Sciences, University of Science and Technology of China, Hefei 230026, China

AIM: An algorithm that carries out finite fault electromagnetic sources inversion. The algorithm is based on the 3-D time-domain finite element method with unstructured edge-based meshes, and an unconditional stable time-stepping method is applied to acquire stable solutions. Several models are provided along with the program. The program is compatible with PC, workstation, and cluster usage with OpenMP and MPI parallelization.

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

./SA -i sph1

# 3 Build your own model

## 3.1 Data preparation

Use TEM to calculate the reference magnetic field, and then copy the folder of all forward simulation results to 'example/forward_data/mod1'

Prepare the observation data in 'example/forward_data/event1/field_true.field'. Every three columns represent the three components of the magnetic field at a measurement point, corresponding to the 'field' file in the forward modeling results.

## 3.2 .inv file

This file is located in the 'example/inv_file'. Here is a sample 'sph1.inv' file that could be used by FETD program:

! Number of faults

fault_num: 1

! Name of the folder where the reference magnetic field is located

forward_data:
mod1

! Name of the folder where the observation data is located

inv_data: 
event1

! Number of sub faults in the strike and dip directions

src_number:
8 22

! The number of azimuth angles

degree_range: 19
110(minimum angle) 10(interval) 290(maximum angle)

! Number of receivers

rec_number: 35

! The number of lines in a 'field' file

result_time_length: 711

! Time interval of observation data

interpolation_time:
0(start time) 1(interval) 250(end time)

! Initial temperature for simulated annealing

initial_temperatrue: 0.1

! Annealing coefficient

annealing_coefficient: 0.985

! Termination temperature

termination_temperature: 1e-07

! The number of iterations at each temperature

move: 200

! Control the range of variation of strength coefficient

intensity_range_factor: 1.6e6

! Weight of model constraint terms

smooth_constraint_weight_a: 0.4
smooth_constraint_weight_t: 0.4
entropy_constraint_weight_a: 0.1
single_constraint_weight: 0.1
sparse_constraint_weight: 0.1

## 3.3 Run and plot

./SA -i sph1

After running the SA program, the inversion result is output to the 'src_inv_f1.txt' file located in the 'inv_desult/sh1_result' folder.

You can use 'result_plot.m' to draw the inversion results.