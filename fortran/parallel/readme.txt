1. Compile with MPI Fortran compiler

module load openmpi

mpif90 parallel.f90 fio.f90 -o img.x

or

make

2. Execute on N = 1, 2, 4, or 8 processes watch execution time

mpirun -n N img.x

3. View the output on Linux

display image.pgm

4. View remotely from Windows machine

convert image.pgm image.png

5. Transfer image.png to local Windows machine

Task: Change the input image with n = 344 and m = 600

