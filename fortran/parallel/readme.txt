# Compile with MPI Fortran compiler

module load openmpi

mpif90 parallel.f90 fio.f90 -o img.x


# Execute on N = 1, 2, 4, or 8 processes watch execution time

mpirun -n N img.x

# View the output on Linux

display image.pgm

# View remotely from Windows machine

convert image.pgm image.png

# Transfer image.png to local Windows machine

# Change the input image with n = 344 and m = 600

