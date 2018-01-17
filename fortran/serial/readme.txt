# Compile with Intel fortran compiler

module load ifort

ifort serial.f90 fio.f90 -o img.x

# Compile with gfortran compiler

gfortran -O3 serial.f90 fio.f90 -o img.x

# Execute

./img.x

# View the output on Linux

display image.pgm

# View remotely from Windows machine

convert image.pgm image.png

# Transfer image.png to local Windows machine

