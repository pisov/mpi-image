1. Compile with Intel fortran compiler

module load ifort

ifort serial.f90 fio.f90 -o img.x

2. Compile with gfortran compiler

or

make

gfortran -O3 serial.f90 fio.f90 -o img.x

3. Execute

./img.x

4. View the output on Linux

display image.pgm

5. View remotely from Windows machine

convert image.pgm image.png

6. Transfer image.png to local Windows machine

