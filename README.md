# About
This program tries to reconstruct the trajectory and position of a fireball given observer's data.

For now, data is defined in `data_values.h`. Soon the program will read it from a file.

Output on current test data looks something like that:
```
Data is initialized

Summary on finding flash position:
    Total square-error (rad):  0.000288
    Mean square-error  (rad):  0.000026
    Standard error     (deg):  0.293309

Summary on finding flash trajectory:
    Total square-error (rad):  0.000225
    Mean square-error  (rad):  0.000013
    Standard error     (deg):  0.208667

Answer: 35.23692 33.46764 17045.49347 27999.82563 37920.18937 -754.43496

                  i     z0        h0        zb        hb         A         t
Processed Answer (1):  50.545 |  13.511 |  82.855 |  40.188 | 235.904 | 1.246094
Processed Answer (2):  51.986 |  39.925 | 129.961 |  54.509 | 287.918 | 1.996094
Processed Answer (3):  34.009 |  12.603 |  78.067 |  48.107 | 240.987 | 1.691406
Processed Answer (4): 340.296 |  14.349 | 359.325 |  68.030 | 202.660 | 1.246094
Processed Answer (5):  83.010 |  14.260 | 113.871 |  39.041 | 236.208 | 1.902344
Processed Answer (6): 104.465 |  17.006 | 124.264 |  36.717 | 226.866 | 1.503906

Ignore: 'altitude end'   for observer 3
Ignore: 'desent_angle'   for observer 5
```
# Compilation
```
$ git clone git@github.com:MaxVerevkin/fireball.git
$ make
```
Optionally, to enable parallelism (OpenMP), make with (remember to `make clean`):
```
$ make paralel
```
