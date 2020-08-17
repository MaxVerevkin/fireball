# About
This program tries to reconstruct the trajectory and position of a fireball given observer's data.

Program reads the data from a text file. The default is 'data.txt'. It can be overriden by passing a file name as an argument.

Output on current test data looks something like that:
```
Data is initialized

Summary on finding flash position:
    Total square-error (rad):  0.000288
    Mean square-error  (rad):  0.000026
    Standard error     (deg):  0.293309

Summary on finding flash trajectory:
    Total square-error (rad):  0.000217
    Mean square-error  (rad):  0.000013
    Standard error     (deg):  0.204691

Answer:
  Location: 35.23692 33.46764 17045.49347
  Velocity: -11978.64533 19317.14237 -26395.93569

                  i     z0        h0        zb        hb         A         t
Processed Answer (1):  50.545 |  13.511 |  82.920 |  40.218 | 235.950 | 1.691119
Processed Answer (2):  51.986 |  39.925 | 129.963 |  54.496 | 287.929 | 2.708481
Processed Answer (3):  34.009 |  12.603 |  78.009 |  48.076 | 240.950 | 2.286615
Processed Answer (4): 340.296 |  14.349 | 359.328 |  67.966 | 202.680 | 1.684388
Processed Answer (5):  83.010 |  14.260 | 113.954 |  39.090 | 236.253 | 2.589706
Processed Answer (6): 104.465 |  17.006 | 124.280 |  36.741 | 226.857 | 2.042091

Ignore: 'altitude end'   for observer 3
Ignore: 'desent_angle'   for observer 5


```
# Compilation
```
$ git clone git@github.com:MaxVerevkin/fireball.git
$ cd fireball && make
```
Optionally, to enable parallelism (OpenMP), make with (remember to `make clean`):
```
$ make paralel
```
