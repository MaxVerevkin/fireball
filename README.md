# About
This program tries to reconstruct the trajectory and position of a fireball given observer's data.

Program reads the data from a text file. The default is 'data.txt'. It can be overriden by passing a file name as an argument.

Output on current test data looks something like that (longitude's sigma and velocity is broken for now):
```
Data is initialized

Summary on finding flash position:
    Total square-error         :  0.002992 rad
    Standard deviation of meaan:  0.572177°

Summary on finding flash trajectory:
    Total square-error:           9.209590 rad
    Standard deviation of meaan: 10.542857°

Answer:
  Location:
    lat =   35.2444 ± 0.0135(°)
    lon =   33.4917 ±    inf(°)
    z   =   17.7726 ± 0.3933(km)
  Velocity:  105.093(km/s)
    Local:
      v_East  =   -37.928(km/s)
      v_North =    30.632(km/s)
      v_z     =   -93.100(km/s)
    Global:
      v_x =  -57.223(km/s)
      v_y =  -83.342(km/s)
      v_z =  -28.708(km/s)

                  i     z0°       h0°       zb°       hb°        A°        t(s)    w(°/s)
Processed Answer (1):  51.173 |  14.624 |  96.655 |  32.743 | 257.371 |   2.649 |  16.986
Processed Answer (2):  53.802 |  39.345 | 113.170 |  36.282 | 292.666 |   2.855 |  16.166
Processed Answer (3):  35.082 |  13.657 |  95.883 |  32.536 | 266.229 |   2.360 |  24.670
Processed Answer (4): 342.532 |  16.220 | 107.915 | -10.831 | 282.243 |   0.712 | 178.111
Processed Answer (5):  82.501 |  14.947 |  97.725 |  31.329 | 223.818 |   2.855 |   7.531
Processed Answer (6): 103.068 |  17.311 | 102.304 |  31.254 | 176.972 |   2.855 |   4.890

Ignore: 'altitude end'   for observer 3
Ignore: 'desent_angle'   for observer 6
Total ingored: 2
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
