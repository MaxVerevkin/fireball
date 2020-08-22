# About
This program tries to reconstruct the trajectory and position of a fireball given observer's data.

Program reads the data from a text file. The default is 'data.txt'. It can be overriden by passing a file name as an argument.

Output on current test data looks something like that:
```
Data is initialized

Summary on finding flash position:
    Total square-error         :  0.000221 rad
    Standard deviation of meaan:  0.081284°

Summary on finding flash trajectory:
    Total square-error:           0.000173 rad
    Standard deviation of meaan:  0.045724°

Answer:
    lat =  35.23706°
    lon =  33.46793°
    z   =  17.06100(km)
    vx  =  -11.95412(km/s)
    vy  =  19.24591(km/s)
    vz  =  -26.26384(km/s)

                  i     z0        h0        zb        hb         A         t   
Processed Answer (1):  50.548 |  13.517 |  82.915 |  40.174 | 235.992 | 1.697540
Processed Answer (2):  51.996 |  39.907 | 129.947 |  54.442 | 287.936 | 2.725168
Processed Answer (3):  34.020 |  12.610 |  78.006 |  48.021 | 240.982 | 2.294969
Processed Answer (4): 340.325 |  14.362 | 359.376 |  67.856 | 202.734 | 1.687220
Processed Answer (5):  82.998 |  14.267 | 113.965 |  39.067 | 236.315 | 2.604759
Processed Answer (6): 104.441 |  17.014 | 124.263 |  36.717 | 226.916 | 2.052218

Ignore: 'altitude end'   for observer 3
Ignore: 'desent_angle'   for observer 5
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
