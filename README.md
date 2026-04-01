Makefile wiil create an executable file named "sitnikov". Run it with ./sitnikov, also try ./sitnikov -h for more information. To plot Poincare surface use ./plot.py or just python plot.py, also try -h flag.
This project uses minIni parser (https://github.com/compuphase/minIni/tree/master) for .ini files.
Runge-Kutta 4th order integrator is used.

Config.ini requirements: 
1. Integration parameters section, where n is the number of points on Poincare map and m is the number of steps in one period 2*PI.
[integration] 
n = 200
m = 1000
2. Eccentricity and initial values. There could be several of them problem sections, just use word "problem" in the beggining of section name. ("problem1", "problem2" are legit). 
[problem...]
e = 0.1
v = 0.1
r = 0.0
*You can also change the maximum amount of problem sections in minIni header files (i suppose), but i dont think you'll need to.


