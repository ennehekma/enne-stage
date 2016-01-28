# text_file = open("../build/example.txt", "r")
# lines = text_file.readlines()
# print lines
# # print len(lines)
# text_file.close()
import math as math
import numpy as np
import math as math
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import pi
from PyKEP import epoch, DAY2SEC, AU, MU_SUN, lambert_problem, epoch_from_string
from PyKEP.planet import jpl_lp
from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler

from numpy import loadtxt
import matplotlib.pyplot as plt
import numpy as np



r1 = np.array([6778136,0,0])
# r2 = np.array([-6978136,0.00007,0])
r2 = np.array([0,6043243.047,3489068])
v0 = np.array([0,7668.558733,0])
# v3 = np.array([0.000000008,-7557.86574,0])
v3 = np.array([-7557.86574,0,0])


r1a = math.sqrt(r1[0]**2+r1[1]**2+r1[2]**2)

MU_EARTH = 3.986004418e14

dt = 3.1415926535* math.sqrt((r1a+100000)**3/MU_EARTH) # seconds
l = lambert_problem(r1, r2, dt*.5, MU_EARTH)
fig2 = plt.figure(2)
axis2 = fig2.gca(projection='3d')
axis2.scatter([0], [0], [0], color='y') 
plot_lambert(l, sol=0, ax=axis2, color='r')
# plot_lambert(l, sol=1, ax=axis2, color='r')
x0 = l.get_x()[0]
plot_kepler(r1, v0, dt*2, MU_EARTH, N=600, units=1, color='b',legend=False, ax=axis2)
plot_kepler(r2, v3, dt*2.2, MU_EARTH, N=600, units=1, color='b',legend=False, ax=axis2)
axis2.set_ylim3d(-1.2*6378000,1.2*6378000)
axis2.set_xlim3d(-1.2*6378000,1.2*6378000)
plt.xlabel('X coordinate [m]')
plt.ylabel('Y coordinate [m]')
# plt.zlabel('Z coordinate [m]')
# plt.ylabel('DV [m/s]')
plt.title('Lambert transfers for 200km altitude increase with different time of flight')
plt.show()

# dt = 3*3.1415926535* math.sqrt((r1a+100000)**3/MU_EARTH) # seconds
l = lambert_problem(r1, r2, dt*2.4, MU_EARTH)
fig2 = plt.figure(2)
axis2 = fig2.gca(projection='3d')
axis2.scatter([0], [0], [0], color='y') 
plot_lambert(l, sol=1, ax=axis2, color='r')
# plot_lambert(l, sol=1, ax=axis2, color='r')
x0 = l.get_x()[0]
plot_kepler(r1, v0, dt*2, MU_EARTH, N=600, units=1, color='b',legend=False, ax=axis2)
plot_kepler(r2, v3, dt*2.2, MU_EARTH, N=600, units=1, color='b',legend=False, ax=axis2)
axis2.set_ylim3d(-1.2*6378000,1.2*6378000)
axis2.set_xlim3d(-1.2*6378000,1.2*6378000)
plt.xlabel('X coordinate [m]')
plt.ylabel('Y coordinate [m]')
# plt.zlabel('Z coordinate [m]')
# plt.ylabel('DV [m/s]')
plt.title('Lambert transfers for 200km altitude increase with different time of flight')
plt.show()


# arr = np.array([[1,2,3,4], [5,6,7,8], [9,10,11,12]])
# 	print arr
# 	arr = np.delete(arr, 1, 0)
# 	print arr