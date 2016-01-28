# text_file = open("../build/example.txt", "r")
# lines = text_file.readlines()
# print lines
# # print len(lines)
# text_file.close()
import math as math
import numpy as np
# import math as math
# import matplotlib
# import matplotlib.cm as cm
# import matplotlib.mlab as mlab
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from scipy.constants import pi
# from PyKEP import epoch, DAY2SEC, AU, MU_SUN, lambert_problem, epoch_from_string
# from PyKEP.planet import jpl_lp
# from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler

from numpy import loadtxt
import matplotlib.pyplot as plt
# import numpy as np
mylegend = ['0 revolutions','1 revolution','2 revolutions','3 revolutions','4 revolutions','5 revolutions']

# test1 = [4,3,2,1,0]#0,5,1)
test1 = np.arange(0,5,1)
for j in (test1):
	filename = "../data/dv_" + str(j) + ".txt"
	# print filename
	counter = loadtxt("../data/i.txt", comments="#", delimiter=" ", unpack=False)
	dV = loadtxt(filename, comments="#", delimiter=" ", unpack=False)
	# for i in np.size(counter):cd
	# 	np.array(i) = dV[i]
	# print dV, np.size(dV)
	# print counter, np.size(counter)
	dvmin = np.amin(dV)
	# print dvmin
	total = np.ndarray(shape=(np.size(dV),2), dtype=float, order='F')
	# total[1,1] = 5
	for i in np.arange(0,np.size(dV)):
		total[i,0] = dV[i]
		total[i,1] = counter[i]
	# itemindex = 1
	itemindex = np.where(total[:,0]<100000)
	indexs  = np.amin(itemindex)
	total = np.delete(total, np.arange(0,indexs,1), 0)
	plt.plot(total[:,1], total[:,0],label=mylegend[j])
# plt.axis([0, np.size(counter)*0.3/10000, 0, 2500])#np.amax(dV)*1.1])
plt.axis([0, np.size(counter)/10000, 0, 14000])#np.amax(dV)*1.1])
plt.xlabel('ToF [orbital periods of lower orbit]')
plt.ylabel('DV [m/s]')
plt.title('dV needed to transfer from 400km to 600km altitude as function of Time of Flight')
plt.legend(loc='lower right')
plt.xticks(np.arange(0.5, 10.5, 1.0))
plt.show()

