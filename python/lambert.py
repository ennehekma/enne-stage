# text_file = open("../build/example.txt", "r")
# lines = text_file.readlines()
# print lines
# # print len(lines)
# text_file.close()

from numpy import loadtxt
import matplotlib.pyplot as plt
import numpy as np

test1 = np.arange(0,5,1)
# test1 = np.arange(5,6,1)
for j in (test1):
	filename = "../data/dv_" + str(j) + ".txt"
	print filename
	counter = loadtxt("../data/i.txt", comments="#", delimiter=" ", unpack=False)
	dV = loadtxt("../data/dv_0.txt", comments="#", delimiter=" ", unpack=False)
	print dV
	dvmin = np.amin(dV)
	print dvmin

	plt.plot(counter, dV)
	plt.axis([0,15,0,20000])
	# plt.axis([np.amin(toff)*.9, np.amax(toff)*1.1, np.amin(dvtotarray)*.9, np.amax(dvtotarray)*1.1])
	plt.xlabel('ToF [orbital periods]')
	plt.ylabel('DV [m/s]')
	plt.title('dV needed as function of Time of Flight')
	plt.show()