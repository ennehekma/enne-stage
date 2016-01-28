import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

data = genfromtxt('../build/file.txt', unpack=True)
dvmin = np.amin(data)
dvminint = math.floor(dvmin)
dvsize = math.sqrt(np.size(data))
for cutoff in [dvmin+1000,dvmin+20,dvmin+1]:
	for i in xrange(0,540):
		for j in xrange(0,540):
			if data[i][j]>cutoff:
				data[i][j]=cutoff


	fig, (ax2) = plt.subplots(1, 1, figsize=[16, 13])
	plt.xlabel('ToF [orbital periods of lower orbit]')
	plt.ylabel('DV [m/s]')
	plt.title('Cutoff dV: '+str(cutoff)+' [m/s]')
	
	axins = inset_axes(ax2,
	                   width="5%",  # width = 10% of parent_bbox width
	                   height="100%",  # height : 50%
	                   loc=3,
	                   bbox_to_anchor=(1.05, 0., 1, 1),
	                   bbox_transform=ax2.transAxes,
	                   borderpad=0,
	                   )
	im = ax2.imshow(data, cmap=plt.get_cmap('hot'))
	# plt.colorbar(im, cax=axins, ticks=np.arange(-5000,25000,5000))
	plt.colorbar(im, cax=axins, ticks=np.arange(dvmin-0.00001,(cutoff-dvminint)*1.01+dvminint,(cutoff-dvminint)/5))
	plt.draw()
	plt.show()


