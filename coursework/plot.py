import math
import pylab
import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

DPI = 200



if __name__ == '__main__':
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	x = numpy.loadtxt("space_tgrid_2.x.txt")
	y = numpy.loadtxt("space_tgrid_2.y.txt")
	z = numpy.loadtxt("space_tgrid_2.data.txt")

	X, Y = np.meshgrid(x, y)

	fig, ax = plt.subplots()
	#locator=ticker.LogLocator(base=math.pow(10, 1/10000))
	cs = ax.contourf(X, Y, z, 55, cmap=cm.coolwarm)
	cbar = fig.colorbar(cs)

	plt.title(r'Graph', fontsize=19)
	plt.xlabel(r'$t_x$', fontsize=10)
	plt.ylabel(r'$t_y$', fontsize=10)
	plt.tick_params(axis='both', labelsize=8)
	plt.grid(alpha=0.5)
	plt.savefig('pic.png', dpi=DPI)
	plt.clf()