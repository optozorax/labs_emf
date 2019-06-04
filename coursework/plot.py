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

def make_image(xpath, ypath, zpath, resultpath, mytitle):
	x = numpy.loadtxt(xpath)
	y = numpy.loadtxt(ypath)
	z = numpy.loadtxt(zpath)	

	X, Y = np.meshgrid(x, y)

	fig, ax = plt.subplots()
	#locator=ticker.LogLocator(base=math.pow(10, 1/10000))
	cs = ax.contourf(X, Y, z, 55, cmap=cm.coolwarm)
	cbar = fig.colorbar(cs)

	plt.title(mytitle, fontsize=19)
	plt.xlabel(r'$t_x$', fontsize=15)
	plt.ylabel(r'$t_y$', fontsize=15)
	plt.tick_params(axis='both', labelsize=10)
	plt.grid(alpha=0.25)
	plt.savefig(resultpath, dpi=DPI)
	plt.clf()

def make_images(i):
	make_image(f"space_tgrid_{i}.x.txt", f"space_tgrid_{i}.y.txt", f"space_tgrid_{i}.integral.txt", f"space_tgrid_{i}_integral.png", r"Integral of functions difference");
	make_image(f"space_tgrid_{i}.x.txt", f"space_tgrid_{i}.y.txt", f"space_tgrid_{i}.norm.txt", f"space_tgrid_{i}_norm.png", r"Norm of $q$ vectors difference");
	make_image(f"space_tgrid_{i}.x.txt", f"space_tgrid_{i}.y.txt", f"space_tgrid_{i}.time.txt", f"space_tgrid_{i}_time.png", r"Solving time");

if __name__ == '__main__':
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	make_images(0)
	make_images(1)
	make_images(2)