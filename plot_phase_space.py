# -*- coding: iso-8859-1 -*-
# Copyright (C) 2014 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

def plot_phase_line(a,b,c,d,e):
	from mpmath import polyroots, sqrt, polyval
	from numpy import linspace
	from pylab import plot, xlim, xticks, yticks, grid, xlabel, ylabel
	assert(a < 0)
	p4roots, err = polyroots([a,b,c,d,e],error = True, maxsteps = 1000)
	real_roots = filter(lambda x: x.imag == 0,p4roots)
	assert(len(real_roots) == 2 or len(real_roots) == 4)
	print(real_roots)
	# Left and right padding for the plot.
	lr_pad = (real_roots[-1] - real_roots[0]) / 10
	# This is the plotting function.
	def func(x):
		retval = sqrt(polyval([a,b,c,d,e],x))
		if retval.imag == 0:
			return retval
		else:
			return float('nan')
	func_m = lambda x: -func(x)
	def plot_lobe(start,end,f):
		delta = (end - start) / 100
		rng = linspace(start + delta,end - delta,1000)
		plot(rng,[f(x) for x in rng],'k-',linewidth=2)
		rng = linspace(start,start + delta,1000)
		plot(rng,[f(x) for x in rng],'k-',linewidth=2)
		rng = linspace(end - delta,end,1000)
		plot(rng,[f(x) for x in rng],'k-',linewidth=2)
	if len(real_roots) == 2:
		plot_lobe(real_roots[0],real_roots[1],func)
		plot_lobe(real_roots[0],real_roots[1],func_m)
	else:
		plot_lobe(real_roots[0],real_roots[1],func)
		plot_lobe(real_roots[0],real_roots[1],func_m)
		plot_lobe(real_roots[2],real_roots[3],func)
		plot_lobe(real_roots[2],real_roots[3],func_m)
	xlim(float(real_roots[0] - lr_pad),float(real_roots[-1] + lr_pad))
	xticks([0],[''])
	yticks([0],[''])
	grid()
	xlabel(r'$H$')
	ylabel(r'$dH/dt$')

def plot_all():
	from pylab import subplot, title, rc
	rc('text', usetex=True)
	subplot(121)
	title(r'(a)')
	plot_phase_line(-2,-5,2,9,4)
	subplot(122)
	title(r'(b)')
	plot_phase_line(-1,1,2,2,4)
