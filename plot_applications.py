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

from mpmath import mp
import spin_gr_theory
float_dps = 180
mp.dps = float_dps

sp = spin_gr_theory.spin_gr_theory()

def plot_pulsar_g():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 1.76e-4
	rot2 = 1256.
	params = {'m2': 1.4 * 1.98892e30,'r2':20000.0,'rot2':rot2,\
		'r1':70000.e3/2,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':600000.e3,'e':.1,'i':0.8,'h':2}
	# Kill spins.
	params['rot2'] = 1E-30
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ein_period = sp.wp_period
	g_time = sp.g_time
	ein_rate = (g_time(ein_period) - g_time(0))/ein_period
	lspace = linspace(0,5*ein_period,250)
	subplot(1,3,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$g\left( ^\circ \right)$')
	plot(lspace,[g_time(t).real*(360/(2*mpmath.pi())) for t in lspace],'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,3,2)
	title('(b)')
	# Restore spins.
	params['rot2'] = rot2
	params['rot1'] = rot1
	sp.parameters = params
	g_time = sp.g_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	g_t_vec = array([g_time(t).real for t in lspace])
	delta_1 = g_t_vec - array([(ein_rate*t).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi())),'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	subplot(1,3,3)
	title('(c)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	rate_spin = (g_time(sp.wp_period) / sp.wp_period).real
	delta_2 = g_t_vec - array([(rate_spin*t).real for t in lspace])
	plot(lspace,delta_2*(360/(2*mpmath.pi())),'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#axial_tilt_time = sp.axial_tilt_time
	#plot(lspace,[sympy.re(axial_tilt_time(t)) for t in lspace])

def plot_pulsar_h():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 1.76e-4
	rot2 = 1256.
	params = {'m2': 1.4 * 1.98892e30,'r2':20000.0,'rot2':rot2,\
		'r1':70000.e3/2,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':600000.e3,'e':.1,'i':0.8,'h':2}
	# Kill the secondary spin.
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	lt_period = sp.wp_period
	h_time = lambda t: sp.hs_time(t) + sp.ht_time(t)
	lt_rate = (h_time(lt_period) - h_time(0))/lt_period
	lspace = linspace(0,5*lt_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$h\left( ^{\prime\prime} \right)$')
	plot(lspace,[h_time(t).real*(360/(2*mpmath.pi()))*3600 for t in lspace],'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	yticks(arange(ylim()[0],ylim()[1],400),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],400)])
	subplot(1,2,2)
	title('(b)')
	# Restore secondary spins.
	params['rot1'] = rot1
	sp.parameters = params
	h_time = lambda t: sp.hs_time(t) + sp.ht_time(t)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta h\left( ^{\prime\prime} \right)$')
	h_t_vec = array([h_time(t).real for t in lspace])
	delta_1 = h_t_vec - array([(lt_rate*t + h_time(0)).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi()))*3600,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])

def plot_pulsar_ht():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 1.76e-4
	rot2 = 1256.
	params = {'m2': 1.4 * 1.98892e30,'r2':20000.0,'rot2':rot2,\
		'r1':70000.e3/2,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':600000.e3,'e':.1,'i':0.8,'h':2}
	# Kill the primary spin.
	params['rot2'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ds_period = sp.wp_period
	ht_time = sp.ht_time
	lspace = linspace(0,5*ds_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\tilde{h}\left( ^\circ \right)$')
	ds_series = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,2,2)
	title('(b)')
	# Restore primary spins.
	params['rot2'] = rot2
	sp.parameters = params
	ht_time = sp.ht_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta \tilde{h}\left( ^\circ \right)$')
	ht_t_vec = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series - ht_t_vec,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])

def plot_pulsar_3d_spins():
	from numpy import linspace, sqrt, sin, cos, array
	from pylab import rc, subplot, title
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import axes3d, Axes3D
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	rot1 = 1.76e-4
	rot2 = 1256.
	params = {'m2': 1.4 * 1.98892e30,'r2':20000.0,'rot2':rot2,\
		'r1':70000.e3/2,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':600000.e3,'e':.1,'i':0.8,'h':2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	H_time = sp.H_time
	hs_time = sp.hs_time
	ht_time = sp.ht_time
	lspace = linspace(0,10*period,500)
	H_series = array([float(H_time(t).real) for t in lspace])
	hs_series = array([float((hs_time(t)).real) for t in lspace])
	ht_series = array([float((ht_time(t)).real) for t in lspace])
	h_series = hs_series + ht_series
	G = sqrt(6.673E-11*1.4*1.98892e30*600000.e3*(1. - .1**2))
	x_series = sqrt(G**2-H_series**2)*sin(h_series)
	y_series = -sqrt(G**2-H_series**2)*cos(h_series)
	z_series = H_series
	ax = fig.add_subplot(121,projection='3d')
	title('(a)')
	ax.plot(x_series,y_series,z_series,'k-',linewidth=2)
	ax = fig.add_subplot(122,projection='3d')
	title('(b)')
	Gt = 2./5. * (70000.e3/2)**2 * rot1
	Hts = Gt * cos(.5) + G * cos(.8)
	Gtxys = sqrt(Gt**2 - (Hts - H_series)**2)
	x_series = Gtxys*sin(ht_series)
	y_series = -Gtxys*cos(ht_series)
	z_series = Hts - H_series
	ax.plot(x_series,y_series,z_series,'k-',linewidth=2.)

def plot_pulsar_obliquity():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	rot1 = 1.76e-4
	rot2 = 1256.
	params = {'m2': 1.4 * 1.98892e30,'r2':20000.0,'rot2':rot2,\
		'r1':70000.e3/2,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':600000.e3,'e':.1,'i':0.8,'h':2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	return
	params = {'m2': 1.4 * 1.98892e30,'r2':20000.0,'rot2':rot2,\
		'r1':15000.e3/2,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':600000.e3,'e':.1,'i':0.8,'h':2}
	# Set parameters.
	sp.parameters = params
	# Period.
	ob = sp.obliquity_time
	ob_series2 = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series-ob_series2,'k-',linewidth=2)

def plot_mercury_obliquity():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	params = {'m2':1.98892e30,'r2':695500000.0,'rot2':2.972e-06,\
		'r1':2440.e3,'rot1':1.24e-06,'i_a':0.0524,'ht': 0.1,\
		'a':57909100.e3,'e':0.206,'i':0.122,'h':0.2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*5E6),[r'$'+str(int(_)*5)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*5E6))]])

def plot_psr_sgr_obliquity():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	params = {'m2':1.98892e30 * 4.3E6,'r2':44E9,'rot2':0.0019,\
		'r1':20E3,'rot1':1256,'i_a':.5,'ht':0.,\
		'a':90000000000000.,'e':.8,'i':0.8,'h':2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*2000))]])

def plot_psr_sgr_g():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 1256
	rot2 = 0.0019
	params = {'m2':1.98892e30 * 4.3E6,'r2':44E9,'rot2':rot2,\
		'r1':20E3,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':90000000000000.,'e':.8,'i':0.8,'h':2}
	# Kill spins.
	params['rot2'] = 1E-30
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ein_period = sp.wp_period
	g_time = sp.g_time
	ein_rate = (g_time(ein_period) - g_time(0))/ein_period
	lspace = linspace(0,5*ein_period,250)
	subplot(1,3,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$g\left( ^\circ \right)$')
	plot(lspace,[g_time(t).real*(360/(2*mpmath.pi())) for t in lspace],'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,3,2)
	title('(b)')
	# Restore spins.
	params['rot2'] = rot2
	params['rot1'] = rot1
	sp.parameters = params
	g_time = sp.g_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	g_t_vec = array([g_time(t).real for t in lspace])
	delta_1 = g_t_vec - array([(ein_rate*t).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	subplot(1,3,3)
	title('(c)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	rate_spin = (g_time(sp.wp_period) / sp.wp_period).real
	delta_2 = g_t_vec - array([(rate_spin*t).real for t in lspace])
	plot(lspace,delta_2*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#axial_tilt_time = sp.axial_tilt_time
	#plot(lspace,[sympy.re(axial_tilt_time(t)) for t in lspace])

def plot_psr_sgr_ht():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 1256
	rot2 = 0.0019
	params = {'m2':1.98892e30 * 4.3E6,'r2':44E9,'rot2':rot2,\
		'r1':20E3,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':90000000000000.,'e':.8,'i':0.8,'h':2}
	# Kill the primary spin.
	params['rot2'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ds_period = sp.wp_period
	ht_time = sp.ht_time
	lspace = linspace(0,5*ds_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\tilde{h}\left( ^\circ \right)$')
	ds_series = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])
	yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,2,2)
	title('(b)')
	# Restore primary spins.
	params['rot2'] = rot2
	sp.parameters = params
	ht_time = sp.ht_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\Delta \tilde{h}\left( ^\circ \right)$')
	ht_t_vec = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series - ht_t_vec,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])

def plot_psr_sgr_3d_spins():
	from numpy import linspace, sqrt, sin, cos, array
	from pylab import rc, subplot, title
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import axes3d, Axes3D
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	rot1 = 1256
	rot2 = 0.0019
	params = {'m2':1.98892e30 * 4.3E6,'r2':44E9,'rot2':rot2,\
		'r1':20E3,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':90000000000000.,'e':.8,'i':0.8,'h':2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	H_time = sp.H_time
	hs_time = sp.hs_time
	ht_time = sp.ht_time
	lspace = linspace(0,10*period,500)
	H_series = array([float(H_time(t).real) for t in lspace])
	hs_series = array([float((hs_time(t)).real) for t in lspace])
	ht_series = array([float((ht_time(t)).real) for t in lspace])
	h_series = hs_series + ht_series
	G = sqrt(6.673E-11*1.98892e30 * 4.3E6 *90000000000000.*(1. - .8**2))
	x_series = sqrt(G**2-H_series**2)*sin(h_series)
	y_series = -sqrt(G**2-H_series**2)*cos(h_series)
	z_series = H_series
	ax = fig.add_subplot(121,projection='3d')
	title('(a)')
	ax.plot(x_series,y_series,z_series,'k-',linewidth=2)
	ax = fig.add_subplot(122,projection='3d')
	title('(b)')
	Gt = 2./5. * (20E3)**2 * rot1
	Hts = Gt * cos(.5) + G * cos(.8)
	Gtxys = sqrt(Gt**2 - (Hts - H_series)**2)
	x_series = Gtxys*sin(ht_series)
	y_series = -Gtxys*cos(ht_series)
	z_series = Hts - H_series
	ax.plot(x_series,y_series,z_series,'k-',linewidth=2.)

def plot_psr_sgr_h():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 1256
	rot2 = 0.0019
	params = {'m2':1.98892e30 * 4.3E6,'r2':44E9,'rot2':rot2,\
		'r1':20E3,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':90000000000000.,'e':.8,'i':0.8,'h':2}
	# Kill the secondary spin.
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	lt_period = sp.wp_period
	h_time = lambda t: sp.hs_time(t) + sp.ht_time(t)
	lt_rate = (h_time(lt_period) - h_time(0))/lt_period
	lspace = linspace(0,5*lt_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$h\left( ^\circ \right)$')
	plot(lspace,[h_time(t).real*(360/(2*mpmath.pi())) for t in lspace],'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])
	#yticks(arange(ylim()[0],ylim()[1],400),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],400)])
	subplot(1,2,2)
	title('(b)')
	# Restore secondary spins.
	params['rot1'] = rot1
	sp.parameters = params
	h_time = lambda t: sp.hs_time(t) + sp.ht_time(t)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\Delta h\left( ^{\prime\prime} \right)$')
	h_t_vec = array([h_time(t).real for t in lspace])
	delta_1 = h_t_vec - array([(lt_rate*t + h_time(0)).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi()))*3600,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])

def psr_sgr_h():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 1256
	rot2 = 0.0019
	params = {'m2':1.98892e30 * 4.3E6,'r2':44E9,'rot2':rot2,\
		'r1':20E3,'rot1':rot1,'i_a':.5,'ht':0.,\
		'a':90000000000000.,'e':.8,'i':0.8,'h':2}
	# Kill the secondary spin.
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	# Restore secondary spins.
	params['rot1'] = rot1
	sp.parameters = params
	return
	lt_period = sp.wp_period
	h_time = lambda t: sp.hs_time(t) + sp.ht_time(t)
	lt_rate = (h_time(lt_period) - h_time(0))/lt_period
	lspace = linspace(0,5*lt_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$h\left( ^\circ \right)$')
	plot(lspace,[h_time(t).real*(360/(2*mpmath.pi())) for t in lspace],'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])
	#yticks(arange(ylim()[0],ylim()[1],400),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],400)])
	subplot(1,2,2)
	title('(b)')
	# Restore secondary spins.
	params['rot1'] = rot1
	sp.parameters = params
	h_time = lambda t: sp.hs_time(t) + sp.ht_time(t)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\Delta h\left( ^{\prime\prime} \right)$')
	h_t_vec = array([h_time(t).real for t in lspace])
	delta_1 = h_t_vec - array([(lt_rate*t + h_time(0)).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi()))*3600,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])

def plot_probe_obliquity():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	params = {'m2':1.98892e30,'r2':695500000.0,'rot2':2.972e-06,\
		'r1':1,'rot1':6.28*100,'i_a':0.5,'ht': 0,\
		'a':57909100.e3 / 2000,'e':0.01,'i':0.8,'h':2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*2000))]])

def plot_probe_jupiter_obliquity():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	#params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
	#	'r1':1,'rot1':6.28*100,'i_a':0.5,'ht': 0,\
	#	'a':421700E3 / 4,'e':0.01,'i':0.8,'h':2}
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
		'r1':0.038,'rot1':6.28*3500/60,'i_a':0.5,'ht': 0,\
		'a':421700E3 / 4,'e':0.01,'i':0.8,'h':2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*10000),[r'$'+str(int(_)*10)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*10000))]])

def plot_probe_jupiter_ht():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 6.28*3500/60
	rot2 = 0.00017585125448028
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':rot2,\
		'r1':0.038,'rot1':rot1,'i_a':0.5,'ht': 0,\
		'a':421700E3 / 4,'e':0.01,'i':0.8,'h':2}
	# Kill the primary spin.
	params['rot2'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ds_period = sp.wp_period
	ht_time = sp.ht_time
	lspace = linspace(0,5*ds_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\tilde{h}\left( ^\circ \right)$')
	ds_series = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series,'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])
	#yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,2,2)
	title('(b)')
	# Restore primary spins.
	params['rot2'] = rot2
	sp.parameters = params
	ht_time = sp.ht_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\Delta \tilde{h}\left( ^\circ \right)$')
	ht_t_vec = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series - ht_t_vec,'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])

def plot_probe_earth_obliquity():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	#params = {'m2':5.97E24,'r2':6371.e3,'rot2':6.28/(24.*3600),\
		#'r1':0.038,'rot1':6.28*3500/60,'i_a':0.5,'ht': 0,\
		#'a':7027E3,'e':0.01,'i':0.8,'h':2}
	params = {'m2':5.97E24,'r2':6371.e3,'rot2':6.28/(24.*3600),\
		'r1':0.038,'rot1':6.28*3400/60,'i_a':90.007 * 2*mpmath.pi()/360,'ht': 0,\
		'a':7013E3,'e':0.0014,'i':mpmath.pi()/2,'h':-mpmath.pi()/2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*2000))]])

def plot_probe_earth_ht():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 6.28*3400/60
	rot2 = 6.28/(24.*3600)
	params = {'m2':5.97E24,'r2':6371.e3,'rot2':rot2,\
		'r1':0.038,'rot1':rot1,'i_a':90.007 * 2*mpmath.pi()/360,'ht': 0,\
		'a':7013E3,'e':0.0014,'i':mpmath.pi()/2,'h':-mpmath.pi()/2}
	# Kill the primary spin.
	params['rot2'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ds_period = sp.wp_period
	ht_time = sp.ht_time
	lspace = linspace(0,ds_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\tilde{h}\left( ^\circ \right)$')
	ds_series = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series,'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])
	#yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,2,2)
	title('(b)')
	# Restore primary spins.
	params['rot2'] = rot2
	sp.parameters = params
	ht_time = sp.ht_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\Delta \tilde{h}\left( ^\circ \right)$')
	ht_t_vec = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ht_t_vec,'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20000),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20000))]])

def plot_earth_moon_ht():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 7.978282018665196e-08
	rot2 = 2.972e-06
	params = {'m2':1.98892e30,'r2':695500000.0,'rot2':rot2,\
		'r1':384748.e3,'rot1':rot1,'i_a':10. * 2*mpmath.pi()/360,'ht': 0,\
		'a':149597870700.,'e':0.017,'i':7. * 2*mpmath.pi()/360,'h':mpmath.pi()/2}
	# Kill the primary spin.
	params['rot2'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ds_period = sp.wp_period
	ht_time = sp.ht_time
	lspace = linspace(0,5*ds_period,250)
	subplot(1,2,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$\tilde{h}\left( ^\circ \right)$')
	ds_series = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*100000000),[r'$'+str(int(_)*100)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*100000000))]])
	yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,2,2)
	title('(b)')
	# Restore primary spins.
	params['rot2'] = rot2
	sp.parameters = params
	ht_time = sp.ht_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$\Delta \tilde{h}\left( ^\circ \right)$')
	ht_t_vec = array([ht_time(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ds_series - ht_t_vec,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*100000000),[r'$'+str(int(_)*100)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*100000000))]])

def plot_earth_moon_obliquity():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	#params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
	#	'r1':1,'rot1':6.28*100,'i_a':0.5,'ht': 0,\
	#	'a':421700E3 / 4,'e':0.01,'i':0.8,'h':2}
	rot1 = 7.978282018665196e-08
	rot2 = 2.972e-06
	params = {'m2':1.98892e30,'r2':695500000.0,'rot2':rot2,\
		'r1':384748.e3,'rot1':rot1,'i_a':10. * 2*mpmath.pi()/360,'ht': 0,\
		'a':149597870700. * 0.1,'e':0.017,'i':7. * 2*mpmath.pi()/360,'h':mpmath.pi()/2}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*100000000),[r'$'+str(int(_)*100)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*100000000))]])

def plot_io_obliquity_1():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	#params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
	#	'r1':1,'rot1':6.28*100,'i_a':0.5,'ht': 0,\
	#	'a':421700E3 / 4,'e':0.01,'i':0.8,'h':2}
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
		'r1':1821E3,'rot1':6.28/(42 * 3600),'i_a':0.017453,'ht': 0,\
		'a':421700E3,'e':0.0041,'i':0.008726646,'h':0.35}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*1000000),[r'$'+str(int(_)*1)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*1000000))]])

def plot_io_obliquity_2():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	#params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
	#	'r1':1,'rot1':6.28*100,'i_a':0.5,'ht': 0,\
	#	'a':421700E3 / 4,'e':0.01,'i':0.8,'h':2}
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
		'r1':1821E3,'rot1':6.28/(21 * 3600),'i_a':0.17453,'ht': 0,\
		'a':421700E3,'e':0.05,'i':0.349,'h':0.35}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*1000000),[r'$'+str(int(_)*1)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*1000000))]])

def plot_io_g_1():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 6.28/(42 * 3600)
	rot2 = 0.00017585125448028
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':rot2,\
		'r1':1821E3,'rot1':rot1,'i_a':0.017453,'ht': 0,\
		'a':421700E3,'e':0.0041,'i':0.008726646,'h':0.35}
	# Kill spins.
	params['rot2'] = 1E-30
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ein_period = sp.wp_period
	g_time = sp.g_time
	ein_rate = (g_time(ein_period) - g_time(0))/ein_period
	lspace = linspace(0,5*ein_period,250)
	subplot(1,3,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$g\left( ^\circ \right)$')
	plot(lspace,[g_time(t).real*(360/(2*mpmath.pi())) for t in lspace],'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,3,2)
	title('(b)')
	# Restore spins.
	params['rot2'] = rot2
	params['rot1'] = rot1
	sp.parameters = params
	g_time = sp.g_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	g_t_vec = array([g_time(t).real for t in lspace])
	delta_1 = g_t_vec - array([(ein_rate*t).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	subplot(1,3,3)
	title('(c)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	rate_spin = (g_time(sp.wp_period) / sp.wp_period).real
	delta_2 = g_t_vec - array([(rate_spin*t).real for t in lspace])
	plot(lspace,delta_2*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#axial_tilt_time = sp.axial_tilt_time
	#plot(lspace,[sympy.re(axial_tilt_time(t)) for t in lspace])

def plot_io_g_2():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 6.28/(21 * 3600)
	rot2 = 0.00017585125448028
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':rot2,\
		'r1':1821E3,'rot1':rot1,'i_a':0.17453,'ht': 0,\
		'a':421700E3,'e':0.05,'i':0.349,'h':0.35}
	# Kill spins.
	params['rot2'] = 1E-30
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ein_period = sp.wp_period
	g_time = sp.g_time
	ein_rate = (g_time(ein_period) - g_time(0))/ein_period
	lspace = linspace(0,5*ein_period,250)
	subplot(1,3,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$g\left( ^\circ \right)$')
	plot(lspace,[g_time(t).real*(360/(2*mpmath.pi())) for t in lspace],'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,3,2)
	title('(b)')
	# Restore spins.
	params['rot2'] = rot2
	params['rot1'] = rot1
	sp.parameters = params
	g_time = sp.g_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	g_t_vec = array([g_time(t).real for t in lspace])
	delta_1 = g_t_vec - array([(ein_rate*t).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	subplot(1,3,3)
	title('(c)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	rate_spin = (g_time(sp.wp_period) / sp.wp_period).real
	delta_2 = g_t_vec - array([(rate_spin*t).real for t in lspace])
	plot(lspace,delta_2*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#axial_tilt_time = sp.axial_tilt_time
	#plot(lspace,[sympy.re(axial_tilt_time(t)) for t in lspace])

def plot_metis_obliquity_1():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
		'r1':20E3,'rot1':6.28/(7 * 3600),'i_a':0.017453,'ht': 0,\
		'a':128000E3,'e':0.0002,'i':0.008726646,'h':0.35}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*10000),[r'$'+str(int(_)*10)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*10000))]])

def plot_metis_obliquity_2():
	from numpy import linspace, array, arange
	from pylab import rc, plot, xlim, xticks, xlabel, ylabel
	import matplotlib.pyplot as plt
	import mpmath
	rc('text', usetex=True)
	fig = plt.figure()
	# Initial params.
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':0.00017585125448028,\
		'r1':20E3,'rot1':6.28/(3.5 * 3600),'i_a':0.17453,'ht': 0,\
		'a':128000E3,'e':0.05,'i':0.349,'h':0.35}
	# Set parameters.
	sp.parameters = params
	# Period.
	period = sp.wp_period
	ob = sp.obliquity_time
	lspace = linspace(0,5*period,250)
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (ka)}$')
	ylabel(r'$\textnormal{Obliquity }\left( ^\circ \right)$')
	ob_series = array([ob(t).real*(360/(2*mpmath.pi())) for t in lspace])
	plot(lspace,ob_series,'k-',linewidth=2)
	xticks(arange(xlim()[0],xlim()[1],365*24*3600*10000),[r'$'+str(int(_)*10)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*10000))]])

def plot_metis_g_1():
	from pylab import plot, subplot, title, rc, xlim, ylim, xlabel, ylabel, xticks, yticks, title
	from numpy import linspace, arange, array
	import mpmath
	rc('text', usetex=True)
	# Initial params.
	rot1 = 6.28/(7 * 3600)
	rot2 = 0.00017585125448028
	params = {'m2':1.9e27,'r2':70000.e3,'rot2':rot2,\
		'r1':20E3,'rot1':rot1,'i_a':0.017453,'ht': 0,\
		'a':128000E3,'e':0.0002,'i':0.008726646,'h':0.35}
	# Kill spins.
	params['rot2'] = 1E-30
	params['rot1'] = 1E-30
	# Set parameters.
	sp.parameters = params
	ein_period = sp.wp_period
	g_time = sp.g_time
	ein_rate = (g_time(ein_period) - g_time(0))/ein_period
	lspace = linspace(0,5*ein_period,250)
	subplot(1,3,1)
	title('(a)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (Ma)}$')
	ylabel(r'$g\left( ^\circ \right)$')
	plot(lspace,[g_time(t).real*(360/(2*mpmath.pi())) for t in lspace],'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#yticks(arange(ylim()[0],ylim()[1],360),[r'$'+str(int(_))+r'$' for _ in arange(ylim()[0],ylim()[1],360)])
	subplot(1,3,2)
	title('(b)')
	# Restore spins.
	params['rot2'] = rot2
	params['rot1'] = rot1
	sp.parameters = params
	g_time = sp.g_time
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	g_t_vec = array([g_time(t).real for t in lspace])
	delta_1 = g_t_vec - array([(ein_rate*t).real for t in lspace])
	plot(lspace,delta_1*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	subplot(1,3,3)
	title('(c)')
	xlim(0,float(lspace[-1]))
	xlabel(r'$\textnormal{Time (years)}$')
	ylabel(r'$\Delta g\left( ^\circ \right)$')
	rate_spin = (g_time(sp.wp_period) / sp.wp_period).real
	delta_2 = g_t_vec - array([(rate_spin*t).real for t in lspace])
	plot(lspace,delta_2*(360/(2*mpmath.pi())),'k-',linewidth=2)
	#xticks(arange(xlim()[0],xlim()[1],365*24*3600*20),[r'$'+str(int(_)*20)+r'$' for _ in [t[0] for t in enumerate(arange(xlim()[0],xlim()[1],365*24*3600*20))]])
	#axial_tilt_time = sp.axial_tilt_time
	#plot(lspace,[sympy.re(axial_tilt_time(t)) for t in lspace])
