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

from pyranha import *

pt = types.poisson_series(types.polynomial(types.rational,types.short))()

# Convenience function to substitute 1/r for its expression via f.
def r_to_f(arg):
	mu = pt(r'\mathcal{G}') * pt('m_2')
	return arg.subs('r',pt('r_1')**-1).subs('r_1',mu*pt('G')**-2*(1+pt('e')*math.cos(pt('f'))))

# Function to register custom derivatives. These are the implicit functions:
# f = f(L,G,l)
# E = E(L,G,l)
# e = e(L,G)
# Gxy = Gxy(G,H)
# Gtxy = Gtxy(Gt,Ht)
# Gtxys = Gtxys(Gt,Hts,H)
# NOTE: r is always replaced by its expression via f, so we can avoid its derivatives
# (provided that r_to_f is always called before invoking partials)
def register_derivatives():
	pt.unregister_all_custom_derivatives()
	G,e,L,H,r,E = [pt(s) for s in 'GeLHrE']
	mu = pt(r'\mathcal{G}') * pt('m_2')
	Gxy = pt('G_{xy}')
	Gtxy = pt(r'\tilde{G}_{xy}')
	Gtxys = pt('\\tilde{G}_{xy\\ast}')
	Gt = pt(r'\tilde{G}')
	Ht = pt(r'\tilde{H}')
	Hts = pt('\\tilde{H}_\\ast')
	e_L = G**2*(e*L**3)**-1
	e_G = -G*(e*L**2)**-1
	f_L = r_to_f(e_L * math.sin(E) * L**3 * (r**2*mu*G)**-1 * (r+G**2*mu**-1))
	f_G = r_to_f(e_G * math.sin(E) * L**3 * (r**2*mu*G)**-1 * (r+G**2*mu**-1))
	f_l = r_to_f(G * L**3 * (mu**2 * r**2)**-1)
	E_L = r_to_f(L**2 * (mu*r)**-1 * e_L * math.sin(E))
	E_G = r_to_f(L**2 * (mu*r)**-1 * e_G * math.sin(E))
	E_l = r_to_f(L**2 * (mu*r)**-1)
	Gxy_G = G * Gxy**-1
	Gxy_H = -H * Gxy**-1
	Gtxy_Gt = Gt * Gtxy**-1
	Gtxy_Ht = -Ht * Gtxy**-1
	Gtxys_H = Gtxys**-1 * (Hts - H)
	Gtxys_Gt = Gtxys**-1 * Gt
	Gtxys_Hts = Gtxys**-1 * (H - Hts)
	pt.register_custom_derivative('L',lambda p: p.partial('L') + p.partial('f')*f_L + p.partial('e')*e_L +\
		p.partial('E')*E_L)
	pt.register_custom_derivative('G',lambda p: p.partial('G') + p.partial('G_{xy}')*Gxy_G +\
		p.partial('f')*f_G + p.partial('E')*E_G + p.partial('e')*e_G)
	pt.register_custom_derivative('H',lambda p: p.partial('H') + p.partial('G_{xy}')*Gxy_H + p.partial('\\tilde{G}_{xy\\ast}')*Gtxys_H)
	pt.register_custom_derivative(r'\tilde{G}',lambda p: p.partial(r'\tilde{G}') + p.partial(r'\tilde{G}_{xy}')*Gtxy_Gt +\
		p.partial('\\tilde{G}_{xy\\ast}')*Gtxys_Gt)
	pt.register_custom_derivative(r'\tilde{H}',lambda p: p.partial(r'\tilde{H}') + p.partial(r'\tilde{G}_{xy}')*Gtxy_Ht)
	pt.register_custom_derivative('\\tilde{H}_\\ast',lambda p: p.partial('\\tilde{H}_\\ast') + p.partial('\\tilde{G}_{xy\\ast}')*Gtxys_Hts)
	pt.register_custom_derivative('l',lambda p: p.partial('l') + p.partial('f')*f_l + p.partial('E')*E_l)

register_derivatives()

import sympy as __sympy
class integral_func(__sympy.Function):
	nargs = 6
	_wp_cache = {}
	@classmethod
	def eval(cls,t,g2,g3,alpha,beta,gamma):
		# Leave unevaluated.
		return

class VV_func(object):
	_wp_cache = {}
	@classmethod
	def eval(cls,t,g2,g3,alpha,beta,gamma):
		import weierstrass_elliptic
		from mpmath import mpc, pi, mpf, mp
		# Check that everything is mpf.
		assert(all([isinstance(x,mpf) for x in [t,g2,g3,alpha,beta,gamma]]))
		invs = (g2,g3)
		if invs in integral_func._wp_cache:
			wp = integral_func._wp_cache[invs]
		else:
			wp = weierstrass_elliptic.weierstrass_elliptic(g2,g3)
			integral_func._wp_cache[invs] = wp
		v = wp.Pinv(beta)
		retval = alpha/wp.Pprime(v) * (mpf(2)*t*wp.zeta(v) + mpc(0,pi()) + wp.ln_sigma(-t+v) - wp.ln_sigma(t+v)) + gamma*t
		return retval

class spin_gr_theory(object):
	# Default parameters: Mercury-like planet around Sun-like star.
	# NOTE: the rotation rate is the equatorial one for the sun. The inclination value for the planet's axis ('i_a')
	# is really kinda picked randomly to be small. Same for the h and ht angles. Everything is in SI units and radians.
	__default_params = {'m2':1.98892e30,'r2':695500000.0,'rot2':2.972e-06,\
		'r1':2440.e3,'rot1':1.24e-06,'i_a':0.0524,'ht': 0.1,\
		'a':57909100.e3,'e':0.206,'i':0.122,'h':0.2}
	@staticmethod
	def __get_g_time(d_eval,t_r,angle_in):
		import mpmath, solutions as s
		assert(all([isinstance(d_eval[_],mpmath.mpf) for _ in d_eval]))
		angle_in = mpmath.mpf(angle_in)
		g2,g3 = d_eval['g_2'],d_eval['g_3']
		def retval(t):
			t = mpmath.mpf(t)
			return angle_in + s.Phig0(d_eval) * t + s.Phig1(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abgg1(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abgg1(d_eval))) \
				+ s.Phig2(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abgg2(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abgg2(d_eval))) \
				+ s.Phig3(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abgg3(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abgg3(d_eval))) \
				+ s.Phig4(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abgg4(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abgg4(d_eval)))
		return retval
	@staticmethod
	def __get_hs_time(d_eval,t_r,angle_in):
		import mpmath, solutions as s
		assert(all([isinstance(d_eval[_],mpmath.mpf) for _ in d_eval]))
		angle_in = mpmath.mpf(angle_in)
		g2,g3 = d_eval['g_2'],d_eval['g_3']
		def retval(t):
			t = mpmath.mpf(t)
			return angle_in + s.Phihs0(d_eval) * t + s.Phihs1(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abghs1(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abghs1(d_eval))) \
				+ s.Phihs2(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abghs2(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abghs2(d_eval))) \
				+ s.Phihs3(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abghs3(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abghs3(d_eval))) \
				+ s.Phihs4(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abghs4(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abghs4(d_eval))) \
				+ s.Phihs5(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abghs5(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abghs5(d_eval)))
		return retval
	@staticmethod
	def __get_ht_time(d_eval,t_r,angle_in):
		import mpmath, solutions as s
		assert(all([isinstance(d_eval[_],mpmath.mpf) for _ in d_eval]))
		angle_in = mpmath.mpf(angle_in)
		g2,g3 = d_eval['g_2'],d_eval['g_3']
		def retval(t):
			t = mpmath.mpf(t)
			return angle_in + s.Phiht0(d_eval) * t + s.Phiht1(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abght1(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abght1(d_eval))) \
				+ s.Phiht2(d_eval) * (VV_func.eval(t-t_r,g2,g3,*s.abght2(d_eval)) - VV_func.eval(-t_r,g2,g3,*s.abght2(d_eval)))
		return retval
	def __verify_solutions(self,eval_dict):
		import solutions as s
		import sympy
		d_eval = dict([(k,sympy.Symbol(k)) for k in eval_dict])
		assert(s.Phig0(d_eval) == self.g_sol[0][-1][1])
		assert(s.Phig1(d_eval) == self.g_sol[0][0][1])
		assert(s.Phig2(d_eval) == self.g_sol[0][1][1])
		assert(s.Phig3(d_eval) == self.g_sol[0][2][1])
		assert(s.Phig4(d_eval) == self.g_sol[0][3][1])
		assert([_.ratsimp() for _ in s.abgg1(d_eval)] == [_.ratsimp() for _ in self.g_sol[1][0][1].args[3:]])
		assert([_.ratsimp() for _ in s.abgg2(d_eval)] == [_.ratsimp() for _ in self.g_sol[1][1][1].args[3:]])
		assert([_.ratsimp() for _ in s.abgg3(d_eval)] == [_.ratsimp() for _ in self.g_sol[1][2][1].args[3:]])
		assert([_.ratsimp() for _ in s.abgg4(d_eval)] == [_.ratsimp() for _ in self.g_sol[1][3][1].args[3:]])
		assert(s.Phihs0(d_eval) == self.hs_sol[0][-1][1])
		assert(s.Phihs1(d_eval) == self.hs_sol[0][0][1])
		assert(s.Phihs2(d_eval) == self.hs_sol[0][1][1])
		assert(s.Phihs3(d_eval) == self.hs_sol[0][2][1])
		assert(s.Phihs4(d_eval) == self.hs_sol[0][3][1])
		assert(s.Phihs5(d_eval) == self.hs_sol[0][4][1])
		assert([_.ratsimp() for _ in s.abghs1(d_eval)] == [_.ratsimp() for _ in self.hs_sol[1][0][1].args[3:]])
		assert([_.ratsimp() for _ in s.abghs2(d_eval)] == [_.ratsimp() for _ in self.hs_sol[1][1][1].args[3:]])
		assert([_.ratsimp() for _ in s.abghs3(d_eval)] == [_.ratsimp() for _ in self.hs_sol[1][2][1].args[3:]])
		assert([_.ratsimp() for _ in s.abghs4(d_eval)] == [_.ratsimp() for _ in self.hs_sol[1][3][1].args[3:]])
		assert([_.ratsimp() for _ in s.abghs5(d_eval)] == [_.ratsimp() for _ in self.hs_sol[1][4][1].args[3:]])
		assert(s.Phiht0(d_eval) == self.ht_sol[0][-1][1])
		assert(s.Phiht1(d_eval) == self.ht_sol[0][0][1])
		assert(s.Phiht2(d_eval) == self.ht_sol[0][1][1])
		assert([_.ratsimp() for _ in s.abght1(d_eval)] == [_.ratsimp() for _ in self.ht_sol[1][0][1].args[3:]])
		assert([_.ratsimp() for _ in s.abght2(d_eval)] == [_.ratsimp() for _ in self.ht_sol[1][1][1].args[3:]])
	def __reachable_root(self,p4roots,H_in):
		# NOTE: this will not work close to equilibirium points, floating-point and/or polyroots accuracy
		# limits might place the roots so that no interval contains the initial value of H.
		from mpmath import mpf, mpc
		# Some checks on the roots. We are expecting either all real roots,
		# or two real roots and two complex.
		# TODO lots of special casing here, multiple roots, etc.
		if all([isinstance(r,mpf) for r in p4roots]):
			case = 0
		elif all([isinstance(r,mpf) for r in p4roots[0:2]]) and all([isinstance(r,mpc) for r in p4roots[2:4]]):
			case = 1
		else:
			raise ValueError('invalid root value(s)')
		if case == 0:
			slist = sorted(p4roots)
			if H_in >= slist[0] and H_in <= slist[1]:
				return slist[0],slist[1],2,0
			if H_in >= slist[2] and H_in <= slist[3]:
				return slist[2],slist[3],2,1
			raise ValueError('initial H is not compatible with real polynomial roots')
		else:
			slist = sorted(p4roots[0:2])
			if H_in >= slist[0] and H_in <= slist[1]:
				return slist[0],slist[1],1,0
			raise ValueError('initial H is not compatible with polynomial roots')
	def __compute_t_r(self,n_lobes,lobe_idx,H_in,Hr,d_eval,p4roots,lead_cf):
		from pyranha import math
		from mpmath import asin, sqrt, ellipf, mpf
		assert(n_lobes == 1 or n_lobes == 2)
		C = -lead_cf
		assert(C > 0)
		# First determine if we are integrating in the upper or lower plane.
		# NOTE: we don't care about eps, we are only interested in the sign.
		if (self.__F1 * math.sin(pt('h_\\ast'))).trim().evaluate(d_eval) > 0:
			sign = 1
		else:
			sign = -1
		if n_lobes == 2:
			assert(lobe_idx == 0 or lobe_idx == 1)
			r0,r1,r2,r3 = p4roots
			# k is the same in both cases.
			k = sqrt(((r3-r2)*(r1-r0))/((r3-r1)*(r2-r0)))
			if lobe_idx == 0:
				assert(Hr == r0)
				phi = asin(sqrt(((r3-r1)*(H_in-r0))/((r1-r0)*(r3-H_in))))
			else:
				assert(Hr == r2)
				phi = asin(sqrt(((r3-r1)*(H_in-r2))/((H_in-r1)*(r3-r2))))
			return -sign * mpf(2) / sqrt(C * (r3 - r1) * (r2 - r0)) * ellipf(phi,k**2)
		else:
			# TODO: single lobe case.
			assert(False)
			pass
	def __set_params(self,d):
		from copy import deepcopy
		from mpmath import mpf, sqrt, polyroots, cos, acos, pi, mp
		from weierstrass_elliptic import weierstrass_elliptic as we
		names = ['m2','r2','rot2','r1','rot1','i_a','ht','a','e','i','h']
		if not all([s in d for s in names]):
			raise ValueError('invalid set of parameters')
		# Convert all the values to mpf and introduce convenience shortcuts.
		m2, r2, rot2, r1, rot1, i_a, ht_in, a, e, i, h_in = [mpf(d[s]) for s in names]
		L_in = sqrt(self.__GG_val * m2 * a)
		G_in = L_in * sqrt(1. - e**2)
		H_in = G_in * cos(i)
		Gt_in = (2 * r1**2 * rot1) / 5
		Ht_in = Gt_in * cos(i_a)
		Hts_in = H_in + Ht_in
		hs_in = h_in - ht_in
		Gxy_in = sqrt(G_in**2 - H_in**2)
		Gtxys_in = sqrt(Gt_in**2 - Ht_in**2)
		J2 = (2 * m2 * r2**2 * rot2) / 5
		II_1 = mpf(5) / (2 * r1**2)
		# Evaluation dictionary.
		eval_names = ['L','G','H','\\tilde{G}','\\tilde{H}_\\ast','h_\\ast','\\tilde{G}_{xy\\ast}','m_2','\\mathcal{G}','J_2','G_{xy}',\
			'\\varepsilon','\\mathcal{I}_1']
		eval_values = [L_in,G_in,H_in,Gt_in,Hts_in,hs_in,Gtxys_in,m2,self.__GG_val,J2,Gxy_in,self.__eps_val,II_1]
		d_eval = dict(zip(eval_names,eval_values))
		# Evaluate Hamiltonian with initial conditions.
		HHp_val = self.__HHp.trim().evaluate(d_eval)
		# Add the value of the Hamiltonian to the eval dictionary.
		d_eval['\\mathcal{H}^\\prime'] = HHp_val
		# Evaluate g2 and g3.
		g2_val, g3_val = self.__g2.trim().evaluate(d_eval), self.__g3.trim().evaluate(d_eval)
		# Create the Weierstrass object.
		wp = we(g2_val,g3_val)
		# Store the period.
		self.__wp_period = wp.periods[0]
		# Now let's find the roots of the quartic polynomial.
		# NOTE: in theory here we could use the exact solution for the quartic.
		p4coeffs = [t[0] * t[1].trim().evaluate(d_eval) for t in zip([1,4,6,4,1],self.__f4_cf)]
		p4roots, err = polyroots(p4coeffs,error = True, maxsteps = 1000)
		# Find a reachable root.
		Hr, H_max, n_lobes, lobe_idx = self.__reachable_root(p4roots,H_in)
		# Determine t_r
		t_r = self.__compute_t_r(n_lobes,lobe_idx,H_in,Hr,d_eval,p4roots,p4coeffs[0])
		# Now evaluate the derivatives of the polynomial. We will need to replace H_in with Hr in the eval dict.
		d_eval['H'] = Hr
		_, f4Hp, f4Hpp, _, _ = self.__f4
		f4p_eval = f4Hp.trim().evaluate(d_eval)
		f4pp_eval = f4Hpp.trim().evaluate(d_eval)
		# Build and store the expression for H(t).
		self.__H_time = lambda t: Hr + f4p_eval / (4 * (wp.P(t - t_r) - f4pp_eval / 24))
		# H will not be needed any more, replace with H_r
		del d_eval['H']
		d_eval['H_r'] = Hr
		# Inject the invariants and the other two constants into the evaluation dictionary.
		d_eval['g_2'] = g2_val
		d_eval['g_3'] = g3_val
		d_eval['A'] = f4p_eval / 4
		d_eval['B'] = f4pp_eval / 24
		# Verify the formulae in solutions.py
		self.__verify_solutions(d_eval)
		# Assuming g = 0 as initial angle.
		self.__g_time = spin_gr_theory.__get_g_time(d_eval,t_r,0.)
		self.__hs_time = spin_gr_theory.__get_hs_time(d_eval,t_r,hs_in)
		self.__ht_time = spin_gr_theory.__get_ht_time(d_eval,t_r,ht_in)
		def obliquity(t):
			from mpmath import acos, cos, sqrt
			H = self.__H_time(t)
			hs = self.__hs_time(t)
			G = d_eval['G']
			Gt = d_eval['\\tilde{G}']
			Hts = d_eval['\\tilde{H}_\\ast']
			Gxy = sqrt(G**2-H**2)
			Gtxys = sqrt(Gt**2-(Hts-H)**2)
			retval = (Gxy*Gtxys*cos(hs)+H*(Hts-H))/(G*Gt)
			return acos(retval)
		self.__obliquity_time = obliquity
		def spin_vector(t):
			import numpy
			ht = self.__ht_time(t)
			H = self.__H_time(t)
			G = d_eval['G']
			Gt = d_eval['\\tilde{G}']
			Hts = d_eval['\\tilde{H}_\\ast']
			Gtxys = sqrt(Gt**2-(Hts-H)**2)
			return numpy.array([x.real for x in [Gtxys*sin(ht),-Gtxys*cos(ht),Hts-H]])
		self.__spin_vector_time = spin_vector
		def orbit_vector(t):
			import numpy
			ht = self.__ht_time(t)
			hs = self.__hs_time(t)
			h = hs + ht
			H = self.__H_time(t)
			G = d_eval['G']
			Gxy = sqrt(G**2-H**2)
			return numpy.array([x.real for x in [Gxy*sin(h),-Gxy*cos(h),H]])
		self.__orbit_vector_time = orbit_vector
		# Store the params of the system.
		self.__params = dict(zip(names,[mpf(d[s]) for s in names]))
		# Final report.
		rad_conv = 360 / (2 * pi())
		print("\x1b[31mAccuracy in the identification of the poly roots:\x1b[0m")
		print(err)
		print("\n\x1b[31mPeriod (years):\x1b[0m")
		print(wp.periods[0] / (3600*24*365))
		print("\n\x1b[31mMin and max orbital inclination (deg):\x1b[0m")
		print(acos(Hr/G_in) * rad_conv,acos(H_max/G_in) * rad_conv)
		print("\n\x1b[31mMin and max axial inclination (deg):\x1b[0m")
		print(acos((Hts_in - Hr)/Gt_in) * rad_conv,acos((Hts_in-H_max)/Gt_in)  * rad_conv)
		print("\n\x1b[31mNumber of lobes:\x1b[0m " + str(n_lobes))
		print("\n\x1b[31mLobe idx:\x1b[0m " + str(lobe_idx))
		# Report the known results for simplified system for comparison.
		H = H_in
		HHp,G,L,GG,eps,m2,Hts,Gt,J2 = [d_eval[s] for s in ['\\mathcal{H}^\\prime','G','L','\\mathcal{G}',\
			'\\varepsilon','m_2','\\tilde{H}_\\ast','\\tilde{G}','J_2']]
		print("\n\x1b[31mEinstein (g):\x1b[0m " + str((3 * eps * GG**4 * m2**4/(G**2*L**3))))
		print("\n\x1b[31mLense-Thirring (g):\x1b[0m " + str(((eps * ((-6*H*J2*GG**4*m2**3)/(G**4*L**3)+3*GG**4*m2**4/(G**2*L**3))))))
		print("\n\x1b[31mLense-Thirring (h):\x1b[0m " + str((2*eps*J2*GG**4*m2**3/(G**3*L**3))))
		# These are the Delta_ constants of quasi-periodicity.
		f_period = self.wp_period
		print("\n\x1b[31mDelta_g:\x1b[0m " + str(self.g_time(f_period) - self.g_time(0)))
		print("\n\x1b[31mg_rate:\x1b[0m " + str((self.g_time(f_period) - self.g_time(0))/f_period))
		Delta_hs = self.hs_time(f_period) - self.hs_time(0)
		print("\n\x1b[31mDelta_hs:\x1b[0m " + str(Delta_hs))
		print("\n\x1b[31mhs_rate:\x1b[0m " + str(Delta_hs/f_period))
		Delta_ht = self.ht_time(f_period) - self.ht_time(0)
		print("\n\x1b[31mDelta_ht:\x1b[0m " + str(Delta_ht))
		print("\n\x1b[31mht_rate:\x1b[0m " + str(Delta_ht/f_period))
		print("\n\x1b[31mDelta_h:\x1b[0m " + str(Delta_ht+Delta_hs))
		print("\n\x1b[31mh_rate:\x1b[0m " + str((Delta_ht+Delta_hs)/f_period))
		print("\n\n")
	@staticmethod
	def __to_sympy(s):
		from copy import deepcopy
		from sympy.parsing.sympy_parser import parse_expr
		from sympy import Symbol
		pyr_syms = ['G_{xy}','J_2','\\mathcal{G}','\\mathcal{I}_1','\\tilde{G}','\\tilde{H}_\\ast','\\varepsilon','m_2','h_\\ast','\\tilde{G}_{xy\\ast}']
		subs_dict = dict([(t[1],'s' + str(t[0])) for t in list(enumerate(pyr_syms))])
		retval = deepcopy(s)
		for k in subs_dict:
			retval = retval.subs(k,type(s)(subs_dict[k]))
		ret = parse_expr(repr(retval))
		for k in subs_dict:
			ret = ret.replace(Symbol(subs_dict[k]),Symbol(k))
		return ret
	@staticmethod
	def __split_derivative(chs,der):
		import sympy
		return [t.cancel().apart(sympy.Symbol('H')) for t in der.replace(sympy.cos(sympy.Symbol('h_\\ast')),chs).expand()\
			.replace(sympy.Symbol('G_{xy}'),sympy.sqrt(sympy.Symbol('G')**2-sympy.Symbol('H')**2))\
			.replace(sympy.Symbol('\\tilde{G}_{xy\\ast}'),sympy.sqrt(sympy.Symbol('\\tilde{G}')**2-(sympy.Symbol('\\tilde{H}_\\ast')-sympy.Symbol('H'))**2)).as_ordered_terms()]
	@staticmethod
	def __collect_forms(fH_l,der_list):
		import sympy
		retval = []
		tmp_expr = sum(der_list)
		for fc in fH_l:
			tmp_d = tmp_expr.collect(fc,evaluate=False,exact=True)
			assert(len(tmp_d)) == 2
			retval.append((fc,tmp_d[fc].ratsimp().expand()))
			tmp_expr = tmp_d[sympy.sympify(1)]
		retval.append((sympy.sympify(1),tmp_expr.ratsimp().expand()))
		assert(not sympy.Symbol('H') in tmp_expr.ratsimp().expand())
		assert((sum([t[0]*t[1] for t in retval]) - sum(der_list)).together().ratsimp() == 0)
		return retval
	@staticmethod
	def __solve_g(chs,dH_dG):
		import sympy
		# dg/dt.
		dg_dt_l = spin_gr_theory.__split_derivative(chs,dH_dG)
		fH_l = [sympy.Symbol('H'),(sympy.Symbol('G')-sympy.Symbol('H'))**-1,(sympy.Symbol('G')+sympy.Symbol('H'))**-1,\
			(sympy.Symbol('G')**2*sympy.Symbol('m_2')-sympy.Symbol('H')*sympy.Symbol('J_2'))**-1]
		retval = spin_gr_theory.__collect_forms(fH_l,dg_dt_l)
		# Now we will integrate the expression by substitution.
		result = []
		for element in list(retval[0:-1]):
			tmp1 = element[0].replace(sympy.Symbol('H'),sympy.Symbol('H_r')+sympy.Symbol('A')/(sympy.Symbol('\\wp')-sympy.Symbol('B')))\
				.apart(sympy.Symbol('\\wp'))
			a,b,c,d,e,f = [sympy.Wild(_,exclude = [sympy.Symbol('\\wp')]) for _ in 'abcdef']
			# NOTE: here the pattern has been crafted to make it work, but it seems a bit brittle (e.g., too many
			# wildcards and the needed minus sign on top).
			match = tmp1.match((-a*b/(c*(d*sympy.Symbol('\\wp')+e))) + f)
			assert(not match is None)
			t,g2,g3 = sympy.symbols('t g_2 g_3')
			tmp2 = integral_func(t,g2,g3,-match[a]*match[b]/(match[c]*match[d]),-match[e]/match[d],match[f])
			result.append((-match[a]*match[b]/(match[c]*(match[d]*sympy.Symbol('\\wp')+match[e]))+match[f],tmp2,element[1]))
		# Last element is the constant one.
		result.append((None,sympy.symbols('t'),retval[-1][1]))
		return (retval,result)
	@property
	def g_sol(self):
		from copy import deepcopy
		return (deepcopy(self.__g_sol[0]),deepcopy(self.__g_sol[1]))
	@staticmethod
	def __solve_hs(chs,dH_dH):
		import sympy
		# dhs/dt.
		dhs_dt_l = spin_gr_theory.__split_derivative(chs,dH_dH)
		fH_l = [(sympy.Symbol('G')-sympy.Symbol('H'))**-1,(sympy.Symbol('G')+sympy.Symbol('H'))**-1,\
			(sympy.Symbol('G')**2*sympy.Symbol('m_2')-sympy.Symbol('H')*sympy.Symbol('J_2'))**-1,\
			(sympy.Symbol('H')+sympy.Symbol('\\tilde{G}')-sympy.Symbol('\\tilde{H}_\\ast'))**-1,(sympy.Symbol('H')-sympy.Symbol('\\tilde{G}')-sympy.Symbol('\\tilde{H}_\\ast'))**-1]
		retval = spin_gr_theory.__collect_forms(fH_l,dhs_dt_l)
		result = []
		for element in list(retval[0:-1]):
			tmp1 = element[0].replace(sympy.Symbol('H'),sympy.Symbol('H_r')+sympy.Symbol('A')/(sympy.Symbol('\\wp')-sympy.Symbol('B'))).apart(sympy.Symbol('\\wp'))
			a,b,c,d,e,f,g = [sympy.Wild(_,exclude = [sympy.Symbol('\\wp')]) for _ in 'abcdefg']
			match = tmp1.match((a*b/(c*(d*sympy.Symbol('\\wp')+g*sympy.Symbol('\\wp')+e))) + f)
			assert(not match is None)
			t,g2,g3 = sympy.symbols('t g_2 g_3')
			tmp2 = integral_func(t,g2,g3,match[a]*match[b]/(match[c]*(match[d]+match[g])),-match[e]/(match[d]+match[g]),match[f])
			result.append((match[a]*match[b]/(match[c]*(match[d]*sympy.Symbol('\\wp')+match[g]*sympy.Symbol('\\wp')+match[e]))+match[f],tmp2,element[1]))
		# Last element is the constant one.
		result.append((None,sympy.symbols('t'),retval[-1][1]))
		return (retval,result)
	@property
	def hs_sol(self):
		from copy import deepcopy
		return (deepcopy(self.__hs_sol[0]),deepcopy(self.__hs_sol[1]))
	@staticmethod
	def __solve_ht(chs,dH_dHts):
		import sympy
		# dht/dt.
		dht_dt_l = spin_gr_theory.__split_derivative(chs,dH_dHts)
		fH_l = [(sympy.Symbol('H')+sympy.Symbol('\\tilde{G}')-sympy.Symbol('\\tilde{H}_\\ast'))**-1,\
			(sympy.Symbol('H')-sympy.Symbol('\\tilde{G}')-sympy.Symbol('\\tilde{H}_\\ast'))**-1]
		retval = spin_gr_theory.__collect_forms(fH_l,dht_dt_l)
		result = []
		for element in list(retval[0:-1]):
			tmp1 = element[0].replace(sympy.Symbol('H'),sympy.Symbol('H_r')+sympy.Symbol('A')/(sympy.Symbol('\\wp')-sympy.Symbol('B'))).apart(sympy.Symbol('\\wp'))
			a,b,c,d,e,f,g = [sympy.Wild(_,exclude = [sympy.Symbol('\\wp')]) for _ in 'abcdefg']
			match = tmp1.match((a*b/(c*(d*sympy.Symbol('\\wp')+g*sympy.Symbol('\\wp')+e))) + f)
			assert(not match is None)
			t,g2,g3 = sympy.symbols('t g_2 g_3')
			tmp2 = integral_func(t,g2,g3,match[a]*match[b]/(match[c]*(match[d]+match[g])),-match[e]/(match[d]+match[g]),match[f])
			result.append((match[a]*match[b]/(match[c]*(match[d]*sympy.Symbol('\\wp')+match[g]*sympy.Symbol('\\wp')+match[e]))+match[f],tmp2,element[1]))
		# Last element is the constant one.
		result.append((None,sympy.symbols('t'),retval[-1][1]))
		return (retval,result)
	@property
	def ht_sol(self):
		from copy import deepcopy
		return (deepcopy(self.__ht_sol[0]),deepcopy(self.__ht_sol[1]))
	def __get_params(self):
		from copy import deepcopy
		return deepcopy(self.__params)
	parameters = property(__get_params,__set_params)
	@property
	def H_time(self):
		from copy import deepcopy
		return deepcopy(self.__H_time)
	@property
	def g_time(self):
		from copy import deepcopy
		return deepcopy(self.__g_time)
	@property
	def hs_time(self):
		from copy import deepcopy
		return deepcopy(self.__hs_time)
	@property
	def ht_time(self):
		from copy import deepcopy
		return deepcopy(self.__ht_time)
	@property
	def obliquity_time(self):
		from copy import deepcopy
		return deepcopy(self.__obliquity_time)
	@property
	def spin_vector_time(self):
		from copy import deepcopy
		return deepcopy(self.__spin_vector_time)
	@property
	def orbit_vector_time(self):
		from copy import deepcopy
		return deepcopy(self.__orbit_vector_time)
	@property
	def HHp(self):
		from copy import deepcopy
		return deepcopy(self.__HHp)
	@property
	def f4(self):
		from copy import deepcopy
		return deepcopy(self.__f4)
	@property
	def f4_coeffs(self):
		from copy import deepcopy
		return deepcopy(self.__f4_coeffs)
	@property
	def invariants(self):
		from copy import deepcopy
		return deepcopy((self.__g2,self.__g3))
	@property
	def wp_period(self):
		from copy import deepcopy
		return deepcopy(self.__wp_period)
	def __init__(self,params = __default_params):
		from numpy import dot
		from mpmath import mpf
		from fractions import Fraction as Frac
		import sympy
		from IPython.parallel import Client
		# Set up constants.
		self.__eps_val = (1./mpf(299792458.))**2
		self.__GG_val = mpf(6.673E-11)
		# Various variables.
		Gt,I1,GG,m2,L,r,a,v2,Gtxy,ht,Ht,Gxy,h,H,J2,g,G,f,e,E,eps,hs,Hts,Gtxys = [pt(name) for name in ['\\tilde{G}','\\mathcal{I}_1',\
			'\\mathcal{G}','m_2','L','r','a','v2','\\tilde{G}_{xy}','\\tilde{h}','\\tilde{H}','G_{xy}','h','H','J_2','g','G','f','e','E',\
			'\\varepsilon','h_\\ast','\\tilde{H}_\\ast','\\tilde{G}_{xy\\ast}']]
		# The unperturbed Hamiltonian.
		H0 = Gt**2 * I1 / 2 - GG**2 * m2**2 * L**-2 / 2
		self.__HH0 = H0
		# Pieces of the perturbed Hamiltonian.
		Gt_vec = [Gtxy * math.sin(ht),-Gtxy * math.cos(ht),Ht]
		J2_vec = [0,0,J2]
		r_vec = dot(celmec.orbitalR([math.cos(g),math.sin(g),H*G**-1,Gxy*G**-1,math.cos(h),math.sin(h)]),[r*math.cos(f),r*math.sin(f),0])
		r_cross_v = [Gxy * math.sin(h),-Gxy * math.cos(h),H]
		H1 = -Frac(1,8) * v2 ** 2 - Frac(3,2) * v2 * GG * m2 * r ** -1 + Frac(1,2) * GG**2 * m2**2 * r**-2 +\
			Frac(3,2) * GG * m2 * r**-3 * dot(Gt_vec,r_cross_v) + 2 * GG * r**-3 * dot(J2_vec,r_cross_v) +\
			GG * r**-3 * (3 * dot(Gt_vec,r_vec) * dot(J2_vec,r_vec) * r**-2 - dot(Gt_vec,J2_vec))
		H1 = H1.subs('v2',GG * m2 * (2 * r**-1 - a ** -1)).subs('a',L ** 2 * (GG * m2)**-1)
		# Verify formula in the paper.
		assert(-Frac(1,8)*GG**4*m2**4*L**-4+r**-1*2*GG**3*m2**3*L**-2-r**-2*3*GG**2*m2**2+GG*r**-3*(\
			2*J2*H+3*J2*Gxy**2*Ht*(G**-2)/2+3*m2*Ht*H/2-J2*Ht+(3*m2/2*Gtxy*Gxy-3*J2/2*H*Gxy*Gtxy*G**-2)*math.cos(ht-h)+\
			3*J2*(-Frac(1,2)*Gxy**2*Ht*G**-2*math.cos(2*f+2*g)-Frac(1,4)*Gxy*Gtxy*G**-1*(1-H*G**-1)*math.cos(2*f+2*g+ht-h)+\
			Frac(1,4)*Gxy*Gtxy*G**-1*(1+H*G**-1)*math.cos(2*f+2*g-ht+h))) == H1)
		# Split the Hamiltonian in parts.
		A0 = H1.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['r']) == 0),t[1])).filter(lambda t: t[1] == pt(1)).subs('r',pt(1))
		A1 = H1.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['r']) == -1),t[1])).filter(lambda t: t[1] == pt(1)).subs('r',pt(1))
		A2 = H1.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['r']) == -2),t[1])).filter(lambda t: t[1] == pt(1)).subs('r',pt(1))
		A3a = H1.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['r']) == -3),t[1])).filter(lambda t: t[1] == pt(1)).subs('r',pt(1))
		A3b = H1.filter(lambda t: t[1] == math.cos(ht - h)).transform(lambda t: (t[0],pt(1))).subs('r',pt(1))
		B0 = H1.filter(lambda t: t[1] == math.cos(2*f + 2*g)).transform(lambda t: (t[0],pt(1))).subs('r',pt(1))
		B1 = H1.filter(lambda t: t[1] == math.cos(2*f + 2*g + ht - h)).transform(lambda t: (t[0],pt(1))).subs('r',pt(1))
		B2 = H1.filter(lambda t: t[1] == math.cos(2*f + 2*g - ht + h)).transform(lambda t: (t[0],pt(1))).subs('r',pt(1))
		# Make sure we got them right.
		assert(A0 + A1 * r**-1 + A2 * r**-2 + r**-3 * (A3a + A3b * math.cos(ht - h)) + r**-3 * (B0 * math.cos(2*f + 2*g) +\
			B1 * math.cos(2*f + 2*g + ht - h) + B2 * math.cos(2*f + 2*g - ht + h)) == H1)
		# This is the integrand in f (without the part that is integrated in E).
		f_int = A2 * r**-2 + r**-3 * (A3a + A3b * math.cos(ht - h)) + r**-3 * (B0 * math.cos(2*f + 2*g) + B1 * math.cos(2*f + 2*g + ht - h) + B2 * math.cos(2*f + 2*g - ht + h))
		# Change the integration variable to f (with the constant parts already taken out of the integral).
		f_int *= r**2
		# Substitute the definition of 1/r in terms of f.
		f_int = f_int.subs('r',pt('rm1')**-1).subs('rm1',GG*m2*G**-2*(1+e*math.cos(f)))
		# This is the integrand in E.
		E_int = A1 * r**-1
		# Change the integration variable to f.
		E_int *= r**2
		# Change the integration variable to E.
		E_int *= L*G*r**-1*GG**-1*m2**-1
		assert(E_int == A1*G*L*GG**-1*m2**-1)
		# K1.
		K1 = GG**2 * m2**2 * G**-1 * L**-3 * (f_int + E_int).filter(lambda t: t[1].t_degree(['f']) == 0)
		# K.
		K = K1 + A0
		# The generator.
		chi = G**-1 * (E_int.integrate('E') + f_int.integrate('f')) - L**3 * GG**-2 * m2**-2 * K1.integrate('l')
		# Verifiy that chi satisfies the homological equation, yielding K.
		assert((math.pbracket(H0,chi,['L','G','H','\\tilde{G}','\\tilde{H}'],['l','g','h','\\tilde{g}','\\tilde{h}']) + H1)\
			.subs('r',pt('rm1')**-1).subs('rm1',GG*m2*G**-2*(1+e*math.cos(f))) == K)
		# This is the complete Hamiltonian, with the two coordinates compressed into a single one.
		HHp = H0 + eps * K.subs('h',ht+hs).subs('\\tilde{H}',Hts-H).subs('\\tilde{G}_{xy}',Gtxys)
		# Record it as a member.
		self.__HHp = HHp
		# F0 and F1.
		F0 = HHp.filter(lambda t: t[1].t_degree(['h_\\ast']) == 0).transform(lambda t: (t[0].filter(lambda t: t[1].degree(['\\varepsilon']) == 1),t[1])).subs('\\varepsilon',pt(1))
		F1 = HHp.filter(lambda t: t[1].t_degree(['h_\\ast']) == 1).transform(lambda t: (t[0].filter(lambda t: t[1].degree(['\\varepsilon']) == 1),t[1])).subs('\\varepsilon',pt(1)).subs('h_\\ast',pt(0))
		assert(H0 + eps * F0 + eps * F1 * math.cos(pt('h_\\ast')) == H0 + eps * K.subs('h',ht+hs).subs('\\tilde{H}',Hts-H).subs('\\tilde{G}_{xy}',Gtxys))
		self.__F0 = F0
		self.__F1 = F1
		# Quartic polynomial.
		f4H = (eps**2 * F1 ** 2 - (pt('\\mathcal{H}^\\prime') - H0 - eps * F0)**2)\
			.ipow_subs('G_{xy}',2,G**2-H**2)\
			.ipow_subs('\\tilde{G}_{xy\\ast}',2,Gt**2-(Hts-H)**2)
		a4 = f4H.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['H']) == 0),t[1]))
		a3 = f4H.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['H']) == 1),t[1])).subs('H',pt(1)) / 4
		a2 = f4H.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['H']) == 2),t[1])).subs('H',pt(1)) / 6
		a1 = f4H.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['H']) == 3),t[1])).subs('H',pt(1)) / 4
		a0 = f4H.transform(lambda t: (t[0].filter(lambda t: t[1].degree(['H']) == 4),t[1])).subs('H',pt(1))
		# NOTE: these are not the polynomial coefficient strictly speaking, they are normalised by 4, 6 and 4
		# as shown above and in the assert below.
		self.__f4_cf = (a0,a1,a2,a3,a4)
		# Check we got them right.
		assert(a4+4*a3*H+6*a2*H**2+4*a1*H**3+a0*H**4 == f4H)
		# Store the coefficients - from high degree to low.
		self.__f4_coeffs = [t[0] * t[1].trim() for t in zip([1,4,6,4,1],self.__f4_cf)]
		# Derivatives of the quartic poly.
		f4Hp = math.partial(f4H,'H')
		f4Hpp = math.partial(f4Hp,'H')
		f4Hppp = math.partial(f4Hpp,'H')
		f4Hpppp = math.partial(f4Hppp,'H')
		# Check the derivatives for consistency.
		assert(f4Hp == 4*a3 + 12*a2*H + 12*a1*H**2 + 4*a0*H**3)
		assert(f4Hpp == 12*a2 + 24*a1*H + 12*a0*H**2)
		assert(f4Hppp == 24*a1 + 24*a0*H)
		assert(f4Hpppp == 24*a0)
		self.__f4 = [f4H, f4Hp, f4Hpp, f4Hppp, f4Hpppp]
		# Invariants for wp.
		g2 = a0*a4 - 4*a1*a3 + 3*a2**2
		g3 = a0*a2*a4 + 2*a1*a2*a3 - (a2**3) - (a0*a3**2) - (a1**2*a4)
		self.__g2 = g2
		self.__g3 = g3
		# Solve the angles.
		# Extract cosine of h_s.
		#chs = sympy.solve((spin_gr_theory.__to_sympy(self.HHp.trim())-sympy.Symbol('\\mathcal{H}^\\prime')).replace(sympy.cos(sympy.Symbol('h_\\ast')),sympy.Symbol('chs')),sympy.Symbol('chs'))[0]
		#ipy_view = Client().load_balanced_view()
		#g_sol = ipy_view.apply_async(spin_gr_theory.__solve_g,chs,spin_gr_theory.__to_sympy(math.partial(self.HHp.trim(),'G')))
		#hs_sol = ipy_view.apply_async(spin_gr_theory.__solve_hs,chs,spin_gr_theory.__to_sympy(math.partial(self.HHp.trim(),'H')))
		#ht_sol = ipy_view.apply_async(spin_gr_theory.__solve_ht,chs,spin_gr_theory.__to_sympy(math.partial(self.HHp.trim(),'\\tilde{H}_\\ast')))
		#self.__g_sol = g_sol.get()
		#self.__hs_sol = hs_sol.get()
		#self.__ht_sol = ht_sol.get()
		import pickle
		self.__g_sol = pickle.load(open('g_sol.pickle','rb'))
		self.__hs_sol = pickle.load(open('hs_sol.pickle','rb'))
		self.__ht_sol = pickle.load(open('ht_sol.pickle','rb'))
		# Set the parameters of the theory.
		self.__set_params(params)
		# Some sanity checks.
		HHp,G,L,H,GG,eps,m2,Hts,Gt,J2,hs,Gxy,Gtxys = [sympy.Symbol(s) for s in ['\\mathcal{H}^\\prime','G','L','H','\\mathcal{G}',\
			'\\varepsilon','m_2','\\tilde{H}_\\ast','\\tilde{G}','J_2','h_\\ast','G_{xy}','\\tilde{G}_{xy\\ast}']]
		# Phi_g^4 going to zero in the equilibrium point.
		assert(self.g_sol[0][3][1].subs(HHp,spin_gr_theory.__to_sympy(self.HHp.ipow_subs('G_{xy}',2,pt('G')**2\
			-pt('H')**2).subs('H',pt('G')**2*pt('m_2')*pt('J_2')**-1))).ratsimp() == 0)
		# Reduction to Einstein precession.
		simpl_HHp_ein = spin_gr_theory.__to_sympy(self.HHp.subs('J_2',0).subs('\\tilde{G}_{xy\\ast}',0).subs('\\tilde{H}_\\ast',pt('H'))\
			.subs('\\tilde{G}',0))
		ein_prec = sum([t[0]*t[1] for t in self.g_sol[0]]).subs('J_2',0).subs(Hts,H).subs(Gt,0)\
			.subs(HHp,simpl_HHp_ein).ratsimp()
		assert(ein_prec == 3 * eps * GG**4 * m2**4/(G**2*L**3))
		# Lense-Thirring precession for g.
		simpl_HHp_lt = spin_gr_theory.__to_sympy(self.HHp.subs('\\tilde{G}_{xy\\ast}',0).subs('\\tilde{H}_\\ast',pt('H'))\
			.subs('\\tilde{G}',0))
		lt_prec = sum([t[0]*t[1] for t in self.g_sol[0]]).subs(Hts,H).subs(Gt,0).subs(HHp,simpl_HHp_lt).ratsimp()
		assert(lt_prec == (eps * ((-6*H*J2*GG**4*m2**3)/(G**4*L**3)+3*GG**4*m2**4/(G**2*L**3))).ratsimp())
		# Geodetic effect on g.
		simpl_HHp_ge = spin_gr_theory.__to_sympy(self.HHp.subs('J_2',0)).subs(sympy.cos(hs),-1).subs(Gxy,sympy.sqrt(G**2-H**2))\
			.subs(Gtxys,sympy.sqrt(Gt**2-(Hts-H)**2)).subs(H,(Hts**2+G**2-Gt**2)/(2*Hts))
		ge_g = sum([t[0]*t[1] for t in self.g_sol[0]]).subs(J2,0).subs(Gxy,sympy.sqrt(G**2-H**2))\
			.subs(Gtxys,sympy.sqrt(Gt**2-(Hts-H)**2)).subs(H,(Hts**2+G**2-Gt**2)/(2*Hts)).subs(HHp,simpl_HHp_ge).ratsimp()
		assert(ge_g == (1 / (4*G**4*L**3) * (15*G**2*GG**4*eps*m2**4+9*GG**4*Gt**2*eps*m2**4-9*GG**4*Hts**2*eps*m2**4)).ratsimp())
		# No precession in h for Einstein case.
		assert((sum([t[0]*t[1] for t in self.hs_sol[0]]) + sum([t[0]*t[1] for t in self.ht_sol[0]])).subs(HHp,simpl_HHp_ein)\
			.subs('J_2',0).subs(Hts,H).ratsimp().subs(Gt,0) == 0)
		# h precession in LT.
		assert((sum([t[0]*t[1] for t in self.hs_sol[0]]) + sum([t[0]*t[1] for t in self.ht_sol[0]])).subs(HHp,simpl_HHp_lt)\
			.subs(Hts,H).ratsimp().subs(Gt,0).ratsimp() == 2*eps*J2*GG**4*m2**3/(G**3*L**3))
		# Geodetic effect for h and ht.
		assert((sum([t[0]*t[1] for t in self.hs_sol[0]]) + sum([t[0]*t[1] for t in self.ht_sol[0]])).subs(HHp,simpl_HHp_ge)\
			.subs(J2,0).subs(H,(Hts**2+G**2-Gt**2)/(2*Hts)).ratsimp() == 3*GG**4*Hts*eps*m2**4/(2*G**3*L**3))
		assert((sum([t[0]*t[1] for t in self.ht_sol[0]])).subs(HHp,simpl_HHp_ge)\
			.subs(J2,0).subs(H,(Hts**2+G**2-Gt**2)/(2*Hts)).ratsimp() == 3*GG**4*Hts*eps*m2**4/(2*G**3*L**3))

class spin_gr_theory_num(object):
	def __set_original_H(self):
		import sympy
		HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,H,hs = [sympy.Symbol(k) for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
			'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H','h_\\ast']]
		H_N = (I1*Gt**2)/2 - GG**2*m2**2/(2*L**2)
		F_0 = J2*GG**4*Hts*m2**3/(2*G**3*L**3)+3*H**3*J2*GG**4*m2**3/(2*G**5*L**3)+15*GG**4*m2**4/(8*L**4)\
			+3*H*J2*GG**4*m2**3/(2*G**3*L**3)-3*H**2*J2*GG**4*Hts*m2**3/(2*G**5*L**3)-3*H**2*GG**4*m2**4/(2*G**3*L**3)\
			+3*H*GG**4*Hts*m2**4/(2*G**3*L**3)-3*GG**4*m2**4/(G*L**3)
		Gxy = sympy.sqrt(G**2-H**2)
		Gtxys = sympy.sqrt(Gt**2-(Hts-H)**2)
		F_1 = -3*Gxy*H*J2*GG**4*Gtxys*m2**3/(2*G**5*L**3)+3*Gxy*GG**4*Gtxys*m2**4/(2*G**3*L**3)
		self.__orig_H = H_N+eps*(F_0+F_1*sympy.cos(hs))
		H_diff = [-self.__orig_H.diff(hs)]
		coord_diff = [self.__orig_H.diff(s) for s in [G,H,Hts]]
		self.__diffs = H_diff + coord_diff
	def __init__(self,params):
		import sympy
		from mpmath import mpf, sqrt, cos
		from copy import deepcopy
		self.__set_original_H()
		# Set up constants.
		self.__eps_val = (1./mpf(299792458.))**2
		self.__GG_val = mpf(6.673E-11)
		# Setup names.
		names = ['m2','r2','rot2','r1','rot1','i_a','ht','a','e','i','h']
		if not all([s in params for s in names]):
			raise ValueError('invalid set of parameters')
		# Convert all the values to mpf and introduce convenience shortcuts.
		m2, r2, rot2, r1, rot1, i_a, ht_in, a, e, i, h_in = [mpf(params[s]) for s in names]
		L_in = sqrt(self.__GG_val * m2 * a)
		G_in = L_in * sqrt(1. - e**2)
		H_in = G_in * cos(i)
		Gt_in = (2 * r1**2 * rot1) / 5
		Ht_in = Gt_in * cos(i_a)
		Hts_in = H_in + Ht_in
		hs_in = h_in - ht_in
		Gxy_in = sqrt(G_in**2 - H_in**2)
		Gtxys_in = sqrt(Gt_in**2 - Ht_in**2)
		J2 = (2 * m2 * r2**2 * rot2) / 5
		II_1 = mpf(5) / (2 * r1**2)
		# Evaluation dictionary.
		eval_names = [sympy.Symbol(s) for s in ['L','G','H','\\tilde{G}','\\tilde{H}_\\ast','h_\\ast','m_2','\\mathcal{G}','J_2',\
			'\\varepsilon','\\mathcal{I}_1','\\tilde{h}','g']]
		eval_values = [sympy.Float(x) for x in [L_in,G_in,H_in,Gt_in,Hts_in,hs_in,m2,self.__GG_val,J2,self.__eps_val,II_1,ht_in,0]]
		d_eval = dict(zip(eval_names,eval_values))
		self.__init_d_eval = deepcopy(d_eval)
		def diff_f(x,_):
			from sympy import Symbol
			H = x[0]
			hs = x[2]
			d_eval[Symbol('H')] = H
			d_eval[Symbol('h_\\ast')] = hs
			return [d.evalf(subs = d_eval) for d in self.__diffs]
		self.__diff_f = diff_f
	def integrate(self,t0,t1,n):
		import sympy
		from numpy import linspace
		from scipy.integrate import odeint
		t = linspace(t0,t1,n)
		init_conditions = [float(self.__init_d_eval[sympy.Symbol(k)]) for k in ['H','g','h_\\ast','\\tilde{h}']]
		return odeint(self.__diff_f,init_conditions,t)
