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

class weierstrass_elliptic(object):
	# TODO:
	# - for testing, use results from Wolfram alpha with all possible branches (Delta positive and negative, g3 positive and negative)
	#   and both real and complex arguments.
	def __cubic_roots(self,a,c,d):
		from mpmath import mpf, mpc, sqrt, cbrt
		assert(all([isinstance(_,mpf) for _ in [a,c,d]]))
		Delta = -4 * a * c*c*c - 27 * a*a * d*d
		self.__Delta = Delta
		# NOTE: this was the original function used for root finding.
		# proots, err = polyroots([a,0,c,d],error=True,maxsteps=5000000)
		# Computation of the cubic roots.
		u_list = [mpf(1),mpc(-1,sqrt(3))/2,mpc(-1,-sqrt(3))/2]
		Delta0 = -3 * a * c
		Delta1 = 27 * a * a * d
		# Here we have two choices for the value from sqrt, positive or negative. Since C
		# is used as a denominator below, we pick the choice with the greatest absolute value.
		# http://en.wikipedia.org/wiki/Cubic_function
		# This should handle the case g2 = 0 gracefully.
		C1 = cbrt((Delta1 + sqrt(Delta1 * Delta1 - 4 * Delta0 * Delta0 * Delta0)) / 2)
		C2 = cbrt((Delta1 - sqrt(Delta1 * Delta1 - 4 * Delta0 * Delta0 * Delta0)) / 2)
		if abs(C1) > abs(C2):
			C = C1
		else:
			C = C2
		proots = [(-1 / (3 * a)) * (u * C + Delta0 / (u * C)) for u in u_list]
		# NOTE: we ignore any residual imaginary part that we know must come from numerical artefacts.
		if Delta < 0:
			# Sort the roots following the convention: complex with negative imaginary, real, complex with positive imaginary.
			# Then assign the roots following the P convention (e2 real root, e1 complex with positive imaginary).
			e3,e2,e1 = sorted(proots,key = lambda c: c.imag)
		else:
			# The convention in this case is to sort in descending order.
			e1,e2,e3 = sorted([_.real for _ in proots],reverse = True)
		return e1,e2,e3
	def __compute_periods(self):
		# A+S 18.9.
		from mpmath import sqrt, ellipk, mpc, pi, mpf
		Delta = self.Delta
		e1, e2, e3 = self.__roots
		if Delta > 0:
			m = (e2 - e3) / (e1 - e3)
			Km = ellipk(m)
			Kpm = ellipk(1 - m)
			om = Km / sqrt(e1 - e3)
			omp = mpc(0,1) * om * Kpm / Km
		elif Delta < 0:
			# NOTE: the expression in the sqrt has to be real and positive, as e1 and e3 are
			# complex conjugate and e2 is real.
			H2 = (sqrt((e2 - e3) * (e2 - e1))).real
			assert(H2 > 0)
			m = mpf(1) / mpf(2) - 3 * e2 / (4 * H2)
			Km = ellipk(m)
			Kpm = ellipk(1 - m)
			om2 = Km / sqrt(H2)
			om2p = mpc(0,1) * Kpm * om2 / Km
			om = (om2 - om2p) / 2
			omp = (om2 + om2p) / 2
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				om = mpf('+inf')
				omp = mpc(0,'+inf')
			else:
				# NOTE: here there is no need for the dichotomy on the sign of g3 because
				# we are already working in a regime in which g3 >= 0 by definition.
				c = e1 / 2
				om = 1 / sqrt(12 * c) * pi()
				omp = mpc(0,'+inf')
		return 2 * om, 2 * omp
	def __compute_user_periods(self):
		from mpmath import mpc
		Delta = self.Delta
		# NOTE: here there is no need to handle Delta == 0 separately,
		# as it falls under the case of om purely real and omp purely imaginary.
		if Delta >= 0:
			T1, T3 = self.__periods
			if self.__ng3:
				return T3.imag, T1.real * mpc(0,1)
			else:
				return T1, T3
		else:
			T1, T3 = (sum(self.__periods)).real, self.__periods[1]
			if self.__ng3:
				return 2 * T3.imag, mpc(-T3.imag,T3.real)
			else:
				return T1, T3
	def __compute_user_invariants(self):
		if self.__ng3:
			return self.__invariants[0],-self.__invariants[1]
		else:
			return self.__invariants
	def __init__(self,g2,g3):
		from mpmath import mpf, mpc
		if isinstance(g2,mpc) and isinstance(g3,mpc):
			if g2.imag != 0 or g3.imag != 0:
				raise ValueError('invalid invariant')
			g2 = g2.real
			g3 = g3.real
		g2 = mpf(g2)
		# Handle negative g3.
		if g3 < 0:
			g3 = -mpf(g3)
			self.__ng3 = True
		else:
			g3 = mpf(g3)
			self.__ng3 = False
		assert(g3 >= 0)
		# Store for future use.
		self.__invariants = (g2,g3)
		self.__user_invariants = self.__compute_user_invariants()
		self.__roots = self.__cubic_roots(mpf(4),-g2,-g3)
		self.__periods = self.__compute_periods()
		self.__user_periods = self.__compute_user_periods()
	@property
	def invariants(self):
		from copy import deepcopy
		return deepcopy(self.__user_invariants)
	@property
	def Delta(self):
		from copy import deepcopy
		return deepcopy(self.__Delta)
	@property
	def periods(self):
		from copy import deepcopy
		return deepcopy(self.__user_periods)
	@property
	def roots(self):
		from copy import deepcopy
		if self.__ng3:
			from functools import cmp_to_key
			# g3 < 0 means that all roots change sign.
			retval = [-x for x in self.__roots]
			# Sort by imaginary part, then real.
			return sorted(retval,key = cmp_to_key(lambda z1, z2: z1.imag - z2.imag if z1.imag != z2.imag else z1.real - z2.real),reverse=True)
		else:
			return deepcopy(self.__roots)
	def __repr__(self):
		retval = 'Invariants:\t' + str(self.invariants) + '\n'
		# NOTE: the Delta does not change for differences in sign of g_3.
		retval += 'Delta:\t\t' + str(self.Delta) + '\n'
		retval += 'Periods:\t' + str(self.periods) + '\n'
		retval += 'Roots:\t\t' + str(self.roots)
		return retval
	def P(self,z):
		# A+S 18.9.
		from mpmath import sqrt, mpc, sin, ellipfun, mpf
		Delta = self.Delta
		e1, e2, e3 = self.__roots
		if self.__ng3:
			z = mpc(0,1) * z
		if Delta > 0:
			zs = sqrt(e1 - e3) * z
			m = (e2 - e3) / (e1 - e3)
			retval = e3 + (e1 - e3) / ellipfun('sn',u=zs,m=m)**2
		elif Delta < 0:
			H2 = (sqrt((e2 - e3) * (e2 - e1))).real
			assert(H2 > 0)
			m = mpf(1) / mpf(2) - 3 * e2 / (4 * H2)
			zp = 2 * z * sqrt(H2)
			retval = e2 + H2 * (1 + ellipfun('cn',u=zp,m=m)) / (1 - ellipfun('cn',u=zp,m=m))
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				retval = 1 / (z**2)
			else:
				c = e1 / 2
				retval = -c + 3 * c / (sin(sqrt(3 * c) * z))**2
		if self.__ng3:
			return -retval
		else:
			return retval
	def Pprime(self,z):
		# A+S 18.9.
		from mpmath import ellipfun, sqrt, cos, sin, mpc, mpf
		Delta = self.Delta
		e1, e2, e3 = self.__roots
		if self.__ng3:
			z = mpc(0,1) * z
		if Delta > 0:
			zs = sqrt(e1 - e3) * z
			m = (e2 - e3) / (e1 - e3)
			retval = -2 * sqrt((e1 - e3)**3) * ellipfun('cn',u=zs,m=m) * ellipfun('dn',u=zs,m=m) / (ellipfun('sn',u=zs,m=m)**3)
		elif Delta < 0:
			H2 = (sqrt((e2 - e3) * (e2 - e1))).real
			assert(H2 > 0)
			m = mpf(1) / mpf(2) - 3 * e2 / (4 * H2)
			zp = 2 * z * sqrt(H2)
			retval = -4 * sqrt(H2**3) * ellipfun('sn',u=zp,m=m) * ellipfun('dn',u=zp,m=m) / ((1 - ellipfun('cn',u=zp,m=m))**2)
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				retval = -2 / (z**3)
			else:
				c = e1 / 2
				A = sqrt(3 * c)
				retval = -6 * c * A * cos(A * z) / (sin(A * z))**3
		if self.__ng3:
			return mpc(0,-1) * retval
		else:
			return retval
	def zeta(self,z):
		# A+S 18.10.
		from mpmath import pi, jtheta, exp, mpc, cot, sqrt
		Delta = self.Delta
		e1, _, _ = self.__roots
		om = self.__periods[0] / 2
		omp = self.__periods[1] / 2
		if self.__ng3:
			z = mpc(0,1) * z
		if Delta > 0:
			tau = omp / om
			# NOTE: here q has to be real.
			q = (exp(mpc(0,1) * pi() * tau)).real
			eta = -(pi()**2 * jtheta(n=1,z=0,q=q,derivative=3)) / (12 * om * jtheta(n=1,z=0,q=q,derivative=1))
			v = (pi() * z) / (2 * om)
			retval = (eta * z) / om + (pi() * jtheta(n=1,z=v,q=q,derivative=1)) / (2 * om * jtheta(n=1,z=v,q=q))
		elif Delta < 0:
			om2 = om + omp
			om2p = omp - om
			tau2 = om2p / (2 * om2)
			# NOTE: here q will be pure imaginary.
			q = mpc(0,(mpc(0,1) * exp(mpc(0,1) * pi() * tau2)).imag)
			eta2 = -(pi()**2 * jtheta(n=1,z=0,q=q,derivative=3)) / (12 * om2 * jtheta(n=1,z=0,q=q,derivative=1))
			v = (pi() * z) / (2 * om2)
			retval = (eta2 * z) / om2 + (pi() * jtheta(n=1,z=v,q=q,derivative=1)) / (2 * om2 * jtheta(n=1,z=v,q=q))
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				retval = 1 / z
			else:
				c = e1 / 2
				A = sqrt(3 * c)
				retval = c*z + A * cot(A*z)
		if self.__ng3:
			return mpc(0,1) * retval
		else:
			return retval
	def sigma(self,z):
		# A+S 18.10.
		from mpmath import pi, jtheta, exp, mpc, sqrt, sin
		Delta = self.Delta
		e1, _, _ = self.__roots
		om = self.__periods[0] / 2
		omp = self.__periods[1] / 2
		if self.__ng3:
			z = mpc(0,1) * z
		if Delta > 0:
			tau = omp / om
			q = (exp(mpc(0,1) * pi() * tau)).real
			eta = -(pi()**2 * jtheta(n=1,z=0,q=q,derivative=3)) / (12 * om * jtheta(n=1,z=0,q=q,derivative=1))
			v = (pi() * z) / (2 * om)
			retval = (2 * om) / pi() * exp((eta * z**2)/(2 * om)) * jtheta(n=1,z=v,q=q)/jtheta(n=1,z=0,q=q,derivative=1)
		elif Delta < 0:
			om2 = om + omp
			om2p = omp - om
			tau2 = om2p / (2 * om2)
			q = mpc(0,(mpc(0,1) * exp(mpc(0,1) * pi() * tau2)).imag)
			eta2 = -(pi()**2 * jtheta(n=1,z=0,q=q,derivative=3)) / (12 * om2 * jtheta(n=1,z=0,q=q,derivative=1))
			v = (pi() * z) / (2 * om2)
			retval = (2 * om2) / pi() * exp((eta2 * z**2)/(2 * om2)) * jtheta(n=1,z=v,q=q)/jtheta(n=1,z=0,q=q,derivative=1)
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				retval = z
			else:
				c = e1 / 2
				A = sqrt(3 * c)
				retval = (1 / A) * sin(A*z) * exp((c*z**2) / 2)
		if self.__ng3:
			return mpc(0,-1) * retval
		else:
			return retval
	def Pinv(self,P):
		from mpmath import ellipf, sqrt, asin, acos, mpc, mpf
		Delta = self.Delta
		e1, e2, e3 = self.__roots
		if self.__ng3:
			P = -P
		if Delta > 0:
			m = (e2 - e3) / (e1 - e3)
			retval = (1 / sqrt(e1 - e3)) * ellipf(asin(sqrt((e1 - e3)/(P - e3))),m=m)
		elif Delta < 0:
			H2 = (sqrt((e2 - e3) * (e2 - e1))).real
			assert(H2 > 0)
			m = mpf(1) / mpf(2) - 3 * e2 / (4 * H2)
			retval = 1 / (2 * sqrt(H2)) * ellipf(acos((e2-P+H2)/(e2-P-H2)),m=m)
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				retval = 1 / sqrt(P)
			else:
				c = e1 / 2
				retval = (1 / sqrt(3 * c)) * asin(sqrt((3 * c)/(P + c)))
		if self.__ng3:
			retval /= mpc(0,1)
		alpha, beta, _, _ = self.reduce_to_fpp(retval)
		T1, T2 = self.periods
		return T1 * alpha + T2 * beta
	def reduce_to_fpp(self,z):
		# TODO: need to check what happens here when periods are infinite.
		from mpmath import floor, matrix, lu_solve
		T1, T2 = self.periods
		R1, R2 = T1.real, T2.real
		I1, I2 = T1.imag, T2.imag
		A = matrix([[R1,R2],[I1,I2]])
		b = matrix([z.real,z.imag])
		x = lu_solve(A,b)
		N = int(floor(x[0]))
		M = int(floor(x[1]))
		alpha = x[0] - N
		beta = x[1] - M
		assert(alpha >= 0 and beta >= 0)
		return alpha,beta,N,M
	@staticmethod
	def __reduce_to_frp(z,om1):
		from mpmath import floor, mpc
		N = int(floor(z.real/(2*om1)))
		xs = z.real - N*2*om1
		return mpc(xs,z.imag),N
	def __ln_expansion(self,u,n_iter = 50):
		from mpmath import ln, pi, sin, exp, mpc
		# TODO: test for pathological cases?
		om1, om3 = self.periods[0]/2, self.periods[1]/2
		assert(u.real < 2*om1 and u.real >= 0)
		assert(u.imag < 2*om3.imag and u.imag >= 0)
		tau = om3/om1
		q = exp(mpc(0,1)*pi()*tau)
		eta1 = self.zeta(om1)
		om1_2 = om1 * 2
		retval = ln(om1_2/pi()) + eta1*u**2/(om1_2) + ln(sin(pi()*u/(om1_2)))
		for r in range(1,n_iter + 1):
			q_2 = q**(2*r)
			retval += q_2/(r * (1 - q_2))*(2. * sin(r*pi()*u/(om1_2)))**2
		return retval
	def ln_sigma(self,u):
		from mpmath import mpc, pi
		# TODO: test for pathological cases?
		om1, om3 = self.periods[0]/2, self.periods[1]/2
		if not (u.imag >= 0 and u.imag < 2*om3.imag):
			raise ValueError('invalid input for ln_sigma')
		further_reduction = False
		# This is used to further reduce, via the properties
		# of homogeneity, the computation to a value farther
		# from the limit where the expansion starts diverging.
		if u.imag / (2 * om3.imag) > .5:
			further_reduction = True
			u = -u + 2*om1 + 2*om3
		u_F,N = weierstrass_elliptic.__reduce_to_frp(u,om1)
		eta1 = self.zeta(om1)
		retval = self.__ln_expansion(u_F) + 2*N*eta1*(u_F + N*om1) - mpc(0,1)*N*pi()
		if further_reduction:
			retval -= (u - om1 - om3) * (2*eta1 + 2*self.zeta(om3))
		return retval
