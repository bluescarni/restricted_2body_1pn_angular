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

def Phig0(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return -5 * HHp / G + 5 * I1 * Gt**2 / (2 * G) - 5 * GG**2 * m2 ** 2 / (2 * G * L**2) + 75 * GG**4 * eps * m2**4 / (8 * G * L**4) \
		- 21 * GG**4 * eps * m2**4 / (2 * G**2 * L**3) - J2 * GG**4 * Hts * eps * m2**3 / (2 * G**4 * L**3)

def Phig1(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return 9 * J2 * GG**4 * eps * m2**3 / (2 * G**4 * L**3)

def Phig2(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return HHp / 2 - I1 * Gt**2 / 4 + GG**2 * m2**2 / (4 * L**2) - 15 * GG**4 * eps * m2**4 / (16 * L**4) + 9 * GG**4 * eps * m2**4 / (4 * G * L**3) \
		- 3 * J2 * GG**4 * eps * m2**3 / (2 * G**2  * L**3)  - 3 * GG**4 * Hts * eps * m2**4 / (4 * G**2 * L**3) + J2 * GG**4 * Hts * eps * m2**3 / (2 * G**3 * L**3)

def Phig3(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return HHp / 2 - I1 * Gt**2 / 4 + GG**2 * m2**2 / (4 * L**2) - 15 * GG**4 * eps * m2**4 / (16 * L**4) + 9 * GG**4 * eps * m2**4 / (4 * G * L**3) \
		+ 3 * J2 * GG**4 * eps * m2**3 / (2 * G**2  * L**3)  + 3 * GG**4 * Hts * eps * m2**4 / (4 * G**2 * L**3) + J2 * GG**4 * Hts * eps * m2**3 / (2 * G**3 * L**3)

def Phig4(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return 2*G*HHp*m2-G*I1*Gt**2*m2+G*GG**2*m2**3/L**2-15*G*GG**4*eps*m2**5/(4*L**4)+3*eps*GG**4*m2**5/L**3-J2*GG**4*Hts*eps*m2**4/(G**2*L**3)

def abgg1(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	alpha = A
	beta = B
	gamma = Hr
	return list((alpha,beta,gamma))

def abgg2(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	alpha = A / (G-Hr)**2
	beta = (A + B*G - B*Hr) / (G - Hr)
	gamma = 1 / (G - Hr)
	return list((alpha,beta,gamma))

def abgg3(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	alpha = -A / (G+Hr)**2
	beta = (-A + B*G + B*Hr) / (G + Hr)
	gamma = 1 / (G + Hr)
	return list((alpha,beta,gamma))

def abgg4(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	alpha = A*J2/(G**2*m2-Hr*J2)**2
	beta = (A*J2+B*G**2*m2-B*Hr*J2)/(G**2*m2-Hr*J2)
	gamma = 1 / (G**2*m2-Hr*J2)
	return list((alpha,beta,gamma))

def Phihs0(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return -9*J2*GG**4*eps*m2**3/(2*G**3*L**3) - 3*J2*GG**4*Gt**2*eps*m2**3/(2*G**5*L**3)

def Phihs1(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return -HHp/2 + I1*Gt**2/4 - GG**2*m2**2/(4*L**2) + 15*GG**4*eps*m2**4/(16*L**4) - 9*GG**4*eps*m2**4/(4*G*L**3) + 3*J2*GG**4*eps*m2**3/(2*G**2*L**3) \
		+ 3*GG**4*Hts*eps*m2**4/(4*G**2*L**3) - J2*GG**4*Hts*eps*m2**3/(2*G**3*L**3)

def Phihs2(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return HHp/2 - I1*Gt**2/4 + GG**2*m2**2/(4*L**2) - 15*GG**4*eps*m2**4/(16*L**4) + 9*GG**4*eps*m2**4/(4*G*L**3) + 3*J2*GG**4*eps*m2**3/(2*G**2*L**3) \
		+ 3*GG**4*Hts*eps*m2**4/(4*G**2*L**3) + J2*GG**4*Hts*eps*m2**3/(2*G**3*L**3)

def Phihs3(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return -J2*HHp+J2*I1*Gt**2/2-J2*GG**2*m2**2/(2*L**2)+15*J2*GG**4*eps*m2**4/(8*L**4)-3*J2*GG**4*eps*m2**4/(2*G*L**3)+J2**2*GG**4*Hts*eps*m2**3/(2*G**3*L**3)

def Phihs4(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return HHp/2-I1*Gt**2/4+GG**2*m2**2/(4*L**2)-15*GG**4*eps*m2**4/(16*L**4)+3*GG**4*eps*m2**4/(2*G*L**3)+3*J2*GG**4*Gt*eps*m2**3/(4*G**3*L**3)-J2*GG**4*Hts*eps*m2**3/(G**3*L**3)+3*GG**4*Gt**2*eps*m2**4/(4*G**3*L**3) \
		-3*GG**4*Gt*Hts*eps*m2**4/(4*G**3*L**3)+3*J2*GG**4*Gt**3*eps*m2**3/(4*G**5*L**3)-3*J2*GG**4*Gt**2*Hts*eps*m2**3/(2*G**5*L**3)+3*J2*GG**4*Gt*Hts**2*eps*m2**3/(4*G**5*L**3)

def Phihs5(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return HHp/2-I1*Gt**2/4+GG**2*m2**2/(4*L**2)-15*GG**4*eps*m2**4/(16*L**4)+3*GG**4*eps*m2**4/(2*G*L**3)-3*J2*GG**4*Gt*eps*m2**3/(4*G**3*L**3)-J2*GG**4*Hts*eps*m2**3/(G**3*L**3)+3*GG**4*Gt**2*eps*m2**4/(4*G**3*L**3) \
		+3*GG**4*Gt*Hts*eps*m2**4/(4*G**3*L**3)-3*J2*GG**4*Gt**3*eps*m2**3/(4*G**5*L**3)-3*J2*GG**4*Gt**2*Hts*eps*m2**3/(2*G**5*L**3)-3*J2*GG**4*Gt*Hts**2*eps*m2**3/(4*G**5*L**3)

def abghs1(eval_dict):
	return abgg2(eval_dict)

def abghs2(eval_dict):
	return abgg3(eval_dict)

def abghs3(eval_dict):
	return abgg4(eval_dict)

def abghs4(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	alpha = -A/(Hr+Gt-Hts)**2
	beta = (-A+B*Hr+B*Gt-B*Hts)/(Hr+Gt-Hts)
	gamma = 1/(Hr+Gt-Hts)
	return list((alpha,beta,gamma))

def abghs5(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	alpha = -A/(Hr-Gt-Hts)**2
	beta = (-A+B*Hr-B*Gt-B*Hts)/(Hr-Gt-Hts)
	gamma = 1/(Hr-Gt-Hts)
	return list((alpha,beta,gamma))

def Phiht0(eval_dict):
	HHp,G,I1,Gt,GG,m2,L,eps,J2,Hts,A,B,Hr = [eval_dict[k] for k in ['\\mathcal{H}^\\prime','G','\\mathcal{I}_1', \
		'\\tilde{G}','\\mathcal{G}','m_2','L','\\varepsilon','J_2','\\tilde{H}_\\ast','A','B','H_r']]
	return 2*J2*GG**4*eps*m2**3/(G**3*L**3) + 3*J2*GG**4*Gt**2*eps*m2**3/(2*G**5*L**3)

def Phiht1(eval_dict):
	return -Phihs4(eval_dict)

def Phiht2(eval_dict):
	return -Phihs5(eval_dict)

def abght1(eval_dict):
	return abghs4(eval_dict)

def abght2(eval_dict):
	return abghs5(eval_dict)
