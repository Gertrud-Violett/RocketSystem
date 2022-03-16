#====================================================
#Flight Calculator Polar 2022 MIT License Trude

import numpy as np
import sympy as sy
from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt
from rocketthrustcalc import ceathrust
from USSTDATM import AmbPres
from USSTDATM import aerodrag

#avoid divide by zero for arctan
def divide(a, b):
    if b == 0:
    	if a== 0:
    		return 0
    	else:
        	return a * 1e13 / abs(a)
    else:
        return a / b

#define vehicle params, alt=altitude,v=airspeed,W=weight
def vehicle(StgNo,Area,alt,vr_a,vr_th,vr_phi,W,Pc,ER,MR,Pamb,Ispeff,Thr,At,EngineQty):
	Pamb,rho,c=AmbPres(alt)
	(Isp, Mdot, Thrust)=ceathrust(Pc, ER, MR, Pamb*1e-6, Ispeff, Thr, At)
	Thrusttot = Thrust*EngineQty
	Mach,v,drag=aerodrag(c,vr_a,vr_th,vr_phi,rho,Area)
	Mdottot = Mdot*EngineQty
	return(Mach,v,drag,Isp,Mdottot,Thrusttot)

#Define Spherical momentum at t=i
def accel_mat(fgrav,Wveh,beta,gamma,Thrusttot,v,drag,rad,theta,phi,vr,dthetat,omega):
	vtot=np.sqrt(vr**2+(rad*dthetat)**2+(rad*omega*np.sin(theta))**2)
	ratio_r=np.arcsin(vr/vtot)
	ratio_th=np.arcsin(rad*dthetat/vtot)
	ratio_ph=np.arcsin(rad*omega*np.sin(theta)/vtot)

	accr=((Thrusttot*np.cos(beta)*np.cos(gamma)-fgrav-ratio_r*drag))/Wveh #R acceleration vector fgrav
	acctheta=(Thrusttot*np.sin(gamma)*np.cos(np.arctan(divide(np.sin(beta),np.sin(gamma))))-ratio_th*drag)/Wveh/rad #theta acceleration vector, use gamma for thrust history Constant orbital plane = 0
	accphi=(Thrusttot*np.sin(beta)*np.sin(np.arctan(abs(divide(np.sin(beta),np.sin(gamma)))))-ratio_ph*drag)/Wveh/rad/np.sin(theta) #phi acceleration vector, use beta for thrust history

	ddrt=accr+rad*dthetat**2+rad*omega**2*np.sin(theta)**2
	ddthetat=(acctheta-2/rad*vr*dthetat+omega**2*np.sin(theta)*np.cos(theta))
	ddphit=(accphi-2/rad*vr*omega*np.sin(theta)+2*dthetat*omega*np.cos(theta)) #domegat
	return(ddrt,ddthetat,ddphit)

#Define polar coordinate momentum equation classic Runge-Kutta method
def RGK_Matrix(dt,Mdottot,Wveh,beta,gamma,Thrusttot,v,drag,rad,theta,phi,vr,dthetat,omega):
	fgrav=Wveh*6.67259e-11*6.02e24/rad**2
	#fgrav=0
	#Acceleration vector for r,theta,phi
	(ddrt1,ddthetat1,ddphit1)=accel_mat(fgrav,Wveh,beta,gamma,Thrusttot,v,drag,rad,theta,phi,vr,dthetat,omega)

   #r,dr,theta,dtheta,phi,dphi
	[kr1,kdr1,kt1,kdt1,kp1,kdp1]=[
	vr,
	ddrt1,
	dthetat,
	ddthetat1,
	omega,
	ddphit1
	]

	(ddrt2,ddthetat2,ddphit2)=accel_mat(fgrav,Wveh,beta,gamma,Thrusttot,v,drag,(rad+kr1/2),(theta+kt1/2),(phi+kp1/2),(vr+kdr1/2),(dthetat+kdt1/2),(omega+kdp1/2))
	[kr2,kdr2,kt2,kdt2,kp2,kdp2]=[
	(vr+kdr1/2),
	ddrt2,
	(dthetat+kdt1/2),
	ddthetat2,
	(omega+kdp1/2),
	ddphit2
	]

	(ddrt3,ddthetat3,ddphit3)=accel_mat(fgrav,Wveh,beta,gamma,Thrusttot,v,drag,(rad+kr2/2),(theta+kt2/2),(phi+kp2/2),(vr+kdr2/2),(dthetat+kdt2/2),(omega+kdp2/2))
	[kr3,kdr3,kt3,kdt3,kp3,kdp3]=[
	(vr+kdr2/2),
	ddrt3,
	(dthetat+kdt2/2),
	ddthetat3,
	(omega+kdp2/2),
	ddphit3
	]

	(ddrt4,ddthetat4,ddphit4)=accel_mat(fgrav,Wveh,beta,gamma,Thrusttot,v,drag,(rad+kr3),(theta+kt3),(phi+kp3),(vr+kdr3),(dthetat+kdt3),(omega+kdp3))
	[kr4,kdr4,kt4,kdt4,kp4,kdp4]=[
	(vr+kdr3),
	ddrt4,
	(dthetat+kdt3),
	ddthetat4,
	(omega+kdp3),
	ddphit4
	]

	[rad,vr,theta,dthetat,phi,omega]=[
	rad+(kr1+2*kr2+2*kr3+kr4)/6*dt,
	vr+(kdr1+2*kdr2+2*kdr3+kdr4)/6*dt,
	theta+(kt1+2*kt2+2*kt3+kt4)/6*dt,
	dthetat+(kdt1+2*kdt2+2*kdt3+kdt4)/6*dt,
	phi+(kp1+2*kp2+2*kp3+kp4)/6*dt,
	omega+(kdp1+2*kdp2+2*kdp3+kdp4)/6*dt
	]

	Wveh-=Mdottot*dt
	return(rad,theta,phi,vr,dthetat,omega,ddrt1,ddthetat1,ddphit1,Wveh)


