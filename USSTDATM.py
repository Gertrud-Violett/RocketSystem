#====================================================
#US Standard Atmosphere and Drag 2022 MIT License makkiblog.com

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#US Standard atm
def AmbPres(alt):
	g0=9.80665
	M=0.0289644
	Rstar=8.3144598
	mw = 0.028964

	if alt <= 11000:
		Pstd=101325
		Tstd=288.15
		Lb=-0.0065
		p=Pstd*((Tstd+Lb*(alt-0))/Tstd)**(-g0*M/Rstar/Lb)
		den=p/(Rstar/mw*Tstd)
		c=(1.4*p/den)**0.5
		return(p,den,c)

	elif alt <= 20000:
		Pstd=22632.1
		Tstd=216.65
		Lb=0.00001
		p=Pstd*((Tstd+Lb*(alt-11000))/Tstd)**(-g0*M/Rstar/Lb)
		den=p/(Rstar/mw*Tstd)
		c=(1.4*p/den)**0.5
		return(p,den,c)


	elif alt <= 32000:
		Pstd=5474.89
		Tstd=216.65
		Lb=0.001
		p=Pstd*((Tstd+Lb*(alt-20000))/Tstd)**(-g0*M/Rstar/Lb)
		den=p/(Rstar/mw*Tstd)
		c=(1.4*p/den)**0.5
		return(p,den,c)


	elif alt <= 47000:
		Pstd=868.019
		Tstd=228.65
		Lb=0.0028
		p=Pstd*((Tstd+Lb*(alt-32000))/Tstd)**(-g0*M/Rstar/Lb)
		den=p/(Rstar/mw*Tstd)
		c=(1.4*p/den)**0.5
		return(p,den,c)


	elif alt <= 51000:
		Pstd=110.906
		Tstd=270.65
		Lb=0.00001
		p=Pstd*((Tstd+Lb*(alt-47000))/Tstd)**(-g0*M/Rstar/Lb)
		den=p/(Rstar/mw*Tstd)
		c=(1.4*p/den)**0.5
		return(p,den,c)


	elif alt <= 71000:
		Pstd=66.9389
		Tstd=270.65
		Lb=-0.0028
		p=Pstd*((Tstd+Lb*(alt-51000))/Tstd)**(-g0*M/Rstar/Lb)
		den=p/(Rstar/mw*Tstd)
		c=(1.4*p/den)**0.5
		return(p,den,c)


	elif alt <= 90000:
		Pstd=3.95642
		Tstd=214.65
		Lb=-0.002
		p=Pstd*((Tstd+Lb*(alt-71000))/Tstd)**(-g0*M/Rstar/Lb)
		den=p/(Rstar/mw*Tstd)
		c=(1.4*p/den)**0.5
		return(p,den,c)


	elif alt <= 110000:  #Outer space above 100km
		p=2.81e-2 #13
		Tstd=184
		den=5.08e-7
		c=(1.4*p/den)**0.5
		return(p,den,c)

	elif alt <= 130000:  #Outer space above 110km
		p=2.17e-3 #13
		Tstd=374
		den=1.80e-8
		c=(1.4*p/den)**0.5
		return(p,den,c)

	elif alt <= 150000:  #Outer space above 130km
		p=7.03e-4 #13
		Tstd=635
		den=3.26e-09
		c=(1.4*p/den)**0.5
		return(p,den,c)

	elif alt <= 250000: 
		p=1e-5 #13
		Tstd=1000
		den=4e-11
		c=(1.4*p/den)**0.5
		return(p,den,c)

	else: 
		p=1e-13
		den=0
		c=1000
		return(p,den,c)



#drag coefficient
def MVspec():
	xarr = np.array([0,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.15,1.2,1.3,1.5,1.6,1.7,1.8,2.0,2.5,3.0,3.5,4.0,4.5])
	yarr = np.array([0.3,0.3,0.31,0.35,0.4,0.5,0.58,0.63,0.68,0.69,0.7,0.69,0.68,0.64,0.61,0.57,0.45,0.41,0.39,0.35,0.33])
	def func(x,a,b,c,d,e):
	 	Y = a*x**4+b*x**3+c*x**2+d*x+e
	 	return Y
	popt, pcov = curve_fit(func,xarr,yarr)
	return(popt[0],popt[1],popt[2],popt[3],popt[4])

#Aero drag==========================================
def CD(Mach):
	x=Mach
	(a,b,c,d,e) = MVspec()
	Y = a*x**4+b*x**3+c*x**2+d*x+e
	return Y

def aerodrag(c,vr_a,vr_th,vr_phi,rho,Area):
	v=np.sqrt(vr_a**2+vr_th**2+vr_phi**2)
	Mach=v/c
	drag=3*0.5*rho*CD(Mach)*v**2*Area  #add multiplier for AoA
	return(Mach,v,drag) #in N