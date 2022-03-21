#====================================================
#Flight Calculator Polar 2022 MIT License makkiblog.com

import numpy as np
import sympy as sy
import sys
from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from rocketthrustcalc import ceathrust
from vehicle import vehicle
from vehicle import RGK_Matrix
import pandas as pd
import toml
#mpl.rcParams['agg.path.chunksize'] = 100000

#Setup (Visual)
import seaborn as sns
sns.set_context("notebook", font_scale=1.0, rc={"lines.linewidth": 1.2})
sns.set_style('whitegrid')
plt.rcParams["figure.figsize"] = (16,9)
plt.rcParams['figure.facecolor'] = 'white'

#Read csv and TOML setting files====================
#1st Isp, 2nd Isp, 1st str, 2nd str, payload, throat area ratio
def runcalc(settingfile,casein,firstIspr,secIspr,firstr,secstr,pay,frem1,frem2,thar1,thar2):
	#settingfile = sys.argv[1] #specify setting file: settings.toml in TOML format

	with open(settingfile) as inputf:
		RSM = toml.load(inputf)

	betafile = RSM['CALC_SETTING']['betafile']
	gammafile = RSM['CALC_SETTING']['gammafile']

	anglecsv_b = np.loadtxt(betafile, delimiter=',', skiprows=1) #phi thrust angle file
	anglecsv_g = np.loadtxt(gammafile, delimiter=',', skiprows=1) # theta thrust angle file

	#=========================================================================================
	#SETTINGS read from setting file=========================================
	#Environment settings
	g0=RSM['ENVIRONMENT']['gravity'] #Gravity constant
	Pambient=RSM['ENVIRONMENT']['barometric_pressure'] #MPa Sea Level

	#VEHICLE Settings
	etaS = [RSM['VEHICLE']['Stage1']['Structural_Efficiency']*firstr,RSM['VEHICLE']['Stage2']['Structural_Efficiency']*secstr] #Structural Efficiency 1st, 2nd stage
	W0 = [RSM['VEHICLE']['Stage1']['Initial_Weight'],RSM['VEHICLE']['Stage2']['Initial_Weight']] #Full weight 1st, 2nd stage without payload
	WPayload = RSM['VEHICLE']['Payload']['Weight']*pay #Payload weight
	Wfuelrem = [RSM['VEHICLE']['Stage1']['Fuel_Remaining']*frem1,RSM['VEHICLE']['Stage2']['Fuel_Remaining']*frem2] #Fuel remaining for 1st, 2nd stage landing 
	Area = [RSM['VEHICLE']['Stage1']['Fuselage_Dia']**2/4*np.pi, RSM['VEHICLE']['Stage2']['Fuselage_Dia']**2/4*np.pi,RSM['VEHICLE']['Payload']['Fuselage_Dia']] #1st ,2nd stage Diameter

	#ENGINE Settings
	At=[RSM['ENGINE']['Stage1']['Throat_Dia']**2/4*np.pi*thar1, RSM['ENGINE']['Stage2']['Throat_Dia']**2/4*np.pi*thar2,0] #Throat dia in m2  [1st,2nd]
	Ispeff=[RSM['ENGINE']['Stage1']['Isp_efficiency']*firstIspr, RSM['ENGINE']['Stage2']['Isp_efficiency']*secIspr,1.0] #Isp efficiency
	Pc=[RSM['ENGINE']['Stage1']['Combustion_pressure'], RSM['ENGINE']['Stage2']['Combustion_pressure'],1.0] #Combustion chamber total pressure in MPa
	ER=[RSM['ENGINE']['Stage1']['Nozzle_ExpansionRatio'], RSM['ENGINE']['Stage2']['Nozzle_ExpansionRatio'],1.0] #Expansion Ratio
	MR=[RSM['ENGINE']['Stage1']['OF_mixture'],RSM['ENGINE']['Stage2']['OF_mixture'],3.8] # Mixture Ratio
	Thrust=[RSM['ENGINE']['Stage1']['Target_thrust'],RSM['ENGINE']['Stage2']['Target_thrust'],1] #kNs
	EngineQty=[RSM['ENGINE']['Stage1']['Engine_qty'],RSM['ENGINE']['Stage2']['Engine_qty'],1] # No of engines

	#VECTOR COORDINATES settings
	Latitude_deg = 90 - RSM['VECTOR']['Latitude'] #Theta degress = Latitude: 90 deg at equator, 0 deg at pole 61.75
	Longitude_deg = RSM['VECTOR']['Longitude'] #Phi degrees in -180~180deg -80.5
	r0 = RSM['VECTOR']['radius'] #Initial radius = default: Radius for Earth
	alt = RSM['VECTOR']['altitude']

	#Initialize velocities
	vr=0 #Initial vertical velocity, =0
	#Latitude (=np.pi/2 for equator, 0 for pole)
	theta0=np.pi/2*Latitude_deg/90
	dthetat0=0
	#Longitude
	phi0=Longitude_deg*np.pi/180 
	vphi0=465*theta0/np.pi*2 #Rotation at equator m/s 465 default
	omega0=vphi0/r0 #Initialize radial velocity

	#struc calcs
	WInit = W0[0] + W0[1] + WPayload
	Wfuel = [W0[0]*etaS[0],W0[1]*etaS[1]] #fuel weight 1st, 2nd stage
	Wfinal1 = WInit - (1-Wfuelrem[0])*Wfuel[0] #Final vehicle weight 1st stage
	W1B=Wfuelrem[0]*Wfuel[0]+W0[0]*(1-etaS[0])
	W1L=W0[0]*(1-etaS[0])
	Wfinal2 = W0[1] + WPayload - (1-Wfuelrem[1])*Wfuel[1] #final vehicle weight before payload release


	#CALCULATION SETTINGS================================
	#Initialize 1st stage================================
	StgNo=0 # 0 for stage 1, 1 for stage 2
	Pamb=Pambient #Initial Pressure
	Wveh = WInit #Initial vehicle weight
	[vr_a,vr_th,vr_phi]=[0,0,0] #Initialize airspeed

	t=0 #start time
	dt=RSM['CALC_SETTING']['timestep'] #timestep
	timelimit = RSM['CALC_SETTING']['timelimit'] #max time to break calculation

	casename = casein
	#casename = RSM['CALC_SETTING']['Casename']
	Case=('Study N_'+str(Latitude_deg)+' E_'+str(Longitude_deg)+'_'+casename)  #Study name

	#SETTINGS END HERE===================================


	#Initialization======================================
	#Theroretical deltaV without loss output
	deltaVth1=342*g0*np.log(WInit/Wfinal1)
	deltaVth2=382*g0*np.log((W0[1]+WPayload)/Wfinal2)
	deltaV1stL=342*g0+np.log(W1B/W1L)

	#MAIN CALCULATION===========================================================================
	#Initialize arrays
	(Mach0,v,drag0,Isp0,Mdottot0,Thrusttot0)=vehicle(StgNo,Area[StgNo],alt,vr_a,vr_th,vr_phi,Wveh,Pc[StgNo], ER[StgNo], MR[StgNo], Pamb, Ispeff[StgNo], Thrust[StgNo],At[StgNo],EngineQty[StgNo])
	timearr= np.array([0])
	Rocketarr = np.array([[Isp0,Mdottot0,Thrusttot0,0]])  #Isp, Mdot, Thrust: Isp0,Mdottot0,Thrusttot0, drag
	attitudearr= np.array([[t,0,0]]) #Beta and gamma angle
	Wveharr= np.array([WInit]) #Vehicle weight Winit

	Rarr=np.array([[r0,theta0,phi0,0,dthetat0,omega0,0,0,0]])   #r, th, phi, and dr, dth, dphi and ddr, ddth, ddphipolar coordinate matrix r0,theta0,0
	velarr=np.array([[0,dthetat0*r0,omega0*r0*np.sin(theta0)]]) #vr, vtheta, vphi
	alt_arr=np.array([[0,vphi0]]) # altitude, ABS Velocity

	[rad,theta,phi,vr,dthetat,omega]=[Rarr[0,0],Rarr[0,1],Rarr[0,2],Rarr[0,3],Rarr[0,4],Rarr[0,5]]

	#Define vehicle calculation each stage
	def vehcalc(t,Wveh,Wfin,rad,theta,phi,vr,dthetat,omega,timearr,Rocketarr,attitudearr,Wveharr,Rarr,velarr,alt_arr):
		while Wveh >  Wfin:
			#Global cooridnate ABS speed
			vel_abs=np.sqrt(vr**2+(rad*dthetat)**2+(rad*omega*np.sin(theta))**2)
			#Air speed relative to ground
			[vr_a,vr_th,vr_phi]=[vr,rad*dthetat,(omega-omega0)*rad*np.sin(theta)]

			#Define beta angle in radians and altitude
			beta=np.radians(np.interp(t,anglecsv_b[:,0],anglecsv_b[:,1]))
			gamma=np.radians(np.interp(t,anglecsv_g[:,0],anglecsv_g[:,1]))
			alt=rad-r0

			#Calc R
			(Mach,v,drag,Isp,Mdottot,Thrusttot)=vehicle(StgNo,Area[StgNo],alt,vr_a,vr_th,vr_phi,Wveh,Pc[StgNo], ER[StgNo], MR[StgNo], Pamb, Ispeff[StgNo], Thrust[StgNo],At[StgNo],EngineQty[StgNo])
			(rad,theta,phi,vr,dthetat,omega,ddrt,ddthetat,ddphit,Wveh)=RGK_Matrix(dt,Mdottot,Wveh,beta,gamma,Thrusttot,v,drag,rad,theta,phi,vr,dthetat,omega)
			accel_abs=np.sqrt(ddrt**2+(rad*ddphit*np.sin(theta))**2+(rad*ddthetat)**2)
			#print(" Accel:",ddrt,ddthetat,ddphit)
			#print(" Stage:%1.0f,Time[s]:%1.1f,Weight[t]:%1.0f, Drag[kN]:%1.2f,Thrust[kN]:%1.0f,Alt[m]:%1.0f, Accel[m/s2]:%1.2f,Velocity[m/s]:%1.0f"%(StgNo+1,t,Wveh*1e-3,drag*1e-3,Thrusttot*1e-3,alt,accel_abs,vel_abs))
			timearr=np.append(timearr,t)
			Rocketarr = np.vstack([Rocketarr,[Isp,Mdottot,Thrusttot,drag]])  #Isp, Mdot, Thrust
			attitudearr=np.vstack([attitudearr,[t,beta,gamma]])
			Wveharr=np.append(Wveharr,Wveh)
			Rarr=np.vstack([Rarr,[rad,theta,phi,vr,dthetat,omega,ddrt,ddthetat,ddphit]])
			velarr=np.vstack([velarr,[vr,dthetat*rad,omega*rad*np.sin(theta)]])
			alt_arr=np.vstack([alt_arr,[alt,vel_abs]])

			t+=dt
			if t>timelimit:
				break
		
		result=(" Downrange[km]:%1.0f\n Weight[t]:%1.0f\n Drag[kN]:%1.2f\n Alt[km]:%1.2f\n Orbit Speed [m/s]%1.2f\n Vehicle ABS speed [m/s]:%1.0f\n Vehicle Airspeed [Mach]:%1.3f\n"%(rad*phi*1e-3,Wveh*1e-3,drag*1e-3,alt*1e-3,rad*omega*np.sin(theta),vel_abs,Mach))
		print(result)
		
		with open('./plot/result_'+Case+'_Stg'+str(StgNo+1)+'.txt',"w",newline='\n') as output:
			print(result, file=output)
		
		return(t,Wveh,rad,theta,phi,vr,dthetat,omega,timearr,Rocketarr,attitudearr,Wveharr,Rarr,velarr,alt_arr)


	#SETTINGS_2=============================================================================
	#RUN CALCULATION _ Set Stages=========================
	#RUN 1st Stage========================================
	(t,Wveh,rad,theta,phi,vr,dthetat,omega,timearr,Rocketarr,attitudearr,Wveharr,Rarr,velarr,alt_arr)=vehcalc(t,Wveh,Wfinal1,rad,theta,phi,vr,dthetat,omega,timearr,Rocketarr,attitudearr,Wveharr,Rarr,velarr,alt_arr)

	#run 2nd stage
	StgNo=1 # 0 for stage 1, 1 for stage 2
	Wveh = W0[1] + WPayload #Initial vehicle weight
	(t2,Wveh2,rad2,theta2,phi2,vr2,dthetat2,omega2,timearr2,Rocketarr2,attitudearr2,Wveharr2,Rarr2,velarr2,alt_arr2)=vehcalc(t,Wveh,Wfinal2,rad,theta,phi,vr,dthetat,omega,timearr,Rocketarr,attitudearr,Wveharr,Rarr,velarr,alt_arr)

	#Calculate Payload
	StgNo=2
	Wveh = WPayload #Payload weight
	(t3,Wveh3,rad3,theta3,phi3,vr3,dthetat3,omega3,timearr3,Rocketarr3,attitudearr3,Wveharr3,Rarr3,velarr3,alt_arr3)=(t2,Wveh,rad2,theta2,phi2,vr2,dthetat2,omega2,timearr2,Rocketarr2,attitudearr2,Wveharr2,Rarr2,velarr2,alt_arr2)

	#CALC SETTINGS END HERE==============================
	#====================================================
	#====================================================


	#PLOTTING==============================================================================
	#Plot results

	with open('./plot/result_'+Case+'_summary.txt',"w",newline='\n') as output:
		altres=(" 1st Time[s],alt[km]: %1.0f,%1.0f\n 2nd Time[s],alt[km]: %1.0f,%1.0f\n Payload Time[s],alt[km]: %1.0f,%1.0f\n"%(t,alt_arr[-1,0]*1e-3,t2,alt_arr2[-1,0]*1e-3,t3,alt_arr3[-1,0]*1e-3))
		deltaVth=("Theta[deg]:%1.1f\n DeltaV_th[m/s]:%1.2f\n DeltaV1st[m/s]:%1.2f\n  Landing DeltaV[m/s]:%1.2f\n" %(theta0/np.pi*2*90,(deltaVth1+deltaVth2),deltaVth1,deltaV1stL)) #for check
		print(altres,deltaVth,file=output)
	print(altres,deltaVth)


	fig1,ax1=plt.subplots(3,3)
	ax1[0,0].set_xlabel('Downrange R x phi [km]')
	ax1[0,0].set_ylabel('Altitude R [km]')
	ax1[0,0].plot(Rarr3[:,2]*Rarr3[:,0]/1e3,alt_arr3[:,0]/1e3)
	ax1[0,0].plot(Rarr2[:,2]*Rarr2[:,0]/1e3,alt_arr2[:,0]/1e3)
	ax1[0,0].plot(Rarr[:,2]*Rarr[:,0]/1e3,alt_arr[:,0]/1e3)

	ax1[0,1].set_xlabel('Downrange velocity VR x phi [m/s]')
	ax1[0,1].set_ylabel('Radial Velocity VR [m/s]')
	ax1[0,1].plot(Rarr3[:,5]*Rarr3[:,0],Rarr3[:,3])
	ax1[0,1].plot(Rarr2[:,5]*Rarr2[:,0],Rarr2[:,3])
	ax1[0,1].plot(Rarr[:,5]*Rarr[:,0],Rarr[:,3])

	ax1[0,2].set_xlabel('Times [s]')
	ax1[0,2].set_ylabel('total Velocity v [m/s]')
	ax1[0,2].plot(timearr3,alt_arr3[:,1])
	ax1[0,2].plot(timearr2,alt_arr2[:,1])
	ax1[0,2].plot(timearr,alt_arr[:,1])

	ax1[1,0].set_xlabel('Times [s]')
	ax1[1,0].set_ylabel('VR x Phi Orbital Velocity [m/s]')
	ax1[1,0].plot(timearr3,velarr3[:,2])
	ax1[1,0].plot(timearr2,velarr2[:,2])
	ax1[1,0].plot(timearr,velarr[:,2])

	ax1[1,1].set_xlabel('Times [s]')
	ax1[1,1].set_ylabel('VR x Theta Velocity [m/s]')
	ax1[1,1].plot(timearr3,velarr3[:,1])
	ax1[1,1].plot(timearr2,velarr2[:,1])
	ax1[1,1].plot(timearr,velarr[:,1])

	ax1[1,2].set_xlabel('Times [s]')
	ax1[1,2].set_ylabel('VR x Radial Velocity [m/s]')
	ax1[1,2].plot(timearr3,Rarr3[:,3])
	ax1[1,2].plot(timearr2,Rarr2[:,3])
	ax1[1,2].plot(timearr,Rarr[:,3])

	ax1[2,0].set_xlabel('Times [s]')
	ax1[2,0].set_ylabel('Phi [km]')
	ax1[2,0].plot(timearr3,Rarr3[:,2]*Rarr3[:,0]/1e3)
	ax1[2,0].plot(timearr2,Rarr2[:,2]*Rarr2[:,0]/1e3)
	ax1[2,0].plot(timearr,Rarr[:,2]*Rarr[:,0]/1e3)

	ax1[2,1].set_xlabel('Times [s]')
	ax1[2,1].set_ylabel('Theta [deg]')
	ax1[2,1].plot(timearr3,Rarr3[:,1]/np.pi*2*360)
	ax1[2,1].plot(timearr2,Rarr2[:,1]/np.pi*2*360)
	ax1[2,1].plot(timearr,Rarr[:,1]/np.pi*2*360)

	ax1[2,2].set_xlabel('Times [s]')
	ax1[2,2].set_ylabel('R [km]')
	ax1[2,2].plot(timearr3,Rarr3[:,0]/1e3)
	ax1[2,2].plot(timearr2,Rarr2[:,0]/1e3)
	ax1[2,2].plot(timearr,Rarr[:,0]/1e3)

	plt.suptitle(Case)
	fig1.tight_layout()
	fig1.savefig('./plot/'+Case+'_plot.jpg', dpi=300)


	plt.clf()
	ax = plt.subplot(111,polar=True)
	ax.set_rlim([0,rad3*1.1])
	ax.set_thetalim([-np.pi, np.pi])
	circle = plt.Circle((0, 0), 6400e3, transform=ax.transData._b,color='royalblue', alpha=0.5)
	ax.add_patch(circle)
	ax.plot(Rarr3[:,2],Rarr3[:,0])
	ax.plot(Rarr2[:,2],Rarr2[:,0])
	ax.plot(Rarr[:,2],Rarr[:,0])
	plt.savefig('./plot/'+Case+'_Trajectory2D.jpg', dpi=300)
	#plt.show()


	#plot 3D
	RRA=Rarr3[:,0]
	RTA=Rarr3[:,1]
	RPA=Rarr3[:,2]
	X=RRA*np.sin(RTA)*np.cos(RPA)*1e-3
	Y=RRA*np.sin(RTA)*np.sin(RPA)*1e-3
	Z=RRA*np.cos(RTA)*1e-3
	u = np.linspace(0, np.pi, 100)
	v = np.linspace(0, 2*np.pi, 100)
	xs = r0 * np.outer(np.cos(u), np.sin(v)) * 1e-3
	ys = r0 * np.outer(np.sin(u), np.sin(v)) * 1e-3
	zs = r0 * np.outer(np.ones(np.size(u)), np.cos(v)) * 1e-3


	#output csv file and plotplt 3D plot
	dfnparr=np.vstack((timearr3,X,Y,Z))
	dfnparr=np.vstack((dfnparr,Rarr3[:,0].transpose()*1e-3,Rarr3[:,1].transpose()*1e-3,Rarr3[:,2].transpose()*1e-3,Rocketarr3.transpose(),Wveharr3.transpose(),alt_arr3.transpose()))
	dfnparr=dfnparr.transpose()
	np.savetxt('./plot/out_rocketarr'+Case+'.csv',Rocketarr3,delimiter=',')
	df_data = pd.DataFrame(dfnparr,columns=['Time[s]','ECI-X[km]','ECI-Y[km]','ECI-Z[km]','ECI-R[m]','ECI-Theta[rad]','ECI-Phi[rad]','Isp[s]','Mdot[kg/s]','Thrust[N]','drag[N]','Weight[kg]','Altitude[m]','Velocity[m/s]'])
	print(df_data)
	df_data.to_csv('./plot/out_XYZ_'+Case+'.csv')


	return(timearr3,Rocketarr3,Rarr3,velarr3,alt_arr3)
