import toml
from FlightCalc_study import runcalc
import numpy as np
import matplotlib.pyplot as plt

#runcalc("baseline.toml","baseline",1.0,1.0,1.0,1.0,1.0,1.0,1.0)
#1st Isp, 2nd Isp, 1st str, 2nd str, payload, fuel remaining, throat area ratio 1st / 2nd
RatioArray=[0.8,0.85,0.9,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.0]
PayloadRatio=[0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
FuelremArray=[0.25,0.5,0.75,1.25,1.5]
ThrustRatio=[0.6,0.7,0.8,0.85,0.9,0.95,1.05,1.1,1.15,1.2]

#Isp
Isp1_arr= np.empty((0,18), float)
for Ratio in RatioArray:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","1st_Isp_"+str(Ratio)+"_",Ratio,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
	Rocketarr[:,0]=Rocketarr[:,0]*Ratio
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Isp1_arr=np.append(Isp1_arr,np.array([apparrL]), axis=0)
print("Isp1")

Isp2_arr= np.empty((0,18), float)
for Ratio in RatioArray:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","2nd_Isp_"+str(Ratio)+"_",1.0,Ratio,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
	Rocketarr[:,0]=Rocketarr[:,0]*Ratio
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Isp2_arr=np.append(Isp2_arr,np.array([apparrL]), axis=0)
print("Isp2")

Str1_arr= np.empty((0,18), float)
for Ratio in RatioArray:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","1st_Str_"+str(Ratio)+"_",1.0,1.0,Ratio,1.0,1.0,1.0,1.0,1.0,1.0)
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Str1_arr=np.append(Str1_arr,np.array([apparrL]), axis=0)
print("Str1")

Str2_arr= np.empty((0,18), float)
for Ratio in RatioArray:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","2nd_Str_"+str(Ratio)+"_",1.0,1.0,1.0,Ratio,1.0,1.0,1.0,1.0,1.0)
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Str2_arr=np.append(Str2_arr,np.array([apparrL]), axis=0)
print("Str2")

Pay_arr= np.empty((0,18), float)
for Ratio in PayloadRatio:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","Payload_"+str(Ratio)+"_",1.0,1.0,1.0,1.0,Ratio,1.0,1.0,1.0,1.0)
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Pay_arr=np.append(Pay_arr,np.array([apparrL]), axis=0)
print("Pay")

Fuelrem_arr1= np.empty((0,18), float)
for Ratio in FuelremArray:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","FuelRem_"+str(Ratio)+"_",1.0,1.0,1.0,1.0,1.0,Ratio,1.0,1.0,1.0)
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Fuelrem_arr1=np.append(Fuelrem_arr1,np.array([apparrL]), axis=0)
print("Fuelrem1")

Fuelrem_arr2= np.empty((0,18), float)
for Ratio in FuelremArray:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","FuelRem_"+str(Ratio)+"_",1.0,1.0,1.0,1.0,1.0,1.0,Ratio,1.0,1.0)
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Fuelrem_arr2=np.append(Fuelrem_arr2,np.array([apparrL]), axis=0)
print("Fuelrem2")

Thrust_arr1= np.empty((0,18), float)
for Ratio in ThrustRatio:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","Thrust_"+str(Ratio)+"_",1.0,1.0,1.0,1.0,1.0,1.0,1.0,Ratio,1.0)
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Thrust_arr1=np.append(Thrust_arr1,np.array([apparrL]), axis=0)
print("Thrust1")

Thrust_arr2= np.empty((0,18), float)
for Ratio in ThrustRatio:
	(timearr,Rocketarr,Rarr,velarr,alt_arr)=runcalc("baseline.toml","Thrust_"+str(Ratio)+"_",1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,Ratio)
	appendarr=np.concatenate([Rocketarr,Rarr,velarr,alt_arr],1)
	apparrL = appendarr[-1,:]
	Thrust_arr2=np.append(Thrust_arr2,np.array([apparrL]), axis=0)
print("Thrust2")


#Plot results
fig1,ax1=plt.subplots(3,3)
ax1[0,0].set_xlabel('Isp(1st) Ratio')
ax1[0,0].set_ylabel('DeltaV [m/s]')
x=RatioArray
y=Isp1_arr[:,17]
res1=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[0,0].plot(x,y, label="Sensitivity"+str(res1))
#ax1[0,0].legend()

ax1[0,1].set_xlabel('Isp(2nd) Ratio')
ax1[0,1].set_ylabel('DeltaV [m/s]')
x=RatioArray
y=Isp2_arr[:,17]
res2=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[0,1].plot(x,y, label="Sensitivity"+str(res2))
#ax1[0,1].legend()

ax1[0,2].set_xlabel('Structural Efficiency(1st)')
ax1[0,2].set_ylabel('DeltaV [m/s]')
x=RatioArray
y=Str1_arr[:,17]
res3=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[0,2].plot(x,y, label="Sensitivity"+str(res3))
#ax1[0,2].legend()

ax1[1,0].set_xlabel('Structural Efficiency(2nd)')
ax1[1,0].set_ylabel('DeltaV [m/s]')
x=RatioArray
y=Str2_arr[:,17]
res21=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[1,0].plot(x,y, label="Sensitivity"+str(res21))
#ax1[0,3].legend()

ax1[1,1].set_xlabel('Fuel Remaining (1st) ')
ax1[1,1].set_ylabel('DeltaV [m/s]')
x=FuelremArray
y=Fuelrem_arr1[:,17]
res22=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[1,1].plot(x,y, label="Sensitivity"+str(res22))
#ax1[1,1].legend()

ax1[1,2].set_xlabel('Fuel Remaining (2nd) ')
ax1[1,2].set_ylabel('DeltaV [m/s]')
x=FuelremArray
y=Fuelrem_arr2[:,17]
res23=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[1,2].plot(x,y, label="Sensitivity"+str(res23))
#ax1[1,1].legend()


ax1[2,0].set_xlabel('Thrust Ratio (1st)')
ax1[2,0].set_ylabel('DeltaV [m/s]')
x=ThrustRatio
y=Thrust_arr1[:,17]
res31=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[2,0].plot(x,y, label="Sensitivity"+str(res31))
#ax1[1,2].legend()

ax1[2,1].set_xlabel('Thrust Ratio (2nd) ')
ax1[2,1].set_ylabel('DeltaV [m/s]')
x=ThrustRatio
y=Thrust_arr2[:,17]
res32=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[2,1].plot(x,y, label="Sensitivity"+str(res32))
#ax1[1,2].legend()


ax1[2,2].set_xlabel('Payload [t]')
ax1[2,2].set_ylabel('DeltaV [m/s]')
x=PayloadRatio
y=Pay_arr[:,17]
res33=np.polyfit(x, y, 1)
#y1=np.poly1d(res1)(x)
ax1[2,2].plot(x,y, label="Sensitivity"+str(res33))
#ax1[1,0].legend()


fig1.tight_layout()
fig1.savefig('./plot/sensitivity_summary.jpg', dpi=300)
