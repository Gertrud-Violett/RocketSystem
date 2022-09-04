import numpy as np
import sympy
from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt

#CEA Full Equilibrim combustion chamber==========================================
#Initial setting ================================================================
g0=9.809 #Gravity constant
At=0.215**2/4*np.pi #Throat dia in m2   

#Syntax: Comp Pressure, Expansion Ratio, Mixture Ratio, Ambient Pressure, Isp Efficiency, Expected Thrust (Can be any, for checking only)
def ceathrust(Pci,ERi,MRi,Pambient,Ispeff,Thrust):
    ceares=CEA_Obj(
        oxName='LOX',
        fuelName='CH4',
        isp_units='sec',
        cstar_units='m/s',
        pressure_units='MPa',
        temperature_units='K',
        sonic_velocity_units='m/s',
        density_units='kg/m^3',
        specific_heat_units='kJ/kg-K'
        )

    (Ispvac, Cstar, Tcomb) = ceares.get_IvacCstrTc(
        Pc=Pci, MR=MRi, eps=ERi,
        frozen=0, frozenAtThroat=0
        )
    (Ispamb, mode) = ceares.estimate_Ambient_Isp(
        Pc=Pci, MR=MRi, eps=ERi,
        Pamb=Pambient, frozen=0, frozenAtThroat=0
        )
    Mach = ceares.get_MachNumber(
        Pc=Pci, MR=MRi, eps=ERi,
        frozen=0, frozenAtThroat=0
        )
    (Cc,Ct,Ce) = ceares.get_SonicVelocities(
        Pc=Pci, MR=MRi, eps=ERi,
        frozen=0, frozenAtThroat=0
        )
    (mw,gam) = ceares.get_Chamber_MolWt_gamma(
        Pc=Pci, MR=MRi, eps=ERi
        )
    Pe = Pci/(ceares.get_PcOvPe(
        Pc=Pci, MR=MRi, eps=ERi,
        frozen=0, frozenAtThroat=0
        ))
    
    Ve=Mach*Ce
    Ispact=Ispamb*Ispeff
    Veq=Ispamb*g0
    Mdot_req = Thrust*10**3/Veq/Ispeff
    Mdot_th = At*Pci*10**6/Tcomb**0.5*np.sqrt(gam/8.31*mw*10**-3)*((gam+1)/2)**(-(gam+1)/2/(gam-1)) #kg/s
    Thrust_th =Mdot_th*Veq*Ispeff
    print(Ispvac, Ispamb, Ispact, Cstar, Tcomb, Pe, Mdot_req, Mdot_th, Thrust_th)
    return(Ispact, Mdot_th, Thrust_th)



#For Debug and Initial At setting===============================================

Ispeff=0.98 #Isp efficiency
Pc=33.0 #Combustion chamber total pressure in MPa
ER=34.0 #Expansion Ratio
MR=3.6 # Mixture Ratio
Ae=At*ER
Pambient=0.1 #MPa
Thrust=2300 #kNs

#Obtain At size
Atarr = np.arange(0.2**2/4*np.pi,0.22**2/4*np.pi, 0.005**2/4*np.pi)
Mdotarr=[]
for Ats in Atarr:
    At = Ats
    (Isp, Mdot, Thrust) = ceathrust(Pc,ER,MR,Pambient,Ispeff,Thrust)
    Mdotarr=np.append(Mdotarr,Mdot)
    #print(Mdot)
plt.plot(Atarr,Mdotarr)
plt.xlabel('Throat Dia. At')
plt.ylabel('Mass flow rate Mdot')
plt.show()


#Check MR and Isp Relation
MRarr = np.arange(3.4, 3.9, 0.02)
Isparr=[]
for MRs in MRarr:
    (Isp, Mdot, Thrust) = ceathrust(Pc,ER,MRs,Pambient,Ispeff,Thrust)
    Isparr=np.append(Isparr,Isp)
plt.plot(MRarr,Isparr)
plt.plot([MR,MR],[Isp-5,Isp+5])
plt.xlabel('Mixture ratio MR')
plt.ylabel('Isp SL')
plt.show()



"""
#OPTION FROZEN Classic Rocket equation Manual Calc mode=========================
#syntax: stage no, atmospheric pressure, Isp efficiency
Pc=30.0*10e6 #Combustion chamber total pressure
At=0.2**2/4*np.pi #Throat dia in m2   
ER=30.0 #Expansion Ratio
EtaIsp=0.9 #Isp efficiency
MR=3.8 # Mixture Ratio
Ae=At*ER

def thrustclassic(presatm,Ispeff):
    #Universal constants
    g0=9.81 #Gravity constant
    gamma=1.2 #heat ratio
    Cp=2.3
    Rconst=(Cp-Cp/gamma)*10e3 #kJ/kgK gas constant
    Tt=2300 #Kelvins
    GAMMA=(gamma+1)/2/(gamma-1)
    
    #Propellant flow rate
    Mdot = At*Pc/Tt**0.5*np.sqrt(gamma/Rconst)*((gamma+1)/2)**(-GAMMA) #kg/s

    #Exhaust Velocity
    x = sympy.Symbol('x')
    Mache = sympy.solve((((gamma+1)/2)**(-GAMMA))*(((1+(gamma-1)/2)*x**2)**(GAMMA))/x - ER, x)
    Mache = [n for n in Mache if n.is_real][0]    

    #Exit T,P,v
    Te=Tt*1/(1+(gamma-1)/2*Mache**2)
    Pe=Pc*(1+(gamma-1)/2*Mache**2)**(-gamma/(gamma-1))
    Ve=Mache*(gamma*Rconst*Te)**0.5

    #Thrust
    Thrust=(Mdot*Ve+(Pe-presatm)*Ae)/10e3*Ispeff  #Thrust in kNs
    Veq=(Ve+(Pe-presatm)*At/Mdot)*Ispeff  #Thrust Vel in m/s
    Isp=Veq/g0  #Isp in s
    return(Thrust,Veq,Isp)
"""


