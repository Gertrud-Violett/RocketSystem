# RocketSystem
Rocket System Calculation Tools

## DIR:Performance Calculator
Simple 3DoF trajctory simulation tool to study Rocket vehicle launch capability.

### RunTrajectoryPlotter_polar.py
Run after FlightCalcPolar to create 3D trajctory plot.
```bash
python RunTrajectoryPlotter_polar.py "result filename.csv"
>>
```
![output1](https://user-images.githubusercontent.com/96937287/158642949-8f863e62-efaf-47c3-8ca8-49ae87976259.png)

### FlightCalcPolar.py
Main Calculation script. Results output to filename.csv based on input from settings.toml 
```bash
python FlightCalcPolar.py
...       ...          ...          ...           ...  ...        ...           ...            ...            ...
 Accel: 0.476672826979712 8.719483392066946e-23 -1.0870483799769239e-07
 Stage:3,Time[s]:499.9,Weight[t]:100, Drag[kN]:0.00,Thrust[kN]:0,Alt[m]:263830, Accel[m/s2]:0.86,Velocity[m/s]:7987
 Downrange[km]:1987
 Weight[t]:100
 Drag[kN]:0.00
 Alt[km]:263.83
 Orbit Speed [m/s]7981.52
 Vehicle ABS speed [m/s]:7987
 Vehicle Airspeed [Mach]:7.503

 1st Time[s],alt[km]: 171,95
 2nd Time[s],alt[km]: 421,242
 Payload Time[s],alt[km]: 500,264
 Theta[deg]:90.0
 DeltaV_th[m/s]:9287.78
 DeltaV1st[m/s]:3796.95
  Landing DeltaV[m/s]:3355.10

      Time[s]    ECI-X[km]    ECI-Y[km]     ECI-Z[km]  ...    drag[N]    Weight[kg]    Altitude[m]  Velocity[m/s]
0         0.0  6371.000000     0.000000  3.901112e-13  ...   0.000000  5.100000e+06       0.000000     465.000000
...       ...          ...          ...           ...  ...        ...           ...            ...            ...

5000    499.9  6339.471606  1957.668819  4.062680e-13  ...   0.000000  1.000000e+05  263830.016153    7987.184431
```
![Study N_90 E_0_plot](https://user-images.githubusercontent.com/96937287/158643120-e0c20ddf-f21f-4745-9fb4-12f9af693b5e.jpg)
![Study N_90 E_0_Trajectory2D](https://user-images.githubusercontent.com/96937287/158643153-4560d1e4-2633-4023-8a36-e4877d8737fa.jpg)


### settings.toml
Setting file for Main Calculation. Create new table for each stage
<span style="color: blue; ">
[ENGINE]
  [ENGINE.Stage1]
    Throat_Dia = float
    Isp_efficiency = float
    Combustion_pressure = float
    Nozzle_ExpansionRatio = float
    OF_mixture = float
    Target_thrust = float
    Engine_qty = int

  [ENGINE.Stage2]

[VEHICLE]
  [VEHICLE.StageX]
    Structural_Efficiency = 0.96
    Initial_Weight = 4000000.0
    Fuel_Remaining = 0.1
    Fuselage_Dia = 9.0

[ENVIRONMENT]
  gravity = 9.80665
  barometric_pressure = 0.1013

[VECTOR]
  Latitude = 0.0
  Longitude = 0.0
  altitude = 0.0
  radius=6371000

[CALC_SETTING]
  Casename = "Demo_run"
  timestep = 0.1
  timelimit = 500
</span>

#### attitude_beta.csv
beta (phi) thrust angle history file. Enter thrust angle relative to longitudinal direction.(0deg=vertical, 90deg=horizontal to ground)

#### attitude_gamma.csv
gamma (theta) thrust angle history file. Enter thrust angle relative to latitudinal direction.(0deg=vertical, 90deg=horizontal to ground)

### USSTDATM.py
FlightCalcPolar subroutine. US Standard Atmosphere and Aero drag specification file.
Default settings M-V aero used. 
https://jaxa.repo.nii.ac.jp/?action=repository_action_common_download&item_id=5551&item_no=1&attribute_id=31&file_no=1

```bash
python USSTDATM.py
>>
```

### rocketthrustcalc.py
FlightCalcPolar subroutine. used to run rocketCEA
```bash
python rocketthrustcalc.py
>>
```

### vehicle.py
FlightCalcPolar subroutine. used to calcualte next time step specs with 4th order approximation.
```bash
python vehicle.py
>>
```

### others
#### PandasHandler.py
plotting subroutine forked from bluedack https://github.com/bluedack-space/TrajectoryPlotter

#### TrajectoryPlotter.py
plotting subroutine forked from bluedack https://github.com/bluedack-space/TrajectoryPlotter

#### earth.xlsx
earth shoreline plotting file. do not edit


### result files:
#### DIR:plot
Save directory for results
