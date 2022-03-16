import dash
import sys
from PandasHandler import *
from TrajectoryPlotter import *

import plotly.express as px
import plotly.graph_objects as go


#SETTINGS===============================================================================================
#Option for csv entry
Casename = sys.argv[1] #Enter path/filename of FCcalc csv output file in input argument e.g. 'out_XYZ_Study N_90 E_0.csv'
df_Data         = pd.read_csv('./plot/'+Casename)

#Option for Direct run FCcalc
#import FlightCalcPolar as FCcalc
#Case = FCcalc.case()
#df_Data = FCcalc.output()


#PLOTS==================================================================================================
#Initial Settings for Earth plot========================================================================
filePathExcel      = 'earth.xlsx'
dfSet              = PandasHandler.readAllSheets_Excel(filePathExcel)
df                 = dfSet[0]

#Plotting
colNameX      = 'Time[min]'
colNameYList1 = ['Density[kg/m3]','Altitude[km]','MagFldInt in East[nT]','MagFldInt in North[nT]','MagFldInt in Up[nT]']

from plotly.subplots import make_subplots

colNames = list(df.columns)
imax     = len(colNames)

fig = make_subplots()

Earth = TrajectoryPlotter.getEarth('#325bff',imax=40)

fig.add_trace( go.Scatter3d( x=df["Trace 0, x"], y=df["Trace 0, y"], z=df["Trace 0, z"], mode='lines', marker=dict(size=1, color='grey'   ), name="", hoverinfo='skip'))
fig.add_trace( go.Scatter3d( x=df["Trace 1, x"], y=df["Trace 1, y"], z=df["Trace 1, z"], mode='lines', marker=dict(size=1, color='grey'   ), name="", hoverinfo='skip'))
fig.add_trace( Earth)
fig.add_trace( go.Scatter3d( x=df["Trace 3, x"], y=df["Trace 3, y"], z=df["Trace 3, z"], mode='lines', marker=dict(size=1, color='white'  ), name="", hoverinfo='skip'))
fig.add_trace( go.Scatter3d( x=df_Data["ECI-X[km]"], y=df_Data["ECI-Y[km]"], z=df_Data["ECI-Z[km]"], mode='lines', marker=dict(size=2, color='red'), line=dict(color='red',width=5), name="Trajectory"))
fig.update_layout(scene=dict(xaxis=dict(title="X-ECI[km]", showgrid=True, showline=True, zeroline=False, backgroundcolor='black', color='grey', gridcolor='grey')))
fig.update_layout(scene=dict(yaxis=dict(title="Y-ECI[km]", showgrid=True, showline=True, zeroline=False, backgroundcolor='black', color='grey', gridcolor='grey')))
fig.update_layout(scene=dict(zaxis=dict(title="Z-ECI[km]", showgrid=True, showline=True, zeroline=False, backgroundcolor='black', color='grey', gridcolor='grey')))
fig.update_layout(scene=dict(bgcolor="black"))
#fig.update_layout(showlegend=False)

fig.show()
fig.write_html('./plot/plotly3d_'+Casename+'.html')