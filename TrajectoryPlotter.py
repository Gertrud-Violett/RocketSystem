import numpy as np
import pymap3d as pm
import plotly.graph_objects as go

class TrajectoryPlotter:

    def __init__(self):
        print("This is Constructor")
        
    def __del__(self):
        print("This is Destructor")

    @staticmethod
    def getEarth(clr,imax=30): 
        latitude   = np.linspace( -90.0, 90.0, imax)
        longitude  = np.linspace(-180.0,180.0, imax)
        lat, lon   = np.meshgrid(latitude, longitude)

        x0  = np.zeros([imax,imax])
        y0  = np.zeros([imax,imax])
        z0  = np.zeros([imax,imax])
        alt = 0.0
        for i1 in range(imax):
            lat = latitude[i1]
            for i2 in range(imax):
                lon = longitude[i2]
                x0[i1,i2],y0[i1,i2],z0[i1,i2] = pm.geodetic2ecef(lat, lon, alt)

        x0 = x0/1000.0
        y0 = y0/1000.0
        z0 = z0/1000.0
    
        # Set up trace
        trace= go.Surface(x=x0, y=y0, z=z0, colorscale=[[0,clr], [1,clr]], hoverinfo='skip')
        trace.update(showscale=False)
        return trace
