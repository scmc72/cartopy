from cartopy import geodesic
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
import numpy as np


g = geodesic.Geodesic()


ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
ax.coastlines()
ax.set_global()

geod = geodesic.Geodesic()

radius_km = 200000
n_samples = 80

geoms = []
for lat in np.linspace(-80, 80, 40):
    for lon in np.linspace(-180, 180, 20, endpoint=False):
        lonlat = geod.circle(lon, lat, radius_km)
        
       # geoms.append(lonlat)
        geoms.append(sgeom.Polygon(zip(lonlat[:,0],lonlat[:,1])))
        
#lonlat = geod.circle(90, 70, radius_km)
#print lonlat
#print geom
#geoms = np.asarray(geoms)
#print geoms.shape
ax.add_geometries(geoms, ccrs.Geodetic(), facecolor='blue', alpha=0.7)
#plt.plot(lonlat[:,0],lonlat[:,1],'ro',transform=ccrs.Geodetic())
plt.show()
