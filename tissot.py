from cartopy import geodesic
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
import numpy as np


g = geodesic.Geodesic()


ax = plt.axes(projection=ccrs.Robinson())
ax.coastlines()

geod = geodesic.Geodesic()

radius_km = 500000
n_samples = 80

geoms = []
for lat in np.linspace(-80, 80, 10):
    for lon in np.linspace(-180, 180, 7, endpoint=False):
        lonlat = geod.circle(lon, lat, radius_km)
        
        geoms.append(sgeom.Polygon(zip(lonlat[:,0],lonlat[:,1])))


ax.add_geometries(geoms, ccrs.Geodetic(), facecolor='blue', alpha=0.7)

plt.show()
