import matplotlib.pyplot as plt
import cartopy.crs as ccrs

ax = plt.axes(projection=ccrs.Mercator())
ax.coastlines()

ax.tissot(lat_num = 8, lon_num = 4, facecolor='green', alpha=0.3)
plt.show()
