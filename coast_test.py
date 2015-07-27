import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

from scipy.spatial import KDTree

ax = plt.axes(projection=ccrs.PlateCarree())

ax.coastlines(resolution = '50m')

ax.set_global()

total = 0

x = []
y = []

for string in cartopy.feature.NaturalEarthFeature('physical', 'coastline', '50m').geometries():
    for line in string:
        points = list(line.coords)
        for point in points:
            #print point
            x.append(point[0])
	    y.append(point[1])

tree = KDTree(zip(x,y))

#plt.plot(x, y, 'ro', transform = ccrs.Geodetic())

locations = [(-3.4951, 50.7287), (10,10),(206.769-360,19.911), (16,5)]

for location in locations:
    plt.plot(location[0], location[1], 'ro', transform = ccrs.Geodetic())

    coast_index = tree.query(location)[1]

    coast_point = tree.data[coast_index]

    plt.plot(coast_point[0], coast_point[1], 'bo', transform = ccrs.Geodetic())

    plt.plot([location[0], coast_point[0]], [location[1], coast_point[1]], transform = ccrs.Geodetic())

print total
#ax.set_global()

ax.add_geometries(cartopy.feature.COASTLINE.geometries().next(), ccrs.Geodetic(), facecolor = 'blue')

plt.show()
