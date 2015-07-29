import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from cartopy import geodesic
from scipy.spatial import cKDTree
import numpy as np

fig = plt.figure()

platcar = ccrs.PlateCarree()

coord_sys = ccrs.InterruptedGoodeHomolosine()

ax = plt.axes(projection=coord_sys)

ax.coastlines(resolution = '110m')

ax.set_global()

total = 0

x = []
y = []

for string in cartopy.feature.NaturalEarthFeature('physical', 'coastline', '110m').geometries():
    for line in string:
        points = list(line.coords)
        for point in points:
            #print point
            x.append(point[0])
	    y.append(point[1])

#tree = cKDTree(zip(x,y))


def brute_neighbours(x,y,location):
    coasts = np.array(zip(x,y))
    print coasts
    emp = np.zeros(coasts.shape)
    location = np.array(location)
    
    emp = emp + location
    print emp
    loc =[]
   
    dist = geodesic.Geodesic().vec_inverse(emp,coasts)
    
             
    return coasts[np.argmin(dist[:,0])]
    
            
            
            
    


#plt.plot(x, y, 'ro', transform = ccrs.Geodetic())

locations = []
def mouse_moved(event):
    
    plt.cla()

    ax.coastlines(resolution = '110m')
    ax.set_global()
    if event.xdata is not None and event.ydata is not None:
        point = platcar.transform_point(event.xdata, event.ydata, coord_sys)
        locations.append(point)
    if len(locations) > 1:
        del locations[0]
    
    for location in locations:
        plt.plot(location[0], location[1], 'ro', transform = platcar)

        #coast_index = tree.query(ccrs.PlateCarree().transform_point(location[0],location[1], platcar))[1]
        
        location = ccrs.PlateCarree().transform_point(location[0],location[1], platcar)

        coast_point = brute_neighbours(x,y,location)

        plt.plot(coast_point[0], coast_point[1], 'bo', transform = platcar)

        plt.plot([location[0], coast_point[0]], [location[1], coast_point[1]], linewidth = 2, transform = platcar)
    
    fig.canvas.draw()

print total
#ax.set_global()

ax.add_geometries(cartopy.feature.COASTLINE.geometries().next(), ccrs.Geodetic(), facecolor = 'blue')



cid = fig.canvas.mpl_connect('motion_notify_event', mouse_moved)

#plt.show()


#locc = [50,60]

#loc = brute_neighbours(x,y,locc)
#print loc
#plt.plot(loc[0],loc[1],'ro',transform = platcar)
#plt.plot(locc[0],locc[1],'bo',transform = platcar)
#plt.plot([locc[0],loc[0]],[locc[1],loc[1]],transform = platcar)
plt.show()



