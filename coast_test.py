import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from cartopy import geodesic
from scipy.spatial import cKDTree
import numpy as np

fig = plt.figure()

platcar = ccrs.PlateCarree()

coord_sys = ccrs.Robinson()

bg_valid = False

ax = plt.axes(projection=coord_sys)

ax.coastlines(resolution = '110m')

ax.set_global()

plt.show(block=False)
plt.draw()
fig.canvas.draw()

background = fig.canvas.copy_from_bbox(ax.bbox)

#plt.figure()

#print dir(background)
#print background.to_string()

#plt.imshow(background)


#plt.show()

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
    #print coasts
    emp = np.zeros(coasts.shape)
    location = np.array(location)
    
    emp = emp + location
    #print emp
    loc =[]
   
    dist = geodesic.Geodesic().vec_inverse(emp,coasts)
    
             
    return coasts[np.argmin(dist[:,0])]
    

#plt.plot(x, y, 'ro', transform = ccrs.Geodetic())

locations = []


def resize(event):
    
    fig.set_figwidth(event.width)
    fig.set_figheight(event.height)
    
    plt.cla()
    plt.show(block=False)
    plt.draw()

    background = fig.canvas.copy_from_bbox(ax.bbox)

def mouse_moved(event):

    fig.canvas.restore_region(background)
    
    
    location = [0,0]
    
    if event.xdata is not None and event.ydata is not None:
        location = platcar.transform_point(event.xdata, event.ydata, coord_sys)

    
    m_pt, = plt.plot(location[0], location[1], 'ro', transform = platcar)

    #coast_index = tree.query(ccrs.PlateCarree().transform_point(location[0],location[1], platcar))[1]
    
    location = ccrs.PlateCarree().transform_point(location[0],location[1], platcar)

    coast_point = brute_neighbours(x,y,location)

    c_pt, = plt.plot(coast_point[0], coast_point[1], 'bo', transform = platcar)

    line, = plt.plot([location[0], coast_point[0]], [location[1], coast_point[1]], linewidth = 2, transform = platcar)
    
    ax.draw_artist(m_pt)
    ax.draw_artist(c_pt)
    ax.draw_artist(line)
    
    fig.canvas.blit(ax.bbox)

print total

ax.add_geometries(cartopy.feature.COASTLINE.geometries().next(), ccrs.Geodetic(), facecolor = 'blue')



cid = fig.canvas.mpl_connect('motion_notify_event', mouse_moved)

fig.canvas.mpl_connect('resize_event', resize)

plt.show()



