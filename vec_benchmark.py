import timeit

setup="""
from cartopy.geodesic import Geodesic
import numpy as np
#random latlon points
points = (np.random.rand(1000,2)-0.5)
points[:,0] *= 180
points[:,1] *= 90

endpoints = (np.random.rand(1000,2)-0.5)
endpoints[:,0] *= 180
endpoints[:,1] *= 90

distance = np.random.rand(1000)*100000

azi = np.random.rand(1000)*360

geod = Geodesic()


def v_d_benchmark(geod, points, azi, distance):
    
    return geod.vec_direct(points,azi,distance)
    
def d_benchmark(geod, points, azi, distance):

    lst = []
    
    for i in range(points.shape[0]):
        
        lst.append(geod.direct(points[i][0], points[i][1], azi[i], distance[i]))
        
    return lst
        
def v_i_benchmark(geod, points, endpoints):
    
    return geod.vec_inverse(points,endpoints)
    
def i_benchmark(geod, points, endpoints):

    lst = []
    
    for i in range(points.shape[0]):
        
        lst.append(geod.direct(points[i][0], points[i][1], endpoints[i][0], endpoints[i][1]))
    
    return lst"""
    
def test_outputs():
    from cartopy.geodesic import Geodesic
    import numpy as np
    #random latlon points
    points = (np.random.rand(5,2)-0.5)
    points[:,0] *= 180
    points[:,1] *= 90
    
    endpoints = (np.random.rand(5,2)-0.5)
    endpoints[:,0] *= 180
    endpoints[:,1] *= 90
    
    print points
    
    
    distance = np.random.rand(5)*100000
    
    azi = np.random.rand(5)*360
    
    geod = Geodesic()
    print "\n TESTING OUTPUTS..."
    print "\nDIRECT\n"
    print "VECTORISED"
    print v_d_benchmark(geod, points, azi, distance)
    print "-"
    print "LOOP"
    print d_benchmark(geod, points, azi, distance)
    print "\nINVERSE\n"
    print "VECTORISED"
    print v_i_benchmark(geod, points, endpoints)
    print "-"
    print "LOOP"
    print i_benchmark(geod, points, endpoints)
    
def v_d_benchmark(geod, points, azi, distance):
    
    return geod.vec_direct(points,azi,distance)
    
def d_benchmark(geod, points, azi, distance):

    lst = []
    
    for i in range(points.shape[0]):
        
        lst.append(geod.direct(points[i][0], points[i][1], azi[i], distance[i]))
        
    return lst
        
def v_i_benchmark(geod, points, endpoints):

    return geod.vec_inverse(points,endpoints)
    
def i_benchmark(geod, points, endpoints):

    lst = []
    
    for i in range(points.shape[0]):
        
        lst.append(geod.inverse(points[i][1], points[i][0], endpoints[i][1], endpoints[i][0]))
    
    return lst
    
test_outputs()

print "\nTESTING SPEEDS... \n\n"

print "VECTORISED DIRECT"

print timeit.timeit("v_d_benchmark(geod, points, azi, distance)", setup=setup, number = 1000)

print "-----------------"

print "BASIC LOOP DIRECT"

print timeit.timeit("d_benchmark(geod, points, azi, distance)", setup=setup, number = 1000)

print "############################"
print "############################"

print "VECTORISED INVERSE"

print timeit.timeit("v_i_benchmark(geod, points, endpoints)", setup=setup, number = 1000)

print "-----------------"

print "BASIC LOOP INVERSE"

print timeit.timeit("i_benchmark(geod, points, endpoints)", setup=setup, number = 1000)
