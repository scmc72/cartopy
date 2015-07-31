# (C) British Crown Copyright 2014, Met Office
#
# This file is part of cartopy.
# 
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.



#a change

"""
This module defines the core CRS class which can interface with Proj.4.
The CRS class is the base-class for all projections defined in :mod:`cartopy.crs`.

"""
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
import numpy as np
cimport numpy as np
from cython.parallel cimport prange


# TODO: Support proj4 <4.9 with spherical approximations?


cdef extern from "geodesic.h":
    cdef struct geod_geodesic:
        pass

    ctypedef geod_geodesic* geodesic_t

    void geod_init(geodesic_t, double, double)
    void geod_direct(geodesic_t, double, double, double, double,
                     double*, double*, double*) nogil
    
    void geod_inverse(geodesic_t, double, double, double, double,
                      double*, double*, double*) nogil


cdef class Geodesic:
    cdef geod_geodesic* geod

    def __cinit__(self, radius=6378137, flattening=1/298.257223563):
        # allocate some memory (filled with random data)
        self.geod = <geod_geodesic*> PyMem_Malloc(sizeof(geod_geodesic))
        if not self.geod:
            raise MemoryError()
        geod_init(self.geod, radius, flattening)

    def direct(self, lon0, lat0, azi0, distance):
        cdef double lat, lon, azi
        geod_direct(self.geod, lat0, lon0, azi0, distance, &lat, &lon, &azi)
        return lon, lat, azi
        
    def vec_direct(self, points, azimuths, distances):
        
        cdef int n_points, i
        cdef double[:,:] pts
        cdef double[:] azims, dists
        
        pts = np.array(points, dtype = np.float64).reshape((-1,2))
        azims = np.array(azimuths, dtype = np.float64).reshape(-1)
        dists = np.array(distances, dtype = np.float64).reshape(-1)
        
        n_points = max(pts.shape[0], azims.size, dists.size)
        
        try:
            tmp = np.zeros((n_points,2))
            tmp[:,0] += pts[:,0]
            tmp[:,1] += pts[:,1]
            
            pts = tmp
            
            azims = np.zeros(n_points) + azims[:]
            
            dists = np.zeros(n_points) + dists[:]
            
        except ValueError:
            raise ValueError("Inputs must have common length n or length one.")
        
        #take in a lat-long array; an azimuth array and a distance array
        cdef double[:,:] return_pts = np.empty((n_points,3))
        
        cdef double[:] lat, lon, azi
        
        lat = np.empty(n_points)
        lon = np.empty(n_points)
        azi = np.empty(n_points)
        
        with nogil:
            for i in prange(n_points):
            
                geod_direct(self.geod, pts[i,1], pts[i,0],
                            azims[i], dists[i], &lat[i], &lon[i], &azi[i])
                return_pts[i,0] = lon[i]
                return_pts[i,1] = lat[i]
                return_pts[i,2] = azi[i]
            
            
        return np.array(return_pts)

    def inverse(self, lon0, lat0, lon1, lat1):
        cdef double dist, azi0, azi1
        geod_inverse(self.geod, lat0, lon0, lat1, lon1, &dist, &azi0, &azi1)
        return dist, azi0, azi1
        
    def vec_inverse(self, points, endpoints):
        
        cdef int n_points, i
        
        cdef double[:,:] pts, epts
        
        pts = np.array(points, dtype = np.float64).reshape((-1,2))
        epts =  np.array(endpoints, dtype = np.float64).reshape((-1,2))
        
        n_points = max(pts.shape[0], epts.shape[0])
        
        try:
            tmp = np.zeros((n_points,2))
            tmp[:,0] += pts[:,0]
            tmp[:,1] += pts[:,1]
            
            pts = tmp
            
            tmp = np.zeros((n_points,2))
            tmp[:,0] += epts[:,0]
            tmp[:,1] += epts[:,1]
            
            epts = tmp
            
        except ValueError:
            raise ValueError("Inputs must have common length n or length one.")
        
        cdef double[:,:] results = np.empty((n_points, 3))
        
        cdef double[:] dist, azi0, azi1
        
        dist = np.empty(n_points)
        azi0 = np.empty(n_points)
        azi1 = np.empty(n_points)
        
        with nogil:
            for i in prange(n_points):
                
                geod_inverse(self.geod, pts[i,0], pts[i,1], epts[i,0],
                             epts[i,1], &dist[i], &azi0[i], &azi1[i])
                
                results[i,0] = dist[i]
                results[i,1] = azi0[i]
                results[i,2] = azi1[i]
        
        return np.array(results)
    
    def circle(self, double lon, double lat, double distance, int n_samples=180, endpoint=False):
        #cdef double lat_o, lon_o, azi_o
        cdef int i

        # Put the input arguments into c-typed values.        
        center = np.array([lon, lat]).reshape((1,2))
        distance_m = np.asarray(distance).reshape(1)
        
        print distance_m.shape

        #result = np.empty([n_samples, 2], dtype=np.double)
        azimuths = np.linspace(360., 0., n_samples, endpoint=endpoint).astype(np.double)
        
        geod = Geodesic()

        return geod.vec_direct(center, azimuths, distance_m)[:, 0:2]
        

    def __dealloc__(self):
        PyMem_Free(self.geod)

def main():
    g = Geodesic()

    #print g.direct(-73.78, 40.64, 45.0, 10e6)
    # 32.64284433 49.01103958
    
    lat1, lon1 = 40.6, -73.8 # JFK Airport
    lat2, lon2 = 51.6, -0.5  # LHR Airport
    
    #print g.inverse(lon1, lat1, lon2, lat2)
    # 5551759.4003186785 KM.
    #print g.circle(0, 0, 1000000, 10)
    
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    
    circle = g.circle(0, 52, 500000, n_samples=360, endpoint=True)
    print circle
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()
    plt.plot(circle[:,0], circle[:,1], transform=ccrs.Geodetic())
    plt.show()
