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


# TODO: Support proj4 <4.9 with spherical approximations?


cdef extern from "geodesic.h":
    cdef struct geod_geodesic:
        pass

    ctypedef geod_geodesic* geodesic_t

    void geod_init(geodesic_t, double, double)
    void geod_direct(geodesic_t, double, double, double, double,
                     double*, double*, double*)
    
    void geod_inverse(geodesic_t, double, double, double, double,
                      double*, double*, double*)


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
        
    def vec_direct(self, lon0lat0, azi0, distance):
        
        if type(lon0lat0) != np.ndarray or type(azi0) != np.ndarray or type(distance) != np.ndarray:
            print 'type error: inputs must be type: "numpy.ndarray"'
      
        
        if lon0lat0.shape[1] != 2:
            print 'ERROR: lon0lat0 must have shape: (npoints,2)'
        
        if lon0lat0.shape[0] != azi0.size or lon0lat0.shape[0] != distance.size or azi0.size != distance.size:
            print 'ERROR: number of points must be the same for all inputs: array sizes do not match!'
            
        
            
        
        
        #take in a lat-long array; an azimuth array and a distance array
        return_pts = np.empty((lon0lat0.shape[0],3), dtype = np.float32)
        
        cdef double lat, lon, azi
        
        for i in range(lon0lat0.shape[0]):
        
            geod_direct(self.geod, lon0lat0[i,1], lon0lat0[i,0], azi0[i], distance[i], &lat, &lon, &azi)
            return_pts[i,0] = lon
            return_pts[i,1] = lat
            return_pts[i,2] = azi
            
            
        return return_pts
        
    def inverse(self, lon0, lat0, lon1, lat1):
        cdef double dist, azi0, azi1
        geod_inverse(self.geod, lat0, lon0, lat1, lon1, &dist, &azi0, &azi1)
        return dist, azi0, azi1
    
    def circle(self, lon, lat, distance, int n_samples=180, endpoint=False):
        cdef double lat_o, lon_o, azi_o
        cdef double center_lon, center_lat, distance_m
        cdef int i

        # Put the input arguments into c-typed values.        
        center_lon, center_lat = lon, lat
        distance_m = distance 

        result = np.empty([n_samples, 2], dtype=np.double)
        azimuths = np.linspace(360., 0., n_samples, endpoint=endpoint).astype(np.double)

        for i in range(n_samples):
            geod_direct(self.geod, center_lat, center_lon, azimuths[i], distance_m, &lat_o, &lon_o, &azi_o)
            result[i, :] = lon_o, lat_o
        return result

    def __dealloc__(self):
        PyMem_Free(self.geod)


def main():
    g = Geodesic()

    #print g.direct(-73.78, 40.64, 45.0, 10e6)
    # 32.64284433 49.01103958
    
    lat1, lon1 = 40.6, -73.8 # JFK Airport
    lat2, lon2 = 51.6, -0.5  # LHR Airport
    
    print g.inverse(lon1, lat1, lon2, lat2)
    # 5551759.4003186785 KM.
    print g.circle(0, 0, 1000000, 10)
    
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    
    x, y = g.circle(0, 52, 500000, n_samples=360, endpoint=True).T
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    plt.plot(x, y, transform=ccrs.Geodetic())
    plt.show()
