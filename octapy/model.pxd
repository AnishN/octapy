from libc.stdlib cimport *
cimport libc.math as c_math
cimport cpython.datetime as dt
from octapy.hycom_data_set cimport HycomDataSet
from octapy.projection cimport Projection
from octapy.particle cimport *
from octapy.libs.netcdf cimport *
from octapy.utils cimport *

cdef class Model:
    cdef HycomDataSet data_set
    cdef Projection projection
    cdef size_t leaf_size
    cdef double power

    cdef int _c_get_dim_id(self, int nc_id, char *name) except -1
    cdef int _c_get_var_id(self, int nc_id, char *name) except -1