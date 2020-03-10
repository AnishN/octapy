cimport libc.math as c_math

cdef double to_degrees(double radians) nogil
cdef double to_radians(double degrees) nogil
cdef double wrap_360(double degrees) nogil
cdef double wrap_180(double degrees) nogil
cdef double wrap_90(double degrees) nogil