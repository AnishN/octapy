cdef double to_degrees(double radians) nogil:
    return radians * 180.0 / c_math.M_PI

cdef double to_radians(double degrees) nogil:
    return degrees * c_math.M_PI / 180.0

cdef double wrap_360(double degrees) nogil:
    if 0 <= degrees < 360:
        return degrees
    return (degrees % 360 + 360) % 360

cdef double wrap_180(double degrees) nogil:
    if -180 <degrees <= 180:
        return degrees
    return (degrees + 540) % 360 - 180

cdef double wrap_90(double degrees) nogil:
    if -90 <= degrees <= 90:
        return degrees
    return c_math.fabs((degrees % 360 + 270) % 360 - 180) - 90 