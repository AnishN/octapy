from octapy.libs.proj cimport *

cdef class Projection:
    cdef:
        projPJ source_crs
        projPJ target_crs

    cdef projPJ c_init_proj(self, str proj_str) except *
    cdef void c_transform(self, double x_source, double y_source, double *x_target, double *y_target) except *
    cdef void c_reverse_transform(self, double x_target, double y_target, double *x_source, double *y_source) except *