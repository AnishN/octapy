from libc.stdio cimport *

cdef class Projection:
    
    def __cinit__(self, str source_crs_str, str target_crs_str):
        self.source_crs = self.c_init_proj(source_crs_str)
        self.target_crs = self.c_init_proj(target_crs_str)
    
    def __dealloc__(self):
        pj_free(self.source_crs)
        pj_free(self.target_crs)
    
    cdef projPJ c_init_proj(self, str proj_str) except *:
        cdef:
            projPJ proj
            int err
            str err_str
        
        proj = pj_init_plus(proj_str.encode("utf-8"))
        if proj == NULL:
            err = pj_get_errno_ref()[0]
            err_str = pj_strerrno(err).decode()
            raise ValueError("Projection: could not initialize proj_str due to {0}".format(err_str))
        return proj

    cdef void c_transform(self, double x_source, double y_source, double *x_target, double *y_target) except *:
        cdef:
            int err
            str err_str
        
        x_target[0] = x_source
        y_target[0] = y_source
        err = pj_transform(self.source_crs, self.target_crs, 1, 1, x_target, y_target, NULL)
        if err != 0:
            err_str = pj_strerrno(err).decode()
            raise ValueError("Projection: could not perform transform due to {0}".format(err_str))

    cdef void c_reverse_transform(self, double x_target, double y_target, double *x_source, double *y_source) except *:
        cdef:
            int err
            str err_str
        
        x_source[0] = x_target
        y_source[0] = y_target
        err = pj_transform(self.target_crs, self.source_crs, 1, 1, x_source, y_source, NULL)
        if err != 0:
            err_str = pj_strerrno(err).decode()
            raise ValueError("Projection: could not perform reverse transform due to {0}".format(err_str))