cimport cpython.datetime as dt

cdef class HycomDataSet:
    cdef:
        str sub_model
        str data_dir_path
        readonly float[::1] lats
        readonly float[::1] lons
        readonly float[::1] depths
        readonly dt.datetime base_date
        readonly dt.datetime start
        readonly dt.datetime end
        readonly dt.datetime start_limit
        readonly dt.datetime end_limit