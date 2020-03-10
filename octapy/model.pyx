import numpy as np
from scipy.spatial.ckdtree import *
import datetime as pydt

point_dtype = [
    ("time", "f8"),
    ("lat", "f4"),
    ("lon", "f4"),
    ("depth", "f4"),
    ("u", "f4"),
    ("v", "f4"),
    ("w", "f4"),
    ("temp", "f4"),
    ("sal", "f4"),
]

cdef class Model:

    def __cinit__(self, HycomDataSet data_set, Projection projection, size_t leaf_size=9, double power=1.0):
        self.data_set = data_set
        self.projection = projection
        self.leaf_size = leaf_size
        self.power = power

    def __dealloc__(self):
        pass
    
    def _create_grid_data(self):
        cdef:
            float[:, ::1] grid_data
            size_t i, j
            size_t index
            double lon, lat
            double x, y
            size_t num_lons = self.data_set.lons.shape[0]
            size_t num_lats = self.data_set.lats.shape[0]

        grid_data = np.zeros((num_lons * num_lats, 2), dtype=np.float32)
        for i in range(num_lons):
            for j in range(num_lats):
                lon = to_radians(self.data_set.lons[i])
                lat = to_radians(self.data_set.lats[j])
                self.projection.c_transform(lon, lat, &x, &y)
                index = i * num_lats + j
                grid_data[index, 0] = <float>x
                grid_data[index, 1] = <float>y
        return grid_data

    cdef int _c_get_dim_id(self, int nc_id, char *name) except -1:
        cdef:
            int err
            int dim_id = -1
        err = nc_inq_dimid(nc_id, name, &dim_id)
        if err:
            raise ValueError("{0} dimension not present".format(name.decode()))
        return dim_id

    cdef int _c_get_var_id(self, int nc_id, char *name) except -1:
        cdef:
            int err
            int var_id = -1
        err = nc_inq_varid(nc_id, name, &var_id)
        if err:
            raise ValueError("{0} variable not present".format(name.decode()))
        return var_id

    def _get_netcdf_data(self, str file_path, float[::1] u_data, float[::1] v_data, 
            float[::1] w_data, float[::1] temp_data, float[::1] sal_data):
        cdef:
            int nc_id
            bytes b_file_path
            int ret
            int depth_id, lat_id, lon_id, mt_id
            int u_id, v_id, w_id, temp_id, sal_id

        b_file_path = <bytes>file_path.encode("utf-8")
        ret = nc_open(b_file_path, NC_NOWRITE, &nc_id)
        if ret != NC_NOERR:
            raise ValueError("Cannot open netcdf file {0}".format(file_path))
        
        #Check if necessary dimensions are present
        #Should really do some size validation, this is skipped here
        depth_id = self._c_get_dim_id(nc_id, b"Depth")
        lat_id = self._c_get_dim_id(nc_id, b"Latitude")
        lon_id = self._c_get_dim_id(nc_id, b"Longitude")
        mt_id = self._c_get_dim_id(nc_id, b"MT")

        #Check if necessary variables are present
        u_id = self._c_get_var_id(nc_id, b"u")
        v_id = self._c_get_var_id(nc_id, b"v")
        w_id = self._c_get_var_id(nc_id, b"w_velocity")
        temp_id = self._c_get_var_id(nc_id, b"temperature")
        sal_id = self._c_get_var_id(nc_id, b"salinity")
        
        #fill memoryviews with variable data
        nc_get_var_float(nc_id, u_id, &u_data[0])
        nc_get_var_float(nc_id, u_id, &v_data[0])
        nc_get_var_float(nc_id, u_id, &w_data[0])
        nc_get_var_float(nc_id, u_id, &temp_data[0])
        nc_get_var_float(nc_id, u_id, &sal_data[0])
        
        nc_close(nc_id)

    def run(self, ParticleC[::1] particles):
        cdef:
            size_t i, j, k
            size_t n = particles.shape[0]
            ParticleC *particle
            size_t num_points
            PointC[::1] trajectory
            #PointC *point
            #PointC *next_point
            PointC *point_a
            PointC *point_b
            double date_time
            object grid
            long index
            float[:, ::1] grid_data
            double[::1] distances
            long[::1] indices
            double[::1] weights
            dt.datetime curr_date
            str curr_path
            size_t data_shape
            double x, y#assumed in meters
            double lat, lon#IN RADIANS, NOT DEGREES

            #Assumes a uniform grid of data for now for each time point
            size_t num_depths = self.data_set.depths.shape[0]
            size_t num_lats = self.data_set.lats.shape[0]
            size_t num_lons = self.data_set.lons.shape[0]
            size_t num_total = num_depths * num_lats * num_lons
            tuple var_shape = (num_total,)
            float[::1] u_data = np.zeros(shape=var_shape, dtype=np.float32)
            float[::1] v_data = np.zeros(shape=var_shape, dtype=np.float32)
            float[::1] w_data = np.zeros(shape=var_shape, dtype=np.float32)
            float[::1] temp_data = np.zeros(shape=var_shape, dtype=np.float32)
            float[::1] sal_data = np.zeros(shape=var_shape, dtype=np.float32)
            float u, v, w, temp, sal
            double weights_sum
            double weight
            float u_w, v_w, w_w, temp_w, sal_w#weighted
            double seconds_per_hour = 3600

        grid_data = self._create_grid_data()
        grid = cKDTree(grid_data, balanced_tree=False)#otherwise build time for kd-tree was MASSIVE!!!

        for i in range(n):
            particle = &particles[i]
            num_points = <size_t>((particle.end_time - particle.start_time) * 24)
            trajectory = np.zeros(shape=(num_points,), dtype=point_dtype)
            
            for j in range(0, num_points - 1):
                point_a = &trajectory[j]
                point_b = &trajectory[j + 1]#this is the one that gets adjusted based off the current data
                
                curr_time = particle.start_time + (j / 24.0)
                curr_date = self.data_set.base_date + pydt.timedelta(days=curr_time)
                curr_path = self.data_set.get_path_from_date_time(curr_date)
                self._get_netcdf_data(curr_path, u_data, v_data, w_data, temp_data, sal_data)

                if j == 0:#set first point data; well, as much as is possible without netcdf file data...
                    point_a.time = curr_time
                    point_a.lat = particle.lat
                    point_a.lon = particle.lon
                    point_a.depth = particle.depth

                lon = to_radians(point_a.lon)
                lat = to_radians(point_a.lat)
                depth = point_a.depth
                self.projection.c_transform(lon, lat, &x, &y)

                distances, indices = grid.query([x, y], k=self.leaf_size)
                weights = 1.0 / np.array(distances) ** self.power
                weights_sum = 0.0
                
                for k in range(indices.shape[0]):
                    index = indices[k]#I bet the index itself is wrong...
                    u = u_data[index]
                    v = v_data[index]
                    w = w_data[index]
                    temp = temp_data[index]
                    sal = sal_data[index]
                    weight = weights[k]
                    u_w = weight * u
                    v_w = weight * v
                    w_w = weight * w
                    temp_w = weight * temp
                    sal_w = weight * sal
                    weights_sum += weight
                u_w /= weights_sum
                v_w /= weights_sum
                w_w /= weights_sum
                temp_w /= weights_sum
                sal_w /= weights_sum

                x += u_w * seconds_per_hour
                y += v_w * seconds_per_hour
                depth += w_w * seconds_per_hour
                self.projection.c_reverse_transform(x, y, &lon, &lat)
                
                point_a.u = u_w
                point_a.v = v_w
                point_a.w = w_w
                point_a.temp = temp_w
                point_a.sal = sal_w

                point_b.time = particle.start_time + ((j + 1) / 24.0)
                point_b.lat = to_degrees(lat)
                point_b.lon = to_degrees(lon)
                point_b.depth += depth
        
        print(np.array(trajectory))