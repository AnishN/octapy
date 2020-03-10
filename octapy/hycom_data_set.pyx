import datetime as pydt
import glob
import numpy as np
import os
import string
import urllib.request
import xml.etree.cElementTree as et

# TODO: need to cache the data_set.xml file. This would avoid re-querying HYCOM!

cdef class HycomDataSet:

    def __cinit__(self, str sub_model, str data_dir_path, dt.datetime start=None, dt.datetime end=None):
        self.base_date = pydt.datetime(year=1900, month=12, day=31, hour=0, minute=0, second=0)
        self.sub_model = sub_model
        self.data_dir_path = data_dir_path
        self._parse_data_set_xml()
        self.start = start if start != None else self.start_limit
        self.end = end if end != None else self.end_limit
        self._download_netcdf_files()
    
    def __dealloc__(self):
        pass

    def _parse_data_set_xml(self):
        #see REST API url description here:
        #https://www.unidata.ucar.edu/software/tds/current/reference/NetcdfSubsetServiceReference.html
        cdef:
            str base_xml_url
            str xml_url
            bytes xml_data
            str axis_name
            tuple axis_shape
            size_t i
            str v

        base_xml_url = "http://ncss.hycom.org/thredds/ncss/{0}/dataset.xml"
        xml_url = base_xml_url.format(self.sub_model)
        xml_req = urllib.request.urlopen(xml_url)
        xml_data = xml_req.read()
        root = et.fromstring(xml_data)
        axes = root.findall("axis")
        for axis in axes:
            i = 0
            ax = axis.attrib
            axis_name = ax["name"]
            axis_shape = (int(ax["shape"]),)
            if axis_name == "Depth":
                self.depths = np.zeros(shape=axis_shape, dtype=np.float32)
                values = axis.find("values")
                for v in values.text.split():
                    self.depths[i] = float(v)
                    i += 1
            elif axis_name == "Latitude":
                self.lats = np.zeros(shape=axis_shape, dtype=np.float32)
                values = axis.find("values")
                for v in values.text.split():
                    self.lats[i] = float(v)
                    i += 1
            elif axis_name == "Longitude":
                self.lons = np.zeros(shape=axis_shape, dtype=np.float32)
                values = axis.find("values")
                for v in values.text.split():
                    self.lons[i] = float(v)
                    i += 1
            elif axis_name == "MT":
                values = axis.find("values")
                self.start_limit = self.base_date + dt.timedelta(days=float(values.attrib["start"]))
                increment = float(values.attrib["increment"])
                self.end_limit = self.start_limit + dt.timedelta(days=increment * axis_shape[0])

    def _download_netcdf_files(self):
        """
        Checks to see if the data directory exists.
        Then validates that all of the hourly netcdf files in the supplied data range are present.
        """
        cdef:
            str base_prefix
            str base_query
            str base_time
            str curr_path
            str curr_time
            dt.timedelta delta
            dt.datetime curr_date_time
        
        base_prefix = "http://ncss.hycom.org/thredds/ncss/{0}/{1}/hrly"
        base_query = "?var=salinity&var=temperature&var=u&var=v&var=w_velocity&time_start={0}&time_end={0}&accept=netcdf"
        base_time = "%Y-%m-%dT%H:%M:%S"
        if not os.path.exists(self.data_dir_path):
            os.makedirs(self.data_dir_path)
        delta = pydt.timedelta(hours=1)
        curr_date_time = self.start
        while curr_date_time < self.end:
            curr_path = self.get_path_from_date_time(curr_date_time)
            if not os.path.exists(curr_path):
                prefix = base_prefix.format(self.sub_model, str(curr_date_time.year))
                curr_time = curr_date_time.strftime(base_time)
                query = base_query.format(curr_time)
                url = prefix + query
                urllib.request.urlretrieve(url, curr_path)
            curr_date_time += delta
    
    def get_path_from_date_time(self, dt.datetime date_time):
        cdef:
            str base_path
            str base_time
            str path

        sub_model = self.sub_model
        for c in self.sub_model:
            if c in string.punctuation:
                sub_model = sub_model.replace(c, "")
        base_path = "{0}/{1}_{2}.nc"
        base_time = "%Y%m%d_%H%M%S"
        path = base_path.format(self.data_dir_path, sub_model, date_time.strftime(base_time))
        return path

    def get_min_lat(self):
        return self.lats[0]

    def get_max_lat(self):
        return self.lats[self.lats.shape[0] - 1]

    def get_central_lat(self):
        cdef:
            double central_lat
            size_t num_lats
            size_t middle

        num_lats = self.lats.shape[0]
        middle = num_lats / 2
        if num_lats % 2:
            central_lat = (self.lats[middle] + self.lats[middle + 1])/2.0
        else:
            central_lat = self.lats[middle + 1]
        return central_lat
    
    def get_min_lon(self):
        return self.lons[0]

    def get_max_lon(self):
        return self.lons[self.lons.shape[0] - 1]
    
    def get_central_lon(self):
        cdef:
            double central_lon
            size_t num_lons
            size_t middle

        num_lons = self.lons.shape[0]
        middle = num_lons / 2
        if num_lons % 2:
            central_lon = (self.lons[middle] + self.lons[middle + 1])/2.0
        else:
            central_lon = self.lons[middle + 1]
        return central_lon