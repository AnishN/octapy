import datetime as dt
from octapy.model import Model
from octapy.hycom_data_set import HycomDataSet
from octapy.projection import Projection
from octapy.particle import *

sub_model = "GOMl0.04/expt_31.0"
data_dir_path = "./data"
start = dt.datetime(year=2009, month=5, day=1)
end = dt.datetime(year=2009, month=5, day=2)

hycom = HycomDataSet(sub_model, data_dir_path, start, end)
central_lat = hycom.get_central_lat()
central_lon = hycom.get_central_lon()
proj = Projection(
    source_crs_str="+proj=longlat +ellps=WGS84",
    target_crs_str="+proj=merc +ellps=WGS84 +lat_0={0} +lon_0={1} +x_0=0.0, +y_0=0.0 +units=m".format(central_lat, central_lon)
)
particles = get_particles_from_release_file("release.csv")
model = Model(hycom, proj)
model.run(particles)