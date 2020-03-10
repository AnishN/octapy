from libc.stdint cimport *

ctypedef packed struct ParticleC:
    uint64_t id
    float lat
    float lon
    float depth
    double start_time
    double curr_time
    double end_time

ctypedef packed struct PointC:
    double time
    float lat
    float lon
    float depth
    float u
    float v
    float w
    float temp
    float sal