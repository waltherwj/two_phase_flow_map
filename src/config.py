"""
Configuration file. Contains the class that is 
used for configuration and the configuration variables
"""


class Config:
    """
    class for a namespace of the configs
    """

    NUMBER_ROUGH_DATAPOINTS = 20
    NUMBER_REFINED_DATAPOINTS = 70
    RESOLUTION = 256
    MIN_ULS = 1e-3
    MAX_ULS = 1e1
    MIN_UGS = 1e-2
    MAX_UGS = 1e2
