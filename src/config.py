"""
Configuration file. Contains the class that is 
used for configuration and the configuration variables
"""
from matplotlib import cm


class Config:
    """
    class for a namespace of the configs
    """

    NUMBER_DATAPOINTS = 300
    MIN_ULS = 1e-3
    MAX_ULS = 1e1
    MIN_UGS = 1e-2
    MAX_UGS = 1e2

    CATEGORIES = {
        "dispersed bubble": 1,
        "stratified": 2,
        "annular": 3,
        "bubbly": 4,
        "elongated bubble": 5,
        "slug": 6,
        "churn": 7,
    }
    CMAP = cm.get_cmap("tab20", lut=len(CATEGORIES))
