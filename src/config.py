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
        "dispersed bubble": 0,
        "stratified": 1,
        "annular": 2,
        "bubbly": 3,
        "elongated bubble": 4,
        "slug": 5,
        "churn": 6,
    }
    CMAP = cm.get_cmap("Dark2", lut=len(CATEGORIES))
