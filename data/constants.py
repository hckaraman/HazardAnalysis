import numpy as np
from pathlib import Path
from socket import gethostname, gethostbyname

HAZARD = ["EARTHQUAKE", "FLOOD", "WIND", "AVALANCHE", "LANDSLIDE", "ROCKFALL", "WILDFIRE"]

TECH_HAZARD = ["TECHFIRE", "TECHBLAST"]

TURKEY_TYPICAL_MAF_TN = [0.0228, 0.01386, 0.002107, 0.000404]

TR = ['68%', '50%', '10%', '2%']

DS = ['DS1', 'DS2', 'DS3', 'DS4']

LIMIT_STATES = 4

INF = float('inf')

ROOT = Path(__file__).resolve().parent.parent

DATA = ROOT / "data"

DROPPABLE = [
        "physical1",    
        "physical2",    
        "physical3",
        "capacity1",
        "capacity2",
        "capacity3",
        "capacity4",
        "value[tot]",
        "value[nstr]",
        "emergency_plan",
        "training",
        "insurance",
        "backup",
        "last_event",
        "advance_notice",
        "ews",
        "offset_loss",
        "coo_internal",
        "agreements",
        "coo_gov",
        "coo_prov",
        "coo_ext_assets",
        "interdependencies"
    ]

ASSET_HEADERS = [
    "id",
    "name",
    "sector",
    "subsector",
    "capacity",
    "category",
    "subcategory",
    "region",
    "city",
    "manager",
    "manager_type",
    "damage",
    "degradation",
    "foundations",
    "fragility",
    "latitude",
    "longitude",
    "length",
    "depth",
    "height",
    "area",
    "geometry",
    "operation_time",
    "osm_type",
    "osmid"
    "owner",
    "owner_type"
    "people",
    "physical",
    "retrofit",
    "statics",
    "material",
    "meta",
    "value",
    "value",
    "seismic_provision",
    "yoc_range"
]

CATEGORY = {
    "ho": "HOSPITAL",
    "ep": "ENERGY_PRODUCTION",
    "st": "STORAGE_PLANT",
    "pl": "PIPELINE",
    "ss": "SUBSTATION",
    "re": "REFINERY",
    "td": "TRANSMISSION_AND_DISTRIBUTION_LINE"
}

SUBCATEGORY = {
    "oil": "THERMAL_POWER_PLANT",
    "gas": "THERMAL_POWER_PLANT",
    "lignite": "THERMAL_POWER_PLANT",
    "coal": "THERMAL_POWER_PLANT",
    "geothermal": "THERMAL_POWER_PLANT",
    "thermal": "THERMAL_POWER_PLANT",
    "hydroelectric": "HYDROELECTRIC_POWER_PLANT",
    "solar": "SOLAR_POWER_PLANT",
    "wind": "WIND_POWER_PLANT",
    "nuclear": "NUCLEAR_POWER_PLANT"
}

ANALYSIS_LIST = {
    "EARTHQUAKE": np.arange(0.01, 2.01, 0.01).tolist(),
    "FLOOD": np.arange(0.1, 5.1, 0.1).tolist(),
    "WIND": np.arange(1, 201, 1).tolist(),
    "AVALANCHE": np.arange(0.5, 100, 0.5).tolist(),
    "LANDSLIDE": np.arange(0.1, 5.1, 0.1).tolist(),
    "ROCKFALL": np.arange(0.1, 5.1, 0.1).tolist(),
    "WILDFIRE": np.arange(3, 6, 0.05).tolist()
}

MATERIALS = {
    "rc": "Reinforced_Concrete",
}

HOSTNAME = gethostname()

ADDRESS = gethostbyname(HOSTNAME)

LOSS_TYPES = ["Deaths", "Direct", "Indirect"]

DEFAULT_ASSET_VALUE = 5_000_000

RESILIENCE_INDICATORS = ["Res1", "Res2", "Res3", "Res_tot"]

# MAPS USED
AltitudeMap = "SRTM_TR_30_shaped.tif"
AvgDischargeMap = "logQavg.tif"
VarDischargeMap = "logQvar.tif"
NationalBorder = "gadm36_TUR_0.shx"
ProvinceBorder = "gadm36_TUR_1.shx"
LandCoverMap = "Landcover_Lat_Lon.tif"
SoilTypeMap = "HWSD_TR_shaped.tif"
WindMap = "TUR_wind-speed_10m.tif"


