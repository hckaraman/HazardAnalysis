#############
import logging
import sys

logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG)  ## CHANGE THIS FOR LOG LEVEL
# logger.debug("Loading module...")
######################


from data.constants import *

# from data.constants import *
from hazard.earthquake import Seismic_hazard, find_pgas, Plot_hazard_curve
from hazard.flood import hazard_flood
from hazard.wind import hazard_wind
from hazard.avalanche import *
from hazard.landslide import hazard_landslide
from hazard.rockfall import hazard_rockfall
from hazard.wildfire import hazard_wildfire


def evaluation_hazard(lon, lat, name='Atlar', sector='Dolasiyor'):
    # logger.debug(f'evaluation_hazard({asset["name"]})')

    hazard_functions = []
    hazard_parameters = []
    result = pd.DataFrame()

    for i, h in enumerate(HAZARD):
        if h == "EARTHQUAKE":
            logger.debug("handling '%s'" % h)
            tmp = Seismic_hazard(
                DATA / "EARTHQUAKE",
                [1, 1, 1],
                [
                    round(k, 4) for k in find_pgas([
                    float(lon),
                    float(lat)
                ])
                ],
                TURKEY_TYPICAL_MAF_TN,
                ANALYSIS_LIST[h],
                sector,
                name,
                draw=False
            )

            hazard_parameters.append(tmp[1])
            hazard_functions.append(tmp[0])

            result['IM'] = ANALYSIS_LIST[h]
            result['MAF'] = tmp[0]
            result['Hazard'] = h
            result['longitude'] = lon
            result['latitude'] = lat
            print("Earh")


        elif h == "FLOOD":
            logger.debug("handling '%s'" % h)
            hazard_parameters.append(
                hazard_flood(lon, lat)
            )

            no_flood_hazard = hazard_parameters[i][0] == "NO_HAZARD" \
                              or hazard_parameters[i][0] == "NaN" \
                              or hazard_parameters[i][0] == "NO HAZARD"

            if no_flood_hazard:
                hazard_functions.append("NO_HAZARD")

                data = pd.DataFrame()
                data['IM'] = ANALYSIS_LIST[h]
                data['MAF'] = 0
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)


            else:
                a = hazard_parameters[i][0]
                k = hazard_parameters[i][1]
                tmp = Plot_hazard_curve(
                    a,
                    k,
                    ANALYSIS_LIST[h],
                    max(ANALYSIS_LIST[h]),
                    h,
                    name,
                    'Exponential',
                    draw=False
                )

                data = {
                    'IM': ANALYSIS_LIST[h],
                    'MAF': tmp,
                }
                data = pd.DataFrame.from_dict(data)
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
            print("Flood")


        elif h == "WIND":
            logger.debug("handling '%s'" % h)
            hazard_parameters.append(
                hazard_wind(
                    lon, lat
                )
            )

            a_wind = hazard_parameters[i][0]
            k_wind = hazard_parameters[i][1]

            no_rockfall_hazard = hazard_parameters[i][0] == "NO_HAZARD" or \
                                 hazard_parameters[i][0] == "NaN" \
                                 or hazard_parameters[i][0] == "NO HAZARD"

            if no_rockfall_hazard:
                hazard_functions.append("NO_HAZARD")
                data = pd.DataFrame()
                data['IM'] = ANALYSIS_LIST[h]
                data['MAF'] = 0
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
            else:
                tmp = Plot_hazard_curve(
                    a_wind,
                    k_wind,
                    ANALYSIS_LIST[h],
                    max(ANALYSIS_LIST[h]),
                    h,
                    name,
                    'Exponential',
                    draw=False
                )
                data = {
                    'IM': ANALYSIS_LIST[h],
                    'MAF': tmp,
                }
                data = pd.DataFrame.from_dict(data)
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
            print("Wind")
        elif h == "AVALANCHE":
            logger.debug("handling '%s'" % h)
            hazard_parameters.append(
                Hazard_avalanche_MAFPowerLaw(lon, lat)
            )

            no_avalanche_hazard = hazard_parameters[i][2] == "NO HAZARD"
            if no_avalanche_hazard:  # This is the term pressure in the def()
                hazard_functions.append("NO_HAZARD")
                data = pd.DataFrame()
                data['IM'] = ANALYSIS_LIST[h]
                data['MAF'] = 0
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
            else:
                a = hazard_parameters[i][0]
                k = hazard_parameters[i][1]
                tmp = Plot_hazard_curve(
                    a,
                    k,
                    ANALYSIS_LIST[h],
                    max(ANALYSIS_LIST[h]),
                    h,
                    name,
                    'Power',
                    draw=False
                )

                data = {
                    'IM': ANALYSIS_LIST[h],
                    'MAF': tmp,
                }
                data = pd.DataFrame.from_dict(data)
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)

        elif h == "LANDSLIDE":
            logger.debug("handling '%s'" % h)
            hazard_parameters.append(
                hazard_landslide(lon, lat)
            )

            if hazard_parameters[i][0] == [None, None]:
                hazard_functions.append("NO_HAZARD")
                data = pd.DataFrame()
                data['IM'] = ANALYSIS_LIST[h]
                data['MAF'] = 0
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
            else:
                a = hazard_parameters[i][0][0]
                k = hazard_parameters[i][0][1]
                tmp = Plot_hazard_curve(
                    a,
                    k,
                    ANALYSIS_LIST[h],
                    max(ANALYSIS_LIST[h]),
                    h,
                    name,
                    'Power',  # Exponential
                    draw=False
                )

                data = {
                    'IM': ANALYSIS_LIST[h],
                    'MAF': tmp,
                }
                data = pd.DataFrame.from_dict(data)
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)

        elif h == "ROCKFALL":
            logger.debug("handling '%s'" % h)
            hazard_parameters.append(
                hazard_rockfall(
                    lon, lat
                )
            )

            if hazard_parameters[i] == ['NO HAZARD', 'NO HAZARD']:
                hazard_functions.append("NO_HAZARD")
                data = pd.DataFrame()
                data['IM'] = ANALYSIS_LIST[h]
                data['MAF'] = 0
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
            else:
                a = hazard_parameters[i][0]
                k = hazard_parameters[i][1]
                tmp = Plot_hazard_curve(
                    a,
                    k,
                    ANALYSIS_LIST[h],
                    max(ANALYSIS_LIST[h]),
                    h,
                    name,
                    'Power',
                    draw=False
                )

                data = {
                    'IM': ANALYSIS_LIST[h],
                    'MAF': tmp,
                }
                data = pd.DataFrame.from_dict(data)
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)

        elif h == "WILDFIRE":
            logger.debug("handling '%s'" % h)
            hazard_parameters.append(
                hazard_wildfire(
                    a_wind, k_wind,
                    lon, lat
                )
            )

            if hazard_parameters[i][0] == 'NO HAZARD':
                hazard_functions.append("NO_HAZARD")
                data = pd.DataFrame()
                data['IM'] = ANALYSIS_LIST[h]
                data['MAF'] = 0
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
            else:
                a = hazard_parameters[i][0]
                k = hazard_parameters[i][1]
                tmp = Plot_hazard_curve(
                    a,
                    k,
                    ANALYSIS_LIST[h],
                    max(ANALYSIS_LIST[h]),
                    h,
                    name,
                    'Exponential',
                    draw=False
                )

                data = {
                    'IM': ANALYSIS_LIST[h],
                    'MAF': tmp,
                }
                data = pd.DataFrame.from_dict(data)
                data['Hazard'] = h
                data['longitude'] = lon
                data['latitude'] = lat
                result = result.append(data)
        else:
            logger.WARN("unhandled hazard '%s'" % h)

    result.to_csv('Result.csv')
    # return hazard_functions, hazard_parameters
    return result


# logger.debug("Module loaded.")

"""
27.37124,39.11272
"""
if __name__ == '__main__':
    # asset = pd.DataFrame()
    # asset = {
    #     'name': 'at',
    #     'longitude': 27.37124,
    #     'latitude': 39.11272,
    #     'sector': 'Energy'
    # }
    # # lon,lat = 29,39
    lon = 27.8593
    # lon = float(sys.argv[1])
    lat = 37.7931
    # lat = float(sys.argv[2])
    #
    res = evaluation_hazard(lon, lat)
    print("a")
    # print("This is the name of the program:", sys.argv[1])
    # print("This is the name of the program:", sys.argv[2])
