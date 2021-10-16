
#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
######################

import rasterio
raster_cache = dict()

def _load_raster(path,lat_image_reverse=False):

    if not path in raster_cache:
      logger.debug("loading '%s' from disk" % path)
      with rasterio.Env():
        with rasterio.open(path) as src:  
            image = src.read(1)  ## first and only band  
      lat = [src.xy(i,0)[1] for i in range(src.height)]
      lon = [src.xy(0,i)[0] for i in range(src.width)]
      if lat_image_reverse:
        lat.reverse() # ordered list from max to min    

      raster_cache[path] = (image, lat, lon)
      
    logger.debug("'%s' cached" % path)
    return raster_cache[path]