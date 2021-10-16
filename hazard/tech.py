# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 10:22:27 2020

@author: RVA04
"""
#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
#logger.debug("Loading module...")
######################

from data.constants import *
import pandas as pd
import os
from fragility import fragility_func_from_parameters

def hazard_tech(asset):
    logger.debug("hazard_tech called")
    ratio_fb = 0.066687
    
    if (asset['category'].split('_')[0] == 'st') or (asset['category'].split('_')[0] == 're'):
        p_ing = 130*10**-6
    elif asset['category'].split('_')[0] == 'ho':
        p_ing = 6.74 * 10 ** -6
    else:
        p_ing = 8 * 10 ** -6

    if asset['category'].split('_')[0] == 'ho':
        area = 32
    elif asset['category'].split('_')[0] == 'st' \
            or asset['category'].split('_')[0] == 're':
        area = 1000
    elif (asset['category'] == 'ss') or (asset['subcategory'] == 'thermal'):
        area = 100
    else:
        area = 0

    h_fire = area * p_ing
    h_blast = h_fire * ratio_fb
        
    return h_fire, h_blast

def fragility_tech(asset):
    logger.debug("fragility_tech called")
    if (asset['category'].split('_')[0] == 'st') or (asset['category'].split('_')[0] == 're') :
        fr = 0.99474267
    elif asset['category'].split('_')[0] == 'ho':
        fr = 0.01860329
    else:
        fr = 0.11110354
    return fr

def risk_tech(asset,Loss_OUTPUTS):  
    logger.debug("fragility_tech called")
    h_fire, h_blast = hazard_tech(asset)
    fr = fragility_tech(asset)
    
    if asset['category'].split('_')[0] == 'ho':
        coef = 0.8
    else: coef =1 
    
    TECHfire_EAL_Deaths = Loss_OUTPUTS[0][3]*h_fire*fr/(Loss_OUTPUTS[1][3]/coef)
    TECHfire_EAL_Direct = Loss_OUTPUTS[1][3]*h_fire*fr/(Loss_OUTPUTS[1][3]/coef)
    TECHfire_EAL_Indirect = Loss_OUTPUTS[2][3]*h_fire*fr/(Loss_OUTPUTS[1][3]/coef)

    TECHblast_EAL_Deaths = Loss_OUTPUTS[0][3]*h_blast*fr/(Loss_OUTPUTS[1][3]/coef)
    TECHblast_EAL_Direct = Loss_OUTPUTS[1][3]*h_blast*fr/(Loss_OUTPUTS[1][3]/coef)
    TECHblast_EAL_Indirect = Loss_OUTPUTS[2][3]*h_blast*fr/(Loss_OUTPUTS[1][3]/coef)
    
    return [TECHfire_EAL_Deaths,TECHfire_EAL_Direct,TECHfire_EAL_Indirect,TECHblast_EAL_Deaths,TECHblast_EAL_Direct,TECHblast_EAL_Indirect]

#logger.debug("Module loaded.")