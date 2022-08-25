"""
###############################################################################
###############################################################################
###############################################################################
###
###  This is a collection of 14th Weather Squadron standardized meteorological
###  functions and SAR created methods for use with Python version 3.
###
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Standard function comment block:

###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  {Function name}
###
### INPUTS:  {List function input variables and units}
###
### OUTPUTS:  {List function output variables and units}
###
### PURPOSE:
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:
###
###############################################################################
###
### END OF PYTHON FUNCTION:  {Function name}
###
###############################################################################
###############################################################################
"""

# Import of modules
import datetime as dt
import math
from math import exp
import pandas as pd
import numpy as np
# import ephem
# import pytz
# from timezonefinder import TimezoneFinder


"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  PRcalc
###
### LEGACY SAS PROGRAM NAME:  prcalc (Pressure Calculation)
###
### INPUTS:
###         ELEV = 'name of elevation_in_meters_column'
###         TEMPC = 'name of tempc_column'
###         ALSTG = 'name of altimeter_column'
###         SLP = 'name of SLP_column'
###         STAPRS = 'name of station_pressure_column'
###
### OUTPUTS:
###         DataFrame with filled Pressure Data
###
###
### BCR TRACKING NUMBER:
###
### PURPOSE: Compute missing pressures, if possible from other pressure
###
### METHOD:  Make sure inputs are within bounds.  If elevation is known
###          and at least one pressure exists, then the other pressures
###          can be derived.
###
### REFERENCES:
###   none
###
### PROCESS NARRATIVE LOCATION:
###
### LEGACY SAS MACROS/SUBROUTINES CALLED:
###   none
###
### FILES/DATABASE/TABLES ACCESSED:
###   none
###
### LEGACY SAS DATA DICTIONARY:
###   alstg       Input/Output: altimeter setting in mb
###   elev        Input: elevation in meters
###   slp         Input/Output: sea level pressure in mb
###   stnpr       Input/Output: station pressure in mb
###   tempc       Input: temperature in Celsius
###   alsin       altimeter setting in inches
###   alstg       altimeter setting in mb
###   elev        elevation in meters
###   slpin       sea level pressure in inches
###   slprs       sea level pressure in mb
###   stnin       station pressure in inches
###   stnpr       station pressure in mb
###   tempk       temperature in Kelvin
###
### REMARKS:
###   derived from the FORTRAN routine PRCALC
###
### UPDATES:
###   08MAR1999  Mr. Kiess/DOC2.  New
###   08AUG2018  SSgt Levi Jones port to Python
###
###############################################################################
"""

def PRcalc(elev,tempc,alstg,slp,staprs):
    if pd.isnull(alstg) == False :
        alsin = alstg / 33.864
    else:
        alsin = np.nan
        
    if pd.isnull(slp) == False :
        slpin = slp * .02953
    else:
        slpin = np.nan
        
    if pd.isnull(staprs) == False :
        stnin = staprs * .02953
    else:
        stnin = np.nan
    #  CHECK FOR MISSING ELEVATION.  IF ELEVATION IS MISSING, THEN
    #  IT IS NOT POSSIBLE TO CALCULATE OTHER VALUES, SO THERE IS NO
    #  NEED TO PROCEED.
    ###### important note pd.isnull handles np.nan & None // np.isnan only handles np.nan 
    if pd.isnull(elev) == False :
        if pd.isnull(tempc) == False :
            #IF INPUT PRESSURE DATA IS OUT OF REALISTIC BOUNDS THEN SET DATA TO Zero.
            if (alsin > 32) or (alsin < 25):
                alsin = np.nan
            if (slp > 1090) or (slp < 860):
                slp = np.nan
            if (staprs > 1090) or (staprs < (900 - elev/10)):
                staprs = np.nan
            tempk = tempc + 273.16

            #If all three are Zero none can be calculated
            if np.isnan(alsin) == True :
                if np.isnan(slpin) == True :
                    if np.isnan(stnin) == True :
                        return alsin,slpin,stnin
                    # Ends Function if all is nan

            # if Station pressure is missing
            if np.isnan(stnin) == True :
                #if altimeter is present
                if np.isnan(alsin) == False :
                    #Calc Station pressure using Altimeter
                    stnin = ((alsin**.19025) - (elev * (29.921**.19025) * (.0065/288.16)))**5.2561

                    # if SLP is missing
                    if np.isnan(slpin) == True :
                        #Calc SLP using Station pressure
                        slpin = stnin / exp((-9.8 * elev) / (287 * tempk))

                #if altimeter is missing but SLP is present
                elif np.isnan(slpin) == False:
                #Calc Station pressure using SLP
                    stnin = slpin * exp((-9.8*elev) / (287 * tempk))
                    # Calc Altimeter using Station Pressure
                    alsin = ((stnin**.19025) + (elev*(29.921**.19025) * (.0065/288.16)))**5.2561

            #if SLP is missing & Station pressure resolved to not null
            elif np.isnan(slpin) == True :
                
                #Calc SLP using Station pressure
                slpin = stnin / exp((-9.8 * elev) / (287 * tempk))
                #if altimeter is missing
                if np.isnan(alsin) == True :
                    alsin = ((stnin**.19025) + (elev*(29.921**.19025) * (.0065/288.16)))**5.2561

            # altimeter has to be missing
            else:
                # Calc Altimeter using Station Pressure
                alsin = ((stnin**.19025) + (elev*(29.921**.19025) * (.0065/288.16)))**5.2561

        elif np.isnan(alsin ) == True :
            #Calc alsin from stnin
            if np.isnan(stnin) == False :
                alsin = ((stnin**.19025) + (elev*(29.921**.19025) * (.0065/288.16)))**5.2561
        else :
            return alstg, slp, staprs
        #Convert back to respective unit and return non-zero value
        if np.isnan(alsin) == False :
            alstgout = round(alsin*33.864,1)
        else : alstgout = np.nan

        if np.isnan(stnin) == False :
            staprsout = round(stnin / .02953,1)
        else : staprsout = np.nan

        if np.isnan(slpin) == False :
            slpout = round(slpin / .02953,1)
        else : slpout = np.nan

        ## Only path to get full output##
        return alstgout, slpout, staprsout   
    else : return alstg, slp, staprs
    
"""
### END OF PYTHON FUNCTION:  PRcalc
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  PRSalt
###
### LEGACY SAS PROGRAM NAME:  prsalt (Pressure Altitude)
###
### BCR TRACKING NUMBER:
###
### PURPOSE: Compute the pressure altitude in meters
###
### METHOD:  Compute station pressure from altimeter setting if necessary
###          Compute pressure altitude in meters
###
### REFERENCES:
###   none
###
### FILES/DATABASE/TABLES ACCESSED:
###   none
###
### LEGACY SAS DATA DICTIONARY:
###   &ainvar      Input: altimeter setting in mb
###   &einvar      Input: elevation in meters
###   &pinvar      Input: station pressure in mb
###   &poutvar     Output: pressure altitude in meters
###   _palt        pressure altitude in meters
###   _stnprs      station pressure in mb
###
### REMARKS:
###   derived from FORTRAN routine PRSALT
###
### UPDATES:
###   02MAR1999  Mr. Kiess/DOC2.  New
###   08AUG2018  SSgt Levi Jones  Port to Python
###
###############################################################################
"""


def PRSalt(staprs):
    if pd.isnull(staprs) == False :
        return round(44330.7216 * (1-(staprs/1013.25)**.19029),2)
    else :
        return np.nan

"""
### END OF PYTHON FUNCTION:  PRSalt
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  DENalt
###
### LEGACY SAS PROGRAM NAME: denalt  (Density Altitude / Relitive Humidity)
###
### BCR TRACKING NUMBER:
###
### PURPOSE: Computes the density altitude in meters & Computes relative
###           humidity from temperature and dew point
###
### METHOD:  Validate input and compute station pressure, if necessary
###          Then compute the virtual temperature and finally density alt.
###         Compute vapor pressure and saturation vapor pressure.
###         Compute relative humidity as the ratio of the mixing ratios.
###
### REFERENCES:
###   none
###
### PROCESS NARRATIVE LOCATION:
###
### FILES/DATABASE/TABLES ACCESSED:
###   none
###
### LEGACY SAS DATA DICTIONARY:
### DENALT---------------------
###   &ainvar      Input: altimeter setting in mb
###   &dinvar      Input: dew point temperature in Celsius
###   &doutvar     Output: density altitude in meters
###   &einvar      Input: elevation in meters
###   &pinvar      Input: station pressure in mb
###   &tinvar      Input: temperature in Celsius
###   _dalt        density altitude in meters
###   _dewp        dew point temperature in Celsius
###   _stnprs      station pressure in mb
###   _stnprsi     station pressure in inches
###   _temp        air temperature in Celsius
###   _tvirt       virtual temperature in Celsius
###   _tvirtr      virtual temperature in Rankine
### RELHUM--------------------
###   &dewp        Input: dew point temperature in user specified units
###   &humrel      Output: relative humidity in percent
###   &stnpr       Input: station pressure in mb
###   &temp        Input: air temperature in user specified units
###   &tunit       Input: specifies input temperature units (C, F, K, R)
###   _dewp        dew point temperature in Kelvin
###   _e           vapor pressure in mb
###   _es          saturation vapor pressure in mb
###   _relhum      relative humidity in percent
###   _stnpr       station pressure in mb#
###   _temp        air temperature in Kelvin
###
### REMARKS:
###   derived from the FORTRAN routine DENALT
###
### UPDATES:
###   05MAR1999  Mr. Kiess/DOC2.  New
###   08AUG2018  SSgt Levi Jones  Port to Python cmbine DENLAT and RELHUM
###                               due to both using vapor macro.
###############################################################################
"""

def DENalt(staprs,tempc,dewpc):
    if pd.isnull(staprs) == False and pd.isnull(tempc) == False and pd.isnull(dewpc) == False :
        staprsin = staprs * .02953
        tempk = tempc + 273.15
        dewpk = dewpc + 273.15
        # vapor pressure in mb GOFF-GRATCH METHOD
        if dewpk > 100 :
            if tempk >= 273.16 :
                vapor_e = 10**(23.832241-(5.02808*math.log10(dewpk))-
                (1.3816*(10**(-7))*(10**(11.334-(0.0303998*dewpk))))+
                (8.1328*(10**(-3))*(10**(3.49149-(1302.8844/dewpk))))-(2949.076/dewpk))
            else :
                vapor_e = 10**((3.56654*math.log10(dewpk))-
                (0.0032098*dewpk)-(2484.956/dewpk)+2.0702294)
        else :
            vapor_e = None
            
        # saturation vapor pressure in mb GOFF-GRATCH METHOD
        if tempk > 100 :
            if tempk >= 273.16 :
                vapor_es = 10**(23.832241-(5.02808*math.log10(tempk))-
                (1.3816*(10**(-7))*(10**(11.334-(0.0303998*tempk))))+
                (8.1328*(10**(-3))*(10**(3.49149-(1302.8844/tempk))))-(2949.076/tempk))
            else :
                vapor_es = 10**((3.56654*math.log10(tempk))-
                (0.0032098*tempk)-(2484.956/tempk)+2.0702294)
        else :
            vapor_es = None
            
        # dalt
        if pd.isnull(vapor_e) == False:
            tvirt = (tempk / (1 - (0.37803*vapor_e/staprs))) - 273.16
            tvirtr = (tvirt * 9.0/5.0 + 32) + 459.69
            dalt = 44330.7216 * (1 - ((17.326 * staprsin /tvirtr)**0.235))
        else :
            dalt = np.nan

        return dalt
    
    else : return np.nan 

"""
### END OF PYTHON FUNCTION:  DENalt
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  RELHUMfunc
###
### LEGACY SAS PROGRAM NAME:  relhum  (Relative Humidity)
###
### BCR TRACKING NUMBER:
###
### PURPOSE:  Computes relative humidity from temperature and dew point
###
### METHOD:  Validate input and convert temperatures to Kelvin.
###          Compute vapor pressure and saturation vapor pressure.
###          Compute relative humidity as the ratio of the mixing ratios.
###
### REFERENCES:
###   none
###
### PROCESS NARRATIVE LOCATION:
###
### LEGACY SAS MACROS/SUBROUTINES CALLED:
###   c2k          converts Celsius to Kelvin
###   f2k          converts Fahrenheit to Kelvin
###   r2k          converts Rankine to Kelvin
###   vapor        calculates vapor pressure in mb
###
### FILES/DATABASE/TABLES ACCESSED:
###   none
###
### LEGACY SAS DATA DICTIONARY:
###   &dewp        Input: dew point temperature in user specified units
###   &humrel      Output: relative humidity in percent
###   &stnpr       Input: station pressure in mb
###   &temp        Input: air temperature in user specified units
###   &tunit       Input: specifies input temperature units (C, F, K, R)
###   _dewp        dew point temperature in Kelvin
###   _e           vapor pressure in mb
###   _es          saturation vapor pressure in mb
###   _relhum      relative humidity in percent
###   _stnpr       station pressure in mb
###   _temp        air temperature in Kelvin
###
### REMARKS:
###   derived from FORTRAN function RHDP2
###
### UPDATES:
###   04MAR1999  Mr. Kiess/DOC2.  New
###   07JAN2019  SSgt Jones Python import
###
###############################################################################
"""

def RELHUMfunc(staprs,tempc,dewpc):
    if pd.isnull(staprs) == False and pd.isnull(tempc) == False and pd.isnull(dewpc) == False :
        staprsin = staprs * .02953
        tempk = tempc + 273.15
        dewpk = dewpc + 273.15
        
        # vapor pressure in mb GOFF-GRATCH METHOD
        if dewpk > 100 :
            if tempk >= 273.16 :
                vapor_e = 10**(23.832241-(5.02808*math.log10(dewpk))-
                (1.3816*(10**(-7))*(10**(11.334-(0.0303998*dewpk))))+
                (8.1328*(10**(-3))*(10**(3.49149-(1302.8844/dewpk))))-(2949.076/dewpk))
            else :
                vapor_e = 10**((3.56654*math.log10(dewpk))-
                (0.0032098*dewpk)-(2484.956/dewpk)+2.0702294)
        else :
            vapor_e = np.nan
            
        # saturation vapor pressure in mb GOFF-GRATCH METHOD
        if tempk > 100 :
            if tempk >= 273.16 :
                vapor_es = 10**(23.832241-(5.02808*math.log10(tempk))-
                (1.3816*(10**(-7))*(10**(11.334-(0.0303998*tempk))))+
                (8.1328*(10**(-3))*(10**(3.49149-(1302.8844/tempk))))-(2949.076/tempk))
            else :
                vapor_es = 10**((3.56654*math.log10(tempk))-
                (0.0032098*tempk)-(2484.956/tempk)+2.0702294)
        else :
            vapor_es = np.nan
            
        # relhum
        if np.isnan(vapor_e) == False and np.isnan(vapor_es) == False :
            relhum = (100 * ( vapor_e * (staprs - vapor_es)) /
            ((staprs - vapor_e) * vapor_es))
        else :
            relhum = np.nan
        
        return relhum

    else : return np.nan

"""
### END OF PYTHON FUNCTION:  RELHUMfunc
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  HEATindex
###
### LEGACY SAS PROGRAM NAME:  heatindex (RH + Temp to Apparent Heat Index)
###
### BCR TRACKING NUMBER: 20038
###
### PURPOSE: To compute an Apparent Heat Index using RH and a given Temperature.
### Macro can use temperatures in Fahrenheit, Celsius or Kelvin.
###
### METHOD: The equation for heat index is a regression to human-perceived
###         apparent temperatures.  It uses the dry bulb temp and the relative
###         humidity to compute them.
###
### REFERENCES:
###   This equation is referenced from Stull, Meteorology for Scientists and
### Engineers available in library QC861.2 .S79 2000, which references it from
### www.zunis.org/16element_heat_index_equation.htm
###
###  LEGACY SAS MACROS/SUBROUTINES CALLED:
###  c2f
###  k2f
###
### FILES/DATABASE/TABLES ACCESSED:
###   none
###
### LEGACY SAS DATA DICTIONARY:
###  t = temperature
###  tunit = temperature units (C,K or F)
###  rh = relative humidity
###  hi = humidity index
###
### REMARKS:
###   There is another heat index equation out there.  Perhaps, some day it
### can be incorporated into this macro as another option.  One can find the code
### in $BIN_DIR/wxito_eucom.sas
###
### UPDATES:
### 10Dec2002  Mr. Kiess/DOMD.  New
### 16Dec2002  Mr. Giese/DOMD.  This macro was tested using 46 different
### combinations of RH and Temps in C, K and F.  Results exactly matched the
### values in Table 3-2, page 60 of Meteorology for Scientists and Engineers,
### Stull (2000).
### 24Oct2003  TSgt Wretlind/DOMD. Added routine to check that RH is within
### normal values (0-100%).  Anything outside of those values returns missing.
### 07Jan2019  SSgt Jones import to Python
###
###############################################################################
"""
def HEATindex(tempf,rh):
    if tempf >= 70 :
        hi = (
            16.923+(1.85212e-1*tempf)+(5.37941*rh)-(1.00254e-1*tempf*rh)
            +(9.41695e-3*tempf**2)+(7.28898e-3*rh**2)+(3.45372e-4*tempf**2*rh)
            -(8.14971e-4*tempf*rh**2)+(1.02102e-5*tempf**2*rh**2)
            -(3.8646e-5*tempf**3)+(2.91583e-5*rh**3)+(1.42721e-6*tempf**3*rh)
            +(1.97483e-7*tempf*rh**3)-(2.18429e-8*tempf**3*rh**2)
            +(8.43296e-10*tempf**2*rh**3)-(4.81975e-11*tempf**3*rh**3)
              )
    else:
        hi = tempf
    if (rh > 100) | (rh < 0):
        hi = np.nan

    return hi

"""
### END OF PYTHON FUNCTION:  HEATindex
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  MPS_to_KTS
###
### INPUTS:  WSPDMPS - wind speed meters per second
###
### OUTPUTS:  WSPDKT - wind speed in knots
###
### PURPOSE: convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""


def MPS_to_KTS(wind):
    if pd.isnull(wind) == False:
        wind = wind * 1.943
    else:
        wind = np.nan
    return wind

"""
### END OF PYTHON FUNCTION:  MPS_to_KTS
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  MM_to_IN
###
### INPUTS:  PRCP - precip in mm
###
### OUTPUTS:  PRCP - precip in inchs
###
### PURPOSE: convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""


def MM_to_IN(prcp):
    if pd.isnull(prcp) == False:
        prcp = prcp/25.4
    else:
        prcp = np.nan
    return prcp


"""
### END OF PYTHON FUNCTION:  MM_to_IN
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  MPS_to_MPH
###
### INPUTS:  WSPDMPS - wind speed meters per second
###
### OUTPUTS:  WSPDMPH - wind speed miles per hour
###
### PURPOSE: convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""

def MPS_to_MPH(wind):
    if pd.isnull(wind) == False:
        wind = wind * 2.237
    else:
        wind = np.nan
    return wind

"""
### END OF PYTHON FUNCTION:  MPS_to_MPH
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  WSPDmax
###
### INPUTS:  WSPDKTS - wind speed in knots
###          GUSTKTS - gusts in knots
###
### OUTPUTS:  WMAX - max of the two values
###
### PURPOSE: determine the worst case value for the observation
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""

def WSPDmax(wspd,gust):
    wmax = max(wspd,gust)
    return wmax

"""
### END OF PYTHON FUNCTION:  WSPDmax
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  DMS_to_DD & DD_to_DMS
###
### INPUTS:  DMS Degrees minutes seconds
###          DD decimal degree
###
### OUTPUTS:  opposte of input
###
### PURPOSE: convert between the 2
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""

## DD = degree, MM = minute, SS = second, d = direction (NESW)

def DMS_to_DD(DMS):
    if pd.isnull(DMS):
        return np.nan
    d = str(DMS)[-1]
    if d in ('W','E'):
        DD = str(DMS)[0:3]
        MM = str(DMS)[3:5]
        SS = str(DMS)[5:7]
    else:
        DD = str(DMS)[0:2]
        MM = str(DMS)[2:4]
        SS = str(DMS)[4:6]
    
    
    if d in ('S','W'):
        neg = -1
    else:
        neg = 1

     

    LL = (int(DD) + (int(MM)/60) + (int(SS)/3600)) * neg

    return LL


###############################################################################
###############################################################################
###############################################################################

#d = decimal degree 

def DD_to_DMS(d):
    if pd.isnull(d):
        return np.nan
    dms = [0,0,0]
    dms[0] = int(d)
    dms[1] = int((d-dms[0])*60)
    dms[2] = (d-dms[0]-(dms[1]/60))*3600
    dms[1] = abs(dms[1])
    dms[2] = abs(dms[2])
    dms = tuple(dms)
    return dms

"""
### END OF PYTHON FUNCTION:  WSPDmax
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  DistBear
###
### NAME: DISTBEAR
###
### PROGRAMMER: CAPT MATT LESKO/DOC
###
### CREATED ON: 21 MAR 97
### AS OF DATE: 04 Jan 99
### BASELINED:  20 FEB 04 for BCR#23677
###
### PURPOSE:  "GREAT CIRCLE DISTANCE" CALCULATOR.
###
###           CALCULATES DISTANCE AND BEARING GIVEN TWO POINTS IN
###           LATITUDE AND LONGITUDE.  THE RESULTING CALCULATIONS TAKE
###           INTO ACCOUNT THE CURVATURE OF THE EARTH.
###           THE BEARING IS CALCULATED FROM THE POINT INITIALLY
###           INDICATED AS "LAT1" AND "LON1"  TO THE POINT INDICATED
###           AS "LAT2" AND "LON2".
###
###           NOTE:  THIS PROGRAM IS ADAPTED FROM A FORTRAN FUNCTION
###           DEVELOPED FROM AFGWC, CALLED "DBEAR".  THE MODELLING
###           SECTION AT AFGWC IS THE INITIAL SOURCE FOR THIS PROGRAM.
###
### METHOD:   1)  Error check that lat & lon are within limits.
###
###           2)  THE QUADRANT CHECK FLAG IS INITIALIZED AND DECLINATIONS
###           ARE CALCULATED.
###
###           3)  THE POLAR ANGLE BETWEEN LONGITUDES AND THE SUBTENDED
###           ANGLE ARE CALCULATED.
###
###           4)  THE DISTANCE IS CALCULATED BASED ON THE MEAN OF THE
###           LATITUDES.
###
###           5)  THE BEARING, FROM POINT 1 TO POINT 2, IS CALCULATED.
###
### SAMPLE MACRO CALLS:
###
###      %DISTBEAR(lat1,lon1, lat2,lon2, dist, bear)
###
### DATASETS USED:
###
###    NONE.  USED WITHIN A DATASTEP.
###
### OUTPUT:  NONE.
###
### HISTORY:
###
###   21 MAR 97 - CREATED BY CAPT MATT LESKO/DOC
###   04 JAN 99 - Changed to have argument calls.  Mr Kiess / DOC2
###   20 FEB 04 - Baseline Code.  Capt Freestrom / DOMD
###   20 MAY 15 - Removed break from prolog comments.  Mr Kiess / WXE
###   08 JAN 19 - Forced BEARING = 1 when BEARING>1 in math.asin(BEARING).
###               Mr Budai / WXC
###
### ------ DATA DICTIONARY ------
###
###    NAME     TYPE   DESCRIPTION
###    ----     ----   -----------
###    BEARING  NUM    BEARING IN DEGREES (POINT 1 TO POINT 2) (OUTPUT)
###    DISTANCE NUM    DISTANCE BETWEEN POINTS IN NAUTICAL MILES (OUTPUT)
###    LAT1     NUM    POINT ONE LATITUDE INDICATOR (USER INPUT)
###    LAT2     NUM    POINT ONE LATITUDE INDICATOR (USER INPUT)
###    LON1     NUM    POINT ONE LONGITUDE INDICATOR (USER INPUT)
###    LON2     NUM    POINT TWO LONGITUDE INDICATOR (USER INPUT)
###    M_BETA   NUM    POLAR ANGLE IN RADIANS BETWEEN LONGITUDES
###    M_COSD   NUM    TEMPORARILY USED IN CALCULATING "D"
###    M_D      NUM    SUBTENDED ANGLE IN RADIANS
###    M_DGTORD NUM    DEGREE TO RADIANS CONVERSION RATIO
###    M_EPSLON NUM    TOLERANCE FOR TESTING REAL VALUES
###    M_EQUAT  NUM    DISTANCE FROM CENTER OF THE EARTH TO THE EQUATOR
###                    IN KILOMETERS
###    M_II     NUM    QUADRANT FLAG FOR DETERMINING BEARING
###    M_LAT    NUM    POINT ONE LATITUDE
###    M_LAT2   NUM    POINT TWO LATITUDE
###    M_LON    NUM    POINT ONE LONGITUDE
###    M_LON2   NUM    POINT TWO LONGITUDE
###    M_PI     NUM    CONSTANT
###    M_PI2    NUM    M_PI DIVIDED BY 2
###    M_POLE   NUM    DISTANCE FROM CENTER OF EARTH TO THE POLE
###                    IN KILOMETERS
###    M_RADE   NUM    RADIUS OF EARTH AT MEAN LATITUDE
###    M_RAD1   NUM    INTERMEDIATE VARIABLE TO OBTAIN RADIUS
###    M_RAD2   NUM    INTERMEDIATE VARIABLE TO OBTAIN RADIUS
###    M_TWOPI  NUM    M_PI MULTIPLIED BY 2
###    M_Y1     NUM    DECLINATION OF POINT 1
###    M_Y2     NUM    DECLINATION OF POINT 2
###    M_YMEAN  NUM    ARITHMETIC MEAN OF LATITUDE 1 AND LATITUDE 2
###############################################################################
"""
def DistBear(LAT,LON,LAT2,LON2,**kwargs):
    if kwargs:
        rtn = kwargs['rtn'].upper()
    else:
        rtn = 'ALL'
    PI = 3.14159265
    DGTORD=PI/180.0
    POLE=(6356912/1853.248)
    EQUAT=(6378388/1853.248)
    TWOPI=(2*PI)
    PI2=PI/2.0
    EPSLON=.0001
    II=-1
    DISTANC = 0
    BEARING = 0
    if (-90>LAT>90) or (-90>LAT2>90) or (-180>LON>180) or (-180>LON2>180):
        print('Check Lat,Lon')
    ##CALCULATE DECLINATIONS (90-LAT1  AND  90-LAT2) IN RADIANS.
    Y1 = PI2-LAT*DGTORD
    Y2 = PI2-LAT2*DGTORD
    ##CALCULATE THE POLAR ANGLE BETWEEN THE LONGITUDES AND TAKE   
    ##THE SMALLEST ANGLE <= PI RADIANS.  IF ANGLE IS > PI          
    ##INITIALLY, SET QUADRANT FLAG TO ZERO.
    BETA = abs(LON-LON2)*DGTORD
    if (BETA > PI):
        BETA = TWOPI - BETA
        II = 0
    ## GET SUBTENDED ANGLE IN RADIANS, "MD".
    COSD = math.cos(Y1)*math.cos(Y2) + math.sin(Y1)*math.sin(Y2)*math.cos(BETA)
    if COSD < -1:
        COSD = -1
    elif COSD > 1:
        COSD = 1
    D = math.acos(COSD)
    if D == 0:
        DISTANC = 0
        BEARING = 0
    else:
        ##CALCULATE THE MEAN LATITUDE AND THE EARTHS RADIUS AT THE  
        ##MEAN LATITUDE.
        YMEAN = (LAT+LAT2)*DGTORD/2.0
        RAD1 = EQUAT*POLE;
        RAD2 = math.sqrt((POLE*math.cos(YMEAN))**2 + (EQUAT*math.sin(YMEAN))**2);
        RADE = RAD1/RAD2;
        ##FROM THE RADIUS AT THE MEAN LATITUDE AND THE ANGLE "D"  
        ##IN RADIANS, CALCULATE THE DISTANCE IN NAUTICAL MILES.    
        DISTANC = D*RADE
        if round(D,8) == PI:
            BEARING = 0
        else:
            BEARING = math.sin(Y2)*math.sin(BETA)/math.sin(D)
            BEAR2 = abs(BEARING)
            if BEAR2 > 1:
                BEAR2 = 1
            elif BEAR2 < 0 :
                BEARING = -1 * BEAR2
            else:
                BEARING = BEAR2
            if BEARING > 1: BEARING = 1
            BEARING = math.asin(BEARING)
            ##CALCULATE LATITUDE OF POINT 2 IF BEARING WAS 90 DEGREES
            TLAT = math.acos(math.cos(Y1)*COSD)
            ##NOW DETERMINE ACTUAL BEARING BASED ON LONGITUDES AND QUADRANT FLAG.
            if (II==-1) & (LON2 < LON):
                if (Y2<TLAT):
                    BEARING = TWOPI-BEARING
                else:
                    BEARING = BEARING+PI
            elif (II==-1) & (Y2>TLAT):
                BEARING = PI-BEARING
            elif (II==0) & (LON<=LON2):
                if (Y2<TLAT):
                    BEARING = TWOPI-BEARING
                else:
                    BEARING = BEARING+PI
            elif (II==0) & (Y2>TLAT):
                BEARING = PI-BEARING
                
            if (abs(LAT-90)<=EPSLON):
                if (LON2<0):
                    BEARING = abs(LON2)*DGTORD
                else:
                    BEARING = (360-LON2)*DGTORD
            elif (abs(LAT2-90) <= EPSLON):
                if LAT2<0 :
                    BEARING = 0

    BEARING = BEARING/DGTORD
    if rtn == 'ALL':
        return DISTANC,BEARING
    elif rtn == 'DIST':
        return DISTANC
    elif rtn == 'BEAR':
        return BEARING 
"""
### END OF PYTHON FUNCTION:  DistBear
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  C_to_F
###
### INPUTS:  TEMPC
###
### OUTPUTS:  TEMPF
###
### PURPOSE: convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""

def C_to_F(tempc):
    if pd.isnull(tempc) == False:
        temp = round((tempc * 9.0/5.0) + 32,2)
    else:
        temp = np.nan
        
    return temp

"""
### END OF PYTHON FUNCTION:  C_to_F
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  F_to_C
###
### INPUTS:  TEMPF
###
### OUTPUTS:  TEMPC
###
### PURPOSE: convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""

def F_to_C(tempf):
    if pd.isnull(tempf) == False:
        temp = round((tempf - 32) * 5/9,2)
    else:
        temp = np.nan
        
    return temp

"""
### END OF PYTHON FUNCTION:  C_to_F
###
###############################################################################
###############################################################################
###############################################################################
"""

"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  C_to_K
###
### INPUTS:  TEMPC
###
### OUTPUTS:  TEMPK
###
### PURPOSE: convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""


def C_to_K(tempc):
    if pd.isnull(tempc) == False:
        temp = round(tempc + 273.15,2)
    else:
        temp = np.nan
        
    return temp

"""
### END OF PYTHON FUNCTION:  C_to_K
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  K_to_C
###
### INPUTS:  TEMPK
###
### OUTPUTS:  TEMPC
###
### PURPOSE: convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""


def K_to_C(tempk):
    if pd.isnull(tempk) == False:
        temp = round(tempk - 273.15,2)
    else:
        temp = np.nan
        
    return temp

"""
### END OF PYTHON FUNCTION:  K_to_C
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  K_to_F
###
### INPUTS:  TEMPK
###
### OUTPUTS:  TEMPF
###
### PURPOSE:convert
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""


def K_to_F(tempk):
    if pd.isnull(tempk) == False:
        tempc = K_to_C(tempk)
        tempf = tempc*1.8 + 32
    else:
        tempf = np.nan
        
    return tempf

"""
### END OF PYTHON FUNCTION:  K_to_F
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  ROLLhour
###
### INPUTS:  obtime - observationtime
###         reptype - reporttypecode (FM-15, FM-12, etc.)
###         kwarg - roll_all = bool
###
### OUTPUTS:  the appropriate time to the 00z of the next hour
###
### PURPOSE: accuratly assign an hourly ob to its hour
###
### METHOD: if roll_all all obs with minute value GT 30 goes to next hour and LT
###          goes to previous
###         else if minute GTE 45 reptype is FM-15 roll to top of the hour
###
### REFERENCES: METARs for the 0100Z go out anywhere between 15min
###              before and 10 min after 0100z
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""

def ROLLhour(obtime,reptype,**kwargs):
    if kwargs and ('roll_all' in kwargs.keys()) and kwargs['roll_all']:
        if pd.isnull(reptype):
            # if reporttypecode is blank return blank
            return pd.Timestamp(obtime)
        obtime = pd.Timestamp(obtime)
        r_obtime = obtime
        if (obtime.minute >= 30) :##and reptype not in ('FM-16','FM-12','SAOSP'):
            r_obtime = obtime + pd.Timedelta(str(abs(obtime.minute - 60))+' minutes')
        else:
            r_obtime = obtime - pd.Timedelta(str(abs(obtime.minute))+' minutes')
        return r_obtime
    else:
        if pd.isnull(reptype):
            return pd.Timestamp(obtime)
        obtime = pd.Timestamp(obtime)
        r_obtime = obtime
        if (obtime.minute >= 45) and reptype not in ('FM-16','FM-12','SAOSP'):
            r_obtime = obtime + pd.Timedelta(str(abs(obtime.minute - 60))+' minutes')
        return r_obtime

"""
### END OF PYTHON FUNCTION:  ROLLhour
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  ChangDate
###
### INPUTS:  observationtime
###          time_conv - time off Zulu (-5 for EST)
###
### OUTPUTS:  local_observationtime
###
### PURPOSE: convert to local time
### 
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones port to Python
###
###############################################################################
"""


def CHANGEdate(obtime,time_conv):
    obtime = pd.Timestamp(obtime)
    L_obtime = obtime + pd.Timedelta(str(time_conv)+' hours')
    return L_obtime

"""
### END OF PYTHON FUNCTION:  CHANGEdate
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  ParsDate
###
### INPUTS:  observationtime
###
### OUTPUTS:  year,month,day,hour,minute
###
### PURPOSE: parse out the time elements in easy function
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07Jan2019  SSgt Jones Created
###
###############################################################################
"""


def ParsDate(obtime):
    obtime = pd.Timestamp(obtime)
    year = obtime.year
    month = obtime.month
    day = obtime.day
    hour = obtime.hour
    minute = obtime.minute
    return year,month,day,hour,minute

"""
### END OF PYTHON FUNCTION:  ParsDate
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  hour_buckets
###
### INPUTS:
### - int that 24 is divisable by
###
### - postion of cardinal hr,('leading','center','trailing')
###
### - ex // hour_buckets(3,'leading') returns hours = [0,3,6...], buckets = [[0,1,2],[3,4,5]....]
###
### - ex // hour_buckets(6,'center') returns hours = [0,6,12,18], buckets = [[22,23,0,1,2,3],[4,5,6,7,8,9]....]
###
### OUTPUTS:  
### - hours = list of hours to be cardnal hours in DF index
### - buckets = list of list containing the actual hours that go into each hour bucket
### PURPOSE:
### - used in Freq to correctly group hours into buckets if required
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:01Jan2020  SSgt Jones Created
###
###############################################################################
"""
def hour_buckets(x,position):
    if (24 % x != 0):
        raise NameError('24 hours is not divisable by selected combinehours')
    hours = [x for x in range(0,23,x)]
    buckets = []
    if position == 'center':
        for i in range(0,23,x):
            l = []
            l.append(i)
            if i==0:
                eb=24
            else:
                eb=i
            ef = i
            if x % 2 == 0:
                for a in range(int((x)/2)):
                    if a+1 == int((x)/2):
                        ef+=1
                        l.append(ef)
                    else:
                        eb-=1
                        ef+=1
                        l.append(eb)
                        l.append(ef)
            else:
                for a in range(int((x-1)/2)):
                    ef+=1
                    eb-=1
                    l.append(eb)
                    l.append(ef)

            buckets.append(l)
            
        return hours,buckets

    elif position == 'leading':
        for i in range(0,23,x):
            l = []
            l.append(i)
            ef = i
            for s in range(x-1):
                ef+=1
                l.append(ef)
            buckets.append(l)
            
        return hours,buckets

    elif position == 'trailing':
        for i in range(0,23,x):
            l = []
            l.append(i)
            if i==0:
                eb=24
            else:
                eb=i
            for s in range(x-1):
                eb-=1
                l.append(eb)
                
            buckets.append(l)
            
        return hours,buckets
             
"""
### END OF PYTHON FUNCTION:  hour_buckets
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  Mean
###
### INPUTS:
### - df is a DataFrame with either 'ROLLED_OBSERVATIONTIME' or
###      'OBSERVATIONTIME'
### - TimeCat is string holding any 'ymdh' for Year Mo Day Hour
### - day =  string max,mean,min // how to do daily res
###                              // used for monthly mean of daily max/min
### - hr = string max,mean,min // how to do hourly res //
###
###
### OUTPUTS:  {List function output variables and units}
### - Output DF with 'ymdh' cols and means of parameter entered
###
### PURPOSE: standard way to get any time parsed mean and daily mean max/min's  
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 1 Jan 2020 SSgt Jones - Created for user ease
###
###############################################################################
"""


def Mean(df,TimeCat,**kwargs):
    TimeCat = TimeCat.lower()
    TCL = len(TimeCat)
    ##########################################################################################
    ########################### Reset Index if needed ########################################
    ##########################################################################################
    if type(df.index)  != pd.core.indexes.range.RangeIndex:
        if ('ROLLED_OBSERVATIONTIME' or 'OBSERVATIONTIME') not in df.index.values:
            df.reset_index(inplace=True,drop=True)
        else:
            df.reset_index(inplace=True)
    ##########################################################################################
    ########################### Raise Error if no OBSERVATIONTIME ############################
    ##########################################################################################
    if ('ROLLED_OBSERVATIONTIME') not in df.columns.values:
        if ('OBSERVATIONTIME') not in df.columns.values:
            raise NameError('ROLLED_OBSERVATIONTIME or OBSERVATIONTIME needs to be in DF' )
    ##########################################################################################
    ############################ Parse Date ##################################################
    ##########################################################################################
    if 'ROLLED_OBSERVATIONTIME' in df.columns:
        df['YEAR'],df['MO'],df['DAY'],df['HR'],df['MINUTE']= np.vectorize(ParsDate)(df['ROLLED_OBSERVATIONTIME'])
        df.drop(['ROLLED_OBSERVATIONTIME','MINUTE'],axis=1,inplace = True)
        df['TRIMO'] = np.where(df['DAY'] < 11,1,0)
        df['TRIMO'] = np.where((df['DAY'] >= 11) & (df['DAY'] < 21),2,df['TRIMO'])
        df['TRIMO'] = np.where(df['DAY'] >= 21,3,df['TRIMO'])
        dfmax = df.groupby(['YEAR','MO','TRIMO','DAY','HR']).max()
        
    else:
        df['YEAR'],df['MO'],df['DAY'],df['HR'],df['MINUTE']= np.vectorize(ParsDate)(df['OBSERVATIONTIME'])
        df.drop(['OBSERVATIONTIME','MINUTE'],axis=1,inplace = True)
        df['TRIMO'] = np.where(df['DAY'] < 11,1,0)
        df['TRIMO'] = np.where((df['DAY'] >= 11) & (df['DAY'] < 21),2,df['TRIMO'])
        df['TRIMO'] = np.where(df['DAY'] >= 21,3,df['TRIMO'])
        dfmax = df.groupby(['YEAR','MO','TRIMO','DAY','HR']).max()

    # getting the overall group by list in order
    GBL = []
    for l in TimeCat:
        if l == 'y':
            GBL.append('YEAR')
        elif l == 'm':
            GBL.append('MO')
        elif l == 't':
            GBL.append('TRIMO')            
        elif l == 'd':
            GBL.append('DAY')
        elif l == 'h':
            GBL.append('HR')

    
    # step 1 resolve dataframe to hour res based on user input or to mean (default) 
    if kwargs and 'hr' in kwargs.keys():
        if kwargs['hr'] == 'max':
            dfmax = df.groupby(['YEAR','MO','TRIMO','DAY','HR']).max()
        elif kwargs['hr'] == 'mean':
            dfmax = df.groupby(['YEAR','MO','TRIMO','DAY','HR']).mean()
        else:
            dfmax = df.groupby(['YEAR','MO','TRIMO','DAY','HR']).min()
    else:
        dfmax = df.groupby(['YEAR','MO','TRIMO','DAY','HR']).mean()

    if 'HR' in GBL:
        dfn = dfmax.groupby(GBL).mean()
        return dfn

    
    # step 2 resolve dataframe to daily res based on user input else move on 
    if kwargs and 'day' in kwargs.keys():
        if 'TRIMO' in GBL:
            if kwargs['day'] == 'max':
                dfmax = dfmax.groupby(['YEAR','MO','TRIMO','DAY']).max()
            elif kwargs['day'] == 'mean':
                dfmax = dfmax.groupby(['YEAR','MO','TRIMO','DAY']).mean()
            else:
                dfmax = dfax.groupby(['YEAR','MO','TRIMO','DAY']).min()
        else:
            if kwargs['day'] == 'max':
                dfmax = dfmax.groupby(['YEAR','MO','DAY']).max()
            elif kwargs['day'] == 'mean':
                dfmax = dfmax.groupby(['YEAR','MO','DAY']).mean()
            else:
                dfmax = dfmax.groupby(['YEAR','MO','DAY']).min()
        if 'DAY' in GBL:
            if (GBL == ['YEAR','MO','DAY']) or (GBL == ['YEAR','MO','TRIMO','DAY']):
                return dfmax
            else:
                dfn = dfmax.groupby(GBL).mean()
                return dfn
    else:
        if 'DAY' in GBL:
            dfn = dfmax.groupby(GBL).mean()
            return dfn
        else:
            if 'TRIMO' in GBL:
                dfmax = dfmax.groupby(['YEAR','MO','TRIMO','DAY']).mean()
            else:
                dfmax = dfmax.groupby(['YEAR','MO','DAY']).mean()
        
    # step 3 resolve dataframe to trimo res based on user input else move to month
    if 'TRIMO' in GBL:
        dfmax = dfmax.groupby(['YEAR','MO','TRIMO']).mean()
    else:
        dfmax = dfmax.groupby(['YEAR','MO']).mean()
        
    if (GBL == ['YEAR','MO']) or (GBL == ['YEAR','MO','TRIMO']):
        return dfmax
    else:
        dfn = dfmax.groupby(GBL).mean()
        return dfn

"""
### END OF PYTHON FUNCTION:  Mean
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  Freq
###
### INPUTS:
### - df is a DataFrame with either 'ROLLED_OBSERVATIONTIME' or
###      'OBSERVATIONTIME' and flaged parameter col should only
###      hold 1 or 0
### - TimeCat is string holding any 'ymdh' for Year Mo Day Hour
###
### - combinehours = int number of hours to combine // will center if odd if
###     even will have one extra in the future // 3 = (23z,00z,01z) // 4 = (23z,00z,01z,02z)
###
###
###
###
###
### OUTPUTS:  
### - Output DF with 'ymdh' cols and % of time parameter occuered
###
### PURPOSE: An easy to use function that correctly handles the
###           Freqency calculation it's also able to bin hours
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 1 Jan 2020 SSgt Jones - Created for user ease
###
###############################################################################
"""


def Freq(df,TimeCat,**kwargs):
    TimeCat = TimeCat.lower()
    df = df.fillna(0)
    TCL = len(TimeCat)
    ##########################################################################################
    ########################### Reset Index if needed ########################################
    ##########################################################################################
    if type(df.index)  != pd.core.indexes.range.RangeIndex:
        if (('ROLLED_OBSERVATIONTIME' not in df.index) or ('OBSERVATIONTIME'not in df.index) or ('LOCAL_OBSERVATIONTIME'not in df.index)) :
            df.reset_index(inplace=True,drop=True)
        else:
            df.reset_index(inplace=True)
    ##########################################################################################
    ########################### Raise Error if no OBSERVATIONTIME ############################
    ##########################################################################################
    if ('ROLLED_OBSERVATIONTIME') not in df.columns.values:
        if ('OBSERVATIONTIME') not in df.columns.values:
            if ('LOCAL_OBSERVATIONTIME') not in df.columns.values:
##                if sum([i in ['YEAR','MO','DAY','HR'] for i in df.columns]) == 4:
##                    df['OBSERVATIONTIME'] = pd.to_datetime(dict(year=df.YEAR,month=df.MO,day=df.DAY,hour=df.HR))
##                elif sum([i in ['YEAR','MO','DAY','HR'] for i in df.columns]) == 3:
##                    df['OBSERVATIONTIME'] = pd.to_datetime(dict(year=df.YEAR,month=df.MO,day=df.DAY))
##                else:
                raise NameError('ROLLED_OBSERVATIONTIME or OBSERVATIONTIME needs to be in DF' )

    ##########################################################################################
    ###################################### Create OBSERVATIONTIME ############################
    ##########################################################################################

    
    
        
    ##########################################################################################
    ############################ Parse Date ##################################################
    ##########################################################################################

    GBL = []
    for l in TimeCat:
        if l == 'y':
            GBL.append('YEAR')
        elif l == 'm':
            GBL.append('MO')
        elif l == 't':
            GBL.append('TRIMO')
        elif l == 'x':
            GBL.append('HEXMO')
        elif l == 'd':
            GBL.append('DAY')
        elif l == 'h':
            GBL.append('HR')

            
    if 'ROLLED_OBSERVATIONTIME' in df.columns:
        df['YEAR'],df['MO'],df['DAY'],df['HR'],df['MINUTE']= np.vectorize(ParsDate)(df['ROLLED_OBSERVATIONTIME'])
        df.drop(['ROLLED_OBSERVATIONTIME','MINUTE'],axis=1,inplace = True)
        df['HEXMO'] = np.vectorize(hexmo)(df['DAY'])
        df['TRIMO'] = np.vectorize(trimo)(df['DAY'])
        dfmax = df.groupby(['YEAR','MO','TRIMO','HEXMO','DAY','HR']).max()

    elif 'LOCAL_OBSERVATIONTIME' in df.columns:
        df['YEAR'],df['MO'],df['DAY'],df['HR'],df['MINUTE']= np.vectorize(ParsDate)(df['LOCAL_OBSERVATIONTIME'])
        df.drop(['LOCAL_OBSERVATIONTIME','MINUTE'],axis=1,inplace = True)
        df['HEXMO'] = np.vectorize(hexmo)(df['DAY'])
        df['TRIMO'] = np.vectorize(trimo)(df['DAY'])
        dfmax = df.groupby(['YEAR','MO','TRIMO','HEXMO','DAY','HR']).max()
        
    else:
        df['YEAR'],df['MO'],df['DAY'],df['HR'],df['MINUTE']= np.vectorize(ParsDate)(df['OBSERVATIONTIME'])
        df.drop(['OBSERVATIONTIME','MINUTE'],axis=1,inplace = True)
        df['HEXMO'] = np.vectorize(hexmo)(df['DAY'])
        df['TRIMO'] = np.vectorize(trimo)(df['DAY'])
        dfmax = df.groupby(['YEAR','MO','TRIMO','HEXMO','DAY','HR']).max()
        
        
    if (kwargs and ('POR' in kwargs.keys()) and (kwargs['POR'] == True)):
        dfmax = dfmax.reset_index()
        dfmax['YEAR'] = 13
        dfmax = dfmax.set_index(['YEAR','MO','TRIMO','HEXMO','DAY','HR'])
        
    if kwargs and ('ANY' in kwargs.keys()) and (kwargs['ANY'] == True):
        #print(dfmax)
        dfmax['ANY'] = dfmax.max(axis=1,skipna=True)
        #print(dfmax)


            
    if kwargs and ('combinehours' in kwargs.keys()) and (kwargs['combinehours']>0) and ('HR' in GBL):
        if 'position' in kwargs.keys():
            position = kwargs['position']
        else:
            position = 'center'
        hours,buckets = hour_buckets(kwargs['combinehours'],position)

        dfmax = dfmax.reset_index()        
        for i in range(len(hours)):
            dfmax.loc[:,'HR'] = np.where(dfmax.loc[:,'HR'].isin(buckets[i]),hours[i],dfmax.loc[:,'HR'])
        dfmax = dfmax.set_index(['YEAR','MO','TRIMO','HEXMO','DAY','HR'])

        dfn = dfmax.groupby(GBL).sum()/dfmax.groupby(GBL).count()
        
    else:
        dfn = dfmax.groupby(GBL).sum()/dfmax.groupby(GBL).count()

    return dfn

"""
### END OF PYTHON FUNCTION:  Freq
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION: U & V to Dir and Speed
###
### INPUTS:  u_component,v_component
###
### OUTPUTS:  degree,speed
###
### PURPOSE:
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:1 Jan 2020 SSgt Jones ported to python
###
###############################################################################
"""

def U_and_V_Calc(u_component,v_component):
    "This function calculates meteorological wind direction from u and v wind components"
    met_degrees=0.0
    if (math.isnan(u_component) or math.isnan(v_component)):
        met_degrees=np.nan
    else:
        math_degrees=math.degrees(math.atan2(v_component,u_component))
        math_degrees=math_degrees+180.0
        if (math_degrees < 0.0):
            math_degrees=math_degrees+360.0
        if (0.0 <= math_degrees <= 90.0 ):
            met_degrees=90.0-math_degrees
        if (math_degrees > 90.0):
            met_degrees=450.0 - math_degrees
        if(u_component==0.0 and v_component==0.0): met_degrees=0.0
    speed=math.sqrt(u_component**2+v_component**2)
    return met_degrees,speed

"""
### END OF PYTHON FUNCTION:  U & V Direction
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION: Absolute Humidity
###
### INPUTS:  RH % (0-100), tempC
###
### OUTPUTS:  Absolute Humidity
###
### PURPOSE:
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07 Jan 2019 SSgt Jones port to Python
###          1 Jan 2020 SSgt Jones added g/m^3
###
###############################################################################
"""
# not best calculation but was the one ported from SAS, not used in object handler
def ABShum(RH,temp):
# RH as a percent, temp in Celsius, abshum returned as kg/m^3
    
    if (temp < -20.0):   # For saturation vapour pressure over ice
        A=6.114742
        N=273.1466
        M=9.778707
    else:                # For saturation vapour pressure over water
        A=6.116441
        N=240.7263
        M=7.591386
    abshum=(2.16679/(temp+273.15))*(RH/100.0)*A*pow(10.0,(M*temp)/(temp+N))

    
    return(abshum)

# best calculation, used in object handler
def CalcABHUM(rh,t):
    ## input vars are RH (%) & TEMPC
    "Calculate absolute humidity in g/m^3 from T(K) and RH(%)"
#    "Calculate dewpoint temperature in degrees Kelvin from T(K) and RH(%)"
    t+=273.15
    tk=t
    if rh==0:
        abhum=np.nan
    else:
        # Saturation vapor pressure over water above 0C
        if (t>=273.15):
            # Convert temperature to degress Celsius for this e_s formula
            t=t-273.15
            e_s=6.112*math.exp( (17.67*t)/(t+243.5) )
        # Saturation vapor pressure over water between -20C and 0C
        # Goff-Gratch equation (Smithsonian Tables, 1984)
        elif ( (t<273.15)&(t>=253.15) ):
            e_s=math.pow(10,-7.90298*(373.15/t-1)+
                            5.02808*math.log10(373.17/t)-
                            1.3816e-7*(math.pow(10,11.344*(1-t/373.16))-1)+
                            8.1328e-3*(math.pow(10,-3.49149*(373.16/t-1))-1)+
                            math.log10(1013.246) )
        # Saturation vapour pressure over ice
        # Goff-Gratch equation (Smithsonian Tables, 1984)
        else:
            e_s=math.pow(10,-9.09718*(273.16/t-1)-
                            3.56654*math.log10(273.16/t)+
                            0.876793*(1-t/273.16)+
                           math.log10(6.1071) )

        e=(rh/100)*e_s
#        dp=(243.5*math.log(e/6.112))/(17.67-math.log(e/6.112))
#        dp=dp+273.15
        abhum=(216.5*e)/tk

    return(abhum)



"""
### END OF PYTHON FUNCTION:  Absolute Humidity
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  FITS
###
### INPUTS:  {List function input variables and units}
###
### OUTPUTS:  {List function output variables and units}
###
### PURPOSE:
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07 Jan 2019 SSgt Jones port to Python
###
###############################################################################
"""

def FITS_C(tempc,wtbulb):
    FitsC = 0.83*wtbulb+.35*tempc+5.08
    return FitsC

def FITS_F(tempf,wtbulb):
    FitsF = 0.83*wtbulb+.35*tempf+9.14
    return FitsF
"""
### END OF PYTHON FUNCTION:  FITS
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  WINDchill
###
### INPUTS:  {List function input variables and units}
###
### OUTPUTS:  {List function output variables and units}
###
### PURPOSE:
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES:07 Jan 2019 SSgt Jones port to Python
###
###############################################################################
"""

def WINDchill(tempf,wspdmph):
    if wspdmph < 3:
        wc = tempf
    else:
        wc = 35.74+(0.6215*tempf)-(35.75*(wspdmph**0.16))+(0.4275*(tempf*(wspdmph**0.16)))
    return wc

"""
### END OF PYTHON FUNCTION:  WINDchill
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  DAYnight
###
### FROM LEGACY SAS PROGRAM:  daynite (Day/Night)
###
### BCR TRACKING NUMBER:
###
### PURPOSE:  To return the elevation angle and a sun flag determining
###            whether it is day or night.
###
### METHOD:  Compute the date and the position of the sun.  From there,
###          compute the hour angles and elevation of the sun.  Default
###          day to be sun elevation angle > -6 degrees.
###
### REFERENCES:
###   none
###
### LEGACY SAS PROCESS NARRATIVE LOCATION:
###
### LEGACY SAS MACROS/SUBROUTINES CALLED:
###   none
###
### FILES/DATABASE/TABLES ACCESSED:
###   none
###
### LEGACY SAS DATA DICTIONARY:
###   &day         Input: day
###   &elvang      Output: the sun elevation angle in degrees
###   &hr          Input: hour
###   &lat         Input: latitude in decimal degrees
###   &lon         Input: longitude in decimal degrees
###   &mo          Input: month
###   &sun         Output: 0 for night and 1 for day
###   &year        Input: year
###   _decl        declination in radians
###   _delta       used in calc jdate (yr-29)
###   _denom       the intermediate denominator for rasc
###   _elvang      elevation angle
###   _gmst        greenwich mean sidereal time
###   _hangl       hour angle in radians
###   _hold1       temp variable used in calc mnlong
###   _hold2       temp variable used in calc mnlong
###   _jdate       the number of days since 00 ut, 1 jan 4713 bc.
###   _jday        julian day of the year ddd
###   _julian      julian date yyddd
###   _lat         latitude of blkstn in decimal degrees
###   _leap        number of leap years
###   _leclp       ecliptic longitude in radians
###   _lon         longitude of blkstn in decimal degrees
###   _lmst        local mean sidereal time
###   _mnanm       mean anomaly of the sun
###   _mnanm1      used in calc of mean anomaly of the sun
###   _mnanm2      used in calc of mean anomaly of the sun
###   _mnlon       mean longitude of the sun
###   _numer       the intermediate numerator for rasc
###   _obecl       obliquity of ecliptic in radians
###   _rasc        right ascension in radians
###   _sasdate     the sas date
###   _sun         sun flag (day = 1, night = 0)
###   _time        time period centered about 12z on 1 jan 2000
###
### REMARKS:
###   The macro is default for -6 degrees twilight. The user
###   has the option to override this in the main program.
###   The macro expects positive longitude east.
###
### UPDATES:
###   05 DEC 97    SRA RICHARD WRIGHT SYS/TEAM2, PORTING PROJECT FROM
###                MAINFRAME TO UNIX
###                        **** NOTE ****
###                PIECED TOGETHER FROM ADHOC CODE
###
###   07 Jan 99    Mr Kiess DOC2:
###                Lat and Lon are now expected to be in decimal degrees
###                Fixed problem of macro arguments not properly used
###                Dropped internal variables before returning
###
###   07 Jan 2019  SSgt Jones ported to Python
#################################################################################
"""


def DAYnight(yr,mo,day,hr,lat,lon):
    #REMARKS:
    #The macro is default for -6 degrees twilight. The user
    #has the option to override this in the main program.
    #The macro expects positive longitude east.
    lon = lon * (-1)

    # Compute Julian Date using .toordinal() finds # of days since 1 Jan 1
    #then add 1721424.5 and hr/24 to get True # of days since, 00 UT 1 JAN 4713 BC--Levi Jones**
    mydate = dt.date(yr,mo,day)
    jdate = (mydate.toordinal() + 1721424.5) + (hr/24.)

    #*  COMPUTE THE JULIAN DATE (JDATE) WHICH IS THE NUMBER OF DAYS
    #SINCE 00 UT, 1 JAN 4713 BC.  THE FIRST TERM IN THE EXPRESSION
    #FOR "RJDATE" IS THE JULIAN DATE FOR MIDNIGHT 0 JANUARY 1929 UT
    #(ACTUALLY 2,400,000 MUST BE ADDED).  IN THE ORIGINAL REFERENCE
    #JULIAN DATE FOR MIDNIGHT 0 JANUARY 1949 UT WAS GIVEN.

    #Had to add +6 hr converion to adjust location--Levi Jones**
    hr_adjust = hr + 6
    if hr_adjust > 23 :
        hr_adjust = hr_adjust - 24
    #jday = int(datetime.date(yr,mo,day).strftime("%j"))
    #delta = yr - 1929
    #leap = int(delta/4)
    #jdate = 25611.5 + (delta * 365) + leap + jday + hr/24

    # Convert Julian DATE TO A TIME PERIOD CENTERED ABOUT 12Z ON 1 JAN 2000.
    time = jdate - 2451545.0
    #CALCULATE THE MEAN LONGITUDE OF THE SUN (MNLON),
    #AND FORCE ITS VALUE TO LIE BETWEEN 0.0 AND PI RADIANS.
    hold1 = 4.89495 + 0.0172028 * time
    hold2 = int(hold1/6.2831853)
    mnlon = hold1 - (hold2 * 6.2831853)
    #CALCULATE THE MEAN ANOMALY OF THE SUN (MNANM) IN RADIANS AND
    #FORCE ITS VALUE TO LIE BETWEEN 0 AND TWO*PI.
    mnanm1 = 6.24004 + 0.017202 * time
    mnanm2 = int(mnanm1/6.2831853)
    mnanm = mnanm1 - (mnanm2 * 6.2831853)
    #COMPUTE ECLIPTIC LONGITUDE (LECLP) AND OBLIQUITY OF THE ECLIPTIC (OBECL) IN RADIANS
    leclp = mnlon + 0.033423 * math.sin(mnanm) + 0.000349 * math.sin(2.0* mnanm)
    obecl = 0.409090 - 0.000000006 *time
    #CALCULATE RIGHT ASCENSION (RASC) IN RADIANS.  THE INTERMEDIATE
    #NUMER AND DENOM ARE CALCULATED IN ORDER TO PERFORM TESTS AND
    #ALTER THE CALCULATED VALUE OF RASC BY PI OR TWO*PI.
    numer = math.cos(obecl) * math.sin(leclp)
    denom = math.cos(leclp)
    rasc  = math.atan(numer/denom)
    if denom < 0 :
        rasc = rasc + 3.1415926
    if numer < 0 :
        rasc = rasc + 6.2831853
    #CALCULATE DECLINATION (DECL) IN RADIANS
    decl = math.asin(math.sin(obecl) * math.sin(leclp))
    #CALCULATE GREENWICH MEAN SIDEREAL TIME (GMST) IN HOURS,
    #AND FORCE IT TO LIE IN THE INTERVAL 0 TO 23 INCLUSIVE.
    gmst = 6.697375 + 0.0657098242 * time + hr_adjust
    gmst = gmst % 24
    if gmst < 0 :
        gmst = gmst + 24
    #CALCULATE LOCAL MEAN SIDEREAL TIME (LMST) IN HOURS, AND FORCE IT
    #TO LIE IN THE INTERVAL 0 TO 23 INCLUSIVE....THEN CONVERT IT TO A RADIAN VALUE.
    lmst = gmst - (lon/15.0)
    lmst = lmst % 24
    if lmst < 0 :
        lmst = lmst + 24
    lmst = lmst * 15 * 0.017453293
    #CALCULATE HOUR ANGLE (HANGL) IN RADIANS BETWEEN +PI AND -PI
    hangl = lmst - rasc
    if hangl < -3.1415926 :
        hangl = hangl + 6.2831853
    if hangl > 3.1415926 :
        hangl = hangl - 6.2831853
    #CALCULATE ELEVATION ANGLE (ELVANG)
    elvang = (math.asin(math.sin(decl) * math.sin(0.017453293 * lat) +
    math.cos(decl) * math.cos(0.017453293 * lat) * math.sin(hangl)))

    elvang = elvang * 180/3.1415926
    sun = 0
    if elvang > -6 :
        sun = 1

    zenith = 1.5708 - (elvang/57.29578)
    return sun , elvang , zenith

"""
### END OF PYTHON FUNCTION:  DAYnight
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  NoD_Normalize
###
### Inputs:
###   df =  DF with ROLLED_OBSERVATIONTIME OR OBSERVATIONTIME and flagged columns only
###   begpor = begpor
###   endpor = endpor
###   mo_r = decimal % of days required to count month .5 = approx 15 days
###   HN = int() Hours needed to count day as hit // Not good for 3 hr reporters
###
### Outputs:
###   DF with ROLLER_OBSERVATIONTIME OR OBSERVATIONTIME
###
### PURPOSE: An easy to use function that correctly handles the
###           Freqency calculation it's also able to bin hours
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 1 Jan 2020 SSgt Jones - Created for user ease
###
###############################################################################
"""
def NoD_Normalize(df,mo_r,**kwargs):
    ###########################################################
    ## Meant to count and normilize # of Days flagged based on days reported
    ## Inputs:
    ##   df =  DF with ROLLED_OBSERVATIONTIME OR OBSERVATIONTIME and flagged columns only
    ##   begpor = begpor
    ##   endpor = endpor 
    ##   mo_r = decimal % of days required to count month .5 = approx 15 days
    ##   HN = int() Hours needed to count day as hit // Not good for 3 hr reporters
    ## Outputs:
    ##   DF with ROLLER_OBSERVATIONTIME OR OBSERVATIONTIME
    ###########################################################
    
    ##########################################################################################
    ########################### Reset Index if needed ########################################
    ##########################################################################################
    if type(df.index)  != pd.core.indexes.range.RangeIndex:
        if (('ROLLED_OBSERVATIONTIME' not in df.index) or ('OBSERVATIONTIME' not in df.index) or ('LOCAL_OBSERVATIONTIME' not in df.index)) :
            df = df.reset_index(drop=True)
        else:
            df = df.reset_index()
    ##########################################################################################
    ########################### Raise Error if no OBSERVATIONTIME ############################
    ##########################################################################################
    if ('ROLLED_OBSERVATIONTIME') not in df.columns.values:
        if ('OBSERVATIONTIME') not in df.columns.values:
            if ('LOCAL_OBSERVATIONTIME') not in df.columns.values:
                raise NameError('ROLLED_OBSERVATIONTIME or OBSERVATIONTIME needs to be in DF' )
        
    ##########################################################################################
    ############################ Parse Date ##################################################
    ##########################################################################################
    if ('ROLLED_OBSERVATIONTIME') in df.columns.values:
        df['YEAR'] = df['ROLLED_OBSERVATIONTIME'].dt.year
        df['MO'] = df['ROLLED_OBSERVATIONTIME'].dt.month
        df['DAY'] = df['ROLLED_OBSERVATIONTIME'].dt.day
        df['HR'] = df['ROLLED_OBSERVATIONTIME'].dt.hour
        begpor = df.ROLLED_OBSERVATIONTIME.min().strftime('%Y%m%d%H%M')
        endpor = df.ROLLED_OBSERVATIONTIME.max().strftime('%Y%m%d%H%M')
        df = df.drop('ROLLED_OBSERVATIONTIME',axis=1)
        
    elif ('LOCAL_OBSERVATIONTIME') in df.columns.values:
        df['YEAR'] = df['LOCAL_OBSERVATIONTIME'].dt.year
        df['MO'] = df['LOCAL_OBSERVATIONTIME'].dt.month
        df['DAY'] = df['LOCAL_OBSERVATIONTIME'].dt.day
        df['HR'] = df['LOCAL_OBSERVATIONTIME'].dt.hour
        begpor = df.LOCAL_OBSERVATIONTIME.min().strftime('%Y%m%d%H%M')
        endpor = df.LOCAL_OBSERVATIONTIME.max().strftime('%Y%m%d%H%M')
        df = df.drop('LOCAL_OBSERVATIONTIME',axis=1)
        
    else:
        df['YEAR'] = df['OBSERVATIONTIME'].dt.year
        df['MO'] = df['OBSERVATIONTIME'].dt.month
        df['DAY'] = df['OBSERVATIONTIME'].dt.day
        df['HR'] = df['OBSERVATIONTIME'].dt.hour
        begpor = df.OBSERVATIONTIME.min().strftime('%Y%m%d%H%M')
        endpor = df.OBSERVATIONTIME.max().strftime('%Y%m%d%H%M')
        df = df.drop('OBSERVATIONTIME',axis=1)
        
        
    ##########################################################################################
    ############################ Max Flags To Daily Res ######################################
    ##########################################################################################
    if kwargs is not None and 'HN' in kwargs.keys():
        df = df.groupby(['YEAR','MO','DAY','HR']).max()
        df = df.groupby(['YEAR','MO','DAY']).sum()
        df = np.where(df < int(kwargs['HN']),0,df)
        df = np.where(df >= int(kwargs['HN']),1,df)
    else:
        df = df.drop('HR',axis=1)
        df = df.groupby(['YEAR','MO','DAY']).max()

    ##########################################################################################
    ############################ YMD Shell ###################################################
    ##########################################################################################
    
    dshelldict = {'DateTime': pd.date_range(str(begpor), str(endpor), freq='D')}
    # shell is data frame to aid in Freq calculations
    dshell = pd.DataFrame(data = dshelldict)
    # Create new Columns
    dshell['YEAR'] = dshell['DateTime'].dt.year
    dshell['MO'] = dshell['DateTime'].dt.month
    dshell['DAY'] = dshell['DateTime'].dt.day
    dshell = dshell.set_index(['YEAR','MO','DAY'])
    
    ##########################################################################################
    ############################ Merge and Purge #############################################
    ##########################################################################################
    
    dfd = dshell.merge(df,how = 'left',left_index = True,right_index=True)
    dfd = dfd.drop('DateTime',axis=1)

    ##########################################################################################
    ############################ Reshape Shell ###############################################
    ##########################################################################################
    
    for i in df.columns:
        dshell[i] = 1
    
    dshell = dshell.drop('DateTime',axis=1)

    ##########################################################################################
    ############################ Counts and Normalize ########################################
    ##########################################################################################
    ## Count hits
    dfdsum = dfd.groupby(['YEAR','MO']).sum()
    ## Count Reporting Days
    dfdcount = dfd.groupby(['YEAR','MO']).count()
    ## Count real days per month 
    dfdtotal = dshell.groupby(['YEAR','MO']).count()
    ## filter based on required % of days // if month dosn't meet standard it is set to np.nan
##    print(dfdsum)
##    print(dfdcount)
##    print(dfdtotal)
    dfdsumout=pd.DataFrame(data = np.where((dfdcount/dfdtotal<mo_r),np.nan,dfdsum),columns=dfdsum.columns,index=dfdsum.index)
    ## Normalize count of hits
    dfsumnorm = dfdsumout*(dfdtotal/dfdcount)
    
    return dfsumnorm


"""
##########################################################################################
############################ End of NoD_Normalize ########################################
##########################################################################################
"""

"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  SOLAR_components
###
### INPUTS:  data frame from SFC_TableLoader object with Rollhour, daynight,
###                 PRcalc and relhum ran 
###
### OUTPUTS:  all the solar components
###
### PURPOSE: a compilation of specific solar calculations
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 01Jan2019 SSgt Jones  - ported to python
###
###############################################################################
"""

def SOLAR_components(df):
    sfcdata = df
    def solarC(mo,day,hr):
        #THE JULIAN DAY IS CALCULATED BELOW.  THE JULIAN DAY IS THE NUMBER
        #OF DAYS ELAPSED SINCE THE DAY 1 JANUARY 4713 BC.  THE FORMULA
        #BELOW IS FROM THE 1987 ALMANAC FOR COMPUTERS, AND IS ACCURATE ONLY
        #FOR THE INTERVAL 1 JANUARY 1801 TO 31 DECEMBER 2099.
        dum = .5 * np.sign(197800 + mo - 190002.5)
        jgdate = (np.round(725926 - int(7 * (1978 + int((mo + 9)/12))/4) +
        int(275 * mo/9) + day + 1721013.5 - hr/24 - dum))
        #THE FOLLOWING SERIES OF EQUATIONS COMPUTES THE EARTH-SUN DISTANCE
        #IN ASTRONOMICAL UNITS.  THE FORMULA APPEARED IN THE 1978 VOLUME
        #OF THE ALMANAC FOR COMPUTERS.
        #_T IS THE TIME IN CENTURIES PAST 1900 AD
        #_G IS THE MEAN ANOMALY OF THE SUN
        t = (jgdate - 2415020.5)/36525
        g = (358.475833 + t*(35999.049750 - 0.000150*t))*.0174533
        #THE VALUE OF G PRODUCED ABOVE IS USED AS AN ANGLE IN THE FORMULA
        #FOR RHO THAT APPEARS BELOW.  THE VALUES OF G RESULTING FROM THE
        #ABOVE EQUATION ARE OUTSIDE OF THE LIMITS THAT THE CURRENT (AS OF
        #7 APRIL 1988) FORTRAN COMPILER CAN HANDLE.  THE FOLLOWING 2 LINES
        #CONVERT THE VALUE OF G TO LIE IN THE INTERVAL 0 TO 2*PI.
        g = g - int(g/6.283185)*6.283185
        if g < 0 :
            g = g + 6.283185
        #RHO IS A NUMBER RELATED TO THE CONCEPT OF AREAL DENSITY.  RHO IS
        #PROPORTIONAL TO THE AREA OF A SUN-CENTERED SPHERE WHOSE RADIUS IS
        #EQUAL TO THE EARTH-SUN DISTANCE (EXPRESSED IN ASTRONOMICAL UNITS)
        rho = 0.000084*t*math.cos(g) - 0.000140*math.cos(2*g) - 0.033503*math.cos(g) + 1.000421
        dist = math.sqrt(rho)
        #THE FOLLOWING PAIR OF EQUATIONS COMPUTES THE RADIATIVE POWER AT
        #THE DISTANCE CALCULATED ABOVE.  THE AVERAGE VALUE OF 1371.0
        #WATTS PER SQUARE METER FOR THE SOLAR CONSTANT WAS OBTAINED FROM
        #ERB-NIMBUS 7 MEASUREMENTS OVER THE PERIOD NOVEMBER 16, 1978
        #THROUGH OCTOBER 31, 1981, AT A DISTANCE OF 1 AU FROM THE SUN.
        solarC = 1371 / (dist**2)
        return solarC

    sfcdata['SOLARC'] = (np.vectorize(solarC,otypes=[np.float64])
    (sfcdata['MO'],sfcdata['DAY'],sfcdata['HR']))

    # Finds the True max value of last 6hrs on snowdeep for solins_albedo func-- Levi Jones**
    sfcdata['SNWDEEP1'] = pd.to_numeric(sfcdata.SNWDEEP1,errors='coerce').fillna(0)
##    def sndpmaxLAST6(df):
##        sddf = df[['YEAR','MO','DAY','HR','SNWDEEP1']]
##        sddf = sddf.groupby(['YEAR','MO','DAY','HR']).max()
##        sddf['SNDP_LAST6'] = sddf['SNWDEEP1'].rolling(min_periods=1,window=6).max()
##        sddf = sddf.reset_index()
##        sddf = sddf[['YEAR','MO','DAY','HR','SNDP_LAST6']]
##        df = df.merge(sddf, how = 'outer', on = ('YEAR','MO','DAY','HR'))
##        return df
##
##    sfcdata = sndpmaxLAST6(sfcdata)
    albedo = float(.35)

    def solins_ALBEDO(ww1,ww2,ww3,ww4,ww5,ww6,ww7,wwa1,wwa2,wwa3,sndp,albedo):#sndp_last6,
        tmpalb = None
        if (ww1 in range(50,60) )|( ww1 in range(80,83) )|( ww1 == 91 )|( ww1 == 92 ):
            tmpalb = 0
        elif (ww1 in range(30,36) )|( ww1 == 98 ):
            tmpalb = .24
        elif (ww1 in (87,89,90,96,99) ):
            tmpalb = .35
        elif (ww1 in (23,68,69,83,84,95,97) ):
            tmpalb = .5
        elif (ww1 in (26,93,94) ):
            tmpalb = .65
        elif (ww1 in range(36,40) ):
            tmpalb = .8
        elif (ww1 in range(70,80) )|( ww1 in (22,85,86) ):
            tmpalb = .85
        elif (ww2 in range(50,60) )|( ww2 in range(80,83))|( ww2 == 91 )|( ww2 == 92 ):
            tmpalb = 0
        elif (ww2 in range(30,36) )|( ww2 == 98 ):
            tmpalb = .24
        elif (ww2 in (87,89,90,96,99) ):
            tmpalb = .35
        elif (ww2 in (23,68,69,83,84,95,97) ):
            tmpalb = .5
        elif (ww2 in (26,93,94) ):
            tmpalb = .65
        elif (ww2 in range(36,40) ):
            tmpalb = .8
        elif (ww2 in range(70,80) )|( ww2 in (22,85,86)):
            tmpalb = .85
        elif (ww3 in range(50,60) )|( ww3 in range(80,83))|( ww3 == 91 )|( ww3 == 92 ):
            tmpalb = 0
        elif (ww3 in range(30,36) )|( ww3 == 98 ):
            tmpalb = .24
        elif (ww3 in (87,89,90,96,99) ):
            tmpalb = .35
        elif (ww3 in (23,68,69,83,84,95,97) ):
            tmpalb = .5
        elif (ww3 in (26,93,94) ):
            tmpalb = .65
        elif (ww3 in range(36,40) ):
            tmpalb = .8
        elif (ww3 in range(70,80) )|( ww3 in (22,85,86) ):
            tmpalb = .85
        elif (ww4 in range(50,60) )|( ww4 in range(80,83))|( ww4 == 91 )|( ww4 == 92 ):
            tmpalb = 0
        elif (ww4 in range(30,36) )|( ww4 == 98 ):
            tmpalb = .24
        elif (ww4 in (87,89,90,96,99) ):
            tmpalb = .35
        elif (ww4 in (23,68,69,83,84,95,97) ):
            tmpalb = .5
        elif (ww4 in (26,93,94) ):
            tmpalb = .65
        elif (ww4 in range(36,40) ):
            tmpalb = .8
        elif (ww4 in range(70,80) )|( ww4 in (22,85,86) ):
            tmpalb = .85
        elif (ww5 in range(50,60) )|( ww5 in range(80,83) )|( ww5 == 91 )|( ww5 == 92 ):
            tmpalb = 0
        elif (ww5 in range(30,36) )|( ww5 == 98 ):
            tmpalb = .24
        elif (ww5 in (87,89,90,96,99) ):
            tmpalb = .35
        elif (ww5 in (23,68,69,83,84,95,97) ):
            tmpalb = .5
        elif (ww5 in (26,93,94) ):
            tmpalb = .65
        elif (ww5 in range(36,40) ):
            tmpalb = .8
        elif (ww5 in range(70,80) )|( ww5 in (22,85,86) ):
            tmpalb = .85
        elif (ww6 in range(50,60) )|( ww6 in range(80,83) )|( ww6 == 91 )|( ww6 == 92 ):
            tmpalb = 0
        elif (ww6 in range(30,36) )|( ww6 == 98 ):
            tmpalb = .24
        elif (ww6 in (87,89,90,96,99) ):
            tmpalb = .35
        elif (ww6 in (23,68,69,83,84,95,97) ):
            tmpalb = .5
        elif (ww6 in (26,93,94) ):
            tmpalb = .65
        elif (ww6 in range(36,40) ):
            tmpalb = .8
        elif (ww6 in range(70,80) )|( ww6 in (22,85,86) ):
            tmpalb = .85
        elif (ww7 in range(50,60) )|( ww7 in range(80,83) )|( ww7 == 91 )|( ww7 == 92 ):
            tmpalb = 0
        elif (ww7 in range(30,36) )|( ww7 == 98 ):
            tmpalb = .24
        elif (ww7 in (87,89,90,96,99) ):
            tmpalb = .35
        elif (ww7 in (23,68,69,83,84,95,97) ):
            tmpalb = .5
        elif (ww7 in (26,93,94) ):
            tmpalb = .65
        elif (ww7 in range(36,40) ):
            tmpalb = .8
        elif (ww7 in range(70,80) )|( ww7 in (22,85,86) ):
            tmpalb = .85
        if pd.isnull(tmpalb) == True:
            if (wwa1 in (21,22,23,80,81,82,83,84,90,91,92,94,95) )|( wwa1 in range(40,69) ):
                tmpalb = 0
            elif (wwa1 in (93,96) ):
                tmpalb = .35
            elif (wwa1 in range(70,80) )|( wwa1 in (85,86,87) ):
                tmpalb = .85
            elif (wwa2 in (21,22,23,80,81,82,83,84,90,91,92,94,95) )|( wwa2 in range(40,69) ):
                tmpalb = 0
            elif (wwa2 in (93,96) ):
                tmpalb = .35
            elif (wwa2 in range(70,80) )|( wwa2 in (85,86,87) ):
                tmpalb = .85
            elif (wwa3 in (21,22,23,80,81,82,83,84,90,91,92,94,95) )|( wwa3 in range(40,69) ):
                tmpalb = 0
            elif (wwa3 in (93,96) ):
                tmpalb = .35
            elif (wwa3 in range(70,80) )|( wwa3 in (85,86,87) ):
                tmpalb = .85
        if sndp > 0 :
            if sndp > 3 :
                tmpalb = .85
##        elif sndp <= 0 :
##            if sndp_last6 :
##                tmpalb = .85
        if pd.isnull(tmpalb) == True:
            return albedo
        else :
            return tmpalb

    sfcdata['ALBEDO'] = (np.vectorize(solins_ALBEDO)(sfcdata['WW1'],
    sfcdata['WW2'],sfcdata['WW3'],sfcdata['WW4'],sfcdata['WW5'],
    sfcdata['WW6'],sfcdata['WW7'],sfcdata['WWA1'],sfcdata['WWA2'],
    sfcdata['WWA3'],sfcdata['SNWDEEP1'],albedo))#,sfcdata['SNDP_LAST6']


    def skyLYR(skylyr1,skylyr2,skylyr3,skylyr4, skybase1,skybase2,skybase3,
    skybase4, cloudtype1,cloudtype2,cloudtype3,cloudtype4, sumcov1,sumcov2,
    sumcov3,sumcov4, sumhgt1,sumhgt2,sumhgt3,sumhgt4):
        lo = 0
        mid = 0
        hi = 0
        cirus = 0
        obscr = 0

        funclist = [(skylyr1,skybase1,cloudtype1,sumcov1,sumhgt1),
                    (skylyr2,skybase2,cloudtype2,sumcov2,sumhgt2),
                    (skylyr3,skybase3,cloudtype3,sumcov3,sumhgt3),
                    (skylyr4,skybase4,cloudtype4,sumcov4,sumhgt4)]
        
        for i,j,k,l,m in funclist:
            if i > 0 :
                if i == 10 :
                    sl = 2
                elif i == 9 :
                    sl = 8
                else :
                    sl = i

                if (k >= 6) & (k <= 9) | (j <= 1800) & (j != 0):
                    lo = sl
                elif (k >= 3) & (k <= 5) | (j > 1800) & (j <= 5400):
                    mid = sl
                elif (k >= 0) & (k <= 2) | (j > 5400):
                    hi = sl
                    if k == 0 :
                        cirus = 1
                else :
                    obscr = 1

            if l > 0 :
                if l == 6 :
                    sl = 2
                elif l == 5 :
                    sl = 8
                else :
                    sl = l*2

                if (m <= 1800) & (m > 0):
                    lo = sl
                elif (m > 1800) & (m <= 5400) :
                    mid = sl
                elif (m > 5400) :
                    hi = sl
                else :
                    obscr = 1
        return lo,mid,hi,cirus,obscr

    sfcdata['lo8'],sfcdata['mid8'],sfcdata['hi8'],sfcdata['CIRUS'],sfcdata['OBSCR'] = (
    np.vectorize(skyLYR)(sfcdata['SKYLYR1'],sfcdata['SKYLYR2'],sfcdata['SKYLYR3'],
    sfcdata['SKYLYR4'],sfcdata['SKYBASE1'],sfcdata['SKYBASE2'],sfcdata['SKYBASE3'],
    sfcdata['SKYBASE4'],sfcdata['CLOUDTYPE1'],sfcdata['CLOUDTYPE2'],sfcdata['CLOUDTYPE3'],
    sfcdata['CLOUDTYPE4'],sfcdata['SUMCOV1'],sfcdata['SUMCOV2'],sfcdata['SUMCOV3'],
    sfcdata['SUMCOV4'],sfcdata['SUMHGT1'],sfcdata['SUMHGT2'],sfcdata['SUMHGT3'],sfcdata['SUMHGT4']))




    def skyLYR2(lo,mid,hi,sky,sun,cavok,ww1,ww2):
        skyWXerror = 0
        if ((lo == 0) | (mid == 0) |
        (hi== 0) | (sun == 1)) :
            # If CAVOK flag 'Y' aways assume 1/8 in hi //
            #Have removed need for user input on how to deal with
            #CAVOK per SSgt Levi Jones & Mr. Zautner Collaboration
            if (sky == 0) & (cavok == 'Y') :
                if hi == 0 :
                    hi = 1
            if (ww1 >= 50) | (ww2 >= 50):
                if sky == 0 :
                    #print(' Has WX code but no Sky Cover')
                    skyWXerror = 1
                else :
                    if lo == 0 :
                        lo = sky
                    if mid == 0 :
                        mid = sky
                    if hi == 0 :
                        hi = sky

            else :
                if (lo == 0) & (sky > 0) :
                    lo = sky

        if (skyWXerror == 0) :
            lo8 = lo/8.
            mid8 = mid/8.
            hi8 = hi/8.
        else:
            lo8 = 0
            mid8 = 0
            hi8 = 0
        return lo8,mid8,hi8

    sfcdata['lo8'],sfcdata['mid8'],sfcdata['hi8'] = (
    np.vectorize(skyLYR2)(sfcdata['lo8'],sfcdata['mid8'],sfcdata['hi8'],
    sfcdata['SKY'],sfcdata['SUN'],sfcdata['CAVOK'],sfcdata['WW1'],sfcdata['WW2']))


    def calcRAD(albedo,lo8,mid8,hi8,cirus,obscr,solarc,zenith):

        def wfunc(cldcvr,cldind,zenith):
        # PURPOSE: Computes the weight of the cloud portion of the atmospheric
        #          transmittance equations.
        #
        # METHOD:
        #        THE WEIGHTS COMPUTED IN THIS ROUTINE ARE DETERMINED FROM AN
        #    EXPRESSION OF SINE TERMS AND CLOUD LAYER FRACTIONS.  THE
        #    COEFFICIENTS OF THESE TERMS WERE EMPIRICALLY DETERMINED FROM THE
        #    SOLMET DATA BASE BY DR RALPH SHAPIRO OF STSC.  THESE COEFFICIENTS
        #    AND THE EQUATION FOR CALCULATING THEM ARE DOCUMENTED IN
        #    AFGL-TR-87-0200 (DTIC NUMBER ADB 114 709).
        #
        # REFERENCES:
        #    AFGL-TR-87-0200 (DTIC NUMBER ADB 114 709)
        # DATA DICTIONARY:
        #   &cldcvr    fractional cloud cover
        #   &cldind    cloud indicator flag - see calcrd for values
        #   &weight    cloud layer weighting factor
        #   &zenith    solar zenith angle measure from vertical in radians
        #   _cldcvr    fractional cloud cover
        #   _cldind    cloud indicator flag - see calcrd for values
            if (cldcvr > .99):
                weight = 1
            elif (cldind == 1):
                weight = cldcvr * (1.512 +
                (-1.176) * math.cos(zenith) +
                (-2.160) * cldcvr +
                ( 1.420) * cldcvr * math.cos(zenith) +
                (-0.032) * ((math.cos(zenith))**2) +
                ( 1.422) * cldcvr**2 )
            elif cldind == 2:
                weight = cldcvr * ( 1.429 +
                (-1.207) * math.cos(zenith) +
                (-2.008) * cldcvr +
                ( 0.853) * cldcvr * math.cos(zenith) +
                ( 0.324) * ((math.cos(zenith))**2) +
                ( 1.582) * cldcvr**2)
            elif cldind == 3:
                weight = cldcvr * ( 1.552 +
                (-1.957) * math.cos(zenith) +
                (-1.762) * cldcvr +
                ( 2.067) * cldcvr * math.cos(zenith) +
                ( 0.448) * ((math.cos(zenith))**2) +
                ( 0.932) * cldcvr**2)
            else :
                weight = cldcvr * ( 0.675 +
                (-3.432) * math.cos(zenith) +
                ( 1.929) * cldcvr +
                ( 0.842) * cldcvr * math.cos(zenith) +
                ( 2.693) * ((math.cos(zenith))**2) +
                (-1.354) * cldcvr**2);
            return weight

        def rfunc(cldind, layind, weight, zenith):
        # PURPOSE: Computes the amount of radiation reflected back toward a surfac
        #          by a layer.
        #
        # METHOD:
        #        THE REFLECTION OF RADIATION BY A LAYER IS CALCULATED BY
        #    COMPUTING THE WEIGHTED SUM OF TWO SINE POWER SERIES EQUATIONS.
        #    ONE OF THESE EQUATIONS REPRESENTS THE EFFECTS OF THE CLOUDINESS OF
        #    A LAYER UPON REFLECTION; THIS POWER SERIES IS MULTIPLIED BY THE
        #    WEIGHT COMPUTED FOR THE LAYER IN THE WFUNC FUNCTION.  THE OTHER
        #    EQUATION REPRESENTS THE EFFECT OF REFLECTION OF RADIATION BY A
        #    CLEAR LAYER; THIS POWER SERIES IS MULTIPLIED BY THE QUANTITY
        #    (1 - W), WHERE W IS THE WEIGHT FOR THE CLOUDY LAYER, COMPUTED IN
        #    FUNCTION WFUNC.
        #
        # REFERENCES:
        #    AFGL-TR-87-0200 (DTIC NUMBER ADB 114 709)
        # DATA DICTIONARY:
        #   &cldind    cloud indicator flag - see calcrd for values
        #   &layind    layer indicator code - see calcrd for values
        #   &reflec    reflected radiation
        #   &weight    cloud layer weighting factor
        #   &zenith    solar zenith angle measure from vertical in radians
        #   _acl0      coefficient for power equation
        #   _acl1      coefficient for power equation
        #   _acl2      coefficient for power equation
        #   _acl3      coefficient for power equation
        #   _aov0      coefficient for power equation
        #   _aov1      coefficient for power equation
        #   _aov2      coefficient for power equation
        #   _aov3      coefficient for power equation
        #   _cldind    cloud indicator flag - see calcrd for values
        #   _coszen    cosine of the zenith angle
        #   _layind    layer indicator code
            coszen = math.cos(zenith)
            if (cldind == 1) :
                aov0 =  0.69143
                aov1 = -0.14419
                aov2 = -0.05100
                aov3 =  0.06682
            elif (cldind == 2) :
                aov0 =  0.61394
                aov1 = -0.01469
                aov2 = -0.17400
                aov3 =  0.14215
            elif (cldind == 3) :
                aov0 =  0.42111
                aov1 = -0.04002
                aov2 = -0.51833
                aov3 =  0.40540
            else:
                aov0 =  0.25674
                aov1 = -0.18077
                aov2 = -0.21961
                aov3 =  0.25272

            if (layind == 1) :
                acl0 =  0.27436
                acl1 = -0.43132
                acl2 =  0.26920
                acl3 = -0.00447
            elif (layind == 2) :
                acl0 =  0.15946
                acl1 = -0.42185
                acl2 =  0.48800
                acl3 = -0.18493
            elif (layind == 3) :
                acl0 =  0.15325
                acl1 = -0.39620
                acl2 =  0.42095
                acl3 = -0.14200
            else :
                acl0 =  0.12395
                acl1 = -0.34765
                acl2 =  0.39478
                acl3 = -0.14627
        #    THE EQUATION BELOW COMPUTES THE FRACTION OF RADIATION REFLECTED
        #    BY A LAYER BY TAKING THE WEIGHTED SUM OF TWO COSINE POWER SERIES
        #    EQUATIONS.  BOTH EQUATIONS ARE OF THE SAME FORM, WITH THE
        #    DIFFERENCE BETWEEN THE TWO BEING A RESULT OF THE COEFFICIENTS USED
        #    AND THE WIEGHTING FACTOR BY WHICH THEY ARE MULTIPLIED.  THESE
        #    EQUATIONS APPEAR ON PG 6 OF AFGL-TR-87-0200.
        #    THE POWER SERIES EQUATION USED IS OF THE FORM
        #        REFLEC = C0 + C1*COS(ZEN) + C2*(COS(ZEN)**2) + C3*(COS(ZEN)**3)
        #    WHERE C0-C3 REPRESENT THE COEFFICIENTS FOR THE LAYER (AOV0-AOV3 OR
        #    ACL0-ACL3) AND ZEN REPRESENTS THE SOLAR ZENITH ANGLE.
        #    THE TOTAL REFLECTION OF RADIATION BY THE LAYER IS THE SUM
        #    REFLEC = WEIGHT*(CLOUDY REFLEC) + (1 - WEIGHT)*(CLEAR REFLEC)
        #    WHERE WEIGHT IS THE WEIGHT ASSIGNED TO THE CLOUD COMPONENT OF
        #    REFLECTION, CLOUDY REFLEC IS THE REFLECTION CALCULATED BY USING
        #    THE COEFFICIENTS AOV0-AOV3 IN THE POWER SERIES EQUATION, AND CLEAR
        #    REFLEC IS THE REFLECTION CALCULATED BY USING THE COEFFICIENTS
        #    ACL0-ACL3 IN THE POWER SERIES EQUATION.
        #    REMEMBER THAT THE REFLECTION OF THE RADIATION IS ASSUMED BY
        #    SHAPIRO TO BE EQUAL TO THE DIFFUSE COMPONENT OF RADIATION
        #    TRANSMITTED THROUGH THE LAYER FOR WHICH THE REFLECTION IS
        #    CALCULATED.
        #    THE CLOUDY COMPONENT LAYER COMPONENT OF THE DIFFUSE RADIATION
        #    APPEARS FIRST IN THE EQUATION BELOW.
            rflec = (weight * (aov0 + aov1 * coszen +
                aov2 * coszen * coszen +
                aov3 * coszen * coszen * coszen) +
                (1 - weight) * (acl0 + acl1 * coszen +
                acl2 * coszen * coszen +
                acl3 * coszen * coszen * coszen))
            return rflec

        def tfunc(cldind, layind, weight, zenith):
        # METHOD:
        #        THE TRANSMITTANCE OF RADIATION THROUGHT A LAYER IS CALCULATED
        #    BY COMPUTING THE WEIGHTED SUM OF TWO SINE POWER SERIES EQUATIONS.
        #    ONE OF THESE EQUATIONS REPRESENTS THE EFFECTS OF THE CLOUDINESS OF
        #    A LAYER UPON TRANSMITTANCE; THIS POWER SERIES IS MULTIPLIED BY
        #    THE WEIGHT COMPUTED FOR THE LAYER IN THE WFUNC FUNCTION.  THE
        #    OTHER EQUATION REPRESENTS THE EFFECT OF TRANSMITTANCE OF RADIATION
        #    THROUGH A CLEAR LAYER; THIS POWER SERIES IS MULTIPLIED BY THE
        #    QUANTITY (1 - W), WHERE W IS THE WEIGHT FOR THE CLOUDY LAYER,
        #    COMPUTED IN FUNCTION WFUNC.
        #
        # REFERENCES:
        #    AFGL-TR-87-0200 (DTIC NUMBER ADB 114 709)
        #
        # DATA DICTIONARY:
        #   &cldind    cloud indicator flag - see calcrd for values
        #   &layind    layer indicator code - see calcrd for values
        #   &weight    cloud layer weighting factor
        #   &trans     transmittance of radiation through the layer
        #   &zenith    solar zenith angle measure from vertical in radians
        #   _bcl0      coefficient for power equation
        #   _bcl1      coefficient for power equation
        #   _bcl2      coefficient for power equation
        #   _bcl3      coefficient for power equation
        #   _bov0      coefficient for power equation
        #   _bov1      coefficient for power equation
        #   _bov2      coefficient for power equation
        #   _bov3      coefficient for power equation
        #   _cldind    cloud indicator flag - see calcrd for values
        #   _coszen    cosine of the zenith angle
        #   _layind    layer indicator code - see calcrd for values
            coszen = math.cos(zenith)
            if(cldind==1):
                bov0=0.15785
                bov1=0.32410
                bov2=-0.14458
                bov3=0.01457
            elif(cldind==2):
                bov0=0.23865
                bov1=0.20143
                bov2=-0.01183
                bov3=-0.07892
            elif(cldind==3):
                bov0=0.43562
                bov1=0.26094
                bov2=0.36428
                bov3=-0.38556
            else:
                bov0=0.63547
                bov1=0.35229
                bov2=0.08709
                bov3=-0.22902

            if(layind==1):
                bcl0=0.55336
                bcl1=0.61511
                bcl2=-0.29816
                bcl3=-0.06663
            elif(layind==2):
                bcl0=0.68679
                bcl1=0.71012
                bcl2=-0.71463
                bcl3=0.22339
            elif(layind==3):
                bcl0=0.69318
                bcl1=0.68227
                bcl2=-0.64289
                bcl3=0.17910
            else:
                bcl0=0.76977
                bcl1=0.49407
                bcl2=-0.44647
                bcl3=0.11558
        #    THE EQUATION BELOW COMPUTES THE FRACTION OF RADIATION TRANSMITTED
        #    THROUGH A LAYER BY TAKING THE WEIGHTED SUM OF TWO COSINE POWER
        #    SERIES EQUATIONS.  BOTH EQUATIONS ARE OF THE SAME FORM, WITH THE
        #    DIFFERENCE BETWEEN THE TWO BEING A RESULT OF THE COEFFICIENTS USED
        #    AND THE WIEGHTING FACTOR BY WHICH THEY ARE MULTIPLIED.  THESE
        #    EQUATIONS APPEAR ON PG 6 OF AFGL-TR-87-0200.
        #    THE POWER SERIES EQUATION USED IS OF THE FORM
        #      TRANS = C0 + C1*COS(ZEN) + C2*(COS(ZEN)**2) + C3*(COS(ZEN)**3)
        #    WHERE C0-C3 REPRESENT THE COEFFICIENTS FOR THE LAYER (BOV0-BOV3 OR
        #    BCL0-BCL3) AND ZEN REPRESENTS THE SOLAR ZENITH ANGLE.
        #    THE TOTAL TRANSMITTANCE OF RADIATION THROUGH THE LAYER IS THE SUM
        #    TRANS = WEIGHT*(CLOUDY TRANS) + (1 - WEIGHT)*(CLEAR TRANS)
        #    WHERE WEIGHT IS THE WEIGHT ASSIGNED TO THE CLOUD COMPONENT OF
        #    DIRECT TRANSMITTANCE, CLOUDY TRANS IS THE TRANSMITTANCE CALCULATED
        #    BY USING THE COEFFICIENTS BOV0-BOV3 IN THE POWER SERIES EQUATION,
        #    AND CLEAR TRANS IS THE TRANSMITTANCE CALCULATED BY USING THE
        #    COEFFICIENTS BCL0-BLC3 IN THE POWER SERIES EQUATION.
        #    THE CLOUDY COMPONENT LAYER COMPONENT OF THE RADIATION APPEARS
        #    FIRST IN THE EQUATION BELOW.
            trans = weight * (bov0 + bov1 * coszen +
                    bov2 * coszen * coszen +
                    bov3 * coszen * coszen * coszen) + \
                    (1 - weight) * (bcl0 + bcl1 * coszen +
                    bcl2 * coszen * coszen +
                    bcl3 * coszen * coszen * coszen)
            return trans

        # METHOD:
        #        THIS ROUTINE IS BASED UPON A SERIES OF EQUATIONS PRINTED IN
        #    AFGL-TR-87-0200 (DTIC NUMBER ADB 114 709). THESE EQUATIONS MODEL
        #    THE TRANSMITTANCE OF SOLAR RADIATION THROUGH THE ATMOSPHERE AS
        #    TRANSMITTANCE THROUGH A NUMBER OF PLANE PARALLEL LAYERS IN THE
        #    ATMOSPHERE, EACH WITH A CONSTANT TRANSMISSIVITY AND REFLECTANCE.
        #    THE RADIATION TRANSMITTED THROUGH ANY LAYER IS A SUM OF THE DIRECT
        #    AND THE DIFFUSE COMPONENTS OF THE RADIATION PASSING THROUGH THE
        #    LAYER.  THE TRANSMISSIVITY OF SOLAR RADIATION THROUGH A LAYER AS
        #    DIFFUSE RADIATION IS ASSUMED TO BE EQUAL TO THE REFLECTIVITY OF
        #    THE LAYER, SO CALCULATION OF THE REFLECTIVITY OF THE LAYER
        #    PROVIDES THE DIFFUSE TRANSMISSIVITY OF THE LAYER.
        #
        #        ONCE THE TRANSMISSIVITIES FOR THE THREE LAYERS HAVE BEEN
        #    DETERMINED (THROUGH CALLS TO WFUNC, RFUNC, AND TFUNC), THE TOTAL
        #    RADIATION AT THE SURFACE AND THE DIRECT RADIATION AT THE SURFACE
        #    ARE CALCULATED.  THE DIRECT RADIATION IS CALCULATED BY MULTIPLYING
        #    THE PRODUCT OF THE TRANSMISSIVITIES OF DIRECT RADIATION FOR EACH
        #    LAYER BY THE RADIATION INCIDENT ON THE TOP OF THE ATMSOPHERE.  THE
        #    TOTAL RADIATION CALCULATION IS MORE COMPLEX, SINCE THE TOTAL
        #    RADIATION PASSING THROUGH THE MIDDLE LAYER IS THE SUM OF THE
        #    DIRECT RADIATION AND THE DIFFUSE RADIATION PASSING THROUGH THE
        #    TOP LAYER, MULTIPLIED BY THE DIFFUSE TRANSMISSIVITY FACTOR FOR THE
        #    MIDDLE LAYER, PLUS THE DIRECT RADIATION TRANSMITTED THROUGH THE
        #    TOP LAYER MULTIPLIED BY THE TRANSMISSIVITY FOR DIRECT RADIATION
        #    THROUGH THE MIDDLE LAYER:
        #
        #    TOTAL = DIFFUSE2#(DIRECT3 + DIFFUSE3) + DIRECT2#DIRECT3.
        #
        #    THE EXPRESSION FOR THE TOTAL RADIATION THROUGHT THE SURFACE LAYER
        #    IS CORRESPONDINGLY MORE COMPLEX.  TO SIMPLIFY THE REPRESENTATION
        #    OF THIS TRANSMITTANCE, IT IS REPRESENTED BY THE EQUATION:
        #
        #    TOTAL = TOTAL1#TOTAL2#TOTAL3#SOLCON/CORRECTION
        #
        #    WHERE TOTAL IS THE TOTAL RADIATION ARRIVING AT THE SURFACE, TOTAL1
        #    IS THE FRACTION OF RADIATION ALLOWED TO PASS THROUGH THE SURFACE
        #    SURFACE LAYER, TOTAL2 IS THE FRACTION OF ALL RADIATION ALLOWED TO
        #    PASS THROUGH THE MIDDLE LAYER, TOTAL3 IS THE FRACTION OF ALL
        #    RADIATION ALLOWED TO PASS THROUGH THE TOP LAYER, SOLCON IS THE
        #    SOLAR CONSTANT AT THE TIME IN QUESTION, AND CORRECTION IS A
        #    CALCULATED TERM USED TO CORRECT THE EQUATION TO PROVIDE THE VALUE
        #    CALCULATED THROUGH THE ACTUAL ALGEBRAIC EXPRESSION FOR THE TOTAL
        #    RADIATION PASSING THROUGH THE SURFACE LAYER.  THE DIFFUSE
        #    RADIATION AT THE SURFACE IS CALCULATED BY TAKING THE DIFFERENCE
        #    BETWEEN THE TOTAL RADIATION AT THE SURFACE AND THE DIRECT
        #    RADIATION AT THE SURFACE.
        #
        # REFERENCES:
        #    AFGL-TR-87-0200 (DTIC NUMBER ADB 114 709)
        #
        # MACROS/SUBROUTINES CALLED:
        #    rfunc   (INTERNAL)  computes radiation reflected back toward a surface
        #    tfunc   (INTERNAL)  computes radiation that passes through a layer
        #    wfunc   (INTERNAL)  computes weight of cloud portion
        #
        # FILES/DATABASE/TABLES ACCESSED:
        #   none
        #
        # DATA DICTIONARY:
        #   &albedo      albedo
        #   &cirus       cirrus flag 0=no or thin, 1=thick cirrus
        #   &difrad      diffuse radiation W/M2
        #   &dirrad      dirret radiation W/M2
        #   &hi8         fraction of high clouds present
        #   &lo8         fraction of low clouds present
        #   &mid8        fraction of mid clouds present
        #   &obscr       obscuration flag 0=no obscuration, 1=obscurations
        #   &solarc      solar constant
        #   &surrad      surface reflect radiation W/M2
        #   &totrad      total radiation W/M2
        #   &zenith      solar zenith angle in radians
        #   _albedo      albedo
        #   _bigd2       corresponds to reference p3 equation, variable D2
        #   _cirus       cirrus flag 0=no or thin, 1=thick cirrus
        #   _difrad      diffuse radiation W/M2
        #   _dirrad      dirret radiation W/M2
        #   _dtrans1     low layer transmittance difference from reflected
        #   _dtrans2     mid layer transmittance difference from reflected
        #   _dtrans3     high layer transmittance difference from reflected
        #   _hi8         fraction of high clouds present
        #   _hicvr       flag specifying high clouds, 3 with cirrus, 4 without or thin
        #   _hilay       flag specifying high layer, always 4
        #   _lild1       corresponds to reference p3 equation, variable d3
        #   _lild2       corresponds to reference p3 equation, variable d2
        #   _lild3       corresponds to reference p3 equation, variable d1
        #   _lo8         fraction of low clouds present
        #   _locvr       flag specifying low clouds, always 1
        #   _lowlay      flag specifying lower layer 1 with obscurants, 2 without
        #   _mid8        fraction of mid clouds present
        #   _midcvr      flag specifying middle clouds, always 2
        #   _midlay      flag specifying middle layer, always 3
        #   _obscr       obscuration flag 0=no obscuration, 1=obscurations
        #   _rflec0      albedo
        #   _rflec1      low layer reflectivity
        #   _rflec2      mid layer reflectivity
        #   _rflec3      high layer reflectivity
        #   _solarc      solar constant
        #   _solrad      corrected solar constant
        #   _surrad      surface reflect radiation W/M2
        #   _totrad      total radiation W/M2
        #   _trans1      low layer transmittance
        #   _trans2      mid layer transmittance
        #   _trans3      high layer transmittance
        #   _weght1      weight of low cloud presence
        #   _weght2      weight of mid cloud presence
        #   _weght3      weight of high cloud presence
        #   _zenith      solar zenith angle in radians

        if ((albedo<0) | (albedo>1) | (lo8<0) | (lo8>1) | (mid8<0) | (mid8>1) |
        (hi8<0) | (hi8>1) | (solarc<1300) | (solarc>1450) | (zenith<0) | (zenith>1.571)):
            difrad=np.nan
            dirrad=np.nan
            surrad=np.nan
            totrad=np.nan
        else:
        # Correct the solar power to represent the incoming radiation
        #  received on a horizontal surface
            solrad = solarc * math.cos(zenith)
        #  DETERMINE THE VALUES OF LOWLAY AND HICVR, ON THE BASIS OF
        #  THE PRESENCE OR ABSENCE OF SURFACE OBSCURATION AND THE
        #  THICKNESS OF ANY CIRRUS CLOUDS THAT MAY BE PRESENT.
            if (obscr == 1) :
                lowlay = 1
            else:
                lowlay = 2

            midlay = 3
            hilay  = 4

            locvr  = 1
            midcvr = 2

            if (cirus == 1) :
                hicvr = 3
            else :
                hicvr = 4

            rflec0 = albedo
        #  determine the weight for each of the three cloud layers
        #  determine the reflectivity for each layer
        #  determine the transmissivity for each layer
            weght1 = wfunc(lo8 , locvr , zenith)
            weght2 = wfunc(mid8, midcvr, zenith)
            weght3 = wfunc(hi8 , hicvr , zenith)
            rflec1 = rfunc(locvr , lowlay, weght1, zenith)
            rflec2 = rfunc(midcvr, midlay, weght2, zenith)
            rflec3 = rfunc(hicvr , hilay , weght3, zenith)
            trans1 = tfunc(locvr , lowlay, weght1, zenith)
            trans2 = tfunc(midcvr, midlay, weght2, zenith)
            trans3 = tfunc(hicvr , hilay , weght3, zenith)
        #    THE AMOUNT OF RADIATION TRANSMITTED THROUGH A LAYER AS DIRECT
        #    RADIATION IS EQUAL TO THE DIFFERENCE OF THE TOTAL RADIATION
        #    TRANSMITTED THROUGH THE LAYER AND THE AMOUNT OF RADIATION
        #    TRANSMITTED THROUGH THE LAYER AS DIFFUSE RADIATION (WHICH IS
        #    ASSUMED TO EQUAL THE REFLECTIVITY OF THE LAYER).  IF THIS
        #    DIFFERENCE IS LESS THAN ZERO, THEN THE DIRECT RADIATION IS SET TO
        #    ZERO, MEANING THAT ALL RADIATION TRANSMITTED THROUGH THE LAYER IS
        #    TRANSMITTED AS DIFFUSE RADIATION.
            dtrans1 = max(0,trans1 - rflec1)
            dtrans2 = max(0,trans2 - rflec2)
            dtrans3 = max(0,trans3 - rflec3)
        #    THE TOTAL RADIATION ARRIVING AT THE SURFACE IS CALCULATED BY A
        #    COMPLEX EQUATION INITIALLY REFERRED TO IN THE EMBEDDED USERS
        #    MANUAL.  THIS EQUATION IS EXPRESSED BELOW AS THE PRODUCT OF THE
        #    TOTAL RADIATION TRANSMITTED THROUGH EACH LEVEL MULTIPLIED BY THE
        #    INCIDENT RADIATION PER UNIT AREA AT THE TOP OF THE ATMOSPHERE,
        #    ALL DIVIDED BY AN ALGEBRAIC CORRECTION FACTOR.  THE CORRECTION
        #    () IS CALCULATED BELOW.       THE SOURCE FOR THESE EQUATIONS IS
        #    P3 OF AFGL-TR-87-0200.
        #
        #    THE LABELS LILD1-LILD3 REFER TO THE LOWERCASE VALUES D1-D3, WHILE
        #    BIGD2 REFERS TO THE UPPERCASE VALUE D2.  WHEN REFERRING TO THE
        #    TECH REPORT, NOTE THAT LILD1 AND LILD3 ARE TRANSPOSED HERE TO
        #    MAINTAIN CONSISTENCY WITH A MODEL ATMOSPHERE WITH THE LOW LAYER
        #    BEING LAYER ONE AND THE HIGH LAYER BEING LAYER THREE.
            lild1 = 1 - rflec1 * rflec0
            lild2 = 1 - rflec2 * rflec1
            lild3 = 1 - rflec3 * rflec2
            bigd2 = (lild1 * (lild3 * lild2 - rflec3 * rflec1 * trans2 * trans2) -
            lild3 * rflec2 * rflec0 * trans1 * trans1 - rflec3 * rflec0 * trans2 *
            trans2 * trans1 * trans1 )
        #    THE TOTAL RADIATION ARRIVING AT THE SURFACE IS THE PRODUCT OF THE
        #    DIRECT RADIATION TRANSMITTED THROUGH EACH LEVEL MULTIPLIED BY THE
        #    INCIDENT RADIATION PER UNIT AREA AT THE TOP OF THE ATMOSPHERE.
        #    THE DIFFUSE RADIATION ARRIVING AT THE SURFACE IS THE DIFFERENCE
        #    BETWEEN THE TOTAL RADIATION ARRIVING AT THE SURFACE & THE DIRECT
        #    RADIATION ARRIVING AT THE SURFACE.  THE SURFACE REFLECTED
        #    RADIATION IS EQUAL TO THE SUM OF THE DIFFUSE RADIATION AND THE
        #    PRODUCT OF THE DIRECT RADIATION AND THE COSINE OF THE ZENITH
        #    ANGLE, MULTIPLIED BY THE ALBEDO OF THE SURFACE.
            dirrad = dtrans1 * dtrans2 * dtrans3 * solrad
            totrad = trans1*trans2*trans3*solrad/bigd2
            difrad = totrad - dirrad
            surrad = (dirrad*math.cos(zenith) + difrad) * albedo
            if pd.isnull(dirrad) == True:
                dirrad = np.nan
            if pd.isnull(totrad) == True:
                totrad = np.nan
            if pd.isnull(difrad) == True:
                difrad = np.nan
            if pd.isnull(surrad) == True:
                surrad = np.nan
        return dirrad,totrad,difrad,surrad

    sfcdata['DIRRAD'],sfcdata['TOTRAD'],sfcdata['DIFRAD'],sfcdata['SURRAD'] = (
    np.vectorize(calcRAD)(sfcdata['ALBEDO'],sfcdata['lo8'],sfcdata['mid8'],
    sfcdata['hi8'],sfcdata['CIRUS'],sfcdata['OBSCR'],sfcdata['SOLARC'],sfcdata['ZENITH']))

    return sfcdata

"""
### END OF PYTHON FUNCTION:  SOLAR_compenents
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  WETBULBglobetemp
###
### PROGRAM NAME: wbgt
###
### TRACKING NUMBER: 34749
###
### Inputs:  dataframe with all of he solar_componants
###
### PURPOSE: Calculate wet bulb globe temperature.
###
### METHOD: WBGT = 0.7Tw + 0.2Tg + 0.1Td
###
### REFERENCES:
###    none
###
### MACROS/SUBROUTINES CALLED:
###    none
###
### DATA DICTIONARY:
###    &wtbulb      Natural wet-bulb temperature (humidity indicator)
###    &blackglobe  Black globe temperature
###    &wbgt        Wet bulb globe temperature
###    &tempk       Temperature (kelvin)
###
### REMARKS: none
###
### UPDATES:
###    29Sept08    Mr. Kiess     New
###
###
### PROGRAM NAME: wtbulb  (Wet Bulb Temperature)
###
### BCR TRACKING NUMBER:
###
### PURPOSE:   Computes the wet bulb temperature from dry bulb and dew
###            point temperatures.
###
### METHOD:
###    If dry bulb temperature, dewpoint temperature and station pressure
###    are present, the wet bulb temperature is calculated iteratively
###    using a bisecting method, which begins with an initial estimate at
###    the wet bulb temperature--the mean of the dry bulb and dew point
###    temperature. The method then uses that value to calculate a wet bulb
###    temperature so that the new estimate is the mean between the old
###    estimate and the wet bulb temperature.  The process is continued
###    until the difference between the new and the old estimate is less
###    than 0.05k or the number of iterations exceeds 50. If 50 iterations
###    is exceeded then a value of missing is returned for wet bulb.
###
### REFERENCES:
###    none
###
### PROCESS NARRATIVE LOCATION:
###
### MACROS/SUBROUTINES CALLED:
###   k2c       converts Kelvin to Celsius
###   k2f       converts Kelvin to Farhenheit
###
### FILES/DATABASE/TABLES ACCESSED:
###   none
###
### DATA DICTIONARY:
###   &dewp       Input: dew point temperature in &tunit units
###   &i          iteration counter of bisecting guesses
###   &iter       iteration counter
###   &j          iteration counter
###   &stnpr      Input: station pressure in mb
###   &temp       Input: dry bulb temperature in &tunit units
###   &tunit      Input: temperature unit designation (C, F, K, or R)
###   &wtbulb     Output: wet bulb temperature in &tunit units
###   _dewp       dew point temperature in Kelvin
###   _dewsat     saturation mixing ratio of dew point temperature
###   _errmin     smallest computed error of wet bulb calculation
###   _estmt      saturation vapor pressure in mb
###   _imin       index pointing to guess value of least error
###   _rcona      constant "a" in saturation vapor pressure computation
###   _rconb      constant "b" in saturation vapor pressure computation
###   _stnpr      station pressure in mb
###   _stop       flag set to 1 when convergence criteria met, otherwise 0
###   _t1-5       narrower array of guesses at wet bulb in Kelvin
###   _temp       dry bulb temperature in Kelvin
###   _tprime1-5  array of guesses at wet bulb in Kelvin
###   _twerr2-4   computed error of wet bulb calculation
###   _wetsat     saturation mixing ration of wet bulb guess
###   _wtbulb     wet bulb temperature in Kelvin
###
### REMARKS:
###   derived from the FORTRAN routine WTBULB
###
### UPDATES:
###   04MAR1999  Mr. Kiess/DOC2.  New
###
###   20JUL1999  TSgt Henderson/SCS - removed "go to" block from end of macro
###              due to problem with multiple calls to the block with reuse of
###              the block name.  Instead put conversion to input unit of
###              measurement after all calculations have been done.
###
###   27JUN2001  Mr. Kiess/DOC2 - Previous "fix" did not work correctly as it
###              removed the jump out of the converging criteria.  Recoded
###              with a do while loop with a flag variable.
###
###   01Jan2019 SSgt Jones  - ported to python
###
###############################################################################
"""

def WETBULBglobetemp(df):
    sfcdata = df
    def wtBULB(temp, dewp, stnpr):
        temp = temp + 273.15
        dewp = dewp + 273.15
        wb_dict = {}
        if (temp >= 216) & (temp <= 334) & (dewp > 216) & (dewp < 334) & (pd.isnull(stnpr) == False):
            if (abs(temp-dewp) < .05) :
                wtbulb = temp - 273.15
                return wtbulb
                
            else:
                wb_dict['tprime1'] = dewp
                wb_dict['tprime5'] = temp
                wb_dict['tprime3'] = wb_dict['tprime1'] + ((wb_dict['tprime5'] - wb_dict['tprime1']) / 2.)
                wb_dict['tprime2'] = wb_dict['tprime1'] + ((wb_dict['tprime3'] - wb_dict['tprime1']) / 2.)
                wb_dict['tprime4'] = wb_dict['tprime3'] + ((wb_dict['tprime5'] - wb_dict['tprime3']) / 2.)
                if (dewp <= 273.16) :
                    rcona = 21.874
                    rconb =  7.66
                else :
                    rcona = 17.269
                    rconb = 35.86
                estmt = 6.11 * math.exp((rcona * (dewp - 273.16)) / (dewp - rconb))
                dewsat = .622 * estmt / stnpr

                for i in range(1,51):
                    errmin = 99999
                    imin = 0

                    for x in range(2,5):
                        if wb_dict['tprime%d'%x] <= 273.16:
                            rcona = 21.874
                            rconb =  7.66
                        else :
                            rcona = 17.269
                            rconb = 35.86
                        estmt = (6.11 * math.exp((rcona * (wb_dict['tprime%d'%x] - 273.16)) /
                        (wb_dict['tprime%d'%x] - rconb)))
                        wetsat = .622 * estmt / stnpr
                        wb_dict['twerr%d'%x] = (abs((wetsat - dewsat) *
                        (2500 / (1.00464 + (dewsat * 1.84603))) - temp + wb_dict['tprime%d'%x]))

                        if wb_dict['twerr%d'%x] < errmin :
                            #print (str(wb_dict['twerr%d'%x]))
                            errmin = wb_dict['twerr%d'%x]
                            if x == 2 :
                                wb_dict['t1'] = wb_dict['tprime1']
                                wb_dict['t3'] = wb_dict['tprime2']
                                wb_dict['t5'] = wb_dict['tprime3']
                                imin = 2
                            elif x == 3 :
                                wb_dict['t1'] = wb_dict['tprime2']
                                wb_dict['t3'] = wb_dict['tprime3']
                                wb_dict['t5'] = wb_dict['tprime4']
                                imin = 3
                            else :
                                wb_dict['t1'] = wb_dict['tprime3']
                                wb_dict['t3'] = wb_dict['tprime4']
                                wb_dict['t5'] = wb_dict['tprime5']
                                imin = 4

                    if errmin > 0.05 :
                        wb_dict['t2'] = (wb_dict['t1'] + wb_dict['t3'])/2.
                        wb_dict['t4'] = (wb_dict['t3'] + wb_dict['t5'])/2.
                        if i == 50 :
                            print('wtbulb cant be found!' + str(errmin))
                        else:
                            for j in range(1,6):
                                wb_dict['tprime%d'%j] = wb_dict['t%d'%j]

                    else :
                        if imin == 2:
                            if wb_dict['tprime2'] > temp :
                                wtbulb = temp
                                return wtbulb
                                break
                            elif wb_dict['tprime2'] < dewp :
                                wtbulb = dewp
                                return wtbulb
                                break
                            else :
                                wtbulb = wb_dict['tprime2'] - 273.15
                                return wtbulb
                                break
                        elif imin == 3:
                            if wb_dict['tprime3'] > temp :
                                wtbulb = temp
                                return wtbulb
                                break
                            elif wb_dict['tprime3'] < dewp :
                                wtbulb = dewp
                                return wtbulb
                                break
                            else :
                                wtbulb = wb_dict['tprime3'] - 273.15
                                return wtbulb
                                break
                        else:
                            if wb_dict['tprime4'] > temp :
                                wtbulb = temp
                                return wtbulb
                                break
                            elif wb_dict['tprime4'] < dewp :
                                wtbulb = dewp
                                return wtbulb
                                break
                            else :
                                wtbulb = wb_dict['tprime4'] - 273.15
                                return wtbulb
                                break

    sfcdata['WTBULBC'] = (np.vectorize(wtBULB,otypes=[np.float64])
    (sfcdata['TEMPC'],sfcdata['DEWPC'],sfcdata['STAPRS']))

    def blackGLOBE(totrad, temp, wspd):

        temp = temp + 273.15
        if totrad > 0 :
            mnradiant_temp = ((totrad/.000000053865) + (temp**4))**0.25
        else:
            mnradiant_temp = temp

        blackglobe = temp
        windvel = math.sqrt(max(1.5,wspd))
        #print(blackglobe)
        for i in range(1,151):
            rad_heat_gain = 6.01 * (mnradiant_temp - blackglobe)
            conv_heat_loss = 6.32 * (0.15**(-0.4)) * windvel * (blackglobe - temp)
            if conv_heat_loss > rad_heat_gain :
                return blackglobe - 273.15
                break
            else:
                blackglobe = blackglobe + 0.5
                #print(str(conv_heat_loss) + '\n' + str(rad_heat_gain))
                if i == 150 :
                    blackglobe == np.nan
                    return blackglobe
                    #print(blackglobe)

    sfcdata['BLACKGLOBEC'] = (np.vectorize(blackGLOBE,otypes=[np.float64])
    (sfcdata['TOTRAD'],sfcdata['TEMPC'],sfcdata['WSPDMPS']))

    def wbgtFUNC(wtbulb, blackglobe, temp):
        # Too Kelvin
        wtbulb = wtbulb + 273.15
        blackglobe = blackglobe + 273.15
        temp = temp + 273.15

        wbgt = 0.7*wtbulb + 0.2*blackglobe + 0.1*temp
        wbgt = wbgt - 273.15
        return  wbgt

    sfcdata['WBGTC'] = (np.vectorize(wbgtFUNC,otypes=[np.float64])
    (sfcdata['WTBULBC'],sfcdata['BLACKGLOBEC'],sfcdata['TEMPC']))

    sfcdata['WTBULBF'] = ((sfcdata['WTBULBC']* 9/5) + 32).round(9)
    sfcdata['BLACKGLOBEF'] = ((sfcdata['BLACKGLOBEC']* 9/5) + 32).round(9)
    sfcdata['WBGTF'] = ((sfcdata['WBGTC']* 9/5) + 32).round(9)

    return sfcdata

"""
### END OF PYTHON FUNCTION:  WETBULBglobetemp
###
###############################################################################
###############################################################################
###############################################################################
"""


"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  UTCIfunc - Universal Thermal Climate Index
###
### Inputs:
###   tempc, dewpc, black_globec, windmps
###   
###
### Outputs:
###   utci 
###
### PURPOSE: An easy to use function that correctly handles the
###           Freqency calculation it's also able to bin hours
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 1 Jan 2020 SSgt Jones - added for a SAR
###
###############################################################################
"""

def UTCIfunc(tempc,dewpc,black_globec,windmps):
    tempk = tempc + 273.15
    dewpk = dewpc + 273.15

    if dewpk > 100 :
        if tempk >= 273.16 :
            vapor_e = 10**(23.832241-(5.02808*math.log10(dewpk))-
            (1.3816*(10**(-7))*(10**(11.334-(0.0303998*dewpk))))+
            (8.1328*(10**(-3))*(10**(3.49149-(1302.8844/dewpk))))-(2949.076/dewpk))
        else :
            vapor_e = 10**((3.56654*math.log10(dewpk))-
            (0.0032098*dewpk)-(2484.956/dewpk)+2.0702294)
    else :
        vapor_e = np.nan
        
    Ta = tempc 
    Pa = vapor_e/10
    D_Tmrt = black_globec - Ta
    va = windmps
    
    if D_Tmrt > 70:
        D_Tmrt = 70
    if D_Tmrt < -30:
        D_Tmrt = -30
    
    UTCI=(Ta+
        ( 6.07562052E-01 )   + 
        ( -2.27712343E-02 ) * Ta + 
        ( 8.06470249E-04 ) * Ta*Ta + 
        ( -1.54271372E-04 ) * Ta*Ta*Ta + 
        ( -3.24651735E-06 ) * Ta*Ta*Ta*Ta + 
        ( 7.32602852E-08 ) * Ta*Ta*Ta*Ta*Ta + 
        ( 1.35959073E-09 ) * Ta*Ta*Ta*Ta*Ta*Ta + 
        ( -2.25836520E+00 ) * va + 
        ( 8.80326035E-02 ) * Ta*va + 
        ( 2.16844454E-03 ) * Ta*Ta*va + 
        ( -1.53347087E-05 ) * Ta*Ta*Ta*va + 
        ( -5.72983704E-07 ) * Ta*Ta*Ta*Ta*va + 
        ( -2.55090145E-09 ) * Ta*Ta*Ta*Ta*Ta*va + 
        ( -7.51269505E-01 ) * va*va + 
        ( -4.08350271E-03 ) * Ta*va*va + 
        ( -5.21670675E-05 ) * Ta*Ta*va*va + 
        ( 1.94544667E-06 ) * Ta*Ta*Ta*va*va + 
        ( 1.14099531E-08 ) * Ta*Ta*Ta*Ta*va*va + 
        ( 1.58137256E-01 ) * va*va*va + 
        ( -6.57263143E-05 ) * Ta*va*va*va + 
        ( 2.22697524E-07 ) * Ta*Ta*va*va*va + 
        ( -4.16117031E-08 ) * Ta*Ta*Ta*va*va*va + 
        ( -1.27762753E-02 ) * va*va*va*va + 
        ( 9.66891875E-06 ) * Ta*va*va*va*va + 
        ( 2.52785852E-09 ) * Ta*Ta*va*va*va*va + 
        ( 4.56306672E-04 ) * va*va*va*va*va + 
        ( -1.74202546E-07 ) * Ta*va*va*va*va*va + 
        ( -5.91491269E-06 ) * va*va*va*va*va*va + 
        ( 3.98374029E-01 ) * D_Tmrt + 
        ( 1.83945314E-04 ) * Ta*D_Tmrt + 
        ( -1.73754510E-04 ) * Ta*Ta*D_Tmrt + 
        ( -7.60781159E-07 ) * Ta*Ta*Ta*D_Tmrt + 
        ( 3.77830287E-08 ) * Ta*Ta*Ta*Ta*D_Tmrt + 
        ( 5.43079673E-10 ) * Ta*Ta*Ta*Ta*Ta*D_Tmrt + 
        ( -2.00518269E-02 ) * va*D_Tmrt + 
        ( 8.92859837E-04 ) * Ta*va*D_Tmrt + 
        ( 3.45433048E-06 ) * Ta*Ta*va*D_Tmrt + 
        ( -3.77925774E-07 ) * Ta*Ta*Ta*va*D_Tmrt + 
        ( -1.69699377E-09 ) * Ta*Ta*Ta*Ta*va*D_Tmrt + 
        ( 1.69992415E-04 ) * va*va*D_Tmrt + 
        ( -4.99204314E-05 ) * Ta*va*va*D_Tmrt + 
        ( 2.47417178E-07 ) * Ta*Ta*va*va*D_Tmrt + 
        ( 1.07596466E-08 ) * Ta*Ta*Ta*va*va*D_Tmrt + 
        ( 8.49242932E-05 ) * va*va*va*D_Tmrt + 
        ( 1.35191328E-06 ) * Ta*va*va*va*D_Tmrt + 
        ( -6.21531254E-09 ) * Ta*Ta*va*va*va*D_Tmrt + 
        ( -4.99410301E-06 ) * va*va*va*va*D_Tmrt + 
        ( -1.89489258E-08 ) * Ta*va*va*va*va*D_Tmrt + 
        ( 8.15300114E-08 ) * va*va*va*va*va*D_Tmrt + 
        ( 7.55043090E-04 ) * D_Tmrt*D_Tmrt + 
        ( -5.65095215E-05 ) * Ta*D_Tmrt*D_Tmrt + 
        ( -4.52166564E-07 ) * Ta*Ta*D_Tmrt*D_Tmrt + 
        ( 2.46688878E-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt + 
        ( 2.42674348E-10 ) * Ta*Ta*Ta*Ta*D_Tmrt*D_Tmrt + 
        ( 1.54547250E-04 ) * va*D_Tmrt*D_Tmrt + 
        ( 5.24110970E-06 ) * Ta*va*D_Tmrt*D_Tmrt + 
        ( -8.75874982E-08 ) * Ta*Ta*va*D_Tmrt*D_Tmrt + 
        ( -1.50743064E-09 ) * Ta*Ta*Ta*va*D_Tmrt*D_Tmrt + 
        ( -1.56236307E-05 ) * va*va*D_Tmrt*D_Tmrt + 
        ( -1.33895614E-07 ) * Ta*va*va*D_Tmrt*D_Tmrt + 
        ( 2.49709824E-09 ) * Ta*Ta*va*va*D_Tmrt*D_Tmrt + 
        ( 6.51711721E-07 ) * va*va*va*D_Tmrt*D_Tmrt + 
        ( 1.94960053E-09 ) * Ta*va*va*va*D_Tmrt*D_Tmrt + 
        ( -1.00361113E-08 ) * va*va*va*va*D_Tmrt*D_Tmrt + 
        ( -1.21206673E-05 ) * D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -2.18203660E-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 7.51269482E-09 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 9.79063848E-11 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 1.25006734E-06 ) * va*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -1.81584736E-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -3.52197671E-10 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -3.36514630E-08 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 1.35908359E-10 ) * Ta*va*va*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 4.17032620E-10 ) * va*va*va*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -1.30369025E-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 4.13908461E-10 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 9.22652254E-12 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -5.08220384E-09 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -2.24730961E-11 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 1.17139133E-10 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 6.62154879E-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 4.03863260E-13 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 1.95087203E-12 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( -4.73602469E-12 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + 
        ( 5.12733497E+00 ) * Pa + 
        ( -3.12788561E-01 ) * Ta*Pa + 
        ( -1.96701861E-02 ) * Ta*Ta*Pa + 
        ( 9.99690870E-04 ) * Ta*Ta*Ta*Pa + 
        ( 9.51738512E-06 ) * Ta*Ta*Ta*Ta*Pa + 
        ( -4.66426341E-07 ) * Ta*Ta*Ta*Ta*Ta*Pa + 
        ( 5.48050612E-01 ) * va*Pa + 
        ( -3.30552823E-03 ) * Ta*va*Pa + 
        ( -1.64119440E-03 ) * Ta*Ta*va*Pa + 
        ( -5.16670694E-06 ) * Ta*Ta*Ta*va*Pa + 
        ( 9.52692432E-07 ) * Ta*Ta*Ta*Ta*va*Pa + 
        ( -4.29223622E-02 ) * va*va*Pa + 
        ( 5.00845667E-03 ) * Ta*va*va*Pa + 
        ( 1.00601257E-06 ) * Ta*Ta*va*va*Pa + 
        ( -1.81748644E-06 ) * Ta*Ta*Ta*va*va*Pa + 
        ( -1.25813502E-03 ) * va*va*va*Pa + 
        ( -1.79330391E-04 ) * Ta*va*va*va*Pa + 
        ( 2.34994441E-06 ) * Ta*Ta*va*va*va*Pa + 
        ( 1.29735808E-04 ) * va*va*va*va*Pa + 
        ( 1.29064870E-06 ) * Ta*va*va*va*va*Pa + 
        ( -2.28558686E-06 ) * va*va*va*va*va*Pa + 
        ( -3.69476348E-02 ) * D_Tmrt*Pa + 
        ( 1.62325322E-03 ) * Ta*D_Tmrt*Pa + 
        ( -3.14279680E-05 ) * Ta*Ta*D_Tmrt*Pa + 
        ( 2.59835559E-06 ) * Ta*Ta*Ta*D_Tmrt*Pa + 
        ( -4.77136523E-08 ) * Ta*Ta*Ta*Ta*D_Tmrt*Pa + 
        ( 8.64203390E-03 ) * va*D_Tmrt*Pa + 
        ( -6.87405181E-04 ) * Ta*va*D_Tmrt*Pa + 
        ( -9.13863872E-06 ) * Ta*Ta*va*D_Tmrt*Pa + 
        ( 5.15916806E-07 ) * Ta*Ta*Ta*va*D_Tmrt*Pa + 
        ( -3.59217476E-05 ) * va*va*D_Tmrt*Pa + 
        ( 3.28696511E-05 ) * Ta*va*va*D_Tmrt*Pa + 
        ( -7.10542454E-07 ) * Ta*Ta*va*va*D_Tmrt*Pa + 
        ( -1.24382300E-05 ) * va*va*va*D_Tmrt*Pa + 
        ( -7.38584400E-09 ) * Ta*va*va*va*D_Tmrt*Pa + 
        ( 2.20609296E-07 ) * va*va*va*va*D_Tmrt*Pa + 
        ( -7.32469180E-04 ) * D_Tmrt*D_Tmrt*Pa + 
        ( -1.87381964E-05 ) * Ta*D_Tmrt*D_Tmrt*Pa + 
        ( 4.80925239E-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa + 
        ( -8.75492040E-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*Pa + 
        ( 2.77862930E-05 ) * va*D_Tmrt*D_Tmrt*Pa + 
        ( -5.06004592E-06 ) * Ta*va*D_Tmrt*D_Tmrt*Pa + 
        ( 1.14325367E-07 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*Pa + 
        ( 2.53016723E-06 ) * va*va*D_Tmrt*D_Tmrt*Pa + 
        ( -1.72857035E-08 ) * Ta*va*va*D_Tmrt*D_Tmrt*Pa + 
        ( -3.95079398E-08 ) * va*va*va*D_Tmrt*D_Tmrt*Pa + 
        ( -3.59413173E-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( 7.04388046E-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( -1.89309167E-08 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( -4.79768731E-07 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( 7.96079978E-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( 1.62897058E-09 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( 3.94367674E-08 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( -1.18566247E-09 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( 3.34678041E-10 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( -1.15606447E-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + 
        ( -2.80626406E+00 ) * Pa*Pa + 
        ( 5.48712484E-01 ) * Ta*Pa*Pa + 
        ( -3.99428410E-03 ) * Ta*Ta*Pa*Pa + 
        ( -9.54009191E-04 ) * Ta*Ta*Ta*Pa*Pa + 
        ( 1.93090978E-05 ) * Ta*Ta*Ta*Ta*Pa*Pa + 
        ( -3.08806365E-01 ) * va*Pa*Pa + 
        ( 1.16952364E-02 ) * Ta*va*Pa*Pa + 
        ( 4.95271903E-04 ) * Ta*Ta*va*Pa*Pa + 
        ( -1.90710882E-05 ) * Ta*Ta*Ta*va*Pa*Pa + 
        ( 2.10787756E-03 ) * va*va*Pa*Pa + 
        ( -6.98445738E-04 ) * Ta*va*va*Pa*Pa + 
        ( 2.30109073E-05 ) * Ta*Ta*va*va*Pa*Pa + 
        ( 4.17856590E-04 ) * va*va*va*Pa*Pa + 
        ( -1.27043871E-05 ) * Ta*va*va*va*Pa*Pa + 
        ( -3.04620472E-06 ) * va*va*va*va*Pa*Pa + 
        ( 5.14507424E-02 ) * D_Tmrt*Pa*Pa + 
        ( -4.32510997E-03 ) * Ta*D_Tmrt*Pa*Pa + 
        ( 8.99281156E-05 ) * Ta*Ta*D_Tmrt*Pa*Pa + 
        ( -7.14663943E-07 ) * Ta*Ta*Ta*D_Tmrt*Pa*Pa + 
        ( -2.66016305E-04 ) * va*D_Tmrt*Pa*Pa + 
        ( 2.63789586E-04 ) * Ta*va*D_Tmrt*Pa*Pa + 
        ( -7.01199003E-06 ) * Ta*Ta*va*D_Tmrt*Pa*Pa + 
        ( -1.06823306E-04 ) * va*va*D_Tmrt*Pa*Pa + 
        ( 3.61341136E-06 ) * Ta*va*va*D_Tmrt*Pa*Pa + 
        ( 2.29748967E-07 ) * va*va*va*D_Tmrt*Pa*Pa + 
        ( 3.04788893E-04 ) * D_Tmrt*D_Tmrt*Pa*Pa + 
        ( -6.42070836E-05 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( 1.16257971E-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( 7.68023384E-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( -5.47446896E-07 ) * Ta*va*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( -3.59937910E-08 ) * va*va*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( -4.36497725E-06 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( 1.68737969E-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( 2.67489271E-08 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( 3.23926897E-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + 
        ( -3.53874123E-02 ) * Pa*Pa*Pa + 
        ( -2.21201190E-01 ) * Ta*Pa*Pa*Pa + 
        ( 1.55126038E-02 ) * Ta*Ta*Pa*Pa*Pa + 
        ( -2.63917279E-04 ) * Ta*Ta*Ta*Pa*Pa*Pa + 
        ( 4.53433455E-02 ) * va*Pa*Pa*Pa + 
        ( -4.32943862E-03 ) * Ta*va*Pa*Pa*Pa + 
        ( 1.45389826E-04 ) * Ta*Ta*va*Pa*Pa*Pa + 
        ( 2.17508610E-04 ) * va*va*Pa*Pa*Pa + 
        ( -6.66724702E-05 ) * Ta*va*va*Pa*Pa*Pa + 
        ( 3.33217140E-05 ) * va*va*va*Pa*Pa*Pa + 
        ( -2.26921615E-03 ) * D_Tmrt*Pa*Pa*Pa + 
        ( 3.80261982E-04 ) * Ta*D_Tmrt*Pa*Pa*Pa + 
        ( -5.45314314E-09 ) * Ta*Ta*D_Tmrt*Pa*Pa*Pa + 
        ( -7.96355448E-04 ) * va*D_Tmrt*Pa*Pa*Pa + 
        ( 2.53458034E-05 ) * Ta*va*D_Tmrt*Pa*Pa*Pa + 
        ( -6.31223658E-06 ) * va*va*D_Tmrt*Pa*Pa*Pa + 
        ( 3.02122035E-04 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa + 
        ( -4.77403547E-06 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa*Pa + 
        ( 1.73825715E-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa*Pa + 
        ( -4.09087898E-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa*Pa + 
        ( 6.14155345E-01 ) * Pa*Pa*Pa*Pa + 
        ( -6.16755931E-02 ) * Ta*Pa*Pa*Pa*Pa + 
        ( 1.33374846E-03 ) * Ta*Ta*Pa*Pa*Pa*Pa + 
        ( 3.55375387E-03 ) * va*Pa*Pa*Pa*Pa + 
        ( -5.13027851E-04 ) * Ta*va*Pa*Pa*Pa*Pa + 
        ( 1.02449757E-04 ) * va*va*Pa*Pa*Pa*Pa + 
        ( -1.48526421E-03 ) * D_Tmrt*Pa*Pa*Pa*Pa + 
        ( -4.11469183E-05 ) * Ta*D_Tmrt*Pa*Pa*Pa*Pa + 
        ( -6.80434415E-06 ) * va*D_Tmrt*Pa*Pa*Pa*Pa + 
        ( -9.77675906E-06 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa*Pa + 
        ( 8.82773108E-02 ) * Pa*Pa*Pa*Pa*Pa + 
        ( -3.01859306E-03 ) * Ta*Pa*Pa*Pa*Pa*Pa + 
        ( 1.04452989E-03 ) * va*Pa*Pa*Pa*Pa*Pa + 
        ( 2.47090539E-04 ) * D_Tmrt*Pa*Pa*Pa*Pa*Pa + 
        ( 1.48348065E-03 ) * Pa*Pa*Pa*Pa*Pa*Pa )

    return UTCI


"""
### END OF PYTHON FUNCTION:  UTCIfunc
###
###############################################################################
###############################################################################
###############################################################################
"""

"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  get_lat_lon
###
### Inputs:
###   lat, lon, dist_nm, bearing
###   
###
### Outputs:
###   lat , lon
###
### PURPOSE: get lat lon of a point x nm's away at a given bearing
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 1 Jan 2020 SSgt Jones - added for a SAR
###
###############################################################################
"""
def get_lat_lon(lat,lon,dist,bear):
    R = 6378.1 #Radius of the Earth
    ### lat , lon in Degree Decimals
    ### distance in NM
    ### bearing in Degrees

    ## convert to radians and km 
    lat = math.radians(lat)
    lon = math.radians(lon)
    bear = math.radians(bear)
    dist = dist*1.852
    
    lat_out = math.asin( math.sin(lat)*math.cos(dist/R) + math.cos(lat)*math.sin(dist/R)*math.cos(bear))
    lon_out = lon + math.atan2(math.sin(bear)*math.sin(dist/R)*math.cos(lat),math.cos(dist/R)-math.sin(lat)*math.sin(lat_out))

    lat_out = math.degrees(lat_out)
    lon_out = math.degrees(lon_out)

    return lat_out,lon_out


"""
### END OF PYTHON FUNCTION:  get_lat_lon
###
###############################################################################
###############################################################################
###############################################################################
"""
"""
###############################################################################
###############################################################################
###############################################################################
###
### PYTHON FUNCTION:  sol_lun
###
### Inputs:
###   SITE_LAT, SITE_LON, SITE_ELEVATION, START_YYYYMMDD, END_YYYYMMDD, SITE_PRESSURE, SITE_TEMPERATURE
###   
###
### Outputs:
###   'Date', 'Nautical Twilight Begin', 'Civil Twilight Begin', 'Sunrise', 'Sunset',
###   'Civil Twilight End', 'Nautical Twilight End', 'Moonrise', 'Moonset', 'Lunar Illumination (%) 00UTC'
###
### PURPOSE: get above columns for a location
###
### METHOD:
###
### REFERENCES:
###
### EXTERNAL ROUTINES CALLED:
###
### FILES/DATABASE/TABLES ACCESSED:
###
### DATA DICTIONARY:  {List function internal variables and units}
###
### REMARKS:
###
### UPDATES: 1 Jan 2020 SSgt Jones - added for a SAR
###
###############################################################################
"""

def sol_lun(SITE_LAT,SITE_LON,SITE_ELEVATION,START_YYYYMMDD,END_YYYYMMDD,SITE_PRESSURE,SITE_TEMPERATURE):
    ephem_base_date=dt.datetime(1899,12,31,12, 0, 0)
    astronominal_epoch=dt.datetime(2020, 1, 1, 0, 0, 0)
    DEFAULT_REFRACTION_ADJUSTMENT=-0.5666667
    USE_DEFAULT_REFRACTION=False
    
    if pd.isnull(SITE_PRESSURE):
        USE_DEFAULT_REFRACTION=True
        SITE_PRESSURE=1013.0
        SITE_TEMPERATURE=15.0
        
    if pd.isnull(SITE_TEMPERATURE):
        USE_DEFAULT_REFRACTION=True
        SITE_PRESSURE=1013.0
        SITE_TEMPERATURE=15.0

    ###########################################################
    #####     Changes not necessary below these lines     #####
     ###                 for standard runs                 ###
      #                                                     #


    ###############
    # Set up site #
    ###############

    site=ephem.Observer()
    site.epoch=astronominal_epoch
    site.lat=str(SITE_LAT)
    site.lon=str(SITE_LON)
    site.elevation=SITE_ELEVATION

    # Refraction settings
    if USE_DEFAULT_REFRACTION:
        refraction_adjustment=DEFAULT_REFRACTION_ADJUSTMENT
        site.pressure=0
    else:
        refraction_adjustment=0.0
        site.pressure=SITE_PRESSURE
        site.temp=SITE_TEMPERATURE

    # Site timezone
    timezone=TimezoneFinder()
    site_timezone=timezone.timezone_at(lng=SITE_LON,lat=SITE_LAT)


    ####################################
    # Echo settings to standard output #
    ####################################

##    print('Site Lat:  '+str(SITE_LAT),flush=True)
##    print('Site Lon:  '+str(SITE_LON),flush=True)
##    print('Site Elevation:  '+str(SITE_ELEVATION)+'m',flush=True)
##    if USE_DEFAULT_REFRACTION:
##        print('Default refraction set to '+str(DEFAULT_REFRACTION_ADJUSTMENT)+' degrees along look-horizons',flush=True)
##    else:
##        print('Dynamic refraction set with...',flush=True)
##        print('    Site Surface Pressure:  '+str(SITE_PRESSURE)+'mb',flush=True)
##        print('    Site Surface Temperature:  '+str(SITE_TEMPERATURE)+'C',flush=True)
##    print('Site Local Timezone:  '+site_timezone,flush=True)
##    print('--------------------------------------------------------------------------',flush=True)


    ####################
    ####################
    ####################
    ### Main program ###
    ####################
    ####################
    ####################

##    print('Calculating event date-times...',flush=True)

    # Initialize output dataframe
    output_columns=['Date','Nautical Twilight Begin','Civil Twilight Begin','Sunrise','Sunset','Civil Twilight End','Nautical Twilight End','Moonrise','Moonset','Lunar Illumination (%) 00UTC']
    solar_lunar_data=pd.DataFrame(columns=output_columns)

    # Loop through days of POR and build output dataframe
    CURRENTDATE=dt.datetime.strptime(START_YYYYMMDD,'%Y%m%d')
    FINALDATE=dt.datetime.strptime(END_YYYYMMDD,'%Y%m%d')

    # Calculate moon phase illumination data for this date at 00UTC
    moon_data=ephem.Moon()
    date_time=dt.datetime(CURRENTDATE.year,CURRENTDATE.month,CURRENTDATE.day,0,0,0,tzinfo=pytz.utc)
    moon_data.compute(date_time)
    moon_illumination=moon_data.phase

    # ...for other calculations, Set base date to solar noon (about 12 local)
    local_tz=pytz.timezone(site_timezone)
    date_time=dt.datetime(CURRENTDATE.year,CURRENTDATE.month,CURRENTDATE.day,12,0,0,tzinfo=local_tz)
    # Convert local noon back to UTC for calculations
    site.date=date_time.astimezone(pytz.utc)
    

    # Calculate moon rise/set data for this date
    moon_data=ephem.Moon()
    moon_data.compute(date_time)

    # Calculate flat-horizon solar/lunar data (results are in days since ephem_epoch)
    site.horizon=str(0.0+refraction_adjustment)
    sunrise=site.previous_rising(ephem.Sun())
    sunset=site.next_setting(ephem.Sun())
    moonrise=site.previous_rising(ephem.Moon())
    moonset=site.next_setting(ephem.Moon())
    
    # Civil rise/set (site.horizon='-6')
    site.horizon=str(-6.0+refraction_adjustment)
    civil_twilight_begin=site.previous_rising(ephem.Sun())
    civil_twilight_end=site.next_setting(ephem.Sun())

    # Nautical rise/set (site.horizon='-12')
    site.horizon=str(-12.0+refraction_adjustment)
    nautical_twilight_begin=site.previous_rising(ephem.Sun())
    nautical_twilight_end=site.next_setting(ephem.Sun())

    # Adjust resultant date-times above from days since ephem_base_date to an actual date-time object
    # (date-times are currently days from ephem_base_date)
    nautical_twilight_begin=ephem_base_date+dt.timedelta(days=nautical_twilight_begin)
    civil_twilight_begin=ephem_base_date+dt.timedelta(days=civil_twilight_begin)
    sunrise=ephem_base_date+dt.timedelta(days=sunrise)
    sunset=ephem_base_date+dt.timedelta(days=sunset)
    civil_twilight_end=ephem_base_date+dt.timedelta(days=civil_twilight_end)
    nautical_twilight_end=ephem_base_date+dt.timedelta(days=nautical_twilight_end)
    moonrise=ephem_base_date+dt.timedelta(days=moonrise)
    moonset=ephem_base_date+dt.timedelta(days=moonset)

    # Add 30 seconds to each event to round to nearest minute, as all event resolutions are minutes
    # and seconds are truncated off
    nautical_twilight_begin=nautical_twilight_begin+dt.timedelta(seconds=30)
    civil_twilight_begin=civil_twilight_begin+dt.timedelta(seconds=30)
    sunrise=sunrise+dt.timedelta(seconds=30)
    sunset=sunset+dt.timedelta(seconds=30)
    civil_twilight_end=civil_twilight_end+dt.timedelta(seconds=30)
    nautical_twilight_end=nautical_twilight_end+dt.timedelta(seconds=30)
    moonrise=moonrise+dt.timedelta(seconds=30)
    moonset=moonset+dt.timedelta(seconds=30)

    # Turn these data into a single row dataframe to append to solar_lunar_data
    return date_time,nautical_twilight_begin.strftime('%Y-%m-%d %H:%M'),civil_twilight_begin.strftime('%Y-%m-%d %H:%M'),sunrise.strftime('%Y-%m-%d %H:%M'),sunset.strftime('%Y-%m-%d %H:%M'),civil_twilight_end.strftime('%Y-%m-%d %H:%M'),nautical_twilight_end.strftime('%Y-%m-%d %H:%M'),moonrise.strftime('%Y-%m-%d %H:%M'),moonset.strftime('%Y-%m-%d %H:%M'),round(moon_illumination,2)

##    tes separate from event times
##    solar_lunar_data['Nautical Twilight Begin Date']=solar_lunar_data['Nautical Twilight Begin'].dt.date
##    solar_lunar_data['Civil Twilight Begin Date']=solar_lunar_data['Civil Twilight Begin'].dt.date
##    solar_lunar_data['Sunrise Date']=solar_lunar_data['Sunrise'].dt.date
##    solar_lunar_data['Sunset Date']=solar_lunar_data['Sunset'].dt.date
##    solar_lunar_data['Nautical Twilight End Date']=solar_lunar_data['Nautical Twilight End'].dt.date
##    solar_lunar_data['Civil Twilight End Date']=solar_lunar_data['Civil Twilight End'].dt.date
##    solar_lunar_data['Moonrise Date']=solar_lunar_data['Moonrise'].dt.date
##    solar_lunar_data['Moonset Date']=solar_lunar_data['Moonset'].dt.date
##    solar_lunar_data['Nautical Twilight Begin Time']=solar_lunar_data['Nautical Twilight Begin'].dt.time
##    solar_lunar_data['Civil Twilight Begin Time']=solar_lunar_data['Civil Twilight Begin'].dt.time
##    solar_lunar_data['Sunrise Time']=solar_lunar_data['Sunrise'].dt.time
##    solar_lunar_data['Sunset Time']=solar_lunar_data['Sunset'].dt.time
##    solar_lunar_data['Nautical Twilight End Time']=solar_lunar_data['Nautical Twilight End'].dt.time
##    solar_lunar_data['Civil Twilight End Time']=solar_lunar_data['Civil Twilight End'].dt.time
##    solar_lunar_data['Moonrise Time']=solar_lunar_data['Moonrise'].dt.time
##    solar_lunar_data['Moonset Time']=solar_lunar_data['Moonset'].dt.time
##    # Generate new table to synchronize YYYY MM DD in events
##    # (On any given row of YYYY MM DD, all of the listed events in that row will be on that YYYY MM DD.)
##    solar_lunar_table=solar_lunar_data.copy()
##    solar_lunar_table.index=solar_lunar_table['Date']
##
##    # The original event columns ('Sunrise', 'Sunset', 'Moonrise'...) will be converted to
##    # hold just the time for the matching YYYY MM DD for that event.  First set to NULL...
##    solar_lunar_table['Nautical Twilight Begin']=''
##    solar_lunar_table['Civil Twilight Begin']=''
##    solar_lunar_table['Sunrise']=''
##    solar_lunar_table['Sunset']=''
##    solar_lunar_table['Nautical Twilight End']=''
##    solar_lunar_table['Civil Twilight End']=''
##    solar_lunar_table['Moonrise']=''
##    solar_lunar_table['Moonset']=''
##
##    for index,row in solar_lunar_table.iterrows():
##        row_date=row.Date
##        for col in ['Nautical Twilight Begin','Civil Twilight Begin','Sunrise','Sunset','Civil Twilight End','Nautical Twilight End','Moonrise','Moonset']:
##            row=solar_lunar_table.iloc[np.where(solar_lunar_table[col+' Date']==row_date)]
##            if len(row)>0: solar_lunar_table.at[index,col]=row[col+' Time'][0].strftime('%H:%M')
##
##    clean_table=solar_lunar_table[['Date','Nautical Twilight Begin','Civil Twilight Begin',
##                                   'Sunrise','Sunset','Civil Twilight End','Nautical Twilight End',
##                                   'Moonrise','Moonset','Lunar Illumination (%) 00UTC']].copy().reset_index(drop=True)
##
##    return clean_table.loc[0,'Date'],clean_table.loc[0,'Nautical Twilight Begin'],clean_table.loc[0,'Civil Twilight Begin'],clean_table.loc[0,'Sunrise'],clean_table.loc[0,'Sunset'],clean_table.loc[0,'Civil Twilight End'],clean_table.loc[0,'Nautical Twilight End'],clean_table.loc[0,'Moonrise'],clean_table.loc[0,'Moonset'],clean_table.loc[0,'Lunar Illumination (%) 00UTC']
##        
    
    
"""
### END OF PYTHON FUNCTION:  sol_lun
###
###############################################################################
###############################################################################
###############################################################################
"""



"""
Formating functions below, along with round to the nearest half and quarter
"""

# used to create a column called trimo by deciding which 10 day bucket the day is in.
def trimo(day):
    if day <= 10:
        out = 1
    elif (day>10) & (day<=20):
        out = 2
    else:
        out = 3
    return out

# used to format the output when trimo is present
def TMformat(tm):
    if tm==1:
        out = '1st-10th'
    elif tm==2:
        out = '11th-20th'
    else:
        out = '21st-EOM'
    return out

# used to create a column called hexmo by deciding which 5 day bucket the day is in.
def hexmo(day):
    if day <= 5:
        out = 1
    elif (day>5) & (day<=10):
        out = 2
    elif (day>10) & (day<=15):
        out = 3
    elif (day>15) & (day<=20):
        out = 4
    elif (day>20) & (day<=25):
        out = 5
    else:
        out = 6
    return out

# used to format the output when hexmo is present
def HMformat(hm):
    if hm==1:
        out = '1st-5th'
    elif hm==2:
        out = '6th-10th'
    elif hm==3:
        out = '11th-15th'
    elif hm==4:
        out = '16th-20th'
    elif hm==5:
        out = '21st-25th'
    else:
        out = '26th-EOM'
    return out 

# round value toe the nearest 1/4 decimal
def q_round(x):
    return round(x*4.0)/4.0

# round value toe the nearest 1/2 decimal
def h_round(x):
    if np.isnan(x):
        return x
    else:
        return round(x*2.0)/2.0



