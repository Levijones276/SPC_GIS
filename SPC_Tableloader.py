import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString
import datetime as dt
import matplotlib.pyplot as plt
##import SFC_TableLoader
import MetFuncs_Raw
##import MetFuncs_Object
##import HeatMap_Plotterelat
##import Logo_plot
idx = pd.IndexSlice
import warnings
warnings.filterwarnings("ignore")




####### https://www.spc.noaa.gov/wcm/#data #######
class SPC_TableLoader:
    def __init__(self,*args,**kwargs):
        self.location = os.getcwd()+'/inputs/'
        if not os.path.exists(self.location):
            os.mkdir(self.location)
        
        if len(args) > 0:
            if args[0].lower() in ['-h','h','help','-help']:

                print("""
SPC_TableLoader:\n
  keywords:
    por = [begpor,endpor]YYYYMMDD
    hail_path = 'path to data if not default'
    torn_path = 'path to data if not default'
    wind_path = 'path to data if not default'
    bbox = [nlat,slat,elon,wlon]
    site = [icao or lat,plat or lon, distance in NM] - still in dev and can't be used with bbox
    """)
                return

            else:
                print("Only non keyword arg is 'help'")
        
        if kwargs:
            if 'hail_path' in kwargs.keys():
                hail_path = kwargs['hail_path']
            else:
                hail_path = self.location+'hail.csv'

            if 'torn_path' in kwargs.keys():
                torn_path = kwargs['torn_path']
            else:
                torn_path = self.location+'torn.csv'

            if 'wind_path' in kwargs.keys():
                wind_path = kwargs['wind_path']
            else:
                wind_path = self.location+'wind.csv'
            if 'update' in kwargs.keys():
                if kwargs['update']:
                    self.Update_CSVs()
                    
        else:
            hail_path = self.location+'hail.csv'
            torn_path = self.location+'torn.csv'
            wind_path = self.location+'wind.csv'

        print('Loading SPC Data....\n..May take a few minutes')

        self.time_type = 'central'
        self.net_plat = False

        try:
            self.process(hail_path,torn_path,wind_path)
            stop = False
        except:
            print('Need to run .Update_CSVs if on mac/linux with wget or \n download SPC data, rename to torn.csv,hail.csv,wind.csv and save in inputs directory')
            stop = True
        if not stop:
            if kwargs:
                if 'por' in kwargs.keys():
                    self.torn['datetime'] = pd.to_datetime({'year':self.torn['yr'].values,'month':self.torn['mo'].values,'day':self.torn['dy'].values,'hour':pd.to_datetime(self.torn['time']).dt.hour})
                    self.hail['datetime'] = pd.to_datetime({'year':self.hail['yr'].values,'month':self.hail['mo'].values,'day':self.hail['dy'].values,'hour':pd.to_datetime(self.hail['time']).dt.hour})
                    self.wind['datetime'] = pd.to_datetime({'year':self.wind['yr'].values,'month':self.wind['mo'].values,'day':self.wind['dy'].values,'hour':pd.to_datetime(self.wind['time']).dt.hour})
                    if type(kwargs['por']) == list:
                        kwargs['por'] = [str(x) for x in kwargs['por']]
                        self.torn = self.torn[(self.torn['datetime']>=dt.datetime.strptime(kwargs['por'][0],'%Y%m%d%H%M'))&
                                            (self.torn['datetime']<=dt.datetime.strptime(kwargs['por'][1],'%Y%m%d%H%M'))]

                        self.hail = self.hail[(self.hail['datetime']>=dt.datetime.strptime(kwargs['por'][0],'%Y%m%d%H%M'))&
                                            (self.hail['datetime']<=dt.datetime.strptime(kwargs['por'][1],'%Y%m%d%H%M'))]

                        self.wind = self.wind[(self.wind['datetime']>=dt.datetime.strptime(kwargs['por'][0],'%Y%m%d%H%M'))&
                                            (self.wind['datetime']<=dt.datetime.strptime(kwargs['por'][1],'%Y%m%d%H%M'))]
                        
                        
                    elif type(kwargs['por']) == tuple:
                        kwargs['por'] = [str(x) for x in kwargs['por']]
                        self.torn = self.torn[(self.torn['datetime']>=dt.datetime.strptime(kwargs['por'][0],'%Y%m%d%H%M'))&
                                            (self.torn['datetime']<=dt.datetime.strptime(kwargs['por'][1],'%Y%m%d%H%M'))]

                        self.hail = self.hail[(self.hail['datetime']>=dt.datetime.strptime(kwargs['por'][0],'%Y%m%d%H%M'))&
                                            (self.hail['datetime']<=dt.datetime.strptime(kwargs['por'][1],'%Y%m%d%H%M'))]

                        self.wind = self.wind[(self.wind['datetime']>=dt.datetime.strptime(kwargs['por'][0],'%Y%m%d%H%M'))&
                                            (self.wind['datetime']<=dt.datetime.strptime(kwargs['por'][1],'%Y%m%d%H%M'))]
                        
                    else:
                        print('por must be a List or Tuple\n[begpor,endpor]YYYMMDD')
                        sys.exit()

                ## bbox and site can't be in the key word dict at the same time ##
                if 'bbox' in kwargs.keys():
                    if 'site' in kwargs.keys():
                        print("Can't pass bbox and site at the same time!!" )
                        sys.exit()
                    if type(kwargs['bbox']) == list:
                        bbox = kwargs['bbox']
                    elif type(kwargs['bbox']) == tuple:
                        bbox = [x for x in kwargs['bbox']]
                    else:
                        print('Bounding must be List or Tuple\n[nlat,slat,elon,wlon]')
                        sys.exit()

                    self.torn = self.torn[((self.torn.slat <= bbox[0]) &
                                        (self.torn.slat >= bbox[1]))|
                                        ((self.torn.elat <= bbox[0]) &
                                        (self.torn.elat >= bbox[1]))]

                    self.torn = self.torn[((self.torn.slon <= bbox[2]) &
                                        (self.torn.slon >= bbox[3]))|
                                        ((self.torn.elon <= bbox[2]) &
                                        (self.torn.elon >= bbox[3]))]
                
                    self.hail = self.hail[((self.hail.slat <= bbox[0]) &
                                        (self.hail.slat >= bbox[1]))|
                                        ((self.hail.elat <= bbox[0]) &
                                        (self.hail.elat >= bbox[1]))]

                    self.hail = self.hail[((self.hail.slon <= bbox[2]) &
                                        (self.hail.slon >= bbox[3]))|
                                        ((self.hail.elon <= bbox[2]) &
                                        (self.hail.elon >= bbox[3]))]

                    self.wind = self.wind[((self.wind.slat <= bbox[0]) &
                                        (self.wind.slat >= bbox[1]))|
                                        ((self.wind.elat <= bbox[0]) &
                                        (self.wind.elat >= bbox[1]))]

                    self.wind = self.wind[((self.wind.slon <= bbox[2]) &
                                        (self.wind.slon >= bbox[3]))|
                                        ((self.wind.elon <= bbox[2]) &
                                        (self.wind.elon >= bbox[3]))]

                    
                    if self.time_type == 'central':
                        self.Time_to_Zulu()

                    self.torn.reset_index(inplace=True,drop=True)
                    self.torn = self.torn[['yr', 'mo', 'dy','hr','time','mag','slat', 'slon', 'elat', 'elon','YEAR','MO','DAY','HR','geometry']]

                    self.hail.reset_index(inplace=True,drop=True)
                    self.hail = self.hail[['yr', 'mo', 'dy','hr','time','mag','slat', 'slon', 'elat', 'elon','YEAR','MO','DAY','HR','geometry']]

                    self.wind.reset_index(inplace=True,drop=True)
                    self.wind = self.wind[['yr', 'mo', 'dy','hr','time','mag','slat', 'slon', 'elat', 'elon','YEAR','MO','DAY','HR','geometry']]
                    

                if 'site' in kwargs.keys():
                    if 'bbox' in kwargs.keys():
                        print("Can't pass bbox and site at the same time!!" )
                        sys.exit()
    ##                print('No Site logic yet')
                    
                    try:
                        float(kwargs['site'][0])
                        self.net_plat = False
                    except ValueError:
                        self.net_plat = True

                    if self.net_plat:
                        net=kwargs['site'][0]
                        plat=kwargs['site'][1]
                        dist=kwargs['site'][2]
                        self.site_df = pd.read_csv(self.location+'stations.csv')
                        self.site_df = self.site_df[(self.site_df.NET == net.upper()) & (self.site_df.PLAT == plat.upper())].reset_index(drop = True)
                        lat = self.site_df.loc[0,'LAT']
                        lon = self.site_df.loc[0,'LON']
                        site_point = Point(lon,lat)
                        
                        
                    else:
                        lat = kwargs['site'][0]
                        lon = kwargs['site'][1]
                        dist= kwargs['site'][2]
                        site_point = Point(lon,lat)
                        
                    


                    self.hail['dist'],self.hail['bear']= np.vectorize(MetFuncs_Raw.DistBear)(lat,lon,self.hail['slat'],self.hail['slon'])
                    self.wind['dist'],self.wind['bear']= np.vectorize(MetFuncs_Raw.DistBear)(lat,lon,self.wind['slat'],self.wind['slon'])
                    self.torn['dist'],self.torn['bear'],self.torn['online']=np.vectorize(self.ClosestPoint)(self.torn['slat'],self.torn['slon'],self.torn['elat'],self.torn['elon'],lat,lon)

                    
                    if self.time_type == 'central':
                        self.Time_to_Zulu()

                    self.torn.reset_index(inplace=True,drop=True)
                    self.torn = self.torn[(self.torn.dist<=dist)&(self.torn.online == True)][['yr', 'mo', 'dy','hr','time','mag','dist','bear','slat', 'slon', 'elat', 'elon','YEAR','MO','DAY','HR','geometry']]
                    self.torn.reset_index(inplace=True,drop=True)

                    self.hail.reset_index(inplace=True,drop=True)
                    self.hail = self.hail[self.hail.dist<=dist][['yr', 'mo', 'dy','hr','time','mag','dist','bear','slat', 'slon', 'elat', 'elon','YEAR','MO','DAY','HR','geometry']]
                    self.hail.reset_index(inplace=True,drop=True)

                    self.wind.reset_index(inplace=True,drop=True)
                    self.wind = self.wind[self.wind.dist<=dist][['yr', 'mo', 'dy','hr','time','mag','dist','bear','slat', 'slon', 'elat', 'elon','YEAR','MO','DAY','HR','geometry']]
                    self.wind.reset_index(inplace=True,drop=True)

                
                if self.time_type == 'central':
                    self.Time_to_Zulu()


        
    def Time_to_Zulu(self):
        
        if self.time_type == 'central':
            #print('Converting Time to Zulu')
            print('Date/Time is in CST')
            self.torn['datetime'] = pd.to_datetime({'year':self.torn['yr'].values,'month':self.torn['mo'].values,'day':self.torn['dy'].values,'hour':pd.to_datetime(self.torn['time']).dt.hour})
            self.hail['datetime'] = pd.to_datetime({'year':self.hail['yr'].values,'month':self.hail['mo'].values,'day':self.hail['dy'].values,'hour':pd.to_datetime(self.hail['time']).dt.hour})
            self.wind['datetime'] = pd.to_datetime({'year':self.wind['yr'].values,'month':self.wind['mo'].values,'day':self.wind['dy'].values,'hour':pd.to_datetime(self.wind['time']).dt.hour})
            
            # to zulu turned off
            # self.torn.datetime = self.torn.datetime+pd.to_timedelta(6,unit='h')
            # self.hail.datetime = self.hail.datetime+pd.to_timedelta(6,unit='h')
            # self.wind.datetime = self.wind.datetime+pd.to_timedelta(6,unit='h')

            self.torn['yr']= self.torn['datetime'].dt.year
            self.torn['mo']= self.torn['datetime'].dt.month
            self.torn['dy']= self.torn['datetime'].dt.day
            self.torn['hr']= self.torn['datetime'].dt.hour

            self.hail['yr']= self.hail['datetime'].dt.year
            self.hail['mo']= self.hail['datetime'].dt.month
            self.hail['dy']= self.hail['datetime'].dt.day
            self.hail['hr']= self.hail['datetime'].dt.hour
 
            self.wind['yr']= self.wind['datetime'].dt.year
            self.wind['mo']= self.wind['datetime'].dt.month
            self.wind['dy']= self.wind['datetime'].dt.day
            self.wind['hr']= self.wind['datetime'].dt.hour
            
            self.time_type = 'zulu'
            to = ['YEAR','MO','DAY','HR']
            fm = ['yr','mo','dy','hr']

            for i in range(4):
                self.torn[to[i]] = self.torn[fm[i]]
                self.hail[to[i]] = self.hail[fm[i]]
                self.wind[to[i]] = self.wind[fm[i]]
                

        else:
            #print('Data already converted to Zulu')
            print('Date/Time is in CST')

    

    def ClosestPoint(self,lat1,lon1,lat2,lon2,clat,clon):
        if (lat1==lat2) & (lon1==lon2):
            return 999,999,False
        k = ((lon2-lon1)*(clat-lat1)-(lat2-lat1)*(clon-lon1))/((lon2-lon1)**2 + (lat2-lat1)**2)
        lat4 = clat-k*(lon2-lon1)
        lon4 = clon+k*(lat2-lat1)
        lat=False
        lon=False
        online=False
        if (lat4 >= min(lat1,lat2))&(lat4 <= max(lat1,lat2)):
            lat=True
        if (lon4 >= min(lon1,lon2))&(lon4 <= max(lon1,lon2)):
            lon=True
        if (lat==True)&(lon==True):
            online = True
            
        dist,bear = MetFuncs_Raw.DistBear(clat,clon,lat4,lon4)
        return dist,bear,online        

    
    def process(self,hail_path,torn_path,wind_path):
        self.torn = pd.read_csv(torn_path)
        self.torn = gpd.GeoDataFrame(self.torn)
        self.torn['spoint'] = self.torn.apply(lambda row: Point([row['slon'], row['slat']]), axis=1)
        self.torn['epoint'] = self.torn.apply(lambda row: Point([row['elon'], row['elat']]), axis=1)
        self.torn['geometry'] = self.torn.apply(lambda row: LineString([row['spoint'], row['epoint']]), axis=1)
        self.hail = pd.read_csv(hail_path) 
        self.hail = gpd.GeoDataFrame(self.hail, geometry=gpd.points_from_xy(self.hail.slon, self.hail.slat))
        self.wind = pd.read_csv(wind_path)
        self.wind = gpd.GeoDataFrame(self.wind, geometry=gpd.points_from_xy(self.wind.slon, self.wind.slat))
        self.Time_to_Zulu()
        print('SPC Data Loaded')


    # def Plot_MH_BE(self,**kwargs):
    #     if kwargs:
    #         if 'show' in kwargs.keys():
    #             show = kwargs['show']
    #         else:
    #             show = False
    #         if 'transparent' in kwargs.keys():
    #             transparent = kwargs['transparent']
    #         else:
    #             transparent = False
            
    #         if 'site_name' in kwargs.keys():
    #             site_name=kwargs['site_name']
    #         else:
    #             site_name=False
    #     else:
    #         show = False
    #         site_name=False
    #         transparent = False

    #     if not(site_name):
    #         if self.net_plat:
    #             site_name= self.data.station_info.loc[0,'STNNAME']
    #         else:
    #             site_name= input('What is the Site name:  ')

        

        
    #     HeatMap_Plotter.HEATmap(1, 'Tornado Counts for\n'+site_name,## charts , Title
    #                             'PoR: '+str(self.torn.yr.min())+'-'+str(self.torn.yr.max())+'\n Produced '+ dt.date.today().strftime('%d %b %Y') +' by 14WS/CXO DSN: 552-8814', ## info
    #                             'SPC_Torn_Counts_{}_{}'.format(str(self.torn.yr.min()),str(self.torn.yr.max())), ## outputname
    #                             df1=self.torn.groupby(['yr','mo','dlocationy','hr'])[['mag']].max().groupby(['mo','hr'])[['mag']].count().unstack(1),df1_round_to=0,##df1_units=(add_unit,unit), ## df1 kwargs ## uni-code ('\u2109','\u2103','\u00b0')
    #                             transparent=transparent,show=show,landscape=True)##,font_size=6,tick_font_size=7) ## chart kwargs

    #     HeatMap_Plotter.HEATmap(1, 'Hail Counts for\n'+site_name,## charts , Title
    #                             'PoR: '+str(self.hail.yr.min())+'-'+str(self.hail.yr.max())+'\n Produced '+ dt.date.today().strftime('%d %b %Y') +' by 14WS/CXO DSN: 552-8814', ## info
    #                             'SPC_Hail_Counts_{}_{}'.format(str(self.hail.yr.min()),str(self.hail.yr.max())), ## outputname
    #                             df1=self.hail.groupby(['yr','mo','dy','hr'])[['mag']].max().groupby(['mo','hr'])[['mag']].count().unstack(1),df1_round_to=0,##df1_units=(add_unit,unit), ## df1 kwargs ## uni-code ('\u2109','\u2103','\u00b0')
    #                             transparent=transparent,show=show,landscape=True)##,font_size=6,tick_font_size=7) ## chart kwargs

    #     HeatMap_Plotter.HEATmap(1, 'SVR Wind Counts for\n'+site_name,## charts , Title
    #                             'PoR: '+str(self.wind.yr.min())+'-'+str(self.wind.yr.max())+'\n Produced '+ dt.date.today().strftime('%d %b %Y') +' by 14WS/CXO DSN: 552-8814', ## info
    #                             'SPC_SVR_Wind_Counts_{}_{}'.format(str(self.wind.yr.min()),str(self.wind.yr.max())), ## outputname
    #                             df1=self.wind.groupby(['yr','mo','dy','hr'])[['mag']].max().groupby(['mo','hr'])[['mag']].count().unstack(1),df1_round_to=0,##df1_units=(add_unit,unit), ## df1 kwargs ## uni-code ('\u2109','\u2103','\u00b0')
    #                             transparent=transparent,show=show,landscape=True)##,font_size=6,tick_font_size=7) ## chart kwargs


    def Update_CSVs(self):
        print('Attenmpting to update CSVs')
        year = int(dt.date.today().strftime('%Y'))
        torn_add = 'https://www.spc.noaa.gov/wcm/data/1950-{}_torn.csv.zip'
        hail_add = 'https://www.spc.noaa.gov/wcm/data/1955-{}_hail.csv.zip'
        wind_add = 'https://www.spc.noaa.gov/wcm/data/1955-{}_wind.csv.zip'

        trys = 0
        tyear = year
        while trys < 5:
            if os.system('wget '+torn_add.format(str(tyear))+' -P '+self.location) != 2048:
                print('Got Tornado Data 1950 to {}'.format(str(tyear)))
                os.system('unzip '+self.location+'1950-{}_torn.csv.zip'.format(str(tyear)))
                os.system('rm '+self.location+'1950-{}_torn.csv.zip'.format(str(tyear)))
                os.system('mv 1950-{}_torn.csv '.format(str(tyear))+self.location+'torn.csv')
                break
            else:
                tyear -= 1
            trys +=1

        trys = 0
        tyear = year
        while trys < 5:
            if os.system('wget '+hail_add.format(str(tyear))+' -P '+self.location) != 2048:
                print('Got Hail Data 1955 to {}'.format(str(tyear)))
                os.system('unzip '+self.location+'1955-{}_hail.csv.zip'.format(str(tyear)))
                os.system('rm '+self.location+'1955-{}_hail.csv.zip'.format(str(tyear)))
                os.system('mv 1955-{}_hail.csv '.format(str(tyear))+self.location+'hail.csv')
                break
            else:
                tyear -= 1
            trys +=1

        trys = 0
        tyear = year
        while trys < 5:
            if os.system('wget '+wind_add.format(str(tyear))+' -P '+self.location) != 2048:
                print('Got Wind Data 1955 to {}'.format(str(tyear)))
                os.system('unzip '+self.location+'1955-{}_wind.csv.zip'.format(str(tyear)))
                os.system('rm '+self.location+'1955-{}_wind.csv.zip'.format(str(tyear)))
                os.system('mv 1955-{}_wind.csv '.format(str(tyear))+self.location+'wind.csv')
                break
            else:
                tyear -= 1
            trys +=1
            
            
        



