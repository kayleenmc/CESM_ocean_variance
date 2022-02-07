# calculate MHT by basin
import netCDF4 as nc
import xarray
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sig
from scipy import stats
import pandas as pd
from mpl_toolkits.basemap import Basemap
from statsmodels.tsa.stattools import acf 
from calculate_MHT import *
from calculate_V import *
import pickle
import timeit
from os import listdir
from os.path import isfile, join

#load data
data=pd.read_csv("/gpfs_common/share01/slarson/ktmcmoni/research/ocn_variance/CESM_ocean_variance/analysis/Indian.csv")
tmp=pd.DataFrame.to_numpy(data) # 74 x 2. column 0 is PW, column 1 is latitude
#tmp[21:24,:]=np.nan
tmp2=np.concatenate((tmp[0:21,:],tmp[25:75,:]))
tmp=tmp2

data=pd.read_csv("/gpfs_common/share01/slarson/ktmcmoni/research/ocn_variance/CESM_ocean_variance/analysis/trenberth_IO_variance.csv")
IO_var=pd.DataFrame.to_numpy(data) # 74 x 2. column 0 is PW, column 1 is latitude

f=nc.Dataset('/gpfs_backup/slarson_data/ktmcmoni/CCSM4/products/MHT_IO_ts.nc')
lat_IO=f['latitude'][:]
MHT_time=f['time_B'][:] # this is last 200 years
MHT_time_MD=f['time_MD'][:]
MHT_time_MDeq=f['time_MDeq'][:]
MHT_IO_B=f['MHT_IO_B'][:]
MHT_IO_MD=f['MHT_IO_MD'][:]
MHT_IO_MDeq=f['MHT_IO_MDeq'][:]
f.close()

MHT_B=MHT_IO_B/100
MHT_MD=MHT_IO_MD/100
MHT_MDeq=MHT_IO_MDeq/100 # because I didn't convert cm/s to m/s before!

f=nc.Dataset('/gpfs_backup/slarson_data/ktmcmoni/CCSM4/products/MHT_IO_ts_anom.nc')
MHT_IO_B_anom=f['MHT_IO_B'][:]
MHT_IO_MD_anom=f['MHT_IO_MD'][:]
MHT_IO_MDeq_anom=f['MHT_IO_MDeq'][:]
f.close()

MHT_B_anom=MHT_IO_B_anom/100
MHT_MD_anom=MHT_IO_MD_anom/100
MHT_MDeq_anom=MHT_IO_MDeq_anom/100 # because I didn't convert cm/s to m/s before!

dpi=100
plt.rc('font',size=36)

fig=plt.figure(num=3,figsize=(2500/dpi,1100/dpi),dpi=dpi)
plt.subplot(1,2,1)
plt.rc('font',size=30)
plt.plot(np.nanmean(MHT_B,axis=0)/(10**15),np.arange(-35,20),'k',label='CTRL')
plt.plot(np.nanmean(MHT_MDeq,axis=0)/(10**15),np.arange(-35,20),'b',label='NoENSO')
# add in observationally based estimates
plt.scatter([-0.87,-0.99,-1.33,-1.1,-1.5,-0.97],[-32,-32,-32,-32,-32,-34],marker='o',c='0.2',label='Obs (direct)')
plt.plot(tmp[:,0],tmp[:,1],'0.6',label='Obs (T&Z19)')
plt.legend(fontsize=20)
plt.ylabel('Latitude',fontsize=36)
plt.xlabel('MHT (PW)',fontsize=36)
plt.ylim(-36,20)
plt.xlim(-2.6,0)
plt.text(-2.5,-34,'a)')
#plt.title('Time mean Indian Ocean MHT',fontsize=30)
#plt.savefig('../analysis/figures/MHT_Indian_mean2.png')

#fig=plt.figure(num=4,figsize=(2500/dpi,2500/dpi),dpi=dpi)
plt.subplot(1,2,2)
plt.plot(np.nanvar(MHT_B_anom/(10**15),axis=0),np.arange(-35,20),'k',label='CTRL')
plt.plot(np.nanvar(MHT_MDeq_anom/(10**15),axis=0),np.arange(-35,20),'b',label='NoENSO')
plt.plot(IO_var[:,0],IO_var[:,1],'0.6',label='Obs (T&Z19)')
#plt.legend(fontsize=20)
plt.xlabel('MHT variance $\mathregular{(PW^2)}$',fontsize=36)
#plt.ylabel('Latitude',fontsize=36)
plt.xticks(np.arange(0.0,0.24,step=0.05))
plt.ylim(-36,20)
plt.xlim(0.0,0.24)
plt.text(0.005,-34,'b)')
#plt.title('Variance of Indian Ocean MHT',fontsize=30)
plt.savefig('../analysis/figures/MHT_Indian_variance2.png')



plt.show()
