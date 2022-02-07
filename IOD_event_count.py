# count IOD events
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from scipy import stats

thresh=0.3 # threshold for defining IOD events

# first, get IOD time series
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0400-dec0721.anom.nc')
IOD_B=f['SST'][:] # 3864, 1, 1
f.close()

IOD_B=IOD_B[-3672:,0,0]

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0685-dec0990.anom.nc')
IOD_MDeq=f['SST'][:] # 3672, 1, 1
f.close()

IOD_MDeq=IOD_MDeq[:,0,0]

# take ASO average of index
# start in jan (0) 7:9 is first ASO
IOD_ASO_B=np.zeros(306)
IOD_ASO_MDeq=np.zeros(306)

for i in range(0,306):
    IOD_ASO_B[i]=np.mean(IOD_B[(i*12)+7:(i*12)+9])
    IOD_ASO_MDeq[i]=np.mean(IOD_MDeq[(i*12)+7:(i*12)+9])

# how many times does index cross threshold?
events_B=np.sum(IOD_ASO_B>thresh) #110 for thresh=0.3, 66 for thresh=0.6
events_MDeq=np.sum(IOD_ASO_MDeq>thresh) # 91, 38



