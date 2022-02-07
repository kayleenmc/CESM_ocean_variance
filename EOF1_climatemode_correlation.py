# look at correlation between EOF1 and various climate modes

import numpy as np
import matplotlib.pyplot as plt
import xarray
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from scipy import stats
import pickle
from statsmodels.tsa.stattools import acf

f=nc.Dataset('/gpfs_backup/slarson_data/ktmcmoni/CCSM4/products/MHT_IO_ts_anom.nc')
lat_IO=f['latitude'][:]
MHT_time=f['time_B'][:] # this is last 110 years
MHT_time_MD=f['time_MD'][:]
MHT_time_MDeq=f['time_MDeq'][:]
MHT_IO_B=f['MHT_IO_B'][:]
MHT_IO_MD=f['MHT_IO_MD'][:]
MHT_IO_MDeq=f['MHT_IO_MDeq'][:]
f.close()

# load SST
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/sst.jan0685-dec0990.anom.nc')
SST_MDeq=f['SST'][:]
lat_SST=f['lat'][:]
lon_SST=f['lon'][:]
f.close()

SST_MDeq=SST_MDeq[-3672:,...] # think size should already be ok? 

MHT_IO_B=MHT_IO_B/(100*10**15)
MHT_IO_MDeq=MHT_IO_MDeq/(100*10**15)
MHT_IO_MD=MHT_IO_MD/(100*10**15)

MHT_MDeq_anom=MHT_IO_MDeq-MHT_IO_MDeq.mean(axis=0)
u,d,vt=np.linalg.svd(MHT_MDeq_anom.T,full_matrices=False)

vecs_svd=vt.T
tvecs_svd=np.dot(u,np.diag(d))

# see how much variance is explained by each mode
norm_sq_d=(d**2).sum()
norm_sq_MHT=(MHT_IO_B**2).sum()

# extract leading EOFs, renormalize
PC1_MDeq=vecs_svd[:,0]/np.std(vecs_svd[:,0])
U1_MDeq=tvecs_svd[:,0]*np.std(vecs_svd[:,0])

# PC1_MDeq is our PC time series
SIOD=np.mean(np.mean(SST_MDeq[:,56:67,44:52],axis=2),axis=1)-np.mean(np.mean(SST_MDeq[:,66:76,72:80],axis=2),axis=1) # difference in SSTA averaged over 55-65E,37-27S and 90-100E, 28-18S
#tmp_SST=SST_MDeq # SSTA averaged over 22-16S, 102-108E and 32-16S, 108-115E
#tmp_SST[:,np.where(lat_SST>-16),:]=0
#tmp_SST[:,:,np.where(lon_SST>115)]=0
#tmp_SST[:,np.where(lat_SST<-32),:]=0
#tmp_SST[:,:,np.where(lon_SST<102)]=0
#tmp_SST[:,np.where(lat_SST<-22&lat_SST>-36),np.where(lon_SST<108)]=0

#ningaloo_nino=
IOB=np.mean(np.mean(SST_MDeq[:,74:117,32:80],axis=2),axis=1)

# make into 2 boxes, average them
NNino1=np.mean(np.mean(SST_MDeq[:,72:79,82:92],axis=2),axis=1)
NNino2=np.mean(np.mean(SST_MDeq[:,57:72,87:92],axis=2),axis=1)
NNino=np.mean(np.column_stack((NNino1,NNino2)),axis=1)

# try more accurate way
NNino_SST=SST_MDeq[:,57:80,81:93]
NNino_lat=lat_SST[57:80]
NNino_lon=lon_SST[81:93]
# "zoomed in" to box we want. shape: 3672, 23, 12

# at every longitude, keep the values at the latitudes we want
NNino_index=np.zeros([3672,])
a=[]
for k in range(0,3672):
    for i in range(12):
        for j in range(23):
            if NNino_lon[i]<108:
                if NNino_lat[j]>-22:
                    a.append(NNino_SST[k,j,i])
            elif NNino_lon[i]<115:
                if NNino_lat[j]>-36:
                    a.append(NNino_SST[k,j,i])
    NNino_index[k]=np.mean(a)
    a=[]

corr_IOB=np.zeros([48,2])
corr_SIOD=np.zeros([48,2])
corr_NNino=np.zeros([48,2])
for i in range(0,48):
    corr_IOB[i,:]=stats.pearsonr(IOB[-3672+24:-24],PC1_MDeq[-3672+i:-48+i])
    corr_SIOD[i,:]=stats.pearsonr(SIOD[-3672+24:-24],PC1_MDeq[-3672+i:-48+i])
    corr_NNino[i,:]=stats.pearsonr(NNino_index[-3672+24:-24],PC1_MDeq[-3672+i:-48+i])
# plot lead/lag correlations of PC1 with these 2 modes

# assess significance
for i,val in enumerate(acf(NNino_index)):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break

dof_NNino=(len(NNino_index)*1)/(2*Te_B)
# 612

for i,val in enumerate(acf(SIOD)):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break

dof_SIOD=(len(SIOD)*1)/(2*Te_B)
# 367

for i,val in enumerate(acf(PC1_MDeq)):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break

dof_PC1=(len(PC1_MDeq)*1)/(2*Te_B)
# 1836

dof_NNino_PC1=612
z_r=np.log((1+corr_NNino[24,0])/(1-corr_NNino[24,0]))/2 # first index is lag
SE_z_r=1/np.sqrt(dof_NNino_PC1-3)
z=z_r/(SE_z_r)
p_NNino=stats.norm.sf(np.abs(z))*2

dof_SIOD_PC1=367
z_r=np.log((1+corr_SIOD[24,0])/(1-corr_SIOD[24,0]))/2
SE_z_r=1/np.sqrt(dof_SIOD_PC1-3)
z=z_r/(SE_z_r)
p_SIOD=stats.norm.sf(np.abs(z))*2

plt.rc('font',size=24)
lags=np.arange(-24,24)

dpi=100
plt.figure(num=1,figsize=(2500/dpi,1100/dpi),dpi=dpi)
plt.plot(lags,corr_IOB[:,0])
max_lag=lags[np.argmax(corr_IOB[:,0])]
c=np.amax(corr_IOB[:,0])
plt.text(-20,0.08,'max corr at lag: ' + str(max_lag))
plt.text(-20,0.05,'max corr: ' + str(round(c,3)))
#plt.text(-100,-0.1,'DOF: '+str(dof))
plt.xlim(-24,24)
plt.xlabel('Lag (months)')
plt.ylabel('Correlation')
plt.title('PC1 correlation with IOB')
plt.savefig('../analysis/figures/PC1_IOB_MDeq.png')

plt.figure(num=2,figsize=(2500/dpi,1100/dpi),dpi=dpi)
plt.plot(lags,corr_SIOD[:,0])
max_lag=lags[np.argmin(corr_SIOD[:,0])]
c=np.amin(corr_SIOD[:,0])
plt.text(-20,0.05,'max corr at lag: ' + str(max_lag))
plt.text(-20,0.02,'corr: ' + str(round(c,3)))
#plt.text(-100,-0.1,'DOF: '+str(dof))
plt.xlim(-24,24)
plt.xlabel('Lag (months)')
plt.ylabel('Correlation')
plt.title('PC1 correlation with SIOD')
plt.savefig('../analysis/figures/PC1_SIOD_MDeq.png')

plt.figure(num=3,figsize=(2500/dpi,1100/dpi),dpi=dpi)
plt.plot(lags,corr_NNino[:,0])
max_lag=lags[np.argmax(corr_NNino[:,0])]
c=np.amax(corr_NNino[:,0])
plt.text(-10,0.15,'max corr at lag: ' + str(max_lag))
plt.text(-10,0.05,'max corr: ' + str(round(c,3)))
#plt.text(-100,-0.1,'DOF: '+str(dof))
plt.xlim(-24,24)
plt.xlabel('Lag (months)')
plt.ylabel('Correlation')
plt.title('PC1 correlation with Ningaloo Nino')
plt.savefig('../analysis/figures/PC1_NNino_MDeq.png')


plt.show()
