# make subplots for the figure showing what does drive MHT in NoENSO

import numpy as np
import matplotlib.pyplot as plt
import xarray
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from scipy import stats
import pickle

f=nc.Dataset('/gpfs_backup/slarson_data/ktmcmoni/CCSM4/products/MHT_IO_ts_anom.nc')
lat_IO=f['latitude'][:]
MHT_time=f['time_B'][:] # this is last 110 years
MHT_time_MD=f['time_MD'][:]
MHT_time_MDeq=f['time_MDeq'][:]
MHT_IO_B=f['MHT_IO_B'][:]
MHT_IO_MD=f['MHT_IO_MD'][:]
MHT_IO_MDeq=f['MHT_IO_MDeq'][:]
f.close()

MHT_IO_B=MHT_IO_B/(100*10**15)
MHT_IO_MDeq=MHT_IO_MDeq/(100*10**15)
MHT_IO_MD=MHT_IO_MD/(100*10**15)

# load SST
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/sst.jan0685-dec0990.anom.nc')
SST_MDeq=f['SST'][:]
lat_SST=f['lat'][:]
lon_SST=f['lon'][:]
f.close()

SST_MDeq=SST_MDeq[-3672:,...] # think size should already be ok? 

# load wind stress
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/TAU_pop/tau.jan0685-dec0990.anom.nc')
TAUX_MDeq=f['TAUX'][:]
TAUY_MDeq=f['TAUY'][:]
f.close()

# load SL
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/PSL/psl.jan0685-dec0990.anom.nc')
SLP_MDeq=f['PSL'][:]
f.close()

# load Z500
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/500mb/upperlevel.500mb.jan0685-dec0990.anom.nc')
Z500_MDeq=f['Z'][:]
f.close()

Z500_MDeq=Z500_MDeq[:,0,...]

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

SLP_MDeq_corr=np.zeros([192,288,5])
SST_MDeq_corr=np.zeros([192,288,5])
TAUX_MDeq_corr=np.zeros([192,288,5])
TAUY_MDeq_corr=np.zeros([192,288,5])
Z500_MDeq_corr=np.zeros([192,288,5])

for i in range(0,192):
    for j in range(0,288):
        SLP_MDeq_corr[i,j,:]=stats.linregress(-1*PC1_MDeq,SLP_MDeq[:,i,j])
        SST_MDeq_corr[i,j,:]=stats.linregress(-1*PC1_MDeq,SST_MDeq[:,i,j])
        TAUX_MDeq_corr[i,j,:]=stats.linregress(-1*PC1_MDeq,TAUX_MDeq[:,i,j])
        TAUY_MDeq_corr[i,j,:]=stats.linregress(-1*PC1_MDeq,TAUY_MDeq[:,i,j])
        Z500_MDeq_corr[i,j,:]=stats.linregress(-1*PC1_MDeq,Z500_MDeq[:,i,j])

# make figures 
plt.rc('font',size=36)

###### Check correlation with IOB
IOB=np.mean(np.mean(SST_MDeq[:,74:117,32:80],axis=2),axis=1)
corr,p=stats.pearsonr(IOB,PC1_MDeq)

[xx,yy]=np.meshgrid(lon_SST,lat_SST)

dpi=300
fig=plt.figure(num=1,figsize=(2500/100,2500/100),dpi=dpi)
plt.subplot(2,1,1)
map=Basemap(projection='cyl',llcrnrlat=-60.,urcrnrlat=65.,resolution='c',llcrnrlon=0.,urcrnrlon=240.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.contourf(xx,yy,SLP_MDeq_corr[:,:,0],levels=np.linspace(-100.,100.,21),cmap='bwr',extend='both')
c=plt.colorbar()
c.set_label('SLP (Pa)')
map.contour(xx,yy,Z500_MDeq_corr[:,:,0],colors='k',linewidth=4.5)
plt.title('c)',loc='left')

plt.subplot(2,1,2)
map=Basemap(projection='cyl',llcrnrlat=-60.,urcrnrlat=65.,resolution='c',llcrnrlon=0.,urcrnrlon=240.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.fillcontinents('lightgrey',lake_color='lightgrey',zorder=1)
map.contourf(xx,yy,SST_MDeq_corr[:,:,0],levels=np.linspace(-0.1,0.1,21),cmap='bwr',extend='both',zorder=0)
c=plt.colorbar()
c.set_label('SST $\mathregular{(^o C)}$')
Q=map.quiver(xx[::8,::8],yy[::8,::8],TAUX_MDeq_corr[::8,::8,0],TAUY_MDeq_corr[::8,::8,0],scale=1,linewidth=3.5,headaxislength=5,zorder=10)
plt.quiverkey(Q,0.6,0.455,0.1,'0.1 dyne $\mathregular{cm^2}$',coordinates='figure')
plt.title('d)',loc='left')
plt.savefig('../analysis/figures/MDeq_reg_subplots.eps')

plt.show()

