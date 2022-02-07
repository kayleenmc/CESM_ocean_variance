# plot the wind stress associated with ASON IOD events, in FC and NoENSO
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

# get wind stress
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/TAU_pop/tau.jan0400-dec0721.anom.nc')
TAUX_B=f['TAUX'][:] # 3864, 192, 288
TAUY_B=f['TAUY'][:]
time_TAU=f['time'][:]
lat_TAU=f['lat'][:]
lon_TAU=f['lon'][:]
f.close()

TAUX_B=TAUX_B[-3672:,...]
TAUY_B=TAUY_B[-3672:,...]

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/TAU_pop/tau.jan0685-dec0990.anom.nc')
TAUX_MDeq=f['TAUX'][:]
TAUY_MDeq=f['TAUY'][:]
f.close()

# get SST
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/sst.jan0400-dec0721.anom.nc')
SST_B=f['SST'][:] # 3864, 192, 288
time_SST=f['time'][:]
lat_SST=f['lat'][:]
lon_SST=f['lon'][:]
f.close()

SST_B=SST_B[-3672:,...]

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/sst.jan0685-dec0990.anom.nc')
SST_MDeq=f['SST'][:]
lat_SST=f['lat'][:]
lon_SST=f['lon'][:]
f.close()

# composites of IOD events that occur in FC and in NoENSO
# identify IOD events

TAUX_B_pos=np.mean(TAUX_B[IOD_B>thresh,:,:],axis=0)
TAUY_B_pos=np.mean(TAUY_B[IOD_B>thresh,:,:],axis=0)
SST_B_pos=np.mean(SST_B[IOD_B>thresh,:,:],axis=0)
TAUX_MDeq_pos=np.mean(TAUX_MDeq[IOD_MDeq>thresh,:,:],axis=0)
TAUY_MDeq_pos=np.mean(TAUY_MDeq[IOD_MDeq>thresh,:,:],axis=0)
SST_MDeq_pos=np.mean(SST_MDeq[IOD_MDeq>thresh,:,:],axis=0)

dpi=300
plt.rc('font',size=24)
[xx,yy]=np.meshgrid(lon_SST,lat_SST)

# full year
fig=plt.figure(num=1,figsize=(2500/dpi,1100/dpi),dpi=dpi)
plt.subplot(1,3,1)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.fillcontinents('k',lake_color='k')
map.contourf(xx,yy,SST_B_pos,levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both')
c=plt.colorbar(shrink=0.25)
#c.set_label('SST $\mathregular{(^o C)}$')
map.quiver(xx[::8,::8],yy[::8,::8],TAUX_B_pos[::8,::8],TAUY_B_pos[::8,::8],scale=1)
plt.title('a) FC',loc='left')

plt.subplot(1,3,2)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.fillcontinents('k',lake_color='k')
map.contourf(xx,yy,SST_MDeq_pos,levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both')
c=plt.colorbar(shrink=0.25)
#c.set_label('SST $\mathregular{(^o C)}$')
map.quiver(xx[::8,::8],yy[::8,::8],TAUX_MDeq_pos[::8,::8],TAUY_MDeq_pos[::8,::8],scale=1)
plt.title('b) NoENSO',loc='left')

plt.subplot(1,3,3)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.fillcontinents('k',lake_color='k')
map.contourf(xx,yy,SST_B_pos-SST_MDeq_pos,levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both')
c=plt.colorbar(shrink=0.25)
c.set_label('SST $\mathregular{(^o C)}$')
map.quiver(xx[::8,::8],yy[::8,::8],TAUX_B_pos[::8,::8]-TAUX_MDeq_pos[::8,::8],TAUY_B_pos[::8,::8]-TAUY_MDeq_pos[::8,::8],scale=1)
plt.title('c) FC-NoENSO',loc='left')
plt.savefig('../analysis/figures/IOD_SST_tau.png')

# do a regression with all months instead
IOD_SST_FC=np.zeros([192,288,5])
IOD_TAUX_FC=np.zeros([192,288,5])
IOD_TAUY_FC=np.zeros([192,288,5])
IOD_SST_MDeq=np.zeros([192,288,5])
IOD_TAUX_MDeq=np.zeros([192,288,5])
IOD_TAUY_MDeq=np.zeros([192,288,5])

for i in range(0,192):
    for j in range(0,288):
        IOD_SST_FC[i,j,:]=stats.linregress(IOD_B,SST_B[:,i,j])
        IOD_TAUX_FC[i,j,:]=stats.linregress(IOD_B,TAUX_B[:,i,j])
        IOD_TAUY_FC[i,j,:]=stats.linregress(IOD_B,TAUY_B[:,i,j])
        IOD_SST_MDeq[i,j,:]=stats.linregress(IOD_MDeq,SST_MDeq[:,i,j])
        IOD_TAUX_MDeq[i,j,:]=stats.linregress(IOD_MDeq,TAUX_MDeq[:,i,j])
        IOD_TAUY_MDeq[i,j,:]=stats.linregress(IOD_MDeq,TAUY_MDeq[:,i,j])

fig=plt.figure(num=2,figsize=(2500/100,1100/100),dpi=dpi)
plt.subplot(3,1,1)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,0])
map.fillcontinents('lightgrey',lake_color='lightgrey',zorder=1)
map.contourf(xx,yy,IOD_SST_FC[:,:,0],levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both',zorder=0)
#c=plt.colorbar(shrink=0.25)
#c.set_label('SST $\mathregular{(^o C)}$')
Q=map.quiver(xx[::8,::8],yy[::8,::8],IOD_TAUX_FC[::8,::8,0],IOD_TAUY_FC[::8,::8,0],scale=2,linewidth=8,headaxislength=6,headlength=6,zorder=10)
#plt.quiverkey(Q,0.9,0.9,1,'m s-1',coordinates='figure')
plt.title('a) CTRL',loc='left')
plt.xticks([])

plt.subplot(3,1,2)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,0])
map.fillcontinents('lightgrey',lake_color='lightgrey',zorder=1)
map.contourf(xx,yy,IOD_SST_MDeq[:,:,0],levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both',zorder=0)
#c=plt.colorbar(shrink=0.25)
#c.set_label('SST $\mathregular{(^o C)}$')
Q=map.quiver(xx[::8,::8],yy[::8,::8],IOD_TAUX_MDeq[::8,::8,0],IOD_TAUY_MDeq[::8,::8,0],scale=2,linewidth=8,headaxislength=6,headlength=6,zorder=10)
#plt.quiverkey(Q,0.9,0.9,1,'m s-1',coordinates='figure')
plt.title('b) NoENSO',loc='left')
plt.xticks([])

plt.subplot(3,1,3)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,0])
map.fillcontinents('lightgrey',lake_color='lightgrey',zorder=1)
m=map.contourf(xx,yy,IOD_SST_FC[:,:,0]-IOD_SST_MDeq[:,:,0],levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both',zorder=0)
#cax=plt.axes([0.9,0.1,0.075,0.8])
#c=plt.colorbar(m,shrink=0.25)
#c.set_label('SST $\mathregular{(^o C)}$')
Q=map.quiver(xx[::8,::8],yy[::8,::8],IOD_TAUX_FC[::8,::8,0]-IOD_TAUX_MDeq[::8,::8,0],IOD_TAUY_FC[::8,::8,0]-IOD_TAUY_MDeq[::8,::8,0],scale=2,linewidth=8,headaxislength=6,headlength=6,zorder=10)
plt.title('c) CTRL-NoENSO',loc='left')
plt.quiverkey(Q,0.52,0.9,0.5,'0.5 dyne $\mathregular{cm^2}$',coordinates='figure')

fig.subplots_adjust(right=0.8)
cbar_ax=fig.add_axes([0.61,0.1,0.01,0.8])
c=fig.colorbar(m,shrink=0.25,cax=cbar_ax)
c.set_label('SST $\mathregular{(^o C)}$')
plt.savefig('../analysis/figures/IOD_SST_tau_reg.eps')

# just ASO: 7,8,9
TAUX_B=np.concatenate((TAUX_B[7::12,:,:],TAUX_B[8::12,:,:],TAUX_B[9::12,:,:]))
TAUY_B=np.concatenate((TAUY_B[7::12,:,:],TAUY_B[8::12,:,:],TAUY_B[9::12,:,:]))
SST_B=np.concatenate((SST_B[7::12,:,:],SST_B[8::12,:,:],SST_B[9::12,:,:]))

TAUX_MDeq=np.concatenate((TAUX_MDeq[7::12,:,:],TAUX_MDeq[8::12,:,:],TAUX_MDeq[9::12,:,:]))
TAUY_MDeq=np.concatenate((TAUY_MDeq[7::12,:,:],TAUY_MDeq[8::12,:,:],TAUY_MDeq[9::12,:,:]))
SST_MDeq=np.concatenate((SST_MDeq[7::12,:,:],SST_MDeq[8::12,:,:],SST_MDeq[9::12,:,:]))

IOD_B=np.concatenate((IOD_B[7::12],IOD_B[8::12],IOD_B[9::12]))
IOD_MDeq=np.concatenate((IOD_MDeq[7::12],IOD_MDeq[8::12],IOD_MDeq[9::12]))

TAUX_B_pos=np.mean(TAUX_B[IOD_B>thresh,:,:],axis=0)
TAUY_B_pos=np.mean(TAUY_B[IOD_B>thresh,:,:],axis=0)
SST_B_pos=np.mean(SST_B[IOD_B>thresh,:,:],axis=0)
TAUX_MDeq_pos=np.mean(TAUX_MDeq[IOD_MDeq>thresh,:,:],axis=0)
TAUY_MDeq_pos=np.mean(TAUY_MDeq[IOD_MDeq>thresh,:,:],axis=0)
SST_MDeq_pos=np.mean(SST_MDeq[IOD_MDeq>thresh,:,:],axis=0)

fig=plt.figure(num=3,figsize=(2500/dpi,1100/dpi),dpi=dpi)
plt.subplot(1,3,1)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.fillcontinents('k',lake_color='k')
map.contourf(xx,yy,SST_B_pos,levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both')
c=plt.colorbar(shrink=0.25)
#c.set_label('SST $\mathregular{(^o C)}$')
map.quiver(xx[::8,::8],yy[::8,::8],TAUX_B_pos[::8,::8],TAUY_B_pos[::8,::8],scale=1)
plt.title('a) FC',loc='left')

plt.subplot(1,3,2)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.fillcontinents('k',lake_color='k')
map.contourf(xx,yy,SST_MDeq_pos,levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both')
c=plt.colorbar(shrink=0.25)
#c.set_label('SST $\mathregular{(^o C)}$')
map.quiver(xx[::8,::8],yy[::8,::8],TAUX_MDeq_pos[::8,::8],TAUY_MDeq_pos[::8,::8],scale=1)
plt.title('b) FC',loc='left')

plt.subplot(1,3,3)
map=Basemap(projection='cyl',llcrnrlat=-40.,urcrnrlat=30.,resolution='c',llcrnrlon=20.,urcrnrlon=180.)
map.drawcoastlines()
parallels=np.arange(-60,90,30)
meridians=np.arange(0,360,60.)
map.drawparallels(parallels,labels=[1,0,0,0])
map.drawmeridians(meridians,labels=[0,0,0,1])
map.fillcontinents('k',lake_color='k')
map.contourf(xx,yy,SST_B_pos-SST_MDeq_pos,levels=np.linspace(-0.75,0.75,31),cmap='bwr',extend='both')
c=plt.colorbar(shrink=0.25)
c.set_label('SST $\mathregular{(^o C)}$')
map.quiver(xx[::8,::8],yy[::8,::8],TAUX_B_pos[::8,::8]-TAUX_MDeq_pos[::8,::8],TAUY_B_pos[::8,::8]-TAUY_MDeq_pos[::8,::8],scale=1)
plt.title('c) FC-NoENSO',loc='left')
plt.savefig('../analysis/figures/IOD_SST_tau_ASO.png')


plt.show()









