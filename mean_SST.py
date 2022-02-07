# plot the time mean SST from ERSST and our 3 model experiments

import netCDF4 as nc
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sig
from scipy import stats
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import cartopy.mpl.ticker as cticker

f=xr.open_dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/sst.jan0400-dec0721.nc')
SST_B=f['SST'].mean('time',keep_attrs=True)

f1=xr.open_dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/sst.jan0685-dec0990.nc')
SST_MDeq=f1['SST'].mean('time',keep_attrs=True)

# MD FILE DOESN'T EXIST!
f2=xr.open_dataset('/gpfs_backup/slarson_data/ktmcmoni/CCSM4/ccsm4_MD_0.9x1.25_gx1v6/postP/SST/sst.jan0700-dec1302.nc')
SST_MD=f2['SST'].mean('time',keep_attrs=True)

# ERSST data
f3=xr.open_dataset('/gpfs_backup/slarson_data/datasets/Obs/ersst.v5/ersst.v5.190001-202012.nc')
ERSST=f3['sst'].mean('time',keep_attrs=True)

# get rid of level dimension
ERSST=xr.DataArray.squeeze(ERSST)

l=np.linspace(10,30,21) # set levels for all plots
s=0.5 # how much to shrink colorbar

dpi=100
fig=plt.figure(figsize=(2500/dpi,1100/dpi),dpi=dpi)
plt.rc('font',size=20)
#gs1=gridspec.GridSpec(2,2)
#gs1.update(wspace=0.01,hspace=0.12)
#ax=fig.add_subplot(gs1[0],projection=ccrs.PlateCarree(central_longitude=180))
ax1=fig.add_axes([0.1,0.5,0.4,0.4],projection=ccrs.PlateCarree(central_longitude=180))
ax1.patch.set_alpha(0.0)
ax1.set_extent([25,300,-35,35],crs=ccrs.PlateCarree())
ERSST.plot.contourf(ax=ax1,levels=l,transform=ccrs.PlateCarree(),cbar_kwargs={'label': SST_B.units,'shrink': s},cmap='coolwarm',extend='both')
ax1.set_title('')
ax1.coastlines()
ax1.set_xticks([60,120,180,-120,-60],crs=ccrs.PlateCarree())
ax1.set_xticklabels([60,120,180,-120,-60])
ax1.set_yticks([-20,0,20],crs=ccrs.PlateCarree())
ax1.set_yticklabels([-20,0,20])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')
ax1.set_ylabel('')
ax1.set_xlabel('')
plt.text(-180,40,'ERSST')

#ax=fig.add_subplot(gs1[1],projection=ccrs.PlateCarree(central_longitude=180))
ax2=fig.add_axes([0.1,0.1,0.4,0.4],projection=ccrs.PlateCarree(central_longitude=180))
ax2.set_extent([25,300,-35,35],crs=ccrs.PlateCarree())
SST_B.plot.contourf(ax=ax2,levels=l,transform=ccrs.PlateCarree(),cbar_kwargs={'label': SST_B.units,'shrink': s},cmap='coolwarm',extend='both')
ax2.coastlines()
ax2.patch.set_alpha(0.0)
ax2.set_xticks([60,120,180,-120,-60],crs=ccrs.PlateCarree())
ax2.set_xticklabels([60,120,180,-120,-60])
ax2.set_yticks([-20,0,20],crs=ccrs.PlateCarree())#
ax2.set_yticklabels([-20,0,20])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')
ax2.set_ylabel('')
ax2.set_xlabel('')
plt.text(-180,40,'FC')

#ax=fig.add_subplot(gs1[2],projection=ccrs.PlateCarree(central_longitude=180))
ax3=fig.add_axes([0.55,0.5,0.4,0.4],projection=ccrs.PlateCarree(central_longitude=180))
ax3.set_extent([25,300,-35,35],crs=ccrs.PlateCarree())
SST_MDeq.plot.contourf(ax=ax3,levels=l,transform=ccrs.PlateCarree(),cbar_kwargs={'label': SST_B.units,'shrink': s},cmap='coolwarm',extend='both')
ax3.coastlines()
ax3.patch.set_alpha(0.0)
ax3.set_xticks([60,120,180,-120,-60],crs=ccrs.PlateCarree())
ax3.set_xticklabels([60,120,180,-120,-60])
ax3.set_yticks([-20,0,20],crs=ccrs.PlateCarree())
ax3.set_yticklabels([-20,0,20])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')
ax3.set_ylabel('')
ax3.set_xlabel('')
plt.text(-180,40,'NoENSO')

#ax=fig.add_subplot(gs1[3],projection=ccrs.PlateCarree(central_longitude=180))
ax4=fig.add_axes([0.55,0.1,0.4,0.4],projection=ccrs.PlateCarree(central_longitude=180))
ax4.set_extent([25,300,-35,35],crs=ccrs.PlateCarree())
SST_MD.plot.contourf(ax=ax4,levels=l,transform=ccrs.PlateCarree(),cbar_kwargs={'label': SST_B.units,'shrink': s},cmap='coolwarm',extend='both')
ax4.coastlines()
ax4.patch.set_alpha(0.0)
ax4.set_xticks([60,120,180,-120,-60],crs=ccrs.PlateCarree())
ax4.set_xticklabels([60,120,180,-120,-60])
ax4.set_yticks([-20,0,20],crs=ccrs.PlateCarree())
ax4.set_yticklabels([-20,0,20])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')
ax4.set_ylabel('')
ax4.set_xlabel('')
plt.text(-180,40,'MD')

plt.savefig('../analysis/figures/mean_SST.png',dpi=dpi,bbox_inches='tight')
plt.show()
