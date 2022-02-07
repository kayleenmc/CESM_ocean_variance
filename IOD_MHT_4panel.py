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
import matplotlib.patheffects as PathEffects
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import ListedColormap

m=20 # set marker size to show significance

f=nc.Dataset('/gpfs_backup/slarson_data/ktmcmoni/CCSM4/products/MHT_IO_ts_anom.nc')
lat_IO=f['latitude'][:]
MHT_time=f['time_B'][:] # this is last 110 years
MHT_IO_B=f['MHT_IO_B'][:]
MHT_IO_MD=f['MHT_IO_MD'][:]
MHT_IO_MDeq=f['MHT_IO_MDeq'][:]
f.close()

MHT_IO_B=MHT_IO_B/100
MHT_IO_MD=MHT_IO_MD/100
MHT_IO_MDeq=MHT_IO_MDeq/100 # because I forgot to change from cm/s to m/s when integrating

#f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0400-dec0721.anom.nc')
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0400-dec0721.anom.nc')
IOD_B=f['SST'][:]
f.close()

#f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0685-dec0990.anom.nc')
f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0685-dec0990.anom.nc')
IOD_MDeq=f['SST'][:]
f.close()

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/nino34.jan0400-dec0721.anom.nc')
ENSO_B=f['SST'][:]
ENSO_B_time=f['time'][:]
f.close()

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/nino34.jan0685-dec0990.anom.nc')
ENSO_MDeq=f['SST'][:]
ENSO_MDeq_time=f['time'][:]
f.close()

IOD_IO_B=np.zeros([55,2])
IOD_IO_MDeq=np.zeros([55,2])
for i in range(0,55):
    IOD_IO_B[i,:]=stats.pearsonr(IOD_B[-3672:,0,0],MHT_IO_B[-3672:,i])
    IOD_IO_MDeq[i,:]=stats.pearsonr(IOD_MDeq[-3672:,0,0],MHT_IO_MDeq[-3672:,i])

# get significanc

for i,val in enumerate(acf(IOD_MDeq)):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break
dof_MDeq_IOD=(len(IOD_MDeq)*1)/(2*Te_B)

for i,val in enumerate(acf(IOD_B)):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break
dof_B_IOD=(len(IOD_B)*1)/(2*Te_B)

for i,val in enumerate(acf(ENSO_B)):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break
dof_B_ENSO=(len(ENSO_B)*1)/(2*Te_B)


# now loop over MHT and correlations


p_MDeq=np.zeros(55)
dof_MHT_MDeq=np.zeros(55)
Te_B=np.zeros(55)

for k in range(0,55):
    for i,val in enumerate(acf(MHT_IO_MDeq[-3672:,k])):
        if val < 1/np.exp(1):
            Te_B[k]=i
            val_Te_B=val
            break
    dof_MHT_MDeq[k]=(len(MHT_IO_MDeq[-3672:,k])*1)/(2*Te_B[k])
    dof_MDeq=np.minimum(dof_MHT_MDeq,dof_MDeq_IOD)
    z_r=np.log((1+IOD_IO_MDeq[k,0])/(1-IOD_IO_MDeq[k,0]))/2
    SE_z_r=1/np.sqrt(dof_MDeq-3)
    z=z_r/(SE_z_r) # See Field (2009) pp.172 for discussion on z score
    t_r=IOD_IO_MDeq[k,0]*np.sqrt(dof_MDeq-2)/np.sqrt(1-(IOD_IO_MDeq[k,0]**2)) # convert to t value
#    p_MDeq[k]=stats.norm.sf(np.abs(z))*2 # this is 2 sided

sig_MDeq=IOD_IO_MDeq[p_MDeq<0.05,0]
sig_MDeq_lat=lat_IO[p_MDeq<0.05]


p_B=np.zeros(55)
dof_MHT_B=np.zeros(55)
Te_B=np.zeros(55)

for k in range(0,55):
    for i,val in enumerate(acf(MHT_IO_B[-3672:,k])):
        if val < 1/np.exp(1):
            Te_B[k]=i
            val_Te_B=val
            break
    dof_MHT_B[k]=(len(MHT_IO_B[-3672:,k])*1)/(2*Te_B[k])
    dof_B=np.minimum(dof_MHT_B,dof_B_IOD)
    z_r=np.log((1+IOD_IO_B[k,0])/(1-IOD_IO_B[k,0]))/2
    SE_z_r=1/np.sqrt(dof_B-3)
    z=z_r/(SE_z_r) # See Field (2009) pp.172 for discussion on z score
    t_r=IOD_IO_B[k,0]*np.sqrt(dof_B-2)/np.sqrt(1-(IOD_IO_B[k,0]**2)) # convert to t value
#    p_B[k]=stats.norm.sf(np.abs(z))*2 # this is 2 sided

sig_B=IOD_IO_B[p_B<0.05,0]
sig_B_lat=lat_IO[p_B<0.05]

dpi=300


# now calculate the same thing for each season!!! ASO is the season of interest
ASO=np.zeros(306*3)
count=0
for i in range(0,306):
    ASO[count]=-3672+7+(12*i)
    ASO[count+1]=-3672+8+(12*i)
    ASO[count+2]=-3672+9+(12*i)
    count=count+3

ASO2=np.zeros(306*3)
count=0
for i in range(0,306):
    ASO2[count]=0+7+(12*i)
    ASO2[count+1]=0+8+(12*i)
    ASO2[count+2]=0+9+(12*i)
    count=count+3

IOD_B_ASO=IOD_B[ASO.astype(int),0,0]
IOD_MDeq_ASO=IOD_MDeq[ASO.astype(int),0,0]

IOD_B_season=np.zeros([306])
IOD_MDeq_season=np.zeros([306])
for i in range(0,306):
    IOD_B_season[i]=np.mean(IOD_B_ASO[i:i+3])
    IOD_MDeq_season[i]=np.mean(IOD_MDeq_ASO[i:i+3])

IOD_IO_B_ASO=np.zeros([55,24,2])
IOD_IO_MDeq_ASO=np.zeros([55,24,2])
for i in range(0,55):
    for j in range(0,12):
        IOD_IO_B_ASO[i,j,:]=stats.pearsonr(IOD_B_season,MHT_IO_B[j::12,i])
        IOD_IO_MDeq_ASO[i,j,:]=stats.pearsonr(IOD_MDeq_season,MHT_IO_MDeq[j::12,i])

for i in range(0,55):
    for j in range(12,24):
        IOD_IO_B_ASO[i,j,:]=stats.pearsonr(IOD_B_season[:-1],MHT_IO_B[j::12,i])
        IOD_IO_MDeq_ASO[i,j,:]=stats.pearsonr(IOD_MDeq_season[:-1],MHT_IO_MDeq[j::12,i])


# how to plot this.... ASO correlation with MHT (not just in ASO)
# every year has 1 ASO IOD. Correlate with Jan MHT, feb MHT, etc
# get significance
IOD_IO_B_sig=np.zeros([55,24])
dof_IOD_IO_B=306
for k in range(0,55):
    for j in range(0,24):
        z_r=np.log((1+IOD_IO_B_ASO[k,j,0])/(1-IOD_IO_B_ASO[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_IOD_IO_B-3)
        z=z_r/(SE_z_r)
        t_r=IOD_IO_B_ASO[k,j,0]*np.sqrt(dof_IOD_IO_B-2)/np.sqrt(1-(IOD_IO_B_ASO[k,j,0]**2))
        IOD_IO_B_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_b_sig=np.zeros([55,24])
panel_b_sig[IOD_IO_B_sig<0.05]=1

IOD_IO_MDeq_sig=np.zeros([55,24])
dof_IOD_IO_MDeq=306
for k in range(0,55):
    for j in range(0,24):
        z_r=np.log((1+IOD_IO_MDeq_ASO[k,j,0])/(1-IOD_IO_MDeq_ASO[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_IOD_IO_MDeq-3)
        z=z_r/(SE_z_r)
        t_r=IOD_IO_MDeq_ASO[k,j,0]*np.sqrt(dof_IOD_IO_MDeq-2)/np.sqrt(1-(IOD_IO_MDeq_ASO[k,j,0]**2))
        IOD_IO_MDeq_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_c_sig=np.zeros([55,24])
panel_c_sig[IOD_IO_MDeq_sig<0.05]=1

[xx,yy]=np.meshgrid(np.arange(0,24),lat_IO)
#labels=['J0','F0','M0','A0','M0','J0','J0','A0','S0','O0','N0','D0','J1','F1','M1','A1','M1','J1','J1','A1','S1','O1','N1','D1']
labels=['J0','M0','M0','J0','S0','N0','J1','M1','M1','J1','S1','N1']

plt.rc('font',size=36)
fig=plt.figure(num=1,figsize=(2500/dpi,2500/dpi),dpi=dpi)
ax=plt.subplot(2,1,1)
plt.contourf(xx,yy,IOD_IO_B_ASO[:,:,0],levels=np.linspace(-0.6,0.6,13),extend='both',cmap='bwr')
cbar=plt.colorbar()
cbar.set_label('r')
plt.contourf(xx,yy,panel_b_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.xlabel('Month')
plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation')
plt.text(0.5,15,'a) FC/IOD')
#plt.text(-23,5,'MHT leads')
#plt.text(12,5,'IOD leads')
#plt.yticks(np.arange(-30,20,step=10))
#plt.xticks(np.arange(-24,24,step=6))
#plt.plot([0,0],[-35,20],'w:')
#plt.ylim(-35,19)
#plt.xlim(-24,23)
#plt.text(0.5,-34,'a)')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(0,24,step=2))
plt.ylim(-35,19)
plt.xlim(0,23)
ax.set_xticklabels(labels)

ax=plt.subplot(2,1,2)
plt.contourf(xx,yy,IOD_IO_MDeq_ASO[:,:,0],levels=np.linspace(-0.6,0.6,13),extend='both',cmap='bwr')
cbar=plt.colorbar()
cbar.set_label('r')
plt.contourf(xx,yy,panel_c_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.xlabel('Month')
plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation, NoENSO')
plt.text(0.5,15,'b) NoENSO/IOD')
#plt.text(-23,5,'MHT leads')
#plt.text(12,5,'IOD leads')
#plt.yticks(np.arange(-30,20,step=10))
#plt.xticks(np.arange(-24,24,step=6))
#plt.plot([0,0],[-35,19],'w:')
#plt.ylim(-35,20)
#plt.xlim(-24,23)
#plt.text(0.5,-34,'b)')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(0,24,step=2))
plt.ylim(-35,19)
plt.xlim(0,23)
ax.set_xticklabels(labels)

plt.savefig('../analysis/figures/IOD_lags_ASO.png')



# significance


for i,val in enumerate(acf(IOD_MDeq[ASO.astype(int),0,0])):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break
dof_MDeq_IOD=(len(IOD_MDeq[ASO.astype(int),0,0])*1)/(2*Te_B)

for i,val in enumerate(acf(IOD_B[ASO.astype(int),0,0])):
    if val < 1/np.exp(1):
        Te_B=i
        val_Te_B=val
        break
dof_B_IOD=(len(IOD_B[ASO.astype(int),0,0])*1)/(2*Te_B)

# now loop over MHT and correlations


p_MDeq=np.zeros(55)


for k in range(0,55):
    for i,val in enumerate(acf(MHT_IO_MDeq[ASO.astype(int),k])):
        if val < 1/np.exp(1):
            Te_B=i
            val_Te_B=val
    dof_MHT_MDeq=(len(MHT_IO_MDeq[ASO.astype(int),k])*1)/(2*Te_B)
    dof_MDeq=np.minimum(dof_MHT_MDeq,dof_MDeq_IOD)
    z_r=np.log((1+IOD_IO_MDeq[k,0])/(1-IOD_IO_MDeq[k,0]))/2
    SE_z_r=1/np.sqrt(dof_MDeq-3)
    z=z_r/(SE_z_r) # See Field (2009) pp.172 for discussion on z score
    t_r=IOD_IO_MDeq[k,0]*np.sqrt(dof_MDeq-2)/np.sqrt(1-(IOD_IO_MDeq[k,0]**2)) # convert to t value
    p_MDeq[k]=stats.norm.sf(np.abs(z))*2 # this is 2 sided

sig_MDeq=IOD_IO_MDeq[p_MDeq<0.05,0]
sig_MDeq_lat=lat_IO[p_MDeq<0.05]


p_B=np.zeros(55)

for k in range(0,55):
    for i,val in enumerate(acf(MHT_IO_B[ASO.astype(int),k])):
        if val < 1/np.exp(1):
            Te_B=i
            val_Te_B=val
    dof_MHT_B=(len(MHT_IO_B[ASO.astype(int),k])*1)/(2*Te_B)
    dof_B=np.minimum(dof_MHT_B,dof_B_IOD)
    z_r=np.log((1+IOD_IO_B[k,0])/(1-IOD_IO_B[k,0]))/2
    SE_z_r=1/np.sqrt(dof_B-3)
    z=z_r/(SE_z_r) # See Field (2009) pp.172 for discussion on z score
    t_r=IOD_IO_B[k,0]*np.sqrt(dof_B-2)/np.sqrt(1-(IOD_IO_B[k,0]**2)) # convert to t value
    p_B[k]=stats.norm.sf(np.abs(z))*2 # this is 2 sided

sig_B=IOD_IO_B[p_B<0.05,0]
sig_B_lat=lat_IO[p_B<0.05]

# make a contourf fig like fig 8 of Trenberth & Zhang
# at every latitude, lag by +/- 24 months
IOD_IO_B=np.zeros([55,48,2])
IOD_IO_MDeq=np.zeros([55,48,2])
for i in range(0,55):
    for j in range(0,48):
        IOD_IO_B[i,j,:]=stats.pearsonr(IOD_B[-3672+24:-24,0,0],MHT_IO_B[-3672+j:-48+j,i])
        IOD_IO_MDeq[i,j,:]=stats.pearsonr(IOD_MDeq[-3672+24:-24,0,0],MHT_IO_MDeq[-3672+j:-48+j,i])

ENSO_IO_B=np.zeros([55,48,2])
ENSO_IO_MDeq=np.zeros([55,48,2])
for i in range(0,55):
    for j in range(0,48):
        ENSO_IO_B[i,j,:]=stats.pearsonr(ENSO_B[-3672+24:-24,0,0],MHT_IO_B[-3672+j:-48+j,i])
        ENSO_IO_MDeq[i,j,:]=stats.pearsonr(ENSO_MDeq[-3672+24:-24,0,0],MHT_IO_MDeq[-3672+j:-48+j,i])

[xx,yy]=np.meshgrid(np.arange(-24,24),lat_IO)

# get significance. just use 276, 386, 459 as DOF
ENSO_IO_B_sig=np.zeros([55,48])
dof_ENSO_IO_B=276
for k in range(0,55):
    for j in range(0,48):
        z_r=np.log((1+ENSO_IO_B[k,j,0])/(1-ENSO_IO_B[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_ENSO_IO_B-3)
        z=z_r/(SE_z_r)
        t_r=ENSO_IO_B[k,j,0]*np.sqrt(dof_ENSO_IO_B-2)/np.sqrt(1-(ENSO_IO_B[k,j,0]**2))
        ENSO_IO_B_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_a_sig=np.zeros([55,48])
panel_a_sig[ENSO_IO_B_sig<0.05]=1

ENSO_IO_MDeq_sig=np.zeros([55,48])
dof_ENSO_IO_MDeq=276 # check this!
for k in range(0,55):
    for j in range(0,48):
        z_r=np.log((1+ENSO_IO_MDeq[k,j,0])/(1-ENSO_IO_MDeq[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_ENSO_IO_MDeq-3)
        z=z_r/(SE_z_r)
        t_r=ENSO_IO_MDeq[k,j,0]*np.sqrt(dof_ENSO_IO_MDeq-2)/np.sqrt(1-(ENSO_IO_MDeq[k,j,0]**2))
        ENSO_IO_MDeq_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_d_sig=np.zeros([55,48])
panel_d_sig[ENSO_IO_MDeq_sig<0.05]=1

IOD_IO_B_sig=np.zeros([55,48])
dof_IOD_IO_B=386
for k in range(0,55):
    for j in range(0,48):
        z_r=np.log((1+IOD_IO_B[k,j,0])/(1-IOD_IO_B[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_IOD_IO_B-3)
        z=z_r/(SE_z_r)
        t_r=IOD_IO_B[k,j,0]*np.sqrt(dof_IOD_IO_B-2)/np.sqrt(1-(IOD_IO_B[k,j,0]**2))
        IOD_IO_B_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_b_sig=np.zeros([55,48])
panel_b_sig[IOD_IO_B_sig<0.05]=1

IOD_IO_MDeq_sig=np.zeros([55,48])
dof_IOD_IO_MDeq=459
for k in range(0,55):
    for j in range(0,48):
        z_r=np.log((1+IOD_IO_MDeq[k,j,0])/(1-IOD_IO_MDeq[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_IOD_IO_MDeq-3)
        z=z_r/(SE_z_r)
        t_r=IOD_IO_MDeq[k,j,0]*np.sqrt(dof_IOD_IO_MDeq-2)/np.sqrt(1-(IOD_IO_MDeq[k,j,0]**2))
        IOD_IO_MDeq_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_c_sig=np.zeros([55,48])
panel_c_sig[IOD_IO_MDeq_sig<0.05]=1

plt.rc('font',size=36)
fig=plt.figure(num=3,figsize=(3500/100,2500/100),dpi=dpi)
plt.subplot(2,2,1)
plt.contourf(xx,yy,ENSO_IO_B[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr',extend='both')
#cbar=plt.colorbar()
#cbar.set_label('r')
plt.contourf(xx,yy,panel_a_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.xlabel('Lag (months)')
plt.ylabel('Latitude')
props=dict(boxstyle='square',facecolor='white',alpha=0.8)
plt.text(-23,15,'c) CTRL',bbox=props)
plt.text(-23,-33,'MHT leads',bbox=props)
plt.text(9,-33,'ENSO leads',bbox=props)
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(-24,24,step=6))
plt.plot([0,0],[-35,20],'w:')
plt.ylim(-35,19)
plt.xlim(-24,23)
#plt.text(-23,-34,'a)')
#plt.title('ENSO, IOD, and MHT correlations')

plt.subplot(2,2,2)
plt.contourf(xx,yy,IOD_IO_B[:,:,0],levels=np.linspace(-0.6,0.6,13),extend='both',cmap='bwr')
#cbar=plt.colorbar()
#cbar.set_label('r')
plt.contourf(xx,yy,panel_b_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.xlabel('Lag (months)')
#plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation')
plt.text(-23,15,'d) CTRL',bbox=props)
plt.text(-23,-33,'MHT leads',bbox=props)
plt.text(9,-33,'IOD leads',bbox=props)
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(-24,24,step=6))
plt.plot([0,0],[-35,20],'w:')
plt.ylim(-35,19)
plt.xlim(-24,23)
#plt.text(-23,-34,'b)')

plt.subplot(2,2,3)
plt.contourf(xx,yy,ENSO_IO_MDeq[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr',extend='both')
#cbar=plt.colorbar()
#cbar.set_label('r')
plt.contourf(xx,yy,panel_d_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.xlabel('Lag (months)')
plt.ylabel('Latitude')
plt.text(-23,15,'e) NoENSO',bbox=props)
plt.text(-23,-33,'MHT leads',bbox=props)
plt.text(9,-33,'ENSO leads',bbox=props)
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(-24,24,step=6))
plt.plot([0,0],[-35,20],'w:')
plt.ylim(-35,19)
plt.xlim(-24,23)

plt.subplot(2,2,4)
#plt.contourf(xx,yy,IOD_IO_MDeq[:,:,0],levels=np.linspace(-0.6,0.6,13),extend='both',cmap='bwr')
#cbar=plt.colorbar()
#cbar.set_label('r')
plt.contourf(xx,yy,panel_c_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'],zorder=2)
plt.xlabel('Lag (months)')
plt.contourf(xx,yy,IOD_IO_MDeq[:,:,0],levels=np.linspace(-0.6,0.6,13),extend='both',cmap='bwr',zorder=1)
#plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation, NoENSO')
plt.text(-23,15,'f) NoENSO',bbox=props)
plt.text(-23,-33,'MHT leads',bbox=props)
plt.text(9,-33,'IOD leads',bbox=props)
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(-24,24,step=6))
plt.plot([0,0],[-35,19],'w:')
plt.ylim(-35,19)
plt.xlim(-24,23)
#plt.text(-23,-34,'c)')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
cbar=plt.colorbar(cax=cbar_ax)
cbar.set_label('r')
plt.savefig('../analysis/figures/IOD_lags_by_lat3_4panel.eps')



# zero lag correlation by month and by latitude
# first, calculate zero lag correlation for each month
IOD_IO_B_month=np.zeros([55,12,2])
IOD_IO_MDeq_month=np.zeros([55,12,2])
ENSO_IO_B_month=np.zeros([55,12,2])
ENSO_IO_MDeq_month=np.zeros([55,12,2])

for i in range(0,55):
    for k in range(0,12): #-3672 is for zero lag
        ENSO_B_month=ENSO_B[-3672+k::12,0,0]
        ENSO_MDeq_month=ENSO_MDeq[-3672+k::12,0,0]
        IOD_B_month=IOD_B[-3672+k::12,0,0]
        IOD_MDeq_month=IOD_MDeq[-3672+k::12,0,0]
        #ENSO_B_month=ENSO_B[-3674+k:-2:12,0,0]
        #ENSO_MDeq_month=ENSO_MDeq[-3674+k:-2:12,0,0]
        #IOD_B_month=IOD_B[-3674+k:-2:12,0,0]
        #IOD_MDeq_month=IOD_MDeq[-3674+k:-2:12,0,0]
        MHT_IO_B_month=MHT_IO_B[-3672+k::12,...]
        MHT_IO_MDeq_month=MHT_IO_MDeq[-3672+k::12,...]
        #MHT_IO_B_month=MHT_IO_B[-3674+k:-2:12,...]
        #MHT_IO_MDeq_month=MHT_IO_MDeq[-3674+k:-2:12,...]
        ENSO_IO_B_month[i,k,:]=stats.pearsonr(ENSO_B_month,MHT_IO_B_month[:,i])
        ENSO_IO_MDeq_month[i,k,:]=stats.pearsonr(ENSO_MDeq_month,MHT_IO_MDeq_month[:,i])
        IOD_IO_B_month[i,k,:]=stats.pearsonr(IOD_B_month,MHT_IO_B_month[:,i])
        IOD_IO_MDeq_month[i,k,:]=stats.pearsonr(IOD_MDeq_month,MHT_IO_MDeq_month[:,i])

[xx,yy]=np.meshgrid(np.arange(0,24),lat_IO)
labels=['J','F','M','A','M','J','J','A','S','O','N','D']

# get significance
ENSO_IO_B_sig=np.zeros([55,12])
dof_ENSO_IO_B=276
for k in range(0,55):
    for j in range(0,12):
        z_r=np.log((1+ENSO_IO_B_month[k,j,0])/(1-ENSO_IO_B_month[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_ENSO_IO_B-3)
        z=z_r/(SE_z_r)
        t_r=ENSO_IO_B_month[k,j,0]*np.sqrt(dof_ENSO_IO_B-2)/np.sqrt(1-(ENSO_IO_B_month[k,j,0]**2))
        ENSO_IO_B_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_a_sig=np.zeros([55,12])
panel_a_sig[:]=np.nan
panel_a_sig[ENSO_IO_B_sig<0.05]=1

ENSO_IO_MDeq_sig=np.zeros([55,12])
dof_ENSO_IO_MDeq=276 # CHECK
for k in range(0,55):
    for j in range(0,12):
        z_r=np.log((1+ENSO_IO_MDeq_month[k,j,0])/(1-ENSO_IO_MDeq_month[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_ENSO_IO_MDeq-3)
        z=z_r/(SE_z_r)
        t_r=ENSO_IO_MDeq_month[k,j,0]*np.sqrt(dof_ENSO_IO_MDeq-2)/np.sqrt(1-(ENSO_IO_MDeq_month[k,j,0]**2))
        ENSO_IO_MDeq_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_d_sig=np.zeros([55,12])
panel_d_sig[:]=np.nan
panel_d_sig[ENSO_IO_MDeq_sig<0.05]=1

IOD_IO_B_sig=np.zeros([55,12])
dof_IOD_IO_B=306
for k in range(0,55):
    for j in range(0,12):
        z_r=np.log((1+IOD_IO_B_month[k,j,0])/(1-IOD_IO_B_month[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_IOD_IO_B-3)
        z=z_r/(SE_z_r)
        t_r=IOD_IO_B_month[k,j,0]*np.sqrt(dof_IOD_IO_B-2)/np.sqrt(1-(IOD_IO_B_month[k,j,0]**2))
        IOD_IO_B_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_b_sig=np.zeros([55,12])
panel_b_sig[:]=np.nan
panel_b_sig[IOD_IO_B_sig<0.05]=1

IOD_IO_MDeq_sig=np.zeros([55,12])
dof_IOD_IO_MDeq=306
for k in range(0,55):
    for j in range(0,12):
        z_r=np.log((1+IOD_IO_MDeq_month[k,j,0])/(1-IOD_IO_MDeq_month[k,j,0]))/2
        SE_z_r=1/np.sqrt(dof_IOD_IO_MDeq-3)
        z=z_r/(SE_z_r)
        t_r=IOD_IO_MDeq_month[k,j,0]*np.sqrt(dof_IOD_IO_MDeq-2)/np.sqrt(1-(IOD_IO_MDeq_month[k,j,0]**2))
        IOD_IO_MDeq_sig[k,j]=stats.norm.sf(np.abs(z))*2

panel_c_sig=np.zeros([55,12])
panel_c_sig[:]=np.nan
panel_c_sig[IOD_IO_MDeq_sig<0.05]=1


fig=plt.figure(num=5,figsize=(3500/100,2500/100),dpi=dpi)
ax=plt.subplot(2,2,1)
levels=np.linspace(-0.6,0.6,13)
cmap=plt.get_cmap('bwr')
norm=BoundaryNorm(levels,ncolors=cmap.N,clip=True)
none_map=ListedColormap(['none'])
#levels=[0.,1,1.5]
#norm2=BoundaryNorm(levels,ncolors=3,clip=True)
plt.pcolormesh(xx,yy,np.concatenate((ENSO_IO_B_month[:,:,0],ENSO_IO_B_month[:,:,0]),1),norm=norm,cmap='bwr')
#cbar=plt.colorbar()
#cbar.set_label('r')
plt.pcolor(xx,yy,np.concatenate((panel_a_sig,panel_a_sig),1),cmap=none_map,hatch='.')
plt.xlabel('Month')
plt.ylabel('Latitude')
txt=plt.text(0.5,15,'a) CTRL with ENSO',color='k',fontsize=36,bbox=props)
#txt.set_path_effects([PathEffects.withStroke(linewidth=1,foreground='k')])
#plt.text(0.5,-33,'a)')
plt.xticks(np.arange(0,12))
plt.yticks(np.arange(-30,20,step=10))
plt.ylim(-35,19)
plt.xlim(0,12)
ax.set_xticklabels(labels,ha='left')
#plt.title('ENSO, IOD, and MHT correlations')

ax=plt.subplot(2,2,2)
plt.pcolormesh(xx,yy,np.concatenate((IOD_IO_B_month[:,:,0],IOD_IO_B_month[:,:,0]),1),norm=norm,cmap='bwr')
#cbar=plt.colorbar()
#cbar.set_label('r')
#plt.contourf(xx,yy,panel_b_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.pcolor(xx,yy,np.concatenate((panel_b_sig,panel_b_sig),1),cmap=none_map,hatch='.')
plt.xlabel('Month')
#plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation')
txt=plt.text(0.5,15,'b) CTRL with IOD',color='k',fontsize=36,bbox=props)
#txt.set_path_effects([PathEffects.withStroke(linewidth=1,foreground='k')])
#plt.text(0.5,-33,'b)')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(0,12))
plt.ylim(-35,19)
plt.xlim(0,12)
ax.set_xticklabels(labels,ha='left')

ax=plt.subplot(2,2,3)
plt.pcolormesh(xx,yy,np.concatenate((ENSO_IO_MDeq_month[:,:,0],ENSO_IO_MDeq_month[:,:,0]),1),norm=norm,cmap='bwr')
#cbar=plt.colorbar()
#cbar.set_label('r')
#plt.contourf(xx,yy,panel_d_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.pcolor(xx,yy,np.concatenate((panel_d_sig,panel_d_sig),1),cmap=none_map,hatch='.')
plt.xlabel('Month')
plt.ylabel('Latitude')
txt=plt.text(0.5,15,'c) NoENSO with ENSO',color='k',fontsize=36,bbox=props)
#txt.set_path_effects([PathEffects.withStroke(linewidth=1,foreground='k')])
#plt.text(0.5,-33,'a)')
plt.xticks(np.arange(0,12))
plt.yticks(np.arange(-30,20,step=10))
plt.ylim(-35,19)
plt.xlim(0,12)
ax.set_xticklabels(labels,ha='left')

ax=plt.subplot(2,2,4)
#plt.pcolormesh(xx,yy,IOD_IO_MDeq_month[:,:,0],norm=norm,cmap='bwr')
#cbar=plt.colorbar()
#cbar.set_label('r')
#plt.contourf(xx,yy,panel_c_sig,levels=[0.,0.5,1],colors='none',hatches=[None,'.'])
plt.pcolor(xx,yy,np.concatenate((panel_c_sig,panel_c_sig),1),cmap=none_map,hatch='.',zorder=2)
plt.pcolormesh(xx,yy,np.concatenate((IOD_IO_MDeq_month[:,:,0],IOD_IO_MDeq_month[:,:,0]),1),norm=norm,cmap='bwr',zorder=1)
plt.xlabel('Month')
#plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation, NoENSO')
txt=plt.text(0.5,15,'d) NoENSO with IOD',color='k',fontsize=36,bbox=props)
#txt.set_path_effects([PathEffects.withStroke(linewidth=1,foreground='k')])
#plt.text(0.5,-33,'c)')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(0,12))
plt.ylim(-35,19)
plt.xlim(0,12)
ax.set_xticklabels(labels,ha='left')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
cbar=plt.colorbar(cax=cbar_ax)
cbar.set_label('r')

plt.savefig('../analysis/figures/IOD_lags_by_month.eps')
# lead means ENSO/IOD leads
# lag means MHT leads

plt.show()
