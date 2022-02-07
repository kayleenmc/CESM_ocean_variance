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

m=20 # set marker size to show significance

f=nc.Dataset('/gpfs_backup/slarson_data/ktmcmoni/CCSM4/products/MHT_IO_ts_anom_fullMD.nc')
lat_IO=f['latitude'][:]
MHT_time=f['time_B'][:] # this is last 110 years
MHT_IO_B=f['MHT_IO_B'][:]
MHT_IO_MD=f['MHT_IO_MD'][:]
MHT_IO_MDeq=f['MHT_IO_MDeq'][:]
f.close()

MHT_IO_B=MHT_IO_B/100
MHT_IO_MD=MHT_IO_MD/100
MHT_IO_MDeq=MHT_IO_MDeq/100 # because I forgot to change from cm/s to m/s when integrating

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0400-dec0721.anom.nc')
IOD_B=f['SST'][:]
f.close()

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0685-dec0990.anom.nc')
IOD_MDeq=f['SST'][:]
f.close()

f=nc.Dataset('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/nino34.jan0400-dec0721.anom.nc')
ENSO_B=f['SST'][:]
ENSO_B_time=f['time'][:]
f.close()

IOD_IO_B=np.zeros([55,2])
IOD_IO_MDeq=np.zeros([55,2])
for i in range(0,55):
    IOD_IO_B[i,:]=stats.pearsonr(IOD_B[-2400:,0,0],MHT_IO_B[-2400:,i])
    IOD_IO_MDeq[i,:]=stats.pearsonr(IOD_MDeq[-2400:,0,0],MHT_IO_MDeq[-2400:,i])

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
    for i,val in enumerate(acf(MHT_IO_MDeq[-2400:,k])):
        if val < 1/np.exp(1):
            Te_B[k]=i
            val_Te_B=val
            break
    dof_MHT_MDeq[k]=(len(MHT_IO_MDeq[-2400:,k])*1)/(2*Te_B[k])
    dof_MDeq=np.minimum(dof_MHT_MDeq,dof_MDeq_IOD)
    z_r=np.log((1+IOD_IO_MDeq[k,0])/(1-IOD_IO_MDeq[k,0]))/2
    SE_z_r=1/np.sqrt(dof_MDeq-3)
    z=z_r/(SE_z_r) # See Field (2009) pp.172 for discussion on z score
    t_r=IOD_IO_MDeq[k,0]*np.sqrt(dof_MDeq-2)/np.sqrt(1-(IOD_IO_MDeq[k,0]**2)) # convert to t value
    p_MDeq[k]=stats.norm.sf(np.abs(z))*2 # this is 2 sided

sig_MDeq=IOD_IO_MDeq[p_MDeq<0.05,0]
sig_MDeq_lat=lat_IO[p_MDeq<0.05]


p_B=np.zeros(55)
dof_MHT_B=np.zeros(55)
Te_B=np.zeros(55)

for k in range(0,55):
    for i,val in enumerate(acf(MHT_IO_B[-2400:,k])):
        if val < 1/np.exp(1):
            Te_B[k]=i
            val_Te_B=val
            break
    dof_MHT_B[k]=(len(MHT_IO_B[-2400:,k])*1)/(2*Te_B[k])
    dof_B=np.minimum(dof_MHT_B,dof_B_IOD)
    z_r=np.log((1+IOD_IO_B[k,0])/(1-IOD_IO_B[k,0]))/2
    SE_z_r=1/np.sqrt(dof_B-3)
    z=z_r/(SE_z_r) # See Field (2009) pp.172 for discussion on z score
    t_r=IOD_IO_B[k,0]*np.sqrt(dof_B-2)/np.sqrt(1-(IOD_IO_B[k,0]**2)) # convert to t value
    p_B[k]=stats.norm.sf(np.abs(z))*2 # this is 2 sided

sig_B=IOD_IO_B[p_B<0.05,0]
sig_B_lat=lat_IO[p_B<0.05]

dpi=100


# now calculate the same thing for each season!!! ASO is the season of interest
ASO=np.zeros(600)
count=0
for i in range(0,200):
    ASO[count]=-2400+7+(12*i)
    ASO[count+1]=-2400+8+(12*i)
    ASO[count+2]=-2400+9+(12*i)
    count=count+3

ASO2=np.zeros(600)
count=0
for i in range(0,200):
    ASO2[count]=0+7+(12*i)
    ASO2[count+1]=0+8+(12*i)
    ASO2[count+2]=0+9+(12*i)
    count=count+3

IOD_IO_B_ASO=np.zeros([55,2])
IOD_IO_MDeq_ASO=np.zeros([55,2])
for i in range(0,55):
    IOD_IO_B_ASO[i,:]=stats.pearsonr(IOD_B[ASO.astype(int),0,0],MHT_IO_B[ASO.astype(int),i])
    IOD_IO_MDeq_ASO[i,:]=stats.pearsonr(IOD_MDeq[ASO.astype(int),0,0],MHT_IO_MDeq[ASO.astype(int),i])

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
        IOD_IO_B[i,j,:]=stats.pearsonr(IOD_B[-2376:-24,0,0],MHT_IO_B[-2400+j:-48+j,i])
        IOD_IO_MDeq[i,j,:]=stats.pearsonr(IOD_MDeq[-2376:-24,0,0],MHT_IO_MDeq[-2400+j:-48+j,i])

ENSO_IO_B=np.zeros([55,48,2])
for i in range(0,55):
    for j in range(0,48):
        ENSO_IO_B[i,j,:]=stats.pearsonr(ENSO_B[-2376:-24,0,0],MHT_IO_B[-2400+j:-48+j,i])

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

tmp=np.zeros([55,48])
tmp[ENSO_IO_B_sig<0.05]=1

plt.rc('font',size=24)
fig=plt.figure(num=3,figsize=(2500/dpi,2500/dpi),dpi=dpi)
plt.subplot(3,1,1)
#plt.contourf(xx,yy,ENSO_IO_B[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr')
plt.contourf(xx,yy,tmp,colors='k',alpha=0.5)
plt.colorbar()
plt.xlabel('Lag (months)')
plt.ylabel('Latitude')
plt.text(-23,-30,'FC/ENSO')
plt.text(-12,-30,'MHT leads')
plt.text(12,-30,'ENSO leads')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(-24,24,step=6))
plt.plot([0,0],[-35,20],'w:')
plt.ylim(-35,19)
plt.xlim(-24,23)
#plt.title('ENSO, IOD, and MHT correlations')

plt.subplot(3,1,2)
plt.contourf(xx,yy,IOD_IO_B[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr')
plt.colorbar()
plt.xlabel('Lag (months)')
plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation')
plt.text(-23,-30,'FC/IOD')
plt.text(-12,-30,'MHT leads')
plt.text(12,-30,'IOD leads')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(-24,24,step=6))
plt.plot([0,0],[-35,20],'w:')
plt.ylim(-35,19)
plt.xlim(-24,23)

plt.subplot(3,1,3)
plt.contourf(xx,yy,IOD_IO_MDeq[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr')
plt.colorbar()
plt.xlabel('Lag (months)')
plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation, NoENSO')
plt.text(-23,-30,'NoENSO/IOD')
plt.text(-9,-30,'MHT leads')
plt.text(12,-30,'IOD leads')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(-24,24,step=6))
plt.plot([0,0],[-35,19],'w:')
plt.ylim(-35,20)
plt.xlim(-24,23)

plt.savefig('../analysis/figures/IOD_lags_by_lat3.png')



# zero lag correlation by month and by latitude
# first, calculate zero lag correlation for each month
IOD_IO_B_month=np.zeros([55,12,2])
IOD_IO_MDeq_month=np.zeros([55,12,2])
ENSO_IO_B_month=np.zeros([55,12,2])

for i in range(0,55):
    for k in range(0,12):
        ENSO_B_month=ENSO_B[-2400+k::12,0,0]
        IOD_B_month=IOD_B[-2400+k::12,0,0]
        IOD_MDeq_month=IOD_MDeq[-2400+k::12,0,0]
        MHT_IO_B_month=MHT_IO_B[-2400+k::12,...]
        MHT_IO_MDeq_month=MHT_IO_MDeq[-2400+k::12,...]
        ENSO_IO_B_month[i,k,:]=stats.pearsonr(ENSO_B_month,MHT_IO_B_month[:,i])
        IOD_IO_B_month[i,k,:]=stats.pearsonr(IOD_B_month,MHT_IO_B_month[:,i])
        IOD_IO_MDeq_month[i,k,:]=stats.pearsonr(IOD_MDeq_month,MHT_IO_MDeq_month[:,i])

[xx,yy]=np.meshgrid(np.arange(0,12),lat_IO)
labels=['J','F','M','A','M','J','J','A','S','O','N','D']

fig=plt.figure(num=5,figsize=(2500/dpi,2500/dpi),dpi=dpi)
ax=plt.subplot(3,1,1)
plt.contourf(xx,yy,ENSO_IO_B_month[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr',extend='both')
plt.colorbar()
plt.xlabel('Month')
plt.ylabel('Latitude')
plt.text(1,-30,'FC/ENSO')
plt.xticks(np.arange(0,12))
plt.yticks(np.arange(-30,20,step=10))
plt.ylim(-35,19)
plt.xlim(0,11)
ax.set_xticklabels(labels)
#plt.title('ENSO, IOD, and MHT correlations')

ax=plt.subplot(3,1,2)
plt.contourf(xx,yy,IOD_IO_B_month[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr',extend='both')
plt.colorbar()
plt.xlabel('Month')
plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation')
plt.text(1,-30,'FC/IOD')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(0,12))
plt.ylim(-35,19)
plt.xlim(0,11)
ax.set_xticklabels(labels)

ax=plt.subplot(3,1,3)
plt.contourf(xx,yy,IOD_IO_MDeq_month[:,:,0],levels=np.linspace(-0.6,0.6,13),cmap='bwr',extend='both')
plt.colorbar()
plt.xlabel('Month')
plt.ylabel('Latitude')
#plt.title('IOD and MHT correlation, NoENSO')
plt.text(1,-30,'NoENSO/IOD')
plt.yticks(np.arange(-30,20,step=10))
plt.xticks(np.arange(0,12))
plt.ylim(-35,19)
plt.xlim(0,11)
ax.set_xticklabels(labels)

plt.savefig('../analysis/figures/IOD_lags_by_month3.png')


plt.show()
