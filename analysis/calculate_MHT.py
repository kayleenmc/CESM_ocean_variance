# calculate MHT from 4D potential temperature and velocity (U or V)
##################################################################

import numpy as np
import pop_tools
from geopy.distance import geodesic

##################################################################
def compute_vertical_integrated_MHT(z, theta, velocity, z_stop: int = None, depth_axis: int = 0) -> float:
    '''
    Integrates the product of rho*Cp*Theta*velocity given numpy arrays for depth (z, 1-dimension), theta, and velocity
    (same size). Units are in W per m, such that when integrated horizontally and divided by 10^15  it yields PW.

    If z_stop is given, the integration stops at z_stop.
    '''

    # raise errors if inputs are incorrect size
    if (z.ndim != 1):
        raise ValueError("z not a rank-1 numpy array")
    if (theta.ndim < 1):
        raise ValueError("theta is not at least a rank-1 numpy array")
    if (velocity.ndim < 1):
        raise ValueError("velocity is not at least a rank-1 numpy array")
    if (depth_axis >= theta.ndim):
        raise ValueError("theta does not have a specified depth axis")
    if (theta.shape[depth_axis] != z.shape[0]):
        raise ValueError("z and depth axis of theta are different lengths")
    if (velocity.shape[depth_axis] != z.shape[0]):
        raise ValueError("z and depth axis of velocity are different lengths")
    if (velocity.shape != theta.shape):
        raise ValueError("velocity and theta are not the same shape")
    # raise error if z is not monotonically increasing
    if (np.min(z[1:]-z[:-1])<0):
        raise ValueError("z is not monotonically increasing")

    # raise error if z_stop is above surface or below deepest level
    if (z_stop is not None):
        if (z_stop <= 0):
            raise ValueError("z_stop is at or above the surface")
        if (z_stop > z[-1]):
            raise ValueError("z_stop is greater than the maximum depth")
    else:
        z_stop=z[-1]

    # set value for rho * Cp
    rho_Cp=4.093*(100**3) # 4.093 J K-1 cm-3 -> multiply by 100^3 to put in J K-1 m-3
    
    # theta and velocity are representative of z_star, halfway between levels above/below it
    z_star=np.mean([z[1:],z[:-1]],axis=0)
    z_star=np.insert(z_star,0,0.0,axis=0)

    # modify based on z_stop
    zs_idx=np.searchsorted(z_star,z_stop)
    if zs_idx == len(z_star):
        z_star=np.append(z_star,z_stop)
    else:
        z_star[zs_idx]=z_stop
        z_star=z_star[:zs_idx+1]

    # define nested function (don't get this?)
    def column_integrate(col_theta):
        return np.nansum(
            [col_theta[i]*(z_star[i+1]-z_star[i])\
                for i in range(len(z_star)-1)]
         )

    # compute integral
    MHT_vertical=np.apply_along_axis(column_integrate,depth_axis,theta*velocity)*rho_Cp

    return MHT_vertical
#################################################################
def compute_regional_integrated_MHT(MHT_vertical, basin, longitude, latitude, gy) -> float:
    '''
    Integrates the vertically integrated MHT in a basin across a given latitude
    Regions are defined by poptools. options are ['Black Sea', 'Baltic Sea', 'Red Sea', 'Southern Ocean', 'Pacific Ocean',
    'Indian Ocean', 'Persian Gulf', 'Atlantic Ocean', 'Mediterranean Sea',
    'Lab. Sea & Baffin Bay', 'GIN Seas', 'Arctic Ocean', 'Hudson Bay']
    Grid must be gx1v6. 
    gy is the latitude of interest (where we calculate MHT)
    MHT_vertical must be dimensions of time x lat x lon
    '''
    # get regions with pop tools
    grid_name='POP_gx1v6'
    ds=pop_tools.get_grid(grid_name)
    TLAT=ds.TLAT
    TLONG=ds.TLONG
    # raise an error if basin is not an option
    
    # raise an error if MHT_vertical is not right shape
    #if (MHT_vertical.ndim < 2):
    #    raise ValueError("MHT_vertical is not at least rank-2")
    #if (MHT_vertical.shape[2] != longitude.shape[0]):
    #    raise ValueError("MHT_vertical axis 1 is not shape of longitude")
    #if (MHT_vertical.shape[1] != latitude.shape[0]):
    #    raise ValueError("MHT_vertical axis 2 is not shape of latitude") 
    # raise an error if latitude of choice is not in the basin

    # select grid. region will be ones, all else will be zeros. can view regions at: https://pop-tools.readthedocs.io/en/latest/examples/re    gion-mask.html#Alternative-region-masks
    mask3d=pop_tools.region_mask_3d(grid_name,mask_name='Pacific-Indian-Atlantic')
    mask2d=mask3d.sel(region=basin)
    #MHT_vertical_masked=np.multiply(MHT_vertical,mask2d)
    #[xx,yy]=np.meshgrid(longitude,latitude)
    #lon_masked=np.multiply(xx,mask2d)
    
    # if basin was Atlantic, make it so that we don't cut down the basin
    if (basin=='Atlantic Ocean'):
        TLATx=np.zeros([384,320])
        TLATx[:,0:160]=TLAT[:,160:320]
        TLATx[:,160:320]=TLAT[:,0:160]
        longitudex=np.zeros(288)
        longitudex[0:144]=longitude[144:288]
        longitudex[144:288]=longitude[0:144]
        mask2dx=np.zeros([384,320])
        mask2dx[:,0:160]=mask2d[:,160:320]
        mask2dx[:,160:320]=mask2d[:,0:160]
        TLAT=TLATx
        mask2d=mask2dx
        longitude=longitudex

    # find closest index to selected latitude
    idy_MHT=np.searchsorted(latitude,gy)
    idy_pop=np.searchsorted(TLAT[:,10],gy) # index doesn't matter here, TLAT[anything] is the same

    # find max and min longitudes
    lon1=np.amin(np.where(mask2d[idy_pop,:]==1))
    lonend=np.amax(np.where(mask2d[idy_pop,:]==1))

    # find first two longitudes of the basin at that latitude

    idx1=np.searchsorted(longitude,TLONG[idy_pop,lon1])
    idxend=np.searchsorted(longitude,TLONG[idy_pop,lonend])

    # find distance across basin
    p1=(latitude[idy_MHT],TLONG[idy_pop,lon1])
    p2=(latitude[idy_MHT],TLONG[idy_pop,lon1+1]) # need to define longitude somehow
    dx=geodesic(p1,p2).km*1000 # turn into meters instead of km

    # integrate across the latitude of choice. need to use cumsum because of nan's
    #tmp=np.cumsum(MHT_vertical[:,idy_MHT,idx1:idxend],axis=1)
    #MHT=tmp[:,-1]
    MHT_vertical=MHT_vertical*dx
    MHT=np.sum(MHT_vertical[:,idy_MHT,idx1:idxend],axis=1)    
    if (basin=='Global'):
        MHT=np.sum(MHT_vertical[:,idy_MHT,:],axis=1)

    return MHT





