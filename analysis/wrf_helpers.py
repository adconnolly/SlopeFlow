import numpy as np
import xarray as xr
from scipy.interpolate import griddata

def get_zgrid_XZcontour(ds, itime=-1, ktop=68, iy=74):

    z_w = ( ds['PHB'][itime,:ktop+1,iy] + ds['PH'][itime,:ktop+1,iy] ) / 9.81
    # z_w_pad = np.concatenate( [z_w.values, z_w.values[:, 0:1]], axis=1)
    
    z_cent = 0.5 * ( z_w[:-1].values + z_w[1:].values )
    z_cent_pad = np.concatenate([ z_cent, z_cent[:, 0:1]], axis=1)
    
    z_u_noBC = 0.5 * ( z_cent_pad[:,:-1] + z_cent_pad[:,1:] )
    z_u = np.concatenate([ z_u_noBC[:, -1:], z_u_noBC ], axis=1)

    return z_cent, z_u, z_w.values

def get_zgrid_meanProfiles(ds, itime=-1, ktop=12, ix=265):

    z_w = ( ds['PHB'][itime,:ktop+1,:,ix-1:ix+1] + ds['PH'][itime,:ktop+1,:,ix-1:ix+1] ).mean(dim="south_north") / 9.81
    z_w = z_w - ds['HGT'][itime,0,ix-1:ix+1]
    
    z_cent = 0.5 * ( z_w[:-1].values + z_w[1:].values )

    z_u = 0.5 * ( z_cent[:,0] + z_cent[:,1] )

    return z_cent, z_u

def get_vel_XZquiver(u_stag, w_stag, x_u=None, x_w=None, z_u=None, z_w=None, X=None, Z=None, mask=None):
    # Destaggers and returns as DataArray if no grid info
    # With grid, will interpolate to new X-Z grid from x-z
    #   useful to feed regular grid for streamplot()

    if (x_u==None).any():
        u = 0.5 * ( u_stag[:, :-1] + u_stag[:, 1:] ).rename({'west_east_stag': 'west_east'})
        w = 0.5 * ( w_stag[:-1, :] + w_stag[1:, :] ).rename({'bottom_top_stag': 'bottom_top'})
        
        return u.assign_coords(west_east=w.coords['west_east']), w.assign_coords(bottom_top=u.coords['bottom_top'])
    
    else:
        u_points = np.column_stack( (x_u.ravel(), z_u.ravel()) )
        w_points = np.column_stack( (x_w.ravel(), z_w.ravel()) )
        
        U = griddata(u_points, u_stag.values.ravel(), (X, Z), method='linear')
        W = griddata(w_points, w_stag.values.ravel(), (X, Z), method='linear')

        U[mask] = np.nan
        W[mask] = np.nan
        
        return U, W

def get_XZgrid_XZquiver(dx = 1, nx = 10, ztop=1000, nz=10):
    # Needs to be regular grid for streamplot()
    dz=ztop/nz
    xv = np.arange(dx/2, (nx-1)*dx, dx)
    zv = np.arange(dz/2, ztop, dz)
    
    return np.meshgrid(xv, zv)