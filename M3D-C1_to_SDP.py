#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 17:16:56 2025
     M3D-C1 to SDP compaitible NetCDF format loader
     
     Needs :
         Te, Ne:  assumed 2D
         rr, zz: assumed to be 1D?
         bb : assumed 2D
         Note: fluctuation 
@author: rianc
"""


# Engaging cluster library connections
from sys import path;path.append('/orcd/home/002/rianc/')
try:
    import C1py
    import fio_py
except:pass


import xarray as xr
    
def get_fields_from_C1(filename, n) :
    
    # Get necessary variables from file .h5
    psi,B1,B0,R,Z,rmagx,zmagx,PsiLCFS,  ne, ne0, te, te0 = \
         get_fields_from_C1(C1file_name,n,False)
    
    # Save in SDP format
    ds = xr.Dataset(data_vars = dict(\
             rr=(['r'],R), zz=(['z'],Z), ne=(['r','z'],ne), \
             te=(['r','z'],te), bb=(['r','z'], B1), ) )
    ds0 = xr.Dataset(data_vars = dict(\
             rr=(['r'],R), zz=(['z'],Z), ne=(['r','z'],ne0), \
             te=(['r','z'],te0), bb=(['r','z'], B0), ) )
    
    ds0.to_netcdf('C1.h5_equ_ufile.cdf')
    ds.to_netcdf('C1.h5_{0:0>4}_ufile.cdf'.format(1))
    
def get_fields_from_C1(filename,n,phi=0,slice=0):
    # Slice acts as timepoint
    # Get psi grid, B grid, R,Z coords
    print('Reading: %s'%filename)
    
    # Get Psi
    psi = C1py.read_field('psi', slice=-1, filename=filename, points=200,
                         rrange=None, zrange=None, iequil=1)
    
    
    b_field = C1py.read_field('bfield', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=None,idiff=False)
    b_field0 = C1py.read_field('bfield', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=True,idiff=False)
    
    ne = C1py.read_field('ne', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=None,idiff=False)
    ne0 = C1py.read_field('bfneield', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=True,idiff=False)
    
    te = C1py.read_field('te', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=None,idiff=False)
    te0 = C1py.read_field('te', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=True,idiff=False)
 
    R=b_field.coords['R'].values;Z=b_field.coords['Z']

    # C1 file data
    isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,filename)
    psi_lcfs = fio_py.eval_series(fio_py.get_series(isrc, fio_py.FIO_LCFS_PSI), 0.)
    rmagx = fio_py.eval_series(fio_py.get_series(isrc, fio_py.FIO_MAGAXIS_R),0)
    zmagx = fio_py.eval_series(fio_py.get_series(isrc, fio_py.FIO_MAGAXIS_Z),0)

    
    return psi.data,b_field.data,b_field0.data,psi.coords['R'].values,\
        psi.coords['Z'].values, rmagx, zmagx, psi_lcfs, ne, ne0, te, te0