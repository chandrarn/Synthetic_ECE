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


from header_SDP import C1py, fio_py, np, xr

def save_C1_SDP_format(filename='/nobackup1/wenhaw42/Linear/01_n1_test_cases/1000_bate1.0_constbz_0_cp0501/C1.h5', n=1) :
    
    # Get necessary variables from file .h5
    psi,B1,B0,R,Z,rmagx,zmagx,PsiLCFS,  ne, ne0, te, te0 = \
         get_fields_from_C1(filename,n,False)
    
    # convert B(R,phi,Z) to ||B||
    B1 = np.sqrt( B1[0]**2 + B1[1]**2 + B1[2]**2 )
    B0 = np.sqrt( B0[0]**2 + B0[1]**2 + B0[2]**2 )
    # Save in SDP format
    ds = xr.Dataset(data_vars = dict(\
             rr=(['r'],R), zz=(['z'],Z), ne=(['r','z'],ne.data.T), \
             te=(['r','z'],te.data.T), bb=(['r','z'], B1.T), ) )
    ds0 = xr.Dataset(data_vars = dict(\
             rr=(['r'],R), zz=(['z'],Z), ne=(['r','z'],ne0.data.T), \
             te=(['r','z'],te0.data.T), bb=(['r','z'], B0.T), ) )
    
    ds0.to_netcdf('C1.h5_equ_ufile.cdf')
    ds.to_netcdf('C1.h5_{0:0>4}_ufile.cdf'.format(1))

######################################################################3
def get_fields_from_C1(filename,n,phi=0,slice=0,debug=True):
    # Slice acts as timepoint
    # Get psi grid, B grid, R,Z coords
    print('Reading: %s'%filename)
    
    # Get Psi
    psi = C1py.read_field('psi', slice=-1, filename=filename, points=200,
                         rrange=None, zrange=None, iequil=1)
    
    
    b_field = C1py.read_field('bfield', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=None,idiff=False)
    b_field0 = C1py.read_field('bfield', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=None,idiff=True)
    # print(b_field.data.shape,b_field.coords['R'].values.shape)
    # raise SyntaxError


    ne = C1py.read_field('ne', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi+1, iequil=None,idiff=False)
    ne0 = C1py.read_field('ne', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=None,idiff=False)
    
    te = C1py.read_field('te', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi+1, iequil=None,idiff=False)
    te0 = C1py.read_field('te', slice=slice, filename=filename, points=200,
                        rrange=None, zrange=None,phi=phi, iequil=None,idiff=False)
 
    R=b_field.coords['R'].values;Z=b_field.coords['Z']
    
    if debug:
        b_field.to_netcdf('b_field.cdf')
        b_field0.to_netcdf('b_field_0.cdf')

    # C1 file data
    isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,filename)
    psi_lcfs = fio_py.eval_series(fio_py.get_series(isrc, fio_py.FIO_LCFS_PSI), 0.)
    rmagx = fio_py.eval_series(fio_py.get_series(isrc, fio_py.FIO_MAGAXIS_R),0)
    zmagx = fio_py.eval_series(fio_py.get_series(isrc, fio_py.FIO_MAGAXIS_Z),0)

    
    return psi.data,b_field.data,b_field0.data,psi.coords['R'].values,\
        psi.coords['Z'].values, rmagx, zmagx, psi_lcfs, ne, ne0, te, te0

###################################################
if __name__ == '__main__': save_C1_SDP_format()
