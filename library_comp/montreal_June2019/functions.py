#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:59:45 2017

@author: rasen
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import astropy.io.ascii as pyascii
import os
import sys
from astropy.convolution import convolve, Gaussian1DKernel


def flux_dens(photfile):
    print('hip77911_phot/'+photfile)
    s_data = pyascii.read('hip77911_phot/'+photfile,data_start=1) #table, row column
    header = s_data.colnames
    filtflux = [] ; filtflux_err = []
    for i in np.arange(len(s_data)):
        row = s_data[i]
        filtflux.append(row[5])
        filtflux_err.append(row[6])
        #if row[1] == -1: #flux given, convert Jy to erg/s/cm2/A
    #        factorconv = 3.*10**(-5)/int(row[3])**2
    #       fluxinerg=row[5]*factorconv
    #        filtflux.append(fluxinerg)
    #    else:

    #        filtflux.append(10**(-0.4*row[1])*row[5]*10**(-9))

    return s_data[header[0]],np.array(s_data[header[3]]),np.array(s_data[header[4]]),np.array(filtflux), np.array(filtflux_err)



def trans_curve(file,filter, wvlmod,fluxmod):
    filtdir = '/Users/rasen/Documents/PhD/filt_tc/'
    filename=file[:file.find('.txt')]
    transcurve=pyascii.read(filtdir+filename+'/'+filter+'.dat',guess=0,
    data_start=1) #wvl and transmission of the filter
    headertc=transcurve.colnames
    wvltc = np.array(transcurve[headertc[0]])
    transtc = np.array(transcurve[headertc[1]])

    #apply Filter
    mod_wvl_filt=wvlmod[np.logical_and(wvlmod>min(wvltc), wvlmod<max(wvltc))]
    mod_flux_filt=fluxmod[np.logical_and(wvlmod>min(wvltc), wvlmod<max(wvltc))]

    #interpolate filter to model's resolution
    ynew=interpolate.interp1d(wvltc,transtc,kind='cubic')
    newtrans=ynew(mod_wvl_filt)

    #multiply model flux and trans Filter and divide by its mean trans
    mod_flux_filt*=newtrans

    return mod_wvl_filt,mod_flux_filt,newtrans

def apply_ifs_model(wvl_M,flx_M,ifswvls):

      #wavelngth and flux of the model given and wvls of the ifs boxes in mic
      #Model same resolution as IFS

      d_lamb_M = np.array([x-y for x,y in zip(wvl_M[1:], wvl_M)],dtype=float)
      d_lamb_ifs = np.median(ifswvls/55.) # R = lambda/delta_lambda =50

      #model sampling compared to IFS
      sampl=np.round(np.median(np.divide(d_lamb_ifs,d_lamb_M)))
      gauss_kernel=Gaussian1DKernel(stddev=int(sampl/2.35))
      smooth_mod = convolve(flx_M,gauss_kernel,boundary='extend')
      #plt.figure()
      #plt.plot(wvl_M,flx_M)
      #plt.plot(wvl_M,smooth_mod)
      #plt.show()

      #go through the bins
      ifsmp = [(y-x)/2 for x,y in zip(ifswvls, ifswvls[1:])]
      ifsmp = [ifsmp[0]] + ifsmp + [ifsmp[-1]]
      m=0
      fluxes=[]
      for ifsb in ifswvls:
            thisbinmax = (ifsb+ifsmp[m+1])
            thisbinmin = (ifsb - ifsmp[m])
            #model in this bin
            ind=np.logical_and(wvl_M>thisbinmin,wvl_M<thisbinmax)
            fluxes.append(np.mean(smooth_mod[ind]))
            m+=1

      return fluxes



def G_value(flx,eflx,flx_M):
      #flx is flux of spec and flx_M flux of model to compare to
      #scaling factor

      C_fac =  np.sum(flx*flx_M/eflx**(2.)) / np.sum(flx_M**(2.)/eflx**(2.))

      Gval = np.sum( (flx - C_fac*flx_M)**(2.) / ((eflx)**(2.) ) )

      return Gval, C_fac


def high_low_res(wvl_M,flx_M,eflx_M,sphere_wvl):

      #take wavelength spacing of model spectrum, compare it to the SPHERE_LSS wavelength array.
      fwhmloc_i = (np.abs(wvl_M-sphere_wvl[0])).argmin()
      fwhmloc_f = (np.abs(wvl_M-sphere_wvl[1])).argmin()
      fwhmloc = [fwhmloc_i, fwhmloc_f]
      if max(fwhmloc) <0:
          fwhmloc_i = (np.abs(wvl_M-sphere_wvl[-2])).argmin()
          fwhmloc_f = (np.abs(wvl_M-sphere_wvl[-1])).argmin()
          fwhmloc = [fwhmloc_i, fwhmloc_f]

      fwhm_fr = float(fwhmloc_f - fwhmloc_i)

      return fwhm_fr


def apply_same_res(wvl_M,flx_M,eflx_M,obj_wvl, obj_flx, obj_eflx, fwhm_fr):

      #The model has higher resolution than the object, smooth the model

    if fwhm_fr > 0:
        std = float(fwhm_fr/2.355)

        gauss_kernel = Gaussian1DKernel(stddev= std) #FWHM = 2.355*sigma

        smooth_mod = convolve(flx_M,gauss_kernel,boundary='extend')
        smooth_emod = convolve(eflx_M,gauss_kernel,boundary='extend')

        #interpolate the smoothed model flux onto the SPHERE grid
        ynew = interpolate.interp1d(wvl_M,smooth_mod,kind='cubic', fill_value="extrapolate")
        smth_mod = ynew(obj_wvl)

        #assume very coarsely that S/N ~ sqrt(S) and that locally the S ~ constant with channel
        smooth_emod = smooth_emod/np.sqrt(fwhm_fr)
        ynew = interpolate.interp1d(wvl_M,smooth_emod,kind='cubic', fill_value="extrapolate")
        smth_emod = ynew(obj_wvl)

        #to return
        mod_wvl_out = obj_wvl
        mod_flx_out = smth_mod
        mod_eflx_out = smth_emod
        obj_wvl_out = obj_wvl
        obj_flx_out = obj_flx
        obj_eflx_out = obj_eflx

    else:
        #Smooth the object
        #take wavelength spacing of object, compare it to the model wavelength array.
        fwhmloc_i = (np.abs(obj_wvl-wvl_M[0])).argmin()
        fwhmloc_f = (np.abs(obj_wvl-wvl_M[1])).argmin()
        fwhmloc = [fwhmloc_i, fwhmloc_f]
        if max(fwhmloc) <0:
             fwhmloc_i = (np.abs(obj_wvl-wvl_M[-2])).argmin()
             fwhmloc_f = (np.abs(obj_wvl-wvl_M[-1])).argmin()
             fwhmloc = [fwhmloc_i, fwhmloc_f]

        fwhm_frr = float(fwhmloc_f - fwhmloc_i)

        std = float(fwhm_frr/2.355)

        gauss_kernel = Gaussian1DKernel(stddev= std) #FWHM = 2.355*sigma

        smooth_obj = convolve(obj_flx,gauss_kernel,boundary='extend')
        smooth_eobj = convolve(obj_eflx,gauss_kernel,boundary='extend')

        #interpolate the smoothed obj flux onto the model grid
        ynew = interpolate.interp1d(obj_wvl,smooth_obj,kind='cubic', fill_value="extrapolate")
        smth_obj = ynew(wvl_M)

        #assume very coarsely that S/N ~ sqrt(S) and that locally the S ~ constant with channel
        smooth_eobj = smooth_eobj/np.sqrt(fwhm_frr)
        ynew = interpolate.interp1d(obj_wvl,smooth_eobj,kind='cubic', fill_value="extrapolate")
        smth_eobj = ynew(wvl_M)

        #to return
        mod_wvl_out = wvl_M
        mod_flx_out = flx_M
        mod_eflx_out = eflx_M
        obj_wvl_out = wvl_M
        obj_flx_out = smth_obj
        obj_eflx_out = smth_eobj


    return mod_wvl_out, mod_flx_out, mod_eflx_out, obj_wvl_out, obj_flx_out, obj_eflx_out




def apply_res_sphere_lss(wvl_M,flx_M,eflx_M,sphere_wvl):

      #take wavelength spacing of model spectrum, compare it to the SPHERE_LSS wavelength array.
      fwhmloc_i = (np.abs(wvl_M-sphere_wvl[0])).argmin()
      fwhmloc_f = (np.abs(wvl_M-sphere_wvl[1])).argmin()
      fwhmloc = [fwhmloc_i, fwhmloc_f]
      if max(fwhmloc) <0:
          fwhmloc_i = (np.abs(wvl_M-sphere_wvl[-2])).argmin()
          fwhmloc_f = (np.abs(wvl_M-sphere_wvl[-1])).argmin()
          fwhmloc = [fwhmloc_i, fwhmloc_f]

      fwhm_fr = float(fwhmloc_f - fwhmloc_i)




      std = float(fwhm_fr/2.355)

      gauss_kernel = Gaussian1DKernel(stddev= std) #FWHM = 2.355*sigma

      smooth_mod = convolve(flx_M,gauss_kernel,boundary='extend')
      smooth_emod = convolve(eflx_M,gauss_kernel,boundary='extend')

      #interpolate the smoothed model flux onto the SPHERE grid
      ynew = interpolate.interp1d(wvl_M,smooth_mod,kind='cubic', fill_value="extrapolate")
      sph_f = ynew(sphere_wvl)

      #assume very coarsely that S/N ~ sqrt(S) and that locally the S ~ constant with channel
      smooth_emod = smooth_emod/np.sqrt(fwhm_fr)
      ynew = interpolate.interp1d(wvl_M,smooth_emod,kind='cubic', fill_value="extrapolate")
      sph_ef = ynew(sphere_wvl)

      return sph_f, sph_ef
