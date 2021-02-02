import numpy as np
from numpy.linalg import inv
import pandas as pd
from scipy import interpolate
from scipy.stats import chisquare
from data_processing import DataProcessing
import sys
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
import numpy as np

class Science:

    '''
    Class container of all science procedures
    '''

    def chi_sq(self, flx, eflx, flx_M, flx_M_err, main_self):
        #flx is flux of object and flx_M flux of model to compare to
        #scaling factor
        C_fac =  np.sum(flx * flx_M / eflx**(2.) ) / np.sum( flx_M**(2.) / eflx**(2.) )

        
        if isinstance(main_self.cov_mat, np.ndarray): #2x2 array 
            ndof=flx.shape[0]
            dimcovar = main_self.cov_mat.shape[0] #assume square arrays
            cmaterr = np.copy(main_self.cov_mat)
            dummy_one = np.array(flx) - C_fac*np.array(flx_M)

            #add the uncorrelated model errors to the covariance matrix
            for j in np.arange(dimcovar):
                cmaterr[j,j]+= (C_fac * flx_M_err[j])**2.

            chi_sq = np.dot(np.dot(np.transpose(dummy_one),inv(cmaterr)),dummy_one)

        else:
            chi_sq = np.sum( (flx - C_fac*flx_M)**(2.) / ( (eflx**2.) + (C_fac*flx_M_err)**2 ) )

        return chi_sq / len(flx) , C_fac


    def overlap(self, object, model):

        ''''
        Function that returns the overlapping data ranges
        of two spectra.
        :param object is a panda dataframe of the first spectrum, with column wavelength named wvl
        :param model is a panda dataframe of the second spectrum, with column wavelength named wvl
        '''
        ini = max(model['wvl'].iloc[0], object['wvl'].iloc[0])
        fin = min(model['wvl'].iloc[-1], object['wvl'].iloc[-1])

        mod_f = model.loc[ (model['wvl'] > ini) & (model['wvl'] < fin) ]
        obj_f = object.loc[ (object['wvl'] > ini) & (object['wvl'] < fin) ]

        return obj_f,  mod_f



    def higher_res(self, main_self, object, model):
        ''''
        Function that returns the name of the higher resolution spectrum and returns the fwhm
        that will be used to convolve it in a different function
        :param object is a panda dataframe of the first spectrum, with column wavelength named wvl
        :param model is a panda dataframe of the second spectrum, with column wavelength named wvl
        ;param main_self is the self in the main program, which is used only to fetch the name of the object
        '''
        #take wavelength spacing of model spectrum, compare it to the object wavelength array.
        fwhmloc_i = (np.abs(model['wvl']- object['wvl'].iloc[0])).idxmin()
        fwhmloc_f = (np.abs(model['wvl']- object['wvl'].iloc[1])).idxmin()
        fwhmloc = [fwhmloc_i, fwhmloc_f]
        if max(fwhmloc) <0:
            fwhmloc_i = (np.abs(model['wvl']- object['wvl'].iloc[-2])).idxmin()
            fwhmloc_f = (np.abs(model['wvl']- object['wvl'].iloc[-1])).idxmin()
            fwhmloc = [fwhmloc_i, fwhmloc_f]

        #The model has higher resolution, will have to be convolved with a fwhm float(fwhmloc_f - fwhmloc_i)
        fwhm_fr = [model['Name'].iloc[0],  float(fwhmloc_f - fwhmloc_i)]

        #if fwhm_fr is 0, it means that the object resolution is higher than that of the object (or the same)
        if fwhm_fr[1] < 1 : #The object will have to be convolved

            #take wavelength spacing of object, compare it to the model wavelength array.
            fwhmloc_i = (np.abs(object['wvl']-model['wvl'].iloc[0])).idxmin()
            fwhmloc_f = (np.abs(object['wvl']-model['wvl'].iloc[1])).idxmin()
            fwhmloc = [fwhmloc_i, fwhmloc_f]
            if max(fwhmloc) <0:
                 fwhmloc_i = (np.abs(object['wvl']-model['wvl'].iloc[-2])).idxmin()
                 fwhmloc_f = (np.abs(object['wvl']-model['wvl'].iloc[-1])).idxmin()
                 fwhmloc = [fwhmloc_i, fwhmloc_f]

            fwhm_fr = [main_self.obj,  float(fwhmloc_f - fwhmloc_i)]

        #If model and spectrum have more or less same resolution,
        #I convolve model to match wvl of object with a fwhm of 1
        if fwhm_fr[1] == 0:
            fwhm_fr = [model['Name'].iloc[0], 1]

        return fwhm_fr


    def smooth_spec(self, obj_ov, this_mod_ov, fwhm_fr):
        ''''
        Function that smooths a spectrum to a lower resolution that matches another input spectrum
        It returns both spectra projected onto the wavelength space of the lower-resolution spectrum

        :param obj_ov is a panda dataframe of the first spectrum, with column wavelength named wvl
        :param this_mod_ov is a panda dataframe of the second spectrum, with column wavelength named wvl
        ;param main_self is the self in the main program, which is used only to fetch the name of the object
        '''
        std = float(fwhm_fr[1]/2.355)
        gauss_kernel = Gaussian1DKernel(stddev= std) #FWHM = 2.355*sigma

        if fwhm_fr[0] == this_mod_ov['Name'].iloc[0]:

            smooth_mod = convolve(this_mod_ov['flux'].values,gauss_kernel,boundary='extend')
            smooth_emod = convolve(this_mod_ov['flux_err'].values,gauss_kernel,boundary='extend')

            #interpolate the smoothed model flux onto the obj grid
            ynew = interpolate.interp1d(this_mod_ov['wvl'].values,smooth_mod,kind='cubic', fill_value="extrapolate")
            new_mod_flux= ynew(obj_ov['wvl'].values)

            #assume very coarsely that S/N ~ sqrt(S) and that locally the S ~ constant with channel
            smooth_emod = smooth_emod/np.sqrt(fwhm_fr[1])
            ynew = interpolate.interp1d(this_mod_ov['wvl'].values,this_mod_ov['flux_err'].values,kind='cubic', fill_value="extrapolate")
            new_mod_fluxerr = ynew(obj_ov['wvl'].values)

            #interchange wavelenths
            new_mod_wvl = obj_ov['wvl'].values

            this_mod_ov = pd.DataFrame( data = {'wvl': new_mod_wvl, 'flux': new_mod_flux, 'flux_err': new_mod_fluxerr})

        else: #object has higher resolution


            smooth_obj = convolve(obj_ov['flux'].values, gauss_kernel,boundary='extend')
            smooth_eobj = convolve(obj_ov['flux_err'].values, gauss_kernel,boundary='extend')

            #interpolate the smoothed obj flux onto the model grid
            ynew = interpolate.interp1d(obj_ov['wvl'].values, smooth_obj,kind='cubic', fill_value="extrapolate")
            new_obj_flux = ynew(this_mod_ov['wvl'].values)

            #assume very coarsely that S/N ~ sqrt(S) and that locally the S ~ constant with channel
            smooth_eobj = smooth_eobj/np.sqrt(fwhm_fr[1])
            ynew = interpolate.interp1d(obj_ov['wvl'].values,smooth_eobj,kind='cubic', fill_value="extrapolate")
            new_obj_fluxerr = ynew(this_mod_ov['wvl'].values)

            #interchange wavelenths
            new_obj_wvl = this_mod_ov['wvl'].values

            obj_ov = pd.DataFrame( data = {'wvl': new_obj_wvl, 'flux': new_obj_flux, 'flux_err': new_obj_fluxerr})

        return obj_ov, this_mod_ov
