import pandas as pd
import numpy as np
from astropy.io import fits
import astropy.io.ascii as pyascii
from os import listdir
import sys

class DataProcessing:
    '''
    Class container of all data processing functions
    '''

    def read_obj(obj):
        '''
        get obj variable from fits file
        :return: pandas dataframe with target information
        '''
        dummy = pyascii.read(obj.dir, data_start=0)
        names=dummy.colnames
        object = pd.DataFrame()
        wvl = dummy[names[0]] ; flx = dummy[names[1]] ; flx_err = dummy[names[2]] #normally it would be names[2]
        
        if obj.mjy == 1: 
            print('---- Input flux of target in mJy converted to W/m2/micron ----')
            c = 3e14 #speed of light in microns
            flx= 1e-26*flx/1000.*c/wvl**2. #mJy to W/m2/micron, as models
            flx_err= 1e-26*flx_err/1000.*c/wvl**2. #mJy to W/m2/micron, as models
        else:  print('---- Input flux of target in W/m2/micron ----')


        for i in np.arange(len(dummy[names[0]])): object = object.append([ [obj.name, wvl[i], flx[i], flx_err[i]] ])
        object.columns = ['Name', 'wvl', 'flux', 'flux_err']

        return object


    def get_models_info(self, main_self):
        '''
        Get the full information about the models from a directory (or library)
        :param library: relative directory where the models are
        :return: pandas dataframe with all models information
        '''

        if 'montreal' not in main_self.library and 'luhman2017' not in main_self.library:
             raise ValueError('Please select a valid library. Either "montreal" or "luhman2017" must be included in the directory')


        if 'montreal' in main_self.library:

            #Read names, URLs and spectral types
            df = pd.read_csv(main_self.library, sep=';')
            n_mont = df['Name'] ; y_mont = df['Young?']
            st_mont = df['SpT'] ; url_mont = df['URL']

            #Clean spectral types and get gravity indices
            st_mont = [i.lstrip('< ') for i in st_mont] #assume latest st
            st_clean = [i.split(' ') for i in st_mont] # [st, grav, smthg else]
            st = [i[0].replace(':', '') for i in st_clean]
            grav = [i[1] if len(i) > 1 else 'old' for i in st_clean]
            grav_str = ['beta', 'gamma', 'delta']
            grav = [i if i in grav_str else '--' for i in grav]

            #Save the model information
            models_info = pd.DataFrame( data = {'Name': n_mont, 'Young?': y_mont, 'SpT': st, 'grav': grav, 'URL': url_mont})
            #models_info.to_csv(main_self.library+'montreal_models_info.csv')

            if main_self.young == 1 : models_info = models_info.loc[ models_info['Young?'] == 'Y']

        if 'luhman2017' in main_self.library:

            fitsluhman = listdir(main_self.library)
            namesluhman = [i[0:-5] for i in fitsluhman if i.endswith('.fits')] #Names are the fit files names, as they're standard templates
            st = [i.capitalize() if '.' not in i else i[0: i.find('.')].capitalize() for i in namesluhman ]
            st = [i[0:2] + '.' + i[2:] if len(i) > 2 else i for i in st ] #add a point to the spectral types that are in between integer types, e.g., M5.5
            y = ['Y' if 'older' not in i else '--' for i in namesluhman ]
            old = ['Y' if 'younger' not in i else '--' for i in namesluhman ]

            models_info = pd.DataFrame( data = {'Name': namesluhman, 'SpT': st, 'Young?': y, 'Old?': old  })

            if main_self.young == 1 :
                models_info = models_info.loc[ models_info['Young?'] == 'Y']
            else:
                models_info = models_info.loc[ models_info['Old?'] == 'Y']

        models_info = models_info.reset_index(drop=True) #reset index because we might have removed some rows (old objects, for instance)

        return models_info


    def get_models_data(self, main_self, models_info):
        '''
        Uses the info dataframe to extract the wavelength and flux of each model
        '''
    
        models_data = pd.DataFrame(columns = ['Name', 'wvl', 'flux', 'flux_err'])
        print('---- Loading the data from the models ----')

        if 'montreal' in main_self.library:

            #########   READ THE MONTREAL OBJECT
            for url, name in zip(models_info['URL'], models_info['Name']):

                hdul = fits.open(url, mode='readonly')
                mont = hdul[0].data.byteswap().newbyteorder()
                models_data = models_data.append( {'Name':name, 'wvl': mont[0], 'flux': mont[1], 'flux_err': mont[2]}, ignore_index=True)
                hdul.close()

        if 'luhman2017' in main_self.library:

            for name in models_info['Name']:
                flux = fits.open(main_self.library+name+'.fits', mode = 'readonly')[0].data.byteswap().newbyteorder()
                flux_err = flux*0.01 #They've got no error bars
                wvl = np.linspace(0.7, 2.5, len(flux)) #They've got no wvl values

                models_data = models_data.append( {'Name': name, 'wvl': wvl, 'flux': flux, 'flux_err': flux_err}, ignore_index=True)

        return models_data #each row constitutes a model (wvl in first column, flux in second, flux_err in third)



    def concat_bands(self, main_self, data):

        #Check that only j,h,k bands are selected
        for i in main_self.bands:
            if i not in ['j', 'h', 'k']: raise Exception('Input bands should be in the [j, h, k] region')

        # J,H,K
        dataj = pd.DataFrame() ; datah = pd.DataFrame() ; datak = pd.DataFrame()

        if 'j' in main_self.bands: dataj = data.loc[ (data['wvl'] > .8) & (data['wvl'] < 1.35) ]
        if 'h' in main_self.bands: datah = data.loc[ (data['wvl'] > 1.5) & (data['wvl'] < 1.8) ]
        if 'k' in main_self.bands: datak = data.loc[ (data['wvl'] > 2.) & (data['wvl'] < 2.5) ]
        bands = [dataj, datah, datak]
        datafinal = pd.concat(bands)

        return datafinal



    def quality_check(self, data):
        '''
        Checks whether this spectrum passes our quality thresholds
        data: pandas dataframe with all spectrum information
        returns a quality flag
        '''
        quality = 1 #Suppose it's good

        if np.median(data['flux_err']) / np.median(data['flux']) > 0.03:
            print(data['Name'].iloc[0] + ' SKIPPED, not good precision')
            quality = 0

        wvl = data['wvl'].values
        if (np.any(wvl[1:] <= wvl[:-1])) == True:
            print(data['Name'].iloc[0] + " SKIPPED because wavelength not increasing monotonically")
            quality = 0

        return quality
