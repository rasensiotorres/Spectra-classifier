from data_processing import DataProcessing
from science import Science
import pandas as pd
from plots import Plots
import sys
import numpy as np
import matplotlib.pyplot as plt

class Main:
    '''
    Main Class, where the script and parameters need to be adjusted by user
    '''
    def __init__(self):

        '''
        Main initial execution
        '''
        #name and .txt file of our object
        self.name = 'Test object'
        self.dir = 'spectrum_inmjy.dat'

        #bands we want to fit and plot
        self.bands = ['j','h','k'] #Only 'j', 'h' and 'k' allowed
        #Use squared Coviariance matrix to compute chi square, if it exists. 0 = do not use
        #self.cov_mat = 0
        cov = np.load('cov_matrix_b.npy')
        bad = [0,5,6,7,8, 13,14,15,16, -1] #index of wavelength bands not corrected by tellurics. 
        cov = np.delete(cov, bad, 1)
        cov = np.delete(cov, bad, 0)
        self.cov_mat = cov

        ########Define the LIBRARY of comparison spectra
        # When we select young population here, we refer to Taurus age, which spectra are equal to USco age outside the M3.5--M6.5 region
        self.library = "library_comp/montreal_June2019/montreal_Spectral_Library.csv"
        #self.library = "/Users/asensio-torres/Documents/codes/Python/spectral_fit/library_comp/luhman2017/"
        #Only with a young flag? (If Montreal)
        self.young = 0
        #flux in mJy? Models are in W/m2/mu
        self.mjy = 1


if __name__ == '__main__':
    #If we run it as main program, 
    
    ''''
    Data Processing
    '''
    ########
    # OBJECT
    ########

    #Read the data of our stellar object
    obj = Main()
    obj_data = DataProcessing.read_obj(obj)

    #Take object data in selected J,H,K bands
    dummy = DataProcessing()
    obj_bands = dummy.concat_bands(obj, obj_data)

    #Plot object
    plot = Plots()

    obj.p_obj, ax = plot.plot_obj(obj, obj_data) #Here I store the figure of the object

    ########
    # MODELS
    ########

    models = DataProcessing()
    models_info = models.get_models_info(obj)
    if 'montreal' in obj.library:  #we don't want to spam their server. This should be updated at some point
        models_info = models_info[:250]

    models_data = models.get_models_data(obj, models_info)

    '''
    Data Analysis
    '''
    dummy_sci = Science()

    #add c_fac and chi_sq columns to models_data
    models_info = models_info.assign(c_fac = np.zeros(len(models_data)) + 1e15 ) #scale factor
    models_info = models_info.assign(chi_sq_red = np.zeros(len(models_data)) + 1e15 ) #reduced chi sq


    #Loop through the models and fit each of them to the object, saving its reduced chi_sq
    for n_m in np.arange(len(models_data)):

        #take a model
        this_mod = models_data.iloc[n_m]
        #Convert it to panda dataframe in columns
        this_mod_pd = pd.DataFrame(data = {'Name': this_mod['Name'], 'wvl': this_mod['wvl'], 'flux': this_mod['flux'], 'flux_err': this_mod['flux_err']} )
        obj.this = this_mod_pd

        #take only the selected bands
        this_mod_bands = dummy.concat_bands(obj, this_mod_pd)

        #Take only good data points of the model
        this_mod_bands = this_mod_bands.loc[ (this_mod_bands['flux'] > 0.) & (this_mod_bands['flux_err'] > 0.) ]

        #Check if model passes our quality thresholds, otherwise skip
        quality = models.quality_check(this_mod_bands)
        if quality == 0: continue

        #Find overlapping region between model and object
        obj_ov, this_mod_ov = dummy_sci.overlap(obj_bands, this_mod_bands)

        #Find which spectrum has a higher resolution and smooth if the model and the object have different resolutions
        fwhm_fr = dummy_sci.higher_res(obj, obj_ov, this_mod_ov)
        sobj_ov, sthis_mod_ov = dummy_sci.smooth_spec(obj_ov, this_mod_ov, fwhm_fr)

        #Get chi_sq and C factor
        chi_sq_red, c_fac = dummy_sci.chi_sq(sobj_ov['flux'].values,  sobj_ov['flux_err'].values, sthis_mod_ov['flux'].values, sthis_mod_ov['flux_err'], obj)

        models_info.at[n_m, 'chi_sq_red'] = chi_sq_red
        models_info.at[n_m, 'c_fac'] = c_fac

        obj.models_info = models_info
        obj.models_data = models_data


    '''
    Plots
    '''

    #A plot of chi_sq vs SpT
    obj.p_chi_st = plot.plot_chisq_st(obj, models_info)
    plt.savefig(obj.name + '_chi_st_luhman.eps' , format='eps')

    #A plot of the best model fit to the object
    obj.p_best = plot.plot_best_mod(obj, models_info, models_data, obj_data)
    plt.savefig(obj.name + '_best_fit_luhman.eps' , format='eps')