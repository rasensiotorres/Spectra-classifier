import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import warnings


class Plots:
    '''
    Class that returns several plots
    '''

    def plot_obj(self, main_self, obj_data):
        '''
        Plot the calibrated spectrum of the object.
        '''
        fig, ax = plt.subplots(figsize = (12, 7))
        ax.set_xlabel(r"Wavelength (microns)", fontsize=22)
        ax.set_ylabel(r"$F\rm_{\lambda} (erg\,s^{-1}\,cm^{-2}\,\mu^{-1})$", fontsize=22)
        ax.tick_params(direction='inout', length=6, labelsize=18, axis='both', which='major')
        ax.tick_params(direction='inout', length=4, labelsize=18, axis='x', which='minor')

        #obj_data is a dataframe containing name, wavelength, flx and err_flx

        if 'j' in main_self.bands:
            obj_j = obj_data.loc[ (obj_data['wvl'] > 1.) & (obj_data['wvl'] < 1.35) ]
            line = ax.plot(obj_j['wvl'] , obj_j['flux'], color = 'darksalmon')
            plt.fill_between(obj_j['wvl'], (obj_j['flux'] - obj_j['flux_err'] ) , (obj_j['flux'] + obj_j['flux_err'] ) , alpha = 1, color = 'khaki')

        if 'h' in main_self.bands:
            obj_h = obj_data.loc[ (obj_data['wvl'] > 1.5) & (obj_data['wvl'] < 1.8) ]
            line = ax.plot(obj_h['wvl'] , obj_h['flux'], color='darksalmon')
            plt.fill_between(obj_h['wvl'], (obj_h['flux'] - obj_h['flux_err'] ) , (obj_h['flux'] + obj_h['flux_err'] ) , alpha = 1, color = 'khaki')

        if 'k' in main_self.bands:
            obj_k = obj_data.loc[ (obj_data['wvl'] > 2.) & (obj_data['wvl'] < 2.5) ]
            line = ax.plot(obj_k['wvl'] , obj_k['flux'], color= 'darksalmon')
            plt.fill_between(obj_k['wvl'], (obj_k['flux'] - obj_k['flux_err'] ) , (obj_k['flux'] + obj_k['flux_err'] ) , alpha = 1, color = 'khaki')

        ax.legend(line, [main_self.name], frameon=False, fontsize=18)

        return fig, ax


    def plot_chisq_st(self, main_self, models_info):
        '''
        Plot chi_sq vs ST for all the fits
        :Only input parameter is the models_info dataset
        '''
        fig, ax = plt.subplots(figsize = (12, 7))
        ax.set_xlabel(r"Spectral Type", fontsize=22)
        ax.set_ylabel(r"$\chi^{2}$ / $\mu$", fontsize=22)
        ax.tick_params(direction='inout', length=6, labelsize=18, axis='both', which='major')
        ax.tick_params(direction='inout', length=4, labelsize=18, axis='x', which='minor')
        ax.set_ylim([0.1,100])
        ax.set_yscale('log')

        st_ticks = [ 'M3.5', 'M4', 'M4.5', 'M5', 'M5.5', 'M6', 'M6.5', 'M7', 'M7.5', 'M8', 'M8.5',\
                'M9', 'M9.5', 'L0', 'L0.5', 'L1', 'L1.5', 'L2', 'L2.5', 'L3', 'L3.5', 'L4', 'L4.5' ] #add later types if needed

        SpT = models_info['SpT']
        chisq = models_info['chi_sq_red']

        if 'montreal' in main_self.library:
            grav = models_info['grav']
            grav_list = ['beta', 'gamma', 'delta', '--']
            grav_leg = ['β', 'γ', 'δ', 'Field']
            color_list = ['tomato', 'forestgreen', 'darkorange', 'silver']

            #Data to plot (Optimize... )
            for g in np.arange(len(grav_list)):
                x_val = [st_ticks.index(SpT[i]) for i in np.arange(len(SpT)) if (SpT[i] in st_ticks and grav[i] == grav_list[g]) ]
                y_val = [chisq[i] for i in np.arange(len(SpT)) if (SpT[i] in st_ticks and grav[i] == grav_list[g]) ]

                plt.scatter(x_val, y_val, c = color_list[g], marker = 'o', label = grav_leg[g] )

        if 'luhman2017' in main_self.library:
            x_val = [st_ticks.index(i) for i in SpT if i in st_ticks]
            y_val = [chisq[i] for i in np.arange(len(SpT)) if SpT[i] in st_ticks ]

            age = 'old'
            if main_self.young ==1: age = 'young'
            plt.scatter(x_val, y_val, c = 'tomato', marker = 'o', s= 34, edgecolors = 'k',  label = main_self.name + '_vs_' + age + '_luhman2017'  )

        plt.xticks(range(len(st_ticks)), st_ticks)
        for label in ax.xaxis.get_ticklabels()[::2]: label.set_visible(False) #Hide every other tick label
        plt.legend(frameon=False, fontsize=18)

        return fig

    def plot_best_mod(self, main_self, models_info, models_data, obj_data):

        '''''
        Plot the scaled best-fit model together with the input object spectrum
        Input parameter is the models_info, models_data and obj_data
        '''''

        fig, ax = plt.subplots(figsize = (14, 6))
        ax.set_xlabel(r"Wavelength (microns)", fontsize=22)
        ax.set_ylabel(r"$F\rm_{\lambda} (erg\,s^{-1}\,cm^{-2}\,\mu^{-1})$", fontsize=22)
        ax.tick_params(direction='inout', length=6, labelsize=18, axis='both', which='major')
        ax.tick_params(direction='inout', length=4, labelsize=18, axis='x', which='minor')

        #Model
        #Find index minimum chi_sq and take data and info o best model
        min_id = models_info['chi_sq_red'].idxmin()
        best_data = models_data.iloc[min_id]
        best_info = models_info.iloc[min_id]

      
        if 'j' in main_self.bands:
            #Object
            obj_j = obj_data.loc[ (obj_data['wvl'] > 0.8) & (obj_data['wvl'] < 1.35) ]
            lo, = ax.plot(obj_j['wvl'] , obj_j['flux'], color = 'gray')
            #plt.fill_between(obj_j['wvl'], (obj_j['flux'] - obj_j['flux_err'] ) , (obj_j['flux'] + obj_j['flux_err'] ) , alpha = 1, color = 'khaki')

            #Model
            best_j_wvl = best_data['wvl'][np.logical_and(best_data['wvl'] > 1., best_data['wvl'] < 1.35)]
            best_j_flx = best_data['flux'][np.logical_and(best_data['wvl'] > 1., best_data['wvl'] < 1.35)]
            lm, = ax.plot(best_j_wvl , best_info['c_fac']*best_j_flx,  color = 'salmon' )


        if 'h' in main_self.bands:
            #Object
            obj_h = obj_data.loc[ (obj_data['wvl'] > 1.5) & (obj_data['wvl'] < 1.8) ]
            #lo, = ax.plot(obj_h['wvl'] , obj_h['flux'], color='gray')
            plt.plot(obj_h['wvl'] , obj_h['flux'], color='gray', label = 'HIP 79098 (AB)b')
            plt.legend(frameon=False, fontsize=18)
            #plt.fill_between(obj_h['wvl'], (obj_h['flux'] - obj_h['flux_err'] ) , (obj_h['flux'] + obj_h['flux_err'] ) , alpha = 1, color = 'khaki')

            #Model
            best_h_wvl = best_data['wvl'][np.logical_and(best_data['wvl'] > 1.5, best_data['wvl'] < 1.8)]
            best_h_flx = best_data['flux'][np.logical_and(best_data['wvl'] > 1.5, best_data['wvl'] < 1.8)]
            lm, = ax.plot(best_h_wvl , best_info['c_fac']*best_h_flx, color = 'salmon')


        if 'k' in main_self.bands:

            #Object
            obj_k = obj_data.loc[ (obj_data['wvl'] > 2.1) & (obj_data['wvl'] < 2.5) ]
            #lo, = ax.plot(obj_k['wvl'] , obj_k['flux'], color= 'gray')
            plt.plot(obj_k['wvl'] , obj_k['flux'], color= 'gray')

            #plt.fill_between(obj_k['wvl'], (obj_k['flux'] - obj_k['flux_err'] ) , (obj_k['flux'] + obj_k['flux_err'] ) , alpha = 1, color = 'khaki')

            #Model
            best_k_wvl = best_data['wvl'][np.logical_and(best_data['wvl'] > 2., best_data['wvl'] < 2.5)]
            best_k_flx = best_data['flux'][np.logical_and(best_data['wvl'] > 2., best_data['wvl'] < 2.5)]
            lm, = ax.plot(best_k_wvl , best_info['c_fac']*best_k_flx, color = 'salmon')

        if  'montreal' in main_self.library: ax.legend( (lo,lm) , ( main_self.name, best_info['Name'] + ' ' + best_info['SpT'] + 'γ'), frameon=False, fontsize=18)

        else: ax.legend( (lo,lm) , ( main_self.name, 'luhman2017' + ' ' + best_info['SpT']), frameon=False, fontsize=18)


        return fig
