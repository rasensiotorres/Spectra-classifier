import numpy as np
import math
import sys
from sympy import DiracDelta
import matplotlib.pyplot as plt
import astropy.io.ascii as pyascii

#CHARIS J+H+K lambdas
lam_min = 1159.561439621683   #values from header keyword
lam_max = 2369.344052089524
dloglam = 0.03402720386800707
ch_res = 1/dloglam #Resolution of the data cube
nlam = 22
ch_lam = (lam_min * np.exp(np.arange(nlam)*dloglam)) /1000.
lam_c = np.median(ch_lam)

#read CHARIS spectrum
spectrum=pyascii.read('../spectrum_c_final.dat')
names=spectrum.colnames
spec_wvl = spectrum[names[0]]
spec_flx = spectrum[names[1]]
spec_eflx = spectrum[names[2]]
c = 3e14 #speed of light in microns
spec_flx= 1e-26*spec_flx/1000.*c/spec_wvl**2. #mJy to W/m2/micron, as models
spec_eflx = 1e-26*spec_eflx/1000.*c/spec_wvl**2.
spec_eflx = spec_flx*0.05
# C companion
rho = 22.5 #units of lambda/D
Arho = 0.12
Alambda = 0.21
Adelta = 0.65
sigma_rho = 2.15
sigma_lambda = 1.69

psi_c = np.zeros((nlam,nlam)) #correlation matrix
speccovar_c = np.zeros((nlam,nlam)) #covariance matrix

for i in np.arange(nlam):
    for j in np.arange(nlam):
        fac = 1 if i-j == 0 else 0
        psi_c[i,j] = Arho * math.exp(-0.5 * ((rho/sigma_rho) * (ch_lam[i] - ch_lam[j])/lam_c)**2) + \
                     Alambda * math.exp(-0.5 * ((1/sigma_lambda) * (ch_lam[i] - ch_lam[j])/lam_c)**2 ) +\
                     Adelta * fac

        #speccovar_c[i,j] = (psi_c[i,j]*np.sqrt(spec_eflx[i]*spec_eflx[j]))**2.
        #speccovar_c[i,j] = psi_c[i,j]*spec_eflx[i]*spec_eflx[j]
        speccovar_c[i,j] = (psi_c[i,j]*np.sqrt(spec_eflx[i]**2.*spec_eflx[j]**2.))

np.save('cov_matrix_c',speccovar_c)

###same for B

#read CHARIS spectrum
spectrum=pyascii.read('../spectrum_b_v1.dat')
names=spectrum.colnames
spec_wvl = spectrum[names[0]]
spec_flx = spectrum[names[1]]
spec_eflx = spectrum[names[2]]
c = 3e14 #speed of light in microns
spec_flx= 1e-26*spec_flx/1000.*c/spec_wvl**2. #mJy to W/m2/micron, as models
spec_eflx = 1e-26*spec_eflx/1000.*c/spec_wvl**2.

rho = 4.25 #units of lambda/D
Arho = 0.42
Alambda = 0.54
Adelta = 0.02
sigma_rho = 0.34
sigma_lambda = 0.57

psi_b = np.zeros((nlam,nlam)) #correlation matrix
speccovar_b = np.zeros((nlam,nlam)) #covariance matrix
speccovar_b_smallerr = np.zeros((nlam,nlam)) #covariance matrix

for i in np.arange(nlam):
    for j in np.arange(nlam):
        fac = 1 if i-j == 0 else 0
        psi_b[i,j] = Arho * math.exp(-0.5 * ((rho/sigma_rho) * (ch_lam[i] - ch_lam[j])/lam_c)**2) + \
                     Alambda * math.exp(-0.5 * ((1/sigma_lambda) * (ch_lam[i] - ch_lam[j])/lam_c)**2 ) +\
                     Adelta * fac

        #speccovar_b[i,j] = (psi_b[i,j]*np.sqrt(spec_eflx[i]*spec_eflx[j]))**2.
        #speccovar_b[i,j] = psi_b[i,j]*np.sqrt(spec_eflx[i]**2.*spec_eflx[j]**2.)
        #speccovar_b_smallerr[i,j] = psi_b[i,j]*spec_eflx[i]*spec_eflx[j]
        speccovar_b[i,j] = psi_b[i,j]*np.sqrt(spec_eflx[i]**2.*spec_eflx[j]**2.)

np.save('cov_matrix_b',speccovar_b)

sys.exit()

ch_lam_labels = [round(i,2) for i in ch_lam]
ch_lam_l = [ch_lam_labels[0], ch_lam_labels[5], ch_lam_labels[10], ch_lam_labels[15] ,ch_lam_labels[20]]
#np.save('cov_matrix_c',speccovar)



f,axs = plt.subplots(1,2,sharey=True,figsize=(16,8))
pcm1 = axs[0].pcolormesh(psi_b, cmap='RdYlBu_r',vmin=0,vmax=1.)
pcm2 = axs[1].pcolormesh(psi_c, cmap='RdYlBu_r')
axs[0].locator_params(nbins=5)
axs[0].set_xticklabels(ch_lam_l)
axs[0].set_yticklabels(ch_lam_l)
axs[0].xaxis.set_tick_params(labelsize=20)
axs[0].yaxis.set_tick_params(labelsize=20)
axs[0].set_title(r'$\Psi\rm_{ij}$ (HIP 79124 B)',fontsize=22)
axs[0].set_xlabel(r'Wavelength ($\mu$m)',fontsize=20)
axs[0].set_ylabel(r'Wavelength ($\mu$m)',fontsize=20)

axs[1].locator_params(nbins=5)
axs[1].set_xticklabels(ch_lam_l)
axs[1].set_yticklabels(ch_lam_l)
axs[1].xaxis.set_tick_params(labelsize=20)
axs[1].yaxis.set_tick_params(labelsize=20)
axs[1].set_title(r'$\Psi\rm_{ij}$ (HIP 79124 C)', fontsize=22)
axs[1].set_xlabel(r'Wavelength ($\mu$m)', fontsize=20)
cbaxes = f.add_axes([0.125, 0.05, 0.78, 0.02])
cb = plt.colorbar(pcm1, cax = cbaxes, orientation='horizontal')
#f.colorbar(pcm1,ax=axs,orientation='horizontal',fraction=.05,pad=0.1)
plt.xticks(fontsize=20)

plt.show()
sys.exit()
