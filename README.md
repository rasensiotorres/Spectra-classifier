# Spectra classifier [In construction] 


<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

[<img src="https://img.shields.io/badge/Python->3.5-yellow.svg?style=flat">](https://www.python.org/)
[<img src="https://img.shields.io/badge/ADS-AsensioTorres (2019)-violet.svg?style=flat">](https://ui.adsabs.harvard.edu/abs/2019A%26A...622A..42A/abstract)
[<img src="https://img.shields.io/badge/Astropy-4.2-green.svg?style=flat">](https://www.astropy.org/)
[<img src="https://img.shields.io/badge/Linkedin-blue.svg?style=flat">](https://www.linkedin.com/in/rub%C3%A9n-asensio-torres-phd-482352169/)


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/rasensiotorres/Spectral-fit/">
    <img src="images/chisq.png" alt="Logo" width="480" height="350">
  </a>

  <h3 align="center">Spectral fit</h3>

  <p align="center">
    A pipeline to fit and classify the infrared spectrum of low-mass stars
    <br />
    <a href="https://github.com/rasensiotorres/Spectral-fit/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/rasensiotorres/Spectral-fit/">View Demo</a>
    ·
    <a href="https://github.com/rasensiotorres/Spectral-fit/issues">Report Bug</a>
    ·
    <a href="https://github.com/othneildrew/rasensiotorres/Spectral-fit/issues">Request Feature</a>
  </p>
</p>


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About The Project
 <img src="images/HIP_79124_B_Luhman.png" alt="Logo" width="340" height="300">
 
This a routine to classify the infrared spectrum of low mass stars. The input spectrum is compared to spectral libraries of young objects via the chi-square goodness of fit statistic, including correlated errors. Youth is mostly indicated in low-resolution near-IR spectra in the triangular H-band continuum shape, which becomes less pronounced as one moves from very low(δ) to low(γ) and intermediate-gravity(β) late M- and L-type dwarfs. Other indicators exist also in the J and K bands, such as FeH absorption (McLean et al. 2003) or the K-band slope.

We can choose between two different spectral libraries::
* Montreal Spectral library (Gagne et al. (2015)): These objects are members of nearby young moving groups (≤120 Myr), with spectral types in the MLT range and δ, γ and β gravities. We consider only high S/N objects, leaving out those with median uncertainties larger than 5% of the median flux value. These spectra were obtained with several instruments, such as Flamingos − 2 and SpeX. Moreover, we include the near-IR Bonnefoy et al. (2014) VLT/SINFONI library of young dwarfs in the M − L transition (M8.5–L4).

* Dereddened near-IR standard spectral templates (Luhman et al. (2017)): Combination of several optical spectra for each subtype in the M spectral region. These
resulting templates are representatives of the Sco-Cen and TWA associations, with an age of about 10 Myrs. 

In the chi-square goodness of fit statistic we incorporate the correlated errors via the covariance matrix C, following Greco & Brandt (2016). This routine has been used in the scientific paper [Asensio-Torres et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...622A..42A/abstract), published in the Journal Astronomy & Astrophysics.


### Built With
This package has been built with Python 3.5, and requires Pandas, Numpy, Scipy, Matplotlib and Astropy 






