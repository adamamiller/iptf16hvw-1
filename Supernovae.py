#import relevant libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import ascii
import json
from IPython.display import display, Image
from specutils import Spectrum1D
from astropy import units
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from astropy.time import Time
from Supernovae import *

#speed of light (km/s)
c = 3e5 


#Define class to hold releveant information for spectra data
class Spectra:
    
    #Initialization function
    def __init__(self, Spectra, epoch, z , MJD_max):
        
        '''
        Spectra (string) - path to JSON formatted spectra file
        epoch (float) - MJD date
        z (float) - redshift of corresponding SN
        MJD_max (float) - date of B band maximum brightness for SN in MJD
        '''
        
        #correct flux for redshift, change wavelength to SN restframe, Normalize flux and store in Spectra
        self.data= Unpack_Spectra(Spectra, z)
        
        #store epoch of obseravation
        self.epoch = float(epoch)
        
        #store phase of observation
        self.phase = float(epoch) - float(MJD_max)



class Lightcurve():
    
    def __init__(self, times, fluxes,  error, band):
        self.band = band
        self.data = pd.DataFrame(list(zip(times, fluxes, error)), columns = ['times', 'flux', 'err'])
        
        
#Create Supernovae class to store Spectral objects
class Supernovae(object):
    
    
    #Initialization function
    def __init__(self, name, redshift, maximum):
        '''
        name (str) - String of SN name
        redshift (float) - redshift of SN
        maximum (float) - date of B band maximum in MJD
        '''
        #Store name of SN
        self.name = name
        
        #Store redshift of SN
        self.redshift = redshift
        
        #Store date of B band maximum brightness
        self.maximum = maximum
        
        #initiate empty list to hold Spectra objects
        self.spectra = []
        
        self.lightcurves = []
    
   
    #define function to return spectra closest to given phase
    def find_spectra(self, phase1):
        '''
        Args:
        phase1 (float )- phase of interest
        
        Returns:
        Spectra object - Spectra object with phase closest to phase1
        '''
        index = np.argmin([ abs(x.phase - phase1) for x in self.spectra])
        return self.spectra[index]
        
    
    
    #define function to store new spectra
    def store_spectra(self, spectra_object):
        ''' 
        Args:
        spectra_object (Spectra) - Spectra object to store
        
        '''
        #Make sure there are no duplicates and that spectra are sorted by date
        if spectra_object  in self.spectra:
            self.spectra.sort(key= lambda x: x.phase)
            print('already exists')
        elif spectra_object.epoch in [x.epoch for x in self.spectra]: 
            self.spectra.sort(key= lambda x: x.phase)
            pass
        
        else:
            self.spectra.append(spectra_object)
            self.spectra.sort(key= lambda x: x.phase)
    
    #define function to store lightcurve
    def store_lightcurve(self, lightcurve_object):
        if lightcurve_object in self.lightcurves:
             print('already exists')
        else:  
            self.lightcurves.append(lightcurve_object)
      
        
    
    
    
#define function that converts wavlengths to restframe and corrects flux for redshift, and normalizes flux
def Unpack_Spectra(Spectra, z, normalization = [5000,6000]):
    '''
    Args:
    Spectra  - one epoch of spectral data in JSON format from OSN
    z (float) - redshift of SN
    normalizationn (list) - 2 item list containing boundaries of region used for normalization
    
    Returns:
    Pandas DataFrame - 2 column dataframe: wavelength and flux
    
    Flux is corrected for redshift and normalized
    
    Wavelength is converted to SN restframe
    
    
    
    '''
    #Extract Wavelengths 
    wavelengths = [float(x[0]) for x in Spectra]
    
    #Extract Fluxes 
    fluxes = [float(x[1]) for x in Spectra]
    
    #correct fluxes for redshift
    fluxes = [correct_flux(flux, z) for flux in fluxes]
    
    #Extract fluxes in normalization range
    rel_flux_range = [x for x in Spectra if (float(x[0])>normalization[0]) & (float(x[0])<normalization[1])]
    
    #Make sure there rel_flux_range isnt empty
    if len(rel_flux_range) == 0:
        #print('No wavelengths in normalization region, not including spectra')
        return None
    
    #Calculate average flux in this range
    flux_sum = 0
    for x in rel_flux_range:
        flux_sum += float(x[1])
    average_flux = flux_sum / float(len(rel_flux_range))
    
    #Normalize flux
    fluxes = [float(flux) / average_flux for flux in fluxes]
    
    #convert wavelength to restframe
    wavelengths = [wavelength / float(1 + z) for wavelength in wavelengths]
    
    #store in pandas dataframe
    df = pd.DataFrame()
    df['Flux'] = fluxes
    df['Wavelength'] =  wavelengths
    return df
    
def correct_flux(flux_obs, z):
    '''
    Args:
    
    flux_obs (int) - observed flux
    z (int) - redshift
    
    Returns:
    int - redshift corrected flux
    '''
    flux_emit = (z * flux_obs) + flux_obs
    return flux_emit

#Define function to get relevant spectra from OSN JSON data file
def create_SN_object(JSON, MJD_max, z):
    '''
    Function to create Supernovae object for given JSON data file from OSN 
    
    Args:
    JSON (str) - path to OSN JSON file of interest
    MJD_max (int) - number of days past maximum brightness
    phase (int) - phase for spectra of interest
    
    Returns:
    Supernovae - Supernovae object with spectra list filled
    '''
    supernovae = Supernovae(str(JSON[0:-5]), z, MJD_max)
    #Load OSN json data
    file = open('../Data/OSN_data/' + str(JSON))
    json_data = json.load(file)
    spectra_data = json_data[JSON[0:-5]]['spectra']
    spectra_data = np.array(spectra_data)
   
    for i in range(len(spectra_data)):
        spectra = Spectra(spectra_data[i]['data'], float(spectra_data[i]['time']) / (1+z), z, MJD_max)
        if spectra.data is None:
            continue
        else:
            supernovae.store_spectra(spectra)
        
    return supernovae

    
#Define function to convert calendar date to MJD
def convert_date_toMJD(date):
    '''
    Args:
    date (str) - string of calendar date (e.g. '2002-8-17')
    
    Returns:
    float - MJD value of given calendar date
    '''
    t = Time(date)
    t.format = 'mjd'
    return t.value
    




#Define function to calculate absorption velocities
def calc_abs_velc(restframe, dopplershifted):
    '''
    Args:
    restframe (float) - restframe wavelength of absorption
    dopplershifted (float) - dopplershifted wavelength of absorption
    
    Returns:
    float - corresponding absorption velocity 
    
    '''
    velocity = ((restframe - dopplershifted) / np.float(restframe))* c
    return velocity

        
    
    
        
