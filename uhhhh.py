# COR2 Data Directories
# https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/P0/a/cor2/
# https://vso3.stanford.edu/cgi-bin/search

# Astropy Documentation
# https://learn.astropy.org/tutorials/FITS-images.html

from astropy.io import fits
from pathlib import Path
import matplotlib.pyplot as plt
import os
from sunpy.time import parse_time, TimeRange
import numpy as np
import pandas as pd
import gc
import requests
import wget
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec #nice alternatie way of doing subplots
#-Scipy Packages-v
from scipy import optimize 
from scipy.signal import find_peaks, find_peaks_cwt #displaying peaks on data
import scipy.special as sc #normalized lower incomplete gamma function is 'sc.gammainc'
import scipy.stats
#-Prieto MTM Packages-v
from multitaper import MTSpec, MTCross, MTSine, SineCross
import multitaper.utils as utils
#-All other packages-v
import numpy as np
import random
from scipy.fft import fft, ifft, irfft, rfft
from scipy.special import gamma
from amtm_signaldetect import *

def DownloadFTSFiles():
#   Secchi data web (other format for other webs)
#   e.g. https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/L0/b/img/cor2/

    date_dir = 20080119
    while date_dir <= 200801123:
        instrument = 'secchi'
        date_dir = str(date_dir)
        camera = 'cor2'
        level = 'L0'
        spacecraft = 'a'
        datatype = 'img'    
        file_path = 'https://stereo-ssc.nascom.nasa.gov/data/ins_data/'+instrument+'/'+level+'/'+spacecraft+'/'+datatype+'/'+camera+'/'+date_dir
        
        
    #   Set local directory
        local_full_path = 'C:/Users/downs/OneDrive/Desktop/intern/STEREO/COR2_files/'
  #      local_full_path = local_path+date_dir
        
    #   Make a directory to download the data in
        if os.path.exists(local_full_path) == False:
            os.makedirs(local_full_path)
        
        print(local_full_path)
           
    #   Get a list of the files available
    #   Request is for htpps webs
        dir_listing = requests.get(file_path)
        html_listing = dir_listing.text
        
    #   Pull each file listing
    #   Splitting with href, because with the print we see the pattern that after href=" comes our file name (this is for this page)
        file_list = []
        seperated_listing = html_listing.split('\n')
        for listing in seperated_listing:
            if 'fts' in listing:
                href_split = listing.split('href="')[1]
                new_file = href_split.split('"')[0]
    #           Add condition to check if file ends with _d4c2b.fts
                if new_file not in file_list and new_file.endswith('_d4c2A.fts'):
                    file_list.append(new_file)
        
        print(file_list)
        
    #   Filter the list to avoid downloading files containing '0800' and '2400'
        soup = BeautifulSoup(html_listing, 'html.parser')
        links = soup.find_all('a')
        file_links = [link.get('href') for link in links if link.get('href') and not ('0800' in link.get('href') or '2400' in link.get('href'))]
        print(file_links) 
        
    #   Download each file
        for file_name in file_list:
    #       Check if file exists already
            if not os.path.isfile(local_full_path+'/'+file_name):
                print("Downloading:",file_name,'to',local_full_path)
                wget.download(file_path+'/'+file_name, local_full_path)
            else:
                print(file_name,'at',local_full_path,'already exists. Skipping')
        date_dir = int(date_dir)
        date_dir += 1
        if date_dir > 20080124:
            break
def PlotAllFtsInDirectory(image_directory):
# Gathers list of all FTS files in a specified directory
    fts_files = [f for f in os.listdir(image_directory) if f.endswith('.fts')]

    for fts_file in fts_files:
# Gets the file path pointed at the FTS file
        file_path = os.path.join(image_directory, fts_file)
# Gets time from the file name
        time = parse_time(str(fts_file)[:-10])
        print(f'Plotting event: {time}')
# Plots image data
        image = fits.getdata(file_path)
        plt.title(f'{time}')
        plt.imshow(image, cmap='gray')
        plt.colorbar()
        plt.show()
        plt.clf()
     
    print('Finished plotting.')
    
    
def PlotRatioImage(image_directory):
# Gathers list of all FTS files in a specified directory

    fts_files = [f for f in os.listdir(image_directory) if f.endswith('.fts')]
# Initializes an empty array to hold all the FTS data arrays - This method only works with arrays that have equivalent shape, size, and length
    array_stack = []

    for fts_file in fts_files:
# for each file, append its data to the array list
        file_path = os.path.join(image_directory, fts_file)
        image_data = fits.getdata(file_path)
        array_stack.append(image_data)
# Gets the first array in the stack
    initial_array = array_stack[0]
    
# Finds the mean values for each index of the arrays in the stack
    array_stack_mean = np.mean(array_stack)
    
# Divides the initial array by the mean array
    ratio_array = initial_array / array_stack_mean
    
# Gets dates and times
    time_start = parse_time(str(fts_files[0])[:-10])
    time_end = parse_time(str(fts_files[-1])[:-10])
    time_range = TimeRange(time_start, time_end)
    
# Plot ratio_array    
    vmin=1.1
    vmax=1.6
    
#    crop_x_start = 0
#    crop_x_end = 1024
#    crop_y_start = 899
#    crop_y_end = 1149
    
    plt.suptitle(f'Ratio Image | vmin: {vmin} | vmax: {vmax} ')
    plt.title(f'{time_start} - {time_end}')
    plt.imshow(ratio_array, vmin=vmin, vmax=vmax)#[crop_y_start:crop_y_end, crop_x_start:crop_x_end]
    plt.gcf().set_dpi(300)
    plt.colorbar()
    plt.show()
    plt.clf()
    print(f'Finished plotting: {time_start} - {time_end}\n')
    print(time_range)


def CheckArrayLengths(image_directory):

    # Initializes an empty dictionary
    array_lengths = {}
    
    for file_name in (f for f in os.listdir(image_directory) if f.endswith('.fts')):
# For each FTS file in a directory, gets the length of the array
            file_path = os.path.join(image_directory, file_name)
            array = fits.getdata(file_path)
            array_lengths[file_name] = len(array)
# Finds the most common length among all of the arrays in the dictionary
    common_length = max(set(array_lengths.values()), key=list(array_lengths.values()).count)

# Iterates through the dictionary and alerts the user if any arrays don't have matching lengths
    for file_name, length in array_lengths.items():
        if length != common_length:
            print(f"\n{file_name} does not match the common length of {common_length}. It has a length of:{len(file_name)}")
        else:
            print('all share common length')
            break
# Change start_i/j to wherever you need to, from beginning = 0,0
def CreatePixelLists(image_directory, start_i=0, start_j=0):
    fts_files = [f for f in os.listdir(image_directory) if f.endswith('.fts')]

    pixels_time_stack = []
    ut_array = []
    
    for fts_file in fts_files:
        file_path = os.path.join(image_directory, fts_file)
        image_data = fits.getdata(file_path)
        event_time = fts_file[:-10]
        # Creates a 2048x2048 3D array, stores the value of the pixels and the time that value occurs for each FTS file
        pixels_time_array = np.zeros((2048, 2048, 1))
        pixels_time_array[:, :, 0] = image_data
#        pixels_time_array[:, :, 1] = event_time
        pixels_time_stack.append(pixels_time_array)
        ut_array.append(event_time)

    for i in range(start_i, 2048):
        for j in range(start_j, 2048):
# For each pixel, take the pixel value from the pixel_time array's z=0 and the times from z=1
            pixel_values = [pixels_time_array[i, j, 0] for pixels_time_array in pixels_time_stack]
            times = [time for time in ut_array]
# Yields instead of returning, otherwise runs out of memory and crashes
            yield (i, j), pixel_values, times


# Change start_i/j to wherever you need to, from beginning = 0,0
def PlotTimeSeries(image_directory, start_i=0, start_j=0):
# Creates a two lists pixel_values and times

# Get time range outside of pixel_data_generator so it doesn't happen every loop
    fts_files = [f for f in os.listdir(image_directory) if f.endswith('.fts')]
    output_directory = Path('C:/Users/jcastig1/Desktop/STEREO/csv/')
    time_range = []
    for fts_file in fts_files:
        event_time = fts_file[:-10]
        time_range.append(event_time)
    time_start = parse_time(str(time_range[0]))
    time_end = parse_time(str(time_range[-1]))
    time_range = f'{time_start} - {time_end}'
    
    pixel_data_generator = CreatePixelLists(image_directory) 
    print('Start plotting')
# For each pixel, create a dataframe where the 0th column is pixel value and 1st is times
    for (i, j), pixel_values, times in pixel_data_generator:
        df = pd.DataFrame({'Pixel_Value': pixel_values, 'Time': times})
# Trim garbage character from string so to_datetime works
        df['Time'] = df['Time'].str[:8] + df["Time"].str[9:]
        
#        csv_file_path = os.path.join(output_directory, f'timeseries_pixel_({i},{j}).csv')
#        df.to_csv(csv_file_path, index=False)
#        print(f'({i}, {j}) csv saved')
        
        df['Time'] = pd.to_datetime(df['Time'])
# Plot timeseries
        fig, ax = plt.subplots()
        ax.plot(df['Time'], df["Pixel_Value"])
        ticks= df['Time'][::5]
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks.dt.strftime('%H'), rotation=0, ha='right')
        plt.xlabel('Time (Hours UT)')
        plt.ylabel('Pixel Value')
        plt.suptitle(f'Pixel Value Time Series for Pixel ({i}, {j})')
        plt.title(f'{time_range}')
        plt.show()
        plt.close('all')
        gc.collect()

    print('Finished plotting.')
    with open('done.txt', 'w') as file:file.write('complete, no crash')
    

def FrequencyMapGenerator(image_directory):
    time_series_generator = CreatePixelLists(image_directory)
    fts_files = [f for f in os.listdir(image_directory) if f.endswith('.fts')]
    time_start = parse_time(str(fts_files[0])[:-10])
    time_end = parse_time(str(fts_files[-1])[:-10]) 
#    frequency_map_array = np.zeros((2048, 2048, 1))
    time_array_hours = np.arange(0, 240, .5)
    for (i, j), pixel_values, times in time_series_generator:
        df = pd.DataFrame({'Pixel_Value': pixel_values, 'Time': times})
        pixel_value_array = df['Pixel_Value'].to_numpy()
        afX = pixel_value_array
        afTime = time_array_hours
        units = 'hrs'
        achFit = 'PL'
        nw_in = 3
        ktpr_in = None
        [afFk, afSk, afAlphak, afFtest] = get_amtm_specs(afX, afTime, nw_in) #NW = 3, Ktpr = 5
        [afFk_bkg, afBk, FITparams] = get_background_psdfit(afTime,afFk, afSk, afAlphak, achFit, nw_in)
        [Fcrit90, Fcrit95, Fcrit99, Ftest_in] = get_ftest_confs(2*nw_in-1, afTime, afFk, afFtest, nw_in, Frange = None)
        [afGamk, Gcrit90, Gcrit95, Gcrit99, Gcrit50] = get_gamtest_confs(afTime, afFk, afSk, afAlphak, afBk, nw_in)    
        [Fpeaks, Gpeaks, FG_pk, FG_pkfreq] = get_gftest_confpeaks(afFk_bkg, afGamk, Ftest_in, Fcrit99, Gcrit99)
    #    print('F99peaks freqs:\n', afFk_bkg[Fpeaks])
    #    print('G99peaks freqs:\n', afFk_bkg[Gpeaks])
    #    print(f'{Fpeaks}\n')
        FGpk_strings = ["%.3f" % x for x in FG_pkfreq]
        achFGfreqpks = ",".join(FGpk_strings)
        achLabFGpks = '+99 FG Peaks:\n[%s] Hz'%(achFGfreqpks)
#    print(f'{achLabFGpks}\n')
        fig = plt.figure(figsize = (25,15))
        gs = GridSpec(2, 2, figure=fig) # (Row x Col) plot
        ax1 = fig.add_subplot(gs[0, 0]) #time series
        ax2 = fig.add_subplot(gs[1, 0]) #aMTM PSD with background fit
        ax3 = fig.add_subplot(gs[0, 1]) #F-test with conf levels
        ax4 = fig.add_subplot(gs[1, 1]) #gamma-test with confidence levels
        plt.suptitle(f'{time_start} - {time_end}')
        fsize = 20 #fontsize of axis labels/ticks
        legsize = 22 #legend fontsize
        tsize = 20 #titlesize
        f_lw = 3
                 #-Plotting Time Series
        ax1.plot(afTime, afX, label = f'Pixel Value Time Series for ({i},{j})')
        ax1.legend(loc = 'upper right', shadow = True, prop={'size': legsize-2})
        #ax1.set_ylabel('Amplitude Y(t)', fontsize = fsize)
        ax1.set_xlabel('time [sec]', fontsize = fsize)
                 #-Plotting PSD with background fit
        ax2.semilogy(afFk, afSk, label = 'aMTM PSD')
        ax2.semilogy(afFk_bkg, afBk, label = 'Background Fit (%s)'%(achFit))#, achPSDfit))
         #ax2.axvline(freq1, color = 'green', ls = '--', lw = f_lw, label = r'$f_0 = %0.3f$Hz'%(freq1))
         #ax2.set_ylabel('PSD S(f)', fontsize = fsize)
        ax2.set_xlabel('Frequency [Hz]', fontsize = fsize)
                 #-Plotting F-test with conf levels
        ax3.plot(afFk_bkg, Ftest_in, label = 'F-test')
        ax3.plot(afFk_bkg[Fpeaks], Ftest_in[Fpeaks], 'o', ms = 10, label = '+99 FPeaks')
        ax3.plot(afFk_bkg[FG_pk], Ftest_in[FG_pk], '*', ms = 12, label = '+99 FG Peaks')
        ax3.axhline(y = Fcrit90,linestyle = 'dotted',lw = f_lw, color = 'orange', label = '90% confidence level')
         #ax3.axhline(y = Fcrit95,linestyle = 'dashed',lw = f_lw, color = 'orange', label = '95% confidence level')
         #ax3.axhline(y = Fcrit99,linestyle = 'dashdot',lw = f_lw, color = 'orange', label = '99% confidence level')
        ax3.set_xlabel('Frequency [Hz]', fontsize = fsize)
                 #-Plotting Gamma-test with conf levels
        ax4.semilogy(afFk_bkg, afGamk, label = r'$\gamma$-Statistic')
        ax4.semilogy(afFk_bkg[Gpeaks], afGamk[Gpeaks], 'o', ms = 10, label = '+99 GPeaks')
        ax4.axhline(y = Gcrit90,linestyle = 'dotted',lw = f_lw, color = 'orange', label = '90% confidence level')
         #ax4.axhline(y = Gcrit95,linestyle = 'dashed',lw = f_lw, color = 'orange', label = '95% confidence level')
         #ax4.axhline(y = Gcrit99,linestyle = 'dashdot',lw = f_lw, color = 'orange', label = '99% confidence level')
         #ax[2].axhline(y = Gcrit50,linestyle = 'dashdot', color = 'green', label = '50% confidence level')
        ax4.semilogy(afFk_bkg[FG_pk], afGamk[FG_pk], '*', ms = 12, label = '%s'%(achLabFGpks))
        ax4.set_xlabel('Frequency [Hz]', fontsize = fsize)
                 #-Legend parameters
        ax1.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
        ax2.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
        ax3.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
        ax4.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
                 #-Setting tick axes labelsize
        ax1.tick_params(axis='both', labelsize=fsize)
        ax2.tick_params(axis='both', labelsize=fsize) 
        ax3.tick_params(axis='both', labelsize=fsize) 
        ax4.tick_params(axis='both', labelsize=fsize) 
        plt.tight_layout()
        plt.show()
        break

# End of Methods
#----------------------------------------------------------------------------------------------------------------------------------
'''FTS File Download'''     
 
'''
# Change date, instrument, location in method at top
DownloadFTSFiles()
'''

#----------------------------------------------------------------------------------------------------------------------------------
'''Ratio Images'''


#image_directory = Path('C:/Users/jcastig1/Desktop/data/fits')
#PlotAllFtsInDirectory(image_directory)
#PlotRatioImage(image_directory)


#----------------------------------------------------------------------------------------------------------------------------------
'''Frequency Maps - Pixel Value Time Series'''

'''
image_directory = Path('C:/Users/jcastig1/Desktop/STEREO/COR2_files/2010paper/')

PlotTimeSeries(image_directory)
'''

#----------------------------------------------------------------------------------------------------------------------------------
'''Frequency Maps - Tests'''

#----------------------------------------------------------------------------------------------------------------------------------
'''Frequency Maps - Spectral Analysis'''

#----------------------------------------------------------------------------------------------------------------------------------

# csv_path = Path('C:/Users/jcastig1/Desktop/STEREO/csv/timeseries_pixel_(0,15).csv/')
# df = pd.read_csv(csv_path)
#DownloadFTSFiles()



image_directory = Path('C:/Users/downs/OneDrive/Desktop/intern/STEREO/COR2_files/')
#CheckArrayLengths(image_directory)

#fts_files = [f for f in os.listdir(image_directory) if f.endswith('.fts')]

#time_start = parse_time(str(fts_files[0])[:-10])
#time_end = parse_time(str(fts_files[-1])[:-10])       
      
FrequencyMapGenerator(image_directory)





'''    
fig = plt.figure(figsize = (25,15))
gs = GridSpec(2, 2, figure=fig) # (Row x Col) plot
ax1 = fig.add_subplot(gs[0, 0]) #time series
ax2 = fig.add_subplot(gs[1, 0]) #aMTM PSD with background fit
ax3 = fig.add_subplot(gs[0, 1]) #F-test with conf levels
ax4 = fig.add_subplot(gs[1, 1]) #gamma-test with confidence levels
fsize = 20 #fontsize of axis labels/ticks
legsize = 22 #legend fontsize
tsize = 20 #titlesize
f_lw = 3
#-Plotting Time Series
ax1.plot(afTime, afX, label = 'Pixel Value Time Series')
ax1.legend(loc = 'upper right', shadow = True, prop={'size': legsize-2})
#ax1.set_ylabel('Amplitude Y(t)', fontsize = fsize)
ax1.set_xlabel('time [sec]', fontsize = fsize)
#-Plotting PSD with background fit
ax2.semilogy(afFk, afSk, label = 'aMTM PSD')
ax2.semilogy(afFk_bkg, afBk, label = 'Background Fit (%s)'%(achFit))#, achPSDfit))
#ax2.axvline(freq1, color = 'green', ls = '--', lw = f_lw, label = r'$f_0 = %0.3f$Hz'%(freq1))
#ax2.set_ylabel('PSD S(f)', fontsize = fsize)
ax2.set_xlabel('Frequency [Hz]', fontsize = fsize)
#-Plotting F-test with conf levels
ax3.plot(afFk_bkg, Ftest_in, label = 'F-test')
ax3.plot(afFk_bkg[Fpeaks], Ftest_in[Fpeaks], 'o', ms = 10, label = '+99 FPeaks')
ax3.plot(afFk_bkg[FG_pk], Ftest_in[FG_pk], '*', ms = 12, label = '+99 FG Peaks')
ax3.axhline(y = Fcrit90,linestyle = 'dotted',lw = f_lw, color = 'orange', label = '90% confidence level')
#ax3.axhline(y = Fcrit95,linestyle = 'dashed',lw = f_lw, color = 'orange', label = '95% confidence level')
#ax3.axhline(y = Fcrit99,linestyle = 'dashdot',lw = f_lw, color = 'orange', label = '99% confidence level')
ax3.set_xlabel('Frequency [Hz]', fontsize = fsize)
#-Plotting Gamma-test with conf levels
ax4.semilogy(afFk_bkg, afGamk, label = r'$\gamma$-Statistic')
ax4.semilogy(afFk_bkg[Gpeaks], afGamk[Gpeaks], 'o', ms = 10, label = '+99 GPeaks')
ax4.axhline(y = Gcrit90,linestyle = 'dotted',lw = f_lw, color = 'orange', label = '90% confidence level')
#ax4.axhline(y = Gcrit95,linestyle = 'dashed',lw = f_lw, color = 'orange', label = '95% confidence level')
#ax4.axhline(y = Gcrit99,linestyle = 'dashdot',lw = f_lw, color = 'orange', label = '99% confidence level')
#ax[2].axhline(y = Gcrit50,linestyle = 'dashdot', color = 'green', label = '50% confidence level')
ax4.semilogy(afFk_bkg[FG_pk], afGamk[FG_pk], '*', ms = 12, label = '%s'%(achLabFGpks))
ax4.set_xlabel('Frequency [Hz]', fontsize = fsize)
#-Legend parameters
ax1.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
ax2.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
ax3.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
ax4.legend(loc = 'upper right', shadow = True,prop={'size': legsize-2})
#-Setting tick axes labelsize
ax1.tick_params(axis='both', labelsize=fsize)
ax2.tick_params(axis='both', labelsize=fsize) 
ax3.tick_params(axis='both', labelsize=fsize) 
ax4.tick_params(axis='both', labelsize=fsize) 
plt.tight_layout()
''' 
