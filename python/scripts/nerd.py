'''Functions for the NeRD method'''
import os

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.interpolate import griddata

from skimage.color import rgb2gray
import rioxarray as rioxr
import xarray as xr

import skimage as sk


def radoNorm(I,theta=np.arange(0, 180,1)):
    '''Normalised radon transform applied to image windows
    '''
    
    ''' resample image to increase number of pixels'''
    # original grid
    x = np.arange(0,I.shape[0]) 
    y = np.arange(0,I.shape[1])
    grid_x_original, grid_y_original = np.meshgrid(x, y)
    points = np.array((grid_x_original.reshape(grid_x_original.size),grid_y_original.reshape(grid_y_original.size))).T
    values = I.reshape(I.size,1) 

    # new grid
    step = 0.5
    x = np.arange(0,I.shape[0]-step,step) 
    y = np.arange(0,I.shape[1]-step,step)
    grid_x_resampled, grid_y_resampled = np.meshgrid(x, y)
    M, N = np.meshgrid(x, y)

    # interpolate image to new grid
    Iresampled = griddata(points, values, (grid_x_resampled, grid_y_resampled), method='linear') 
    Iresampled = Iresampled.squeeze() # dimensions from [X,Y,1] to [X,Y]

    ''' Initialise output structure '''
    # max length of projection axis rho (equal to window's diagonal):
    rho_full_length = np.ceil( np.sqrt(2)*M.max() ) + 1 # round to integer value; +1 to ensure indexing includes end-value

    # set-up ouptut matrices (filled with NaN)
    radon_norm = np.nan*np.ones((int(rho_full_length),int(len(theta))))
    n_pixels  = np.nan*np.ones((int(rho_full_length),int(len(theta))))
    
    ''' Perform radon transform for every projection angle'''
    for i in theta: 

        ''' Rotate image'''
        
        # rotate frame
        Mrot =  np.cos(np.radians(i)) * (M) + np.sin(np.radians(i))*(N)
        Nrot = -np.sin(np.radians(i)) * (M) + np.cos(np.radians(i))*(N)
        Mnew = Mrot + np.abs( Mrot.min() )
        Nnew = Nrot + np.abs( Nrot.min() )

        # Sort rotated frame to ascending x-axis (Nnew)
    #     Nsort, idx_sort = sort(Nnew); 
        Nnew_array = Nnew.reshape(Nnew.size)
        Nsort = np.sort(Nnew_array)
        idx_sort = np.argsort(Nnew_array) # indices of sorted values

        # Extract Image values that belong to the rotated, sorted frame
        Iresampled_array = Iresampled.reshape(Iresampled.size)
        Inew = Iresampled_array[idx_sort]

        ''' Line integration ''' 

        # Round coordinates to neares integer and bin x-values to each integer value)
        # NB: Nrounded == index_integer_groups_in_N! That was not true in matlab because of indexing starting at 1 vs 0
        Nrounded = np.round(Nsort)
        Nunique, index_integer_groups_in_N = np.unique(Nrounded, return_inverse=True)

        # create dataframe
        df = pd.DataFrame(data={'Nround':Nrounded,'Inew':Inew,'pixel_count':np.ones(Nrounded.size)})

        # Group (/bin) by N-values, sum Image value for each group (resulting in line integration)
        grouped = df.groupby('Nround')
        line_integration_df = grouped.sum()
        normalised_integral = (line_integration_df['Inew']/line_integration_df['pixel_count']).to_numpy()

        ''' Structuring of Output '''

        # Calculate required padding to fill full length projection axis. Set padding equal to left&right 
        n_bins_to_fill = 0;
        n_bins_to_fill = rho_full_length - len(Nunique)
        if n_bins_to_fill > 0: # padding required
            if n_bins_to_fill % 2 == 0: # equal padding left/right
                pad_left  = int(n_bins_to_fill/2);
                pad_right = pad_left
            else: # unequal padding
                pad_left  = int(np.ceil(n_bins_to_fill/2))
                pad_right = int(np.floor(n_bins_to_fill/2))

            radon_norm_theta = np.pad(normalised_integral,(pad_left, pad_right),'constant', constant_values=np.nan )   
            n_pix  = np.pad(line_integration_df['pixel_count'].to_numpy(),(pad_left, pad_right),'constant', constant_values=np.nan )   
            
        else: # no padding
            radon_norm_theta = normalised_integral
            n_pix = line_integration_df['pixel_count'].to_numpy()

        # gather output for each theta in one matrix
        radon_norm[:,i] = radon_norm_theta    
        n_pixels[:,i] = n_pix

#     return radon_norm, n_pixels
    return radon_norm





def radonIce(window,im_filter=False):#, levels):
    '''Function to apply radoNorm function and extract damage signal from output.
    '''
    
    ##
    ''' Apply Radon Transform '''
    
    I = np.squeeze(window)
        
    # filter image with laplacian filter
    if im_filter:
        ksize = 3
        I = sk.filters.laplace(I, ksize=ksize, mask=None)

    radon_norm = radoNorm(I) # Radon projection for level X , shape (rho,180) with rho depending on window size
    # radoNorm returns NaN when window has uniform values
        
    ''' Extract response value '''
    
    stdRad = np.nanstd(radon_norm,axis=0) # std for each projection, shape [180,]
    medRad = median_filter(stdRad,size=3,mode='wrap') # filter signal with median filter (+- 1pixel around center) to remove noise

    # find peaks and their location (= projection angle, theta)
    peaks_loc, _ = find_peaks(medRad, distance=10) # returns location of peaks.
    peaks = medRad[peaks_loc]

    # sort peaks highest > lowest
    idx_pksort = np.argsort(peaks)[::-1] # idx of lowest to highest peaks, switched back to front with [::-1] to get highest > lowest order
    loc_sort = peaks_loc[idx_pksort]
    peaks_sort = peaks[idx_pksort];
    
    
    ''' Structuring of Output '''
    
    # initialise output for this image-windwo as (1,1,8) array
    out = np.zeros((1,1,8)) * np.nan

    if len(peaks_sort) > 4:
        pair_array = np.vstack((loc_sort[:4] , peaks_sort[:4])).reshape(-1,order='F') # reshape to shape -1 means 1D shape with size inferring from N elements
        out[:,:,:] = pair_array
    elif len(peaks_sort) <= 4: 
        # select all significant peaks; stack in one array [2, N_peaks] 
        # Reshape [2, N_peaks] to  1D array [ pk_angle_1 , pk_value_1, ... pk_angle_N , pk_value_N]
        pair_array = np.vstack((loc_sort[:len(peaks_sort)] , 
                                peaks_sort[:len(peaks_sort)])).reshape(-1,order='F') # reshape to shape -1 means 1D shape with size inferring from N elements

        out[:,:,:len(pair_array)] = pair_array

    return out




# def myloopfunc(rolling_da):
def process_img_windows(rolling_da, nan_value=-999):
    '''Function to loop image-windows, suitable for multiprocessing
    Input windows is a xarray DataArray, rolling window construct.
    Structure of rolling_da is either shape (x, y, x_win, y_win) or (x_win, y_win, sample) where (x) and (y) are stacked '''
#     print(rolling_da.shape)

    if len(rolling_da.shape) > 3: # window_construct in (x,y,x_win, y_win) shape:

        # initialise 
        result = np.zeros( (rolling_da.shape[0] , rolling_da.shape[1], 8 ))
        
        # loop
        for x in range(0,rolling_da.shape[0],1):
            for y in range(0,rolling_da.shape[1],1):

                # select window
                window = rolling_da[ x,y,:,: ]

                # apply function
                out = radonIce(np.array(window)) # shape (1,1,8)

                # save output in preallocated matrix
                result[x,y,:] = out
                    
    elif len(rolling_da.shape) == 3: # window_construct in (x_win, y_win, sample) shape:
        
        # initialise 
        result = np.zeros( (rolling_da.shape[2], 8 )) 
        
        # loop
        for n in range(0,rolling_da.shape[2],1):
            # select window
            # window = rolling_da[ :,:,n ]
            window = rolling_da.isel(sample=n)

            # check if window contains NaN values (then skip)
            if (window == nan_value).any():
                result[n,:] = np.nan*np.ones((1,8))
            else:
                # apply function
                # out = radonIce(np.array(window)) # shape (1,1,8)
                out = radonIce(window.values) # shape (1,1,8)

                # save output in preallocated matrix
                result[n,:] = out[0,:,:] # shape (1,8) instead of (1,1,8)

            
    return result



def mask_nan_containing_windows(xarray_to_mask, windows_1D, min_count=None):
    # windows_1D should have dims (x_win,y_win,sample)
    windows_1D_nan_mask = windows_1D.where(windows_1D != -999)  # replace all values equal to -999 with np.nan
    
    if min_count is None:
        min_count = windows_1D_nan_mask.attrs['window_size']

    w_mask_unstack = windows_1D_nan_mask.sum(dim='x_win',skipna=True, min_count=wsize) \
                        .sum(dim='y_win',skipna=True,min_count=wsize) \
                        .unstack(dim='sample').transpose("y","x")
    
    xarray_masked = xarray_to_mask.where(~w_mask_unstack.isnull())
    return xarray_masked




def read_img_to_grayscale(imPath, imName,dbmin=-15, dbmax=5):
     # open image and convert RGB to Grayscale
        
    img = rioxr.open_rasterio(os.path.join(imPath , imName))
    img_xyz = img.transpose("y","x","band")

    if img_xyz.shape[2] == 3: # RGB band
        if img_xyz.min().values < 0 or img_xyz.max().values>255:
            raise Exception("img min max of [{:.1f},{:.1f}] but expected [0 255]".format(img_xyz.min().values, img_xyz.max().values)) 
            
        img_gray = xr.DataArray(data= rgb2gray(img_rgb), 
                                coords=(img["y"], img["x"] ), 
                                dims=("y","x"), name="gray_image", 
                                attrs=img.attrs, indexes=img.indexes)#, fastpath=False)
        
    elif img_xyz.shape[2] == 1: # single band
        
        if img_xyz.min().values < 0 or img_xyz.max().values>1:
            print('Found values of img outside bounds [0,1]: [{:.1f} {:.1f}]'
                  .format(img_xyz.min().values, img_xyz.max().values) )
            if img_xyz.quantile(0.9).values < 1 and img_xyz.quantile(0.1).values > 0:
                print(' --> 10 and 90th quantile are between [0,1] --> clip values to [0, 1]')
                img_gray = img_xyz.isel(band=0).clip(0,1)
            else:
                print('--> clip and normalise img using min,max: [{},{}]'.format(dbmin, dbmax) )
                img_clipped = img_xyz.isel(band=0).values
                img_clipped[img_clipped < dbmin] = dbmin
                img_clipped[img_clipped > dbmax] = dbmax
                img_norm = (img_clipped - dbmin) / (dbmax - dbmin)

                img_gray = xr.DataArray(data= img_norm, 
                                coords=(img["y"], img["x"] ), 
                                dims=("y","x"), name="gray_image", 
                                attrs=img.attrs, indexes=img.indexes)#, fastpath=False)

        else:    
            img_gray = img_xyz.isel(band=0)
    
    img_gray.attrs["long_name"] = "grayscale"
    
    # test: fillNA
    # nan_mask = img_gray.isnull()
    img_gray = img_gray.fillna(-999)
    # nan_mask_2 = img_gray.where(img_gray == -999)  # replace all values equal to -9999 with np.nan
    
    return img_gray

def cut_img_to_windows(img_gray,wsize):
    
    # replace nan values with -999 (needed because .dropna() is used later but we want to keep nan-containing windows
    img_gray = img_gray.fillna(-999)
    
    # Set up window construct
    r = img_gray.rolling(x=wsize, y=wsize, min_periods=None)
    rolling_da = r.construct(x="x_win", y="y_win", stride=wsize) # (y x ywin xwin)
    # print(rolling_da.shape, rolling_da.dims)

    # issue: cutouts along the edge of the original raster are filled with NaN's: drop NaN-containing windows
    windows = rolling_da.stack(sample=["x", "y"],)
    windows = windows.dropna(dim="sample",how="any")
    
    print(windows.dims,windows.shape)

    windows_1D = windows # (x,y,sample)
    windows_1D.attrs["window_size"] = wsize
    # windows_2D = windows.unstack('sample').transpose("y","x","x_win","y_win")
    # windows_2D.shape
    return windows_1D#, nan_mask
