#!/usr/bin/env python3.9

''' Main script to apply NeRD method to an image '''

import sys
import configparser
import os
import time
import numpy as np
import json

from pathos.multiprocessing import ProcessingPool as Pool

import nerd
import xarray as xr

xr.set_options(keep_attrs=True)

def main(configFile,imageFile):

    ''' -------
    Configuration 
    -----------'''
    
    if configFile is None:
        raise NameError('No config file specified. Run script as "python main.py /path/to/config_file.ini /path/to/image.tif"')
    else:
        config = configparser.ConfigParser()
        config.read(os.path.join(configFile))

    if imageFile is None:
        raise NameError('No image file specified. Run script as "python main.py /path/to/config_file.ini /path/to/image.tif"')
    else:
        print('img ' , imageFile)
        imPath,imName = os.path.split(imageFile) #  split file from path
        if not imPath:
            raise NameError('No path to image file provided. Provide file as ./path/to/image.tif')
    
    path2threshold = config['PATHS']['path2files']
    threshold_fname= config['PATHS']['threshold_fname']
    try:
        tau_type=config['PATHS']['tau_type'] # 'mean' or 'pct095'
        if tau_type != 'mean' or tau_type != 'pct095':
            raise ValueError('Expected threshold type ''mean'' or ''pct095'', got {}'.format(tau_type))
    except:
        tau_type = None

    img_res = int(config['DATA']['imRes'])
    source = config['DATA']['source']
    dbmin = int(config['DATA']['img_bounds_min'])
    dbmax = int(config['DATA']['img_bounds_max'])

    wsize = int(config['NERD']['window_size_px'])
    cores= int(config['NERD']['cores'])
    window_range = wsize*img_res


    # -- settings following from config

    outPath = os.path.join(imPath,'damage_detection')
    if not os.path.isdir(outPath):
        os.mkdir(outPath)

    path2save = os.path.join(outPath, 'geotiffs_python/')
    if not os.path.isdir(path2save):
        os.mkdir(path2save)


    ''' -------
    Load Image 
    -----------'''


    fname_out = imName[:-4] + '_output_' + str(wsize) + 'px' # [imName_output_Npx]

    if os.path.exists(os.path.join(outPath,fname_out + '.nc')):
        print('Output already exists for {} -- stop'.format(fname_out))

    else:

        print('---- \n  Read img ' , imName)


        # -- read img to grayscale
        img = nerd.read_img_to_grayscale(imPath, imName, dbmin, dbmax) # (y,x)

        # -- check img resolution
        dx = np.unique(img['x'].diff(dim='x'))
        if len(dx) == 1:
            if not dx[0] == img_res:
                raise Exception("Configured img resolution ({}) does not equal grid resolution ({})".format(img_res,dx[0]))
        else:
            raise Exception("Inconsistent grid spacing; dx values are {}".format(dx)) 


        # -- cut windows
        print('.. processing img of res {}m on {}px windows'.format(img_res,wsize))
        windows_df = nerd.cut_img_to_windows(img,wsize=wsize) # (x_win, y_win, sample)



        ''' -------
        MultiProcess NeRD method with Multi Processing
        -----------'''

        start_time = time.time()

        # split data in parts (equal to N cores) --> converts to list
        windows_split = np.array_split(windows_df, cores, axis=2) 
        print( 'Split img windows into {} parts of {} windows for multi-processing '.format(cores, [part.shape[2] for part in windows_split] ) )

        # Map function over each split data-part, and stack the parts back together
        with Pool(cores) as pool: # create the multiprocessing pool 
            print('started pool.map on df_split')
            pool_out = pool.map( nerd.process_img_windows, windows_split) # list of arrays with shape (n_samples_split,8)
            df_out = np.concatenate( pool_out ) # array with (samples,8)


        pool.close()
        pool.join()
        pool.clear()    


        print('.. done with pool.map after: {:.2f}hour'.format((time.time() - start_time)/(60*60)) )


        ''' -------
        Assemble results back to data array
        -----------'''


        results = df_out # array with (samples,8)

        # put back to xarray dataArray to convert back to 2D
        try:
            da_result = xr.DataArray(results,
                                    dims=("sample","out"), 
                                    coords=(windows_df["sample"], range(8)), 
                                    name="output", 
                                    attrs=img.attrs, indexes=img.indexes) 
        except ValueError: # Explicitly passing indexes via the `indexes` argument is not supported when `fastpath=False`. Use the `coords` argument instead.
            da_result = xr.DataArray(results,
                                    dims=("sample","out"), 
                                    coords=(windows_df["sample"], range(8)), 
                                    name="output", 
                                    attrs=img.attrs ) 

        da_result.attrs['long_name'] = 'Output_NeRD'
        da_result.attrs['descriptions'] = '[theta_1,signal_1, theta_2,signal_2, theta_3,signal_3, theta_4,signal_4]'


        # back to 2D
        da_result = da_result.unstack('sample').transpose("y","x","out")
        # print('..Reshaped output {} back to (x,y,out): {}'.format(df_out.shape,da_result.shape))
        # da_result.attrs

        da_result.attrs['img_res'] = img_res
        da_result.attrs['window_size(px)'] = wsize
        da_result.attrs['window_range(m)'] = window_range
        da_result.attrs['crs']='EPSG:3031'


        ''' -------
        Save to netCDF
        -----------'''
        da_result.to_netcdf(os.path.join(outPath,fname_out+'.nc'))

        print('.. output data saved to {}{}------ '.format(outPath,fname_out))

        ''' -------
        Processing of output 
        - Calculate dmg from crevasse signal value
        - Radon angle to crevasse aangle
        - calculate delta_alpha and delta_thetap?
        -----------'''

        alpha_c = da_result.isel(out=0) - 90 # alpha_c = theta_NERD - 90, such that alpha_c [-90 to 90]
        crevSig = da_result.isel(out=1)

        # convert crevSig to dmg
        if tau_type is not None:
            dmg, threshold = nerd.crevsig_to_dmg(crevSig, os.path.join(path2threshold,threshold_fname), 
                                                 source, img_res,wsize,tau_type=tau_type)
        else:
            dmg, threshold = nerd.crevsig_to_dmg(crevSig, os.path.join(path2threshold,threshold_fname), 
                                                 source, img_res,wsize)
        
        ''' -------
        Save to geotiffs
        -----------'''

        fname_crevSig = fname_out + '_crevSig'
        fname_alpha_c = fname_out + '_alphaC'

        fname_dmg = fname_out + '_dmg'

        # export a single band to Cloud Optimzed geotiff
        alpha_c.rio.to_raster( os.path.join(path2save, fname_alpha_c + '.tif'),driver="COG")
        crevSig.rio.to_raster( os.path.join(path2save, fname_crevSig + '.tif'),driver="COG")

        if threshold is not None:
            dmg.where(dmg>0).rio.to_raster( os.path.join(path2save, fname_dmg + '.tif'),driver="COG")


        print(' ----- geotiffs saved to {}------\n '.format(path2save))

    print(' Done \n')
    

if __name__ == '__main__':
    #  Run script as "python path/to/main.py /path/to/config_file.ini /path/to/image.tif"
        
    # retrieve config and image file name from command line
    config = sys.argv[1] if len(sys.argv) > 1 else None
    imageFile = sys.argv[2] if len(sys.argv) > 2 else None
    # run script
    main(config,imageFile)