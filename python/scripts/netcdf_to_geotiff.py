import sys
import os
import nerd
import xarray as xr
import json

xr.set_options(keep_attrs=True)


# nc_name='S1_recentComposite_PiG_2019-12-1_2020-3-1_leeFilter_output_10px.nc'
# # fname_out = nc_name[:-3]
# path_to_nc_file = os.path.join(homedir,'PhD/CrevasseDetection/GoogleEarth/imageExport/fixedFrame/compositeExamples/damage_detection/',nc_name)
# path_to_threshold_file = os.path.join(homedir,'PhD/CrevasseDetection/NormalisedRadonTransform/files/dmg_threshold_dictionary.json')


def check_and_save_geotiffs(path_to_nc_file, path_to_threshold_file , path2tiff=None):
    ''' -------
    Load info / settings

    The netCDF file convention is: [imName_output_Npx.nc]
    And geotiff file convention i: [imName_output_Npx_varName.tif]

    It is assumed that the netCDF file is in a folder ./damage_detection/this_netcdf_file.nc
    and that the geotiffs are in a subfolder: ./damage_detection/geotiffs_python/
    -----------'''
    # get netcdf filename
    path_netcdf, netcdf_fname = os.path.split(path_to_nc_file)

    fname_out = netcdf_fname.split('.')[0]
    if path2tiff is None:
        path2tiff = os.path.join(path_netcdf,'geotiffs_python')

    # define variables that should be saved
    fname_crevSig = fname_out + '_crevSig'
    fname_alpha_c = fname_out + '_alphaC'
    fname_dmg = fname_out + '_dmg'

    ''' -------
    Read netcdf
    -----------'''

    if os.path.exists(os.path.join(path2tiff,fname_crevSig + '.tif')) and os.path.exists(os.path.join(path2tiff,fname_alpha_c + '.tif')) and os.path.exists(os.path.join(path2tiff,fname_dmg + '.tif')):
        print('..All geotiffs exist for {}'.format(netcdf_fname) )

    else:
        print('Read {}'.format(netcdf_fname))
        
        da_result = xr.open_dataarray(path_to_nc_file)

        # infer settings 
        source=netcdf_fname[0:2] # 'S1'
        if source in ['S1','S2','L7','L8']:
            print('hoi')
        elif netcdf_fname.startswith('relorb'):
            source='S1'
        else:
            raise('Error: cannot infer data source from filename')
        
        wsize = da_result.attrs['window_size(px)'] # e.g. 10px 
        wrange= da_result.attrs['window_range(m)'] # e.g. 300m
        try: # check if runs with error
            img_res = da_result.attrs['img_res']
        except KeyError: # img_res is not an attribute in netcdf (which is the case for old versions)
            img_res = int(wrange/wsize) # calculate from window range and Npix


    ''' -------
    Save to geotiffs (if they do not exist yet)
    -- export a single band to Cloud Optimzed geotiff
    -----------'''

    # crevsig
    if not os.path.exists(os.path.join(path2tiff,fname_crevSig + '.tif')):
        print('..Saving {}'.format(fname_crevSig))
        crevSig = da_result.isel(out=1)
        crevSig.rio.to_raster( os.path.join(path2tiff, fname_crevSig + '.tif'),driver="COG")

    # alpha_c
    if not os.path.exists(os.path.join(path2tiff,fname_alpha_c + '.tif')):
        print('..Saving {}'.format(fname_alpha_c))
        alpha_c = da_result.isel(out=0) - 90 # alpha_c = theta_NERD - 90, such that alpha_c [-90 to 90]
        alpha_c.rio.to_raster( os.path.join(path2tiff, fname_alpha_c + '.tif'),driver="COG")

    # dmg
    if not os.path.exists(os.path.join(path2tiff,fname_dmg + '.tif')):
        print('..Saving {}'.format(fname_dmg))    
        crevSig = da_result.isel(out=1)
        # convert crevSig to dmg
        dmg, threshold = nerd.crevsig_to_dmg(crevSig, path_to_threshold_file, source, img_res,wsize)
        if dmg is not None:
            dmg.where(dmg>0).rio.to_raster( os.path.join(path2tiff, fname_dmg + '.tif'),driver="COG")


def main(netcdf_dir, path_to_threshold_file, path2tiff=None):
    
    if netcdf_dir is None:
        raise NameError('No directory specified. Run script as "python this_script.py /path/to/netcdf_dir/ /path/to/threshold_dict.json"')

    if path_to_threshold_file is None:
        raise NameError('No threshold file specified. Run script as "python this_script.py /path/to/netcdf_dir/ /path/to/threshold_dict.json"')
            
    netcdf_filelist = [file for file in os.listdir(netcdf_dir) if file.endswith('.nc')]
    
    # check all netcdf files in directory for corresponding geotiffs
    print('Checking {} netcdf files for corresponding geotiffs'.format(len(netcdf_filelist)))
    
    for nc_file in netcdf_filelist:
        path_to_nc_file = os.path.join(netcdf_dir, nc_file)
        check_and_save_geotiffs(path_to_nc_file, path_to_threshold_file )

if __name__ == '__main__':
    #  Run script as "python path/to/netcdf_to_geotiff.py /path/to/netcdf_dir/ /path/to/threshold/dictionary.json
    # e.g: python ./scripts/netcdf_to_geotiff.py ../Data/S1_relorbs/damage_detection/ ./files/threshold_dict.json
        
    # retrieve paths from command line
    netcdf_dir = sys.argv[1] if len(sys.argv) > 1 else None
    path_to_threshold_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    # run script
    main(netcdf_dir,path_to_threshold_file)        