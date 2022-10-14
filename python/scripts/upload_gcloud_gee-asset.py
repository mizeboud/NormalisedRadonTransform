#!/usr/bin/env python3.9
'''
upload files from gcloud bucket to GEE asset

'''

import sys
import configparser
import os
import subprocess
import ee
# import rioxarray as rioxr
# import geopandas as gpd

from datetime import datetime

ee.Initialize()

''' ----------------------------------------------------------
         Functions
-------------------------------------------------------------- '''
def is_datetime(datetime_str,date_format='%Y%m%dT%H%M%S'):
    try:
        return datetime.strptime(datetime_str,date_format )
    except ValueError:
        # skip parts that are not datetime
        pass

''' ----------------------------------------------------------
         Config
-------------------------------------------------------------- '''

def main(configFile):

    ''' -------
    Configuration 
    -----------'''
    
    if configFile is None:
        raise NameError('No config file specified. Run script as "python this_script.py /path/to/config_file.ini"')
    else:
        config = configparser.ConfigParser(allow_no_value=True)
        config.read(os.path.join(configFile))

    gcloud_dir = config['ASSET']['gcloud_bucket']
    asset_ID = config['ASSET']['gee_assetID'] 
    descript = asset_ID # 'S1_PiG_40m_toAsset'
    include_vars = [varName.strip() for varName in config['ASSET']['include_vars'].split(',')] # ['crevSig','alphaC','dmg'] 
    imRes = int(config['DATA']['imRes'])
    Npix = int(config['DATA']['nPix'])
    CRS = config['DATA']['CRS']
    start_upload = True if config['DATA']['start_upload'] == 'True' else False

    # -- settings following from config

    # pattern = '*'+str(Npix)+'px_'+bandname+'*' # e.g.: '*10px_dmg*'
    pattern = '*'+str(imRes)+'m*'+str(Npix)+'px_'+'*' # e.g.: '*40m*10px*'
    export_scale = imRes*Npix
    include_vars.remove('crevSig')

    ''' ----------------------------------------------------------
            Find all files with pattern in gcloud
    -------------------------------------------------------------- '''

    content_clouddir = subprocess.run('gsutil ls ' + gcloud_dir + pattern ,shell=True , check=True, capture_output=True, text=True).stdout.splitlines()
    # all imgs
    gcloud_tif_list = [file for file in content_clouddir if file.endswith('.tif')]
    # remove duplicates for multiple vars (image_alphaC; image_dmg; image_crevSig)
    gcloud_tif_crevSig = [file for file in gcloud_tif_list if 'crevSig' in file]
    # gcloud_tif_alphaC = [file for file in gcloud_tif_list if 'alphaC' in file]
    # gcloud_tif_dmg = [file for file in gcloud_tif_list if 'dmg' in file]

    n_files = len(gcloud_tif_crevSig)
    # if len(gcloud_tif_crevSig) != len(gcloud_tif_alphaC) != len(gcloud_tif_dmg):
    #     raise('Not all variables (crevSig, alphaC, dmg) available for every img')

    print('.. Found {} files in gcloud {}'.format(n_files,gcloud_dir) )
    print('.. Uploading to Asset: {}'.format(asset_ID))

    ''' ----------------------------------------------------------
            Export all images to Asset
    -------------------------------------------------------------- '''

    # files_to_upload = gcloud_tif_crevSig # maybe filter this list for files that already exist? -- if img aalready exists in imCol the task yields an error
    counter = 1
    for gcFile in gcloud_tif_crevSig:
        # print('upload {}/{}'.format(counter, n_files ) )

        # -- load crevSig
        cloudImage = ee.Image.loadGeoTIFF(gcFile).rename('crevSig')

        # -- load other variables
        for outvar in include_vars:
            cloudImage_band = ee.Image.loadGeoTIFF(gcFile.replace('crevSig',outvar)).rename(outvar)
            cloudImage = cloudImage.addBands(cloudImage_band)

        # print('bands in img:', cloudImage.bandNames().getInfo())
        
        # -- add timestamp as metadata to img
        img_name = os.path.basename(gcFile)
        img_date = [is_datetime(part).strftime("%Y-%m-%d") for part in img_name.split('_') if is_datetime(part) is not None][0] # date in 'YYYY-mm-dd'
        relorb_id = img_name.replace('relorb_','').split('_'+str(imRes)+'m_output')[0] # remove prefix [relorb_] and suffix [_40m_output_Npix_var.tif]
        
        cloudImage = cloudImage.set({'system:time_start':img_date, 
                                    'relorb_id':relorb_id})
        
        fileName = img_name.replace('_crevSig','').split('.')[0] # filename without file extention
        # print('.. asset name: ', fileName) 
        
        # TO DO: handle already uploaded imgs in imCol
        task = ee.batch.Export.image.toAsset(**{
            'image': cloudImage,
            'description': descript+str(counter),
            'assetId': asset_ID + '/' + fileName , # if img already exists in imCol, the task yields an error
            'scale': export_scale,
            'crs': CRS,
            'region':cloudImage.geometry()
        })
        
        if start_upload:
            print('.. started upload: ', fileName)
            task.start()
        counter=counter+1
        
    print('Done')


if __name__ == '__main__':
    #  Run script as "python path/to/script.py /path/to/config_file.ini"
        
    # retrieve config filename from command line
    config = sys.argv[1] if len(sys.argv) > 1 else None

    # run script
    main(config)    