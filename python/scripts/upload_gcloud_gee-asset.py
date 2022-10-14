#!/usr/bin/env python3.9
'''
upload files from gcloud bucket to GEE asset

'''

# import sys
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

# -- RAMP
# gcloud_dir = 'gs://ee-export_s1_relorbs/RAMP/damage_detection/geotiffs/'
# pattern = '*crevSig*'
# asset_ID = 'users/izeboudmaaike/damage_detection/RAMP_100m_10px_crevSig'
# descript = 'RAMPtoAsset'
# imRes = 100
# Npix = 10 

# -- S1 relorbs
# # gcloud_dir = 'gs://ee-export_s1_relorbs/r100m/damage_detection/geotiffs/'
# # asset_ID = 'users/izeboudmaaike/damage_detection/S1_relorbs_100m_10px_crevSig'
# # descript = 'S1_r100m_toAsset'
# gcloud_dir = 'gs://ee-export_s1_relorbs/r40m/damage_detection/geotiffs/'
# asset_ID = 'users/izeboudmaaike/damage_detection/S1_relorbs_40m_25px_crevSig'
# descript = 'S1_r40m_toAsset'
# pattern = '*25px_crevSig*'
# bandname = 'crevSig'
# imRes = 40
# Npix = 25 

# -- S1 --> to convert to config.ini
gcloud_dir = 'gs://ee-export_s1_relorbs/pineIsland/damage_detection/geotiffs/'
asset_ID = 'users/izeboudmaaike/damage_detection/S1_relorbs_IW_40m-10px_dmg'
descript = 'S1_PiG_40m_toAsset'
# bandname = 'dmg'
include_vars = ['crevSig','alphaC','dmg']
imRes = 40
Npix = 10
CRS = 'EPSG:3031'
start_upload = True # if config['DATA']['start_export'] == 'True' else False

# -- general (interp of .ini)
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