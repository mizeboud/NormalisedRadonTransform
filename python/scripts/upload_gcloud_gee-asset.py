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

ee.Initialize()

''' settings '''

# -- RAMP
# gcloud_dir = 'gs://ee-export_s1_relorbs/RAMP/damage_detection/geotiffs/'
# pattern = '*crevSig*'
# asset_ID = 'users/izeboudmaaike/damage_detection/RAMP_100m_10px_crevSig'
# descript = 'RAMPtoAsset'
# imRes = 100
# Npix = 10 

# -- S1 relorbs
gcloud_dir = 'gs://ee-export_s1_relorbs/r100m/damage_detection/geotiffs/'
pattern = '*crevSig*'
asset_ID = 'users/izeboudmaaike/damage_detection/S1_relorbs_100m_10px_crevSig'
descript = 'S1_r100m_toAsset'
bandname = 'crevSig'
imRes = 100
Npix = 10 

# -- general
export_scale = imRes*Npix
CRS = 'EPSG:3031'


''' Find all files with pattern in gcloud '''

content_clouddir = subprocess.run('gsutil ls ' + gcloud_dir + pattern ,shell=True , check=True, capture_output=True, text=True).stdout.splitlines()
gcloud_tif_list = [file for file in content_clouddir if file.endswith('.tif')]

n_files = len(gcloud_tif_list)
print('.. Found {} files in gcloud {}'.format(n_files,gcloud_dir) )
print('.. Uploading to {}'.format(asset_ID))

''' Export all images to Asset '''

files_to_upload = gcloud_tif_list # maybe filter this list for files that already exist? -- if img aalready exists in imCol the task yields an error
counter = 1
for gcFile in files_to_upload:
    
    cloudImage = ee.Image.loadGeoTIFF(gcFile);
    # cloudImage.rename(bandname)
    
    # print('upload {}/{}'.format(counter, n_files ) )
    
    fileName = os.path.split(gcFile)[1].split('.')[0] # without file extention
    # print('.. asset name: ', fileName) 
    
    task = ee.batch.Export.image.toAsset(**{
        'image': cloudImage.rename(bandname),
        'description': descript+str(counter),
        'assetId': asset_ID + '/' + fileName , # if img already exists in imCol, the task yields an error
        'scale': export_scale,
        'crs': CRS,
        'region':cloudImage.geometry()
    })
    
    task.start()
    counter=counter+1
    
print('Done')