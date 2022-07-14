import ee
import relorbs
import configparser
import os
import time # not needed here

ee.Initialize()

''' ---------------------------------------------------------------------------
        Configuration
-------------------------------------------------------------- '''

configPath = '/Users/tud500158/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - TUD500158/PhD/CrevasseDetection/radonTransform/config_files'
# configPath = '/net/labdata/maaike/NERD/python/config_files'
configFile = 'config_GEE_export_S1.ini'
config = configparser.ConfigParser()
config.read(os.path.join(configPath,configFile))

path2files = config['PATHS']['path2files']
my_bucket = config['PATHS']['gcloud_bucket']

mode = config['DATA']['mode']
t_strt = config['DATA']['t_strt']
t_end = config['DATA']['t_end']
vismin = int(config['DATA']['img_bounds_min'])
vismax = int(config['DATA']['img_bounds_max'])
bnds = config['DATA']['bnds'] 
CRS = config['DATA']['CRS']
scale = int(config['DATA']['imRes']) # test scale
clip_coast = True if config['DATA']['clip_coast'] == 'True' else False

try:
    tileNums = config['DATA']['tileNums']
    tileNums = [int(t) for t in tileNums.split()]
except:
    # no tileNums specified
    tileNums = None


''' ---------------------------------------------------------------------------
        Select all images
-------------------------------------------------------------- '''

fCol_relorbs_list = relorbs.get_S1_relorb_ids_list(t_strt, t_end, bnds=bnds, mode=mode, filter_tileNums=tileNums ) # filter_tileNums=None

# if ROI defined
if tileNums is not None:
    print('Selected {} relorbs for period {} to {}, and tileNums {}'.format(len(fCol_relorbs_list),t_strt,t_end,tileNums))
    relorb_list_fname = 'List_relorbs_S1_'+mode+'_'+bnds+'_'+t_strt+'_'+t_end+'_tileNums'+str(tileNums)+'.txt'
else:
    print('Selected {} relorbs for period {} to {}'.format(len(fCol_relorbs_list),t_strt,t_end))
    relorb_list_fname = 'List_relorbs_S1_'+mode+'_'+bnds+'_'+t_strt+'_'+t_end+'.txt'

# -- save list of relorb img names 
with open(os.path.join(path2files,relorb_list_fname), 'w') as f:
    for img_id in fCol_relorbs_list:
        f.write(str(img_id) + '\n')
    print('.. Written relorbs to {}'.format(relorb_list_fname))

    
''' ---------------------------------------------------------------------------
        Export images to Cloud Bucket
        TO DO: add something that checks if file already exists!!
-------------------------------------------------------------- '''


with open(os.path.join(path2files,relorb_list_fname), 'r') as f:
    img_id_load = [line.rstrip('\n') for line in f]

    
im_task_list = []
print('.. Export to gcloud bucket ', my_bucket)
for i in range(0,len(fCol_relorbs_list)):
    # get img
    imName = img_id_load[i] # select single img id from orbits_to_use
    eeImg = ee.Image('COPERNICUS/S1_GRD/' + imName) 
    eeImg_meta = eeImg.getInfo() # reads metadata
    
    img_geom = eeImg.geometry()
    
    # buffer coastline
    if clip_coast:
        print('.. Clipping img to coastline')
        coastline_buffer = relorbs.get_coastline_buffer(size=20e3,err=1e3) 
        eeImg = eeImg.clipToCollection(coastline_buffer)
        
    
    # export filename
    file_name = 'relorb_'+ imName + '_' + str(scale) +'m'
    
    im_task = relorbs.export_img_to_GCS(eeImg.select(bnds),
                                        bucket=my_bucket,
                                        file_name=file_name,
                                        vismin=vismin,vismax=vismax,scale=scale,CRS=CRS,
                                        export_geometry=img_geom,
                                        start_task=True)
    
    # store img tasks in list to check status later
    im_task_list.append(im_task)
    
# Check status of exports
relorbs.status_task_list(im_task_list)    
    
print('Done')



''' ---------------------------------------------------------------------------
        Images are uploading to gcloud bucket.
        This might take a while.
-------------------------------------------------------------- '''