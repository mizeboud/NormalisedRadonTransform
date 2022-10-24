import ee
import relorbs
import configparser
import os
import sys
import time # not needed here
import json

# ee.Authenticate()
ee.Initialize()


def main(configFile):
    ''' ---------------------------------------------------------------------------
            Configuration
    -------------------------------------------------------------- '''

    # configPath = '/Users/tud500158/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - TUD500158/PhD/CrevasseDetection/NormalisedRadonTransform/config_files'
    # # configPath = '/net/labdata/maaike/NERD/python/config_files'
    # configFile = 'config_GEE_export_S1.ini'
    # config = configparser.ConfigParser(allow_no_value=True)
    # config.read(os.path.join(configPath,configFile))

    if configFile is None:
        raise NameError('No config file specified. Run script as "python this_script.py /path/to/config_file.ini"')
    else:
        config = configparser.ConfigParser(allow_no_value=True)
        config.read(os.path.join(configFile))

    path2files = config['PATHS']['path2relorbList']
    my_bucket = config['PATHS']['gcloud_bucket']        # full path:   e.g. ee-export_s1_relorbs/path/to/dir
    bucket_base   = my_bucket.split('/')[0]             # main bucket: e.g. ee-export_s1_relorbs
    bucket_subdir = my_bucket.replace(bucket_base,'').lstrip('/')    # e.g. path/to/dir/
    if bucket_subdir: # if string is not empty, make sure subdir ends with a trailing "/"
        if bucket_subdir[-1] != '/':
            bucket_subdir += '/' # add trailing "/" if not present

    mode = config['DATA']['mode']
    t_strt = config['DATA']['t_strt']
    t_end = config['DATA']['t_end']
    vismin = int(config['DATA']['img_bounds_min'])
    vismax = int(config['DATA']['img_bounds_max'])
    bnds = config['DATA']['bnds'] 
    CRS = config['DATA']['CRS']
    scale = int(config['DATA']['imRes']) # test scale
    clip_coast = True if config['DATA']['clip_coast'] == 'True' else False
    start_export = True if config['DATA']['start_export'] == 'True' else False
    filter_geometry = None if config['DATA']['AOI'] == 'None' else json.loads(config.get("DATA","AOI"))
    tileNums=config['DATA']['tileNums']
    if tileNums is not None:
        tileNums = config['DATA']['tileNums']
        tileNums = [int(t) for t in tileNums.split()]

    try:
        relorb_list_fname = config['PATHS']['fileRelorbList']
    except:
        # no file specified
        relorb_list_fname = None

    # -- Print some settings for information
    AOI=filter_geometry or tileNums
    relorb_list_rw='Read list' if relorb_list_fname is not None else 'Create & Write list'
    print('Loaded settings: \n \
        AOI:       {}\n \
        clipCoast: {}\n \
        relorbs:   {}\n \
        bucket:    {}\
        '.format(AOI,clip_coast,relorb_list_rw,bucket_subdir))

    ''' ---------------------------------------------------------------------------
            Select all images
    -------------------------------------------------------------- '''

    if relorb_list_fname is None:
        fCol_relorbs_list = relorbs.get_S1_relorb_ids_list(t_strt, t_end, bnds=bnds, mode=mode, filter_tileNums=tileNums ,filterGeom=filter_geometry) # filter_tileNums=None

        # if ROI defined
        if tileNums is not None:
            print('Selected {} relorbs for period {} to {}, and tileNums {}'.format(len(fCol_relorbs_list),t_strt,t_end,tileNums))
            relorb_list_fname = 'List_relorbs_S1_'+mode+'_'+bnds+'_'+t_strt+'_'+t_end+'_tileNums'+str(tileNums)+'.txt'
            # relorb_list_fname = 'List_relorbs_S1_'+mode+'_'+bnds+'_'+t_strt+'_'+t_end+'_tileNums-Ross.txt'
        elif filter_geometry is not None:
            print('Selected {} relorbs for period {} to {}, in the specified AOI'.format(len(fCol_relorbs_list),t_strt,t_end))
        else:
            print('Selected {} relorbs for period {} to {}, for all ice shelves'.format(len(fCol_relorbs_list),t_strt,t_end))
            relorb_list_fname = 'List_relorbs_S1_'+mode+'_'+bnds+'_'+t_strt+'_'+t_end+ '_'+ str(scale) + 'm.txt'

        # -- save list of relorb img names 
        with open(os.path.join(path2files,relorb_list_fname), 'w') as f:
            for img_id in fCol_relorbs_list:
                # f.write(str(img_id) + '_' + str(scale) +'m\n')
                f.write(str(img_id) +'\n')
            print('.. Written relorbs to {}'.format(relorb_list_fname))        

        
    ''' ---------------------------------------------------------------------------
            Export images to Cloud Bucket
    -------------------------------------------------------------- '''


    with open(os.path.join(path2files,relorb_list_fname), 'r') as f:
        img_id_load = [line.rstrip('\n') for line in f]

        
    im_task_list = []

    print('.. Export to gcloud bucket ', my_bucket)
    for i in range(0,len(img_id_load)):
        # get img
        imName = img_id_load[i] # select single img id from orbits_to_use
        eeImg = ee.Image('COPERNICUS/S1_GRD/' + imName) 
        # eeImg_meta = eeImg.getInfo() # reads metadata
        
        # -- Erode Edges: buffer img geometry (to be a bit smaller)
        export_geom = eeImg.geometry().buffer(-5e3,1e3)
            
        
        # buffer coastline
        if clip_coast:
            print('.. Clipping img to coastline')
            coastline_buffer = relorbs.get_coastline_buffer(size=20e3,err=1e3) 
            eeImg = eeImg.clipToCollection(coastline_buffer)
            
        
        # export filename
        file_name = 'relorb_'+ imName + '_' + str(scale) +'m'
        
        im_task = relorbs.export_img_to_GCS(eeImg.select(bnds),
                                            bucket=bucket_base,
                                            file_name= bucket_subdir + file_name,
                                            vismin=vismin,vismax=vismax,scale=scale,CRS=CRS,
                                            export_geometry=export_geom,
                                            start_task=start_export)
        
        # store img tasks in list to check status later
        im_task_list.append(im_task)
        
    # Check status of exports
    # if start_export:
    #     relorbs.status_task_list(im_task_list)    
        
    print('Done')


    ''' ---------------------------------------------------------------------------
            Images are uploading to gcloud bucket.
            This might take a while.
    -------------------------------------------------------------- '''

if __name__ == '__main__':
    #  Run script as "python path/to/script.py /path/to/config_file.ini"
        
    # retrieve config filename from command line
    config = sys.argv[1] if len(sys.argv) > 1 else None

    # run script
    main(config)   