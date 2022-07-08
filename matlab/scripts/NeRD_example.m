%%The following script produces a damage map of the example image
% following the NormalisEd Radon transform Damage detection (NeRD) 
% method as proposed by Izeboud et al (2022, in revision).
%
% M. Izeboud - 2022

%% admin
clc;
fprintf(['For application please cite:\n Izeboud et al. 2022.\n', ...
    '''Damage Detection on Antarctic Ice Shelves using the Normalised Radon Transform'' \n', ...
    'Remote Sensing of Environment.\n\n']);

clear all
close all


%% Read imgs and accompanying velocity-imgs (if supplied)

imPath = './example/data/';
imName = 'example_S2_median_2020-12-1_2021-3-1.tif';
imName_vx = 'example_G0240_vx.tif';
imName_vy = 'example_G0240_vy.tif';

source = 'S2';

% -- set path to save results:
outPath = ['./example/damage_detection/'];
path2save = [outPath 'geotiffs/'];


%% Set radon transform details

% -- Radon-Transform: set details
imRes = 30; % export resolution of input image 
Npix = 10;  % choose number of pixels to process NERD on
blockSiz = Npix*imRes;
I_bounds = [0 255]; % min max values of image, needed for converting to grayscale
levels = 1;

% -- save output?
save_mat = true;

% -- define which vars to export as geotiff

save_crevSig = true;
save_alpha_c = true;
save_D = true;
save_dmg = true;
save_Dalpha = true;
save_Dtheta = true;

%% Process img and store output in .mat and geotiffs

% -- define name of output file
outputfile = ['damage_detection_' imName '_' num2str(blockSiz) 'm.mat'];

if ~exist([outPath outputfile])
    fprintf([' ------------- \n'])
    fprintf(['Applying Damage detection on img ' imName  ' \n'])
    disp(outputfile)

    [ I  ,R  ] = readgeoraster([imPath imName ]);
    [ Ivx, Rvx ] = readgeoraster([imPath imName_vx]); % reads accompanying VX image
    [ Ivy, Rvy ] = readgeoraster([imPath imName_vy]); % reads accompanying VY image

    info = geotiffinfo([imPath imName]);
    key = info.GeoTIFFTags.GeoKeyDirectoryTag;

    %% Create NaN-mask; convert to graayscale 

    INAN = mask_img_ocean_from_values(I,I_bounds); 
    I = INAN;
    

    %% Radon-transform
    tic
    fprintf('processing Radon Transform over imagery \n');
    fprintf(['using square blocks of ' num2str(blockSiz) ' meter wide \n\n']); 

    warning off
    % -- process image blockwise. Apply function 'fun' to each block.
    fun = @(block_struct) radonIce_norm(block_struct.data,levels);%,prefilter
    Theta = blockproc(I, [Npix Npix], fun, 'UseParallel', true,...
                'PadPartialBlocks',true,'PadMethod','replicate'); 
    Theta = relaxLabel_multiLevel_output(Theta); % - get the most smooth transition between levels
    warning on
    clear fun
    toc

    % convert angle to be wrt img axis:
    primeDir = Theta(:,:,1);
    alpha_c  = primeDir - 90;

    % -- crev-strength 
    crevSig1 = Score1;
    
     %% remove img border with timestamp mask
    timestamp_file = [ strrep(imName,'__','_') '_timestamp.tif'];
    if exist([imPath timestamp_file]) 
        crevSig1 = mask_img_border_timestamp( crevSig1 , timestamp_file);
        alpha_c  = mask_img_border_timestamp( alpha_c , timestamp_file);
    end
    
    
    %% Crevasse angle relative to local velocity field / principal strain
    [delta_alpha, R_resz] = calculate_delta_alpha(alpha_c, Ivx , Ivy, Rvx);
    [delta_theta, ~       ] = calculate_delta_theta(alpha_c, Ivx , Ivy, Rvy);
    

    %% Calculate damage
    threshold = select_threshold_value(source,blockSiz);
    
    if isnan(threshold)
        save_dmg = false;
        disp(' --> no tresold --> skip dmg calculation')
    else
        % -- classify damage
        dmg = crevSig1-threshold; dmg(dmg<0) = 0; % damage continuoous
        dmg(isnan(dmg)) = 0;
        D = ceil(dmg);   % damage classified binary
        D(isnan(D)) = 0; % set NaN back to 0  
    end
    
   

    %% Save result as .mat structure

    output = struct([]);
    output(1).crevSig1 = crevSig1;
    output(1).alpha_c = alpha_c;
    output(1).delta_alpha = delta_alpha;
    output(1).delta_theta = delta_theta;
    output(1).R_resz = R_resz;
    output(1).blockSize = blockSiz;
    output(1).processedImage = imName;

    fprintf([' ------------- \n'])
    fprintf(['Done with img ' imName  '; \n'])
    if save_mat
        save([outPath outputfile],'output')
        fprintf(['output saved.\n'])
    end
    
    fprintf([' ------------- \n'])

    %% Saving as GeoTiff

    R_resz.RasterSize = size(crevSig1);

    % -- save crevasse signal
    if save_crevSig
        filename = [imName '_' num2str(blockSiz)  'm_crevSig1.tif'];
        geotiffwrite([path2save filename],crevSig1,R_resz,'CoordRefSysCode','EPSG:3031','GeoKeyDirectoryTag',key)
    end

    % -- save crevasse angle
    if save_alpha_c
        filename = [imName '_' num2str(blockSiz)  'm_alpha_c.tif'];
        geotiffwrite([path2save filename],alpha_c,R_resz,'CoordRefSysCode','EPSG:3031','GeoKeyDirectoryTag',key)
    end
    
    % -- save damage map: binary
    if save_D
        filename = [imName '_' num2str(blockSiz)  'm_D.tif'];
        geotiffwrite([path2save filename],D,R_resz,'CoordRefSysCode','EPSG:3031','GeoKeyDirectoryTag',key)
    end

    % -- save damage map: continuous
    if save_dmg
        filename = [imName '_' num2str(blockSiz)  'm_dmg.tif'];
        geotiffwrite([path2save filename],dmg,R_resz,'CoordRefSysCode','EPSG:3031','GeoKeyDirectoryTag',key)
    end
    
    % -- save delta_alpha
    if save_Dalpha
        filename = [imName '_' num2str(blockSiz)  'm_delta_alpha.tif'];
        geotiffwrite([path2save filename],delta_alpha,R_resz,'CoordRefSysCode','EPSG:3031','GeoKeyDirectoryTag',key)
    end
    % -- save delta_theta
    if save_Dtheta
        filename = [imName '_' num2str(blockSiz)  'm_delta_theta.tif'];
        geotiffwrite([path2save filename],delta_theta,R_resz,'CoordRefSysCode','EPSG:3031','GeoKeyDirectoryTag',key)
    end
%         
else % load existing results
    fprintf(['---- \n Already processed ' imName ' \n']) 

end

fprintf('---- \n Finished.\n');



       
       


