function img_gray_masked = mask_img_ocean_from_values(I,I_bounds)
% Mask satellite image to filter out open ocean
% and convert to grayscale

% Set pixel values that have RGB [0 0 0] to NaN; these are pixels 
% that are not in satellite image (if image doesnt fill complete
% frame) or are open ocean. Apply moving-average filter over the 
% indexes to exclude pixels that validly have value [0 0 0] 
       
%         
% - filter 0-data areas and replace with nan
% -- convert to grayscale, 
% -- apply NAN mask

if size(I,3)>1 % RGB img
    if size(I,3)>3 % select first three bands 
        I = I(:,:,1:3);
    end
    idx_nan = all( I<50 , 3); 
    idx_nan = movmean(idx_nan,[20 20],'omitnan'); 
    idx_nan( idx_nan < 1) = 0; 
    I = rgb2gray(I); 

else % single band img
    
    idx_nan = I==0;
    idx_nan = movmean(idx_nan,[20 20],'omitnan'); 
    
    ct = 0.4; % 8/20 
    idx_nan( idx_nan <= ct ) = 0;  % cutoff of how much of window should be nans  
    idx_nan( idx_nan > ct ) = 1;   
end

% mask idx_nan
img_masked = double(I);
img_masked(logical(idx_nan)) = NaN; 

% -- CONVERT TO GRAYSCALE scale intensity to [0,1]  (nans will become 0)
img_gray_masked = mat2gray(img_masked,I_bounds);  % (extremes will be mapped onto boundary)
img_gray_masked(logical(idx_nan)) = NaN; 

            