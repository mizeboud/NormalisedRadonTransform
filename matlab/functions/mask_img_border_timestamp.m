function [image2mask] = mask_img_border_timestamp( image2mask, timestamp_file)

[ datum, R ] = readgeoraster([imPath 'timestamps/' timestamp_file]); % original resolution

ratio = size(datum,1)/size(image2mask,1);

[datum_resz, ~] = mapresize(datum,R,1/ratio,'nearest');
datum_resz = fix_matrix_size(datum_resz,image2mask);

% -- Method A: mask where date=0
data_nan = datum_resz == 0;

% -- expand border of data-mask by applying mask
idx_nan = movmean(data_nan,[20 20],'omitnan'); 
ct = 0.4; % 8/20
idx_nan( idx_nan <= ct ) = 0;  % cutoff of how much of window should be nans  
idx_nan( idx_nan > ct ) = 1;   

timestamp_mask = logical(idx_nan); % 1 for nodata

% -- apply mask to data:
image2mask(timestamp_mask) = NaN;
    
end