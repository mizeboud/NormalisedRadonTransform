function [I_resample, R_resample] = downsample_img(image,R,kernel_size,resample_method)
% Resample image 

if nargin < 4
    resample_method = 'nearest';
end

[I_resample,R_resample] = mapresize(image,R,1/kernel_size,resample_method);
end
