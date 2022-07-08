%% Example: downscale an img
% Requires function downsample_img.m
% M. Izeboud, 2022


imPath = '../../example/data/';
imName = 'example_S2_median_2020-12-1_2021-3-1.tif';
resample_method = 'nearest';
kernel_size = 10;

[ I  ,R  ] = readgeoraster([imPath imName]);

[I_resample, R_resample] = downsample_img(I,R,kernel_size,resample_method);

%% show difference
figure

subplot(121);
imagesc(I); axis equal, axis off
ax1 = gca;
title('Original')

subplot(122)
imagesc(I_resample); axis equal, axis off
ax2 = gca;
title('Downsampled')
