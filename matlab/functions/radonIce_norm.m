function [ theta ] = radonIce_norm( I, levels, prefilter)
%RADONICE: function that processes the normalised Radon transform for each image wndow cut-out.
%   The function looks for the projection angle where the Normalised Radon
%   transform returns the highest variance.
%   It returns the signal strenght and the accompanying projection
%   orientation.
% 
%   [I] is a NxN window of a grayscale image and should be betweene [0,1]
%   (float). [I] may contain NaN values.
%   [LEVELS] is the desired number of detector outputs. If LEVELS=1 then
%   only the dominant (primary) signal and its orientation is returned. If
%   LEVELS=2 the primary and secondary signals are returned. LEVELS should
%   be <= 4 and defaults to 1. 
%   PREFILTER provides the option to filter the image with a contrast
%   enhancement method. 
%
%   Output THETA is in shape (1,1,8) 
%   Such that for every NxN window 8 values are returned, representing:
%   [angle_level_1,signal_level_1, ... , angle_level_4, signal_level_4]

method = 'median'; 

if nargin < 3
    prefilter = false;
    if nargin < 2
        levels = 1; %2
    end
end
domain = 1:1:180; % angles
subdomain = .1:.1:180;

%% check NaN's in window; discard if too many
prct = 20; % percentage of pixels that is allowed to be NaN before removal

N_NaN = sum(isnan(I),'all');
% N_NaN = sum(isnan(I(:))); % for matlab <2017
N_pix = size(I,1)*size(I,2);
prct_nan = N_NaN/N_pix * 100;

if prct_nan > prct
    theta =NaN.*ones(1,1,8);
else

    %% create bi-nomial filters
    pascalTri = fliplr(pascal(2.*levels)); % triangle of Pascal

    Rad = [];
    for i = 1:2:2.*levels
        if prefilter % filter the image

            % binomial filter (1D)
            binom = diag(pascalTri,i); % takes main diagonal of matrix: level1= [1; 2 ; 1], level2= [1]
            % laplacian filter (1D)
            lap = conv(binom,[1 -2 1]);
            lap = repmat(lap(:),1,length(lap));
            % laplacian filter (2D)
            lap2 = lap'+lap;
            % high pass filtering
            E = imfilter(I, lap2, 'symmetric'); % [-1 -2 -1; -2 12 -2; -1 -2 -1] -- maaike: nee, het is [2 -1 2; -1 -4 -1; 2 -1 2];

        else % Apply radon without pre-processing the image with a filter 

            E = I; 
        end

        % randon transform; normalised
        RadonLevel = radoNorm(E,domain);

        % stack results from multiple levels 
        Rad = cat(1, Rad, RadonLevel); 
    end


    %% statistics
    stdRad = std(Rad,1,'omitnan'); % standard deviation
    medRad = median(cat(1, ... % running median filter
        circshift(stdRad,[0 -1]), stdRad, circshift(stdRad,[0 1])),1);
    % estimate peak through values
    [pks,loc] = findpeaks(medRad, 'minpeakdistance', 10, 'sortstr', 'descend');

    %% gather output

    if ~isempty(pks)
        n_peaks = length(pks);
        if n_peaks>=4 % max 4 signal peaks
            theta = permute(reshape([loc(1:4); pks(1:4)], 1, []), [1 3 2]);
        else 
            theta = NaN.*ones(4,2);
            theta(1:n_peaks,1:2) = [loc(1:n_peaks); pks(1:n_peaks)]'; % all the angles & pks(=signalstrength) that voldoen aan SNR>alpha are passed to output (=theta). 
            theta = permute(reshape(theta', 1, []), [1 3 2]);
        end
    else
        theta =NaN.*ones(1,1,8);
    end    

end

end

