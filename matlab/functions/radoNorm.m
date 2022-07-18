function [ Rnorm ] = radoNorm( I, theta )
%RADONORM Compute Normalised Radon Transform of image I for projection angles THETA
%   RNORM = radoNorm(I) calculates the radon transform for THETA=1:180
%   degree
%
%   RNORM = radoNorm(I,theta) calculates the radon transform for angles THETA
%
%   I should be a 2D matrix and may contain NaN values.
%   THETA should be a scalar or an array
%
%   Author: M. Izeboud, TU Delft, 2021
%   See also: Izeboud et. al (2022, in review) "Damage Detection on 
%   Antarctic Ice Shelves using the Normalised Radon Transform", 
%   Remote Sensing of Environment
    
    
if nargin<2, theta = 1:180; 
end

    % irregular interpolation
    [m,n,b] = size(I); 
    [Mf,Nf] = ndgrid(1:m,1:n); % full numbers
    Out = isnan(I); % find NaN's in img
        
    I = double(I); % Make sure datatype is 'double'

    % Interpolate any remaining NaN values by fitting surface through
    % I(x,y) on grid Mf and Nf
    % In case no NaN in I, F==I.
    F = scatteredInterpolant(Mf(~Out),Nf(~Out),I(~Out), ...
        'linear', ...   % interpolation method
        'none');        % extrapolation method

    [M,N] = ndgrid(1:.5:m,1:.5:n); % expand grid resolution: 2x as many pixels.

    % resample image on expanded grid
    Iresampled = F(M,N);

    % calculate centerpoint of grid; used to rotate grid
    center = [size(M,1)./2 size(N,2)./2]; 

    % pre-define output arrays
    rho_max = floor(sqrt(2)*max(M(:)));
    Npix_all = nan(rho_max,length(theta));
    R = nan(rho_max,length(theta));
    Rnorm = nan(rho_max,length(theta));

    % calculate Normalised Radon Transform for all theta
    for i = theta

        % -- rotate frame 
        Mrot = (cosd(i).*(M-center(1)) + sind(i).*(N-center(2)));
        Nrot = (-sind(i).*(M-center(1)) + cosd(i).*(N-center(2)));
        project_x = Nrot + abs( min(Nrot(:)) ); % projection axis rho (virtual x-axis)
        project_y = Mrot + abs( min(Mrot(:)) ); % virtual y-axis; to be integrated

        % -- reorder projection axis and create bins from unique values
        [proj_x_sort,idx_x_sort] = sort(project_x(:)); 
        [proj_x_unique,~,idx_in_bin] = unique(round(proj_x_sort)); % Unique (integer) values of projection axis

        % Get img values in rotated frame using idx_sort of projection axis
        I_rot = Iresampled(idx_x_sort);

        % -- Determine the required padding of projection axis for current theta 
        count_nan = rho_max - length(proj_x_unique); 
        if count_nan > 0 % padding required
            if mod(count_nan, 2) > 0 % unequal padding left/right
                pad_left = ceil(count_nan/2);
                pad_right = floor(count_nan/2);
            else % equal padding
                pad_left = count_nan/2;
                pad_right = pad_left;
            end
        else
            pad_left = 1; 
            pad_rigth = 1;
        end

        % Calculate normalised line integral
        if numel(I_rot)==0
            R(:,i) = 0;
        else

            try
                % Number of elements (pixels) that is integrated for each bin
                Npix = groupcounts(idx_in_bin); 

                % Calculate Normalised Line Integral: sum all img values in 
                % each bin along the projection axis, and 
                % normalise to the number of elements in the bin.
                % The splitapply function performs this step for each bin.
                % Apply padding to projection axis.

                Rnorm( pad_left:pad_left+length(Npix)-1, i) = splitapply(@nansum,I_rot,idx_in_bin)./Npix; % normalised value

                % update use of nansum since matlab2022a
%                     Rnorm( pad_left:pad_left+length(Npix)-1, i) = splitapply(@sum,Inew,'omitnan',i_Nsort)./Npix; % normalised value
            catch
                keyboard;
            end
        end
    end

end

