function [ Rnorm ] = radoNorm( I, theta )
%RADONAN Summary of this function goes here
%   Detailed explanation goes here
    
    
if nargin<2, theta = 1:180; 
end

% if any(isnan(I(:)))
%     tic;    
    % irregular interpolation
    [m,n,b] = size(I); 
    [Mf,Nf] = ndgrid(1:m,1:n); % full numbers
    
    [M,N] = ndgrid(1:.5:m,1:.5:n); % expand grid resolution: 2x as many pixels.
    Out = isnan(I); % exclude NaN's from the analysis
    if sum(Out(:)) > numel(I)-4 % to many Out: return 0 as output
%         Rnorm = zeros(length(1:.5:m),length(theta));
%         
%         rho_max = floor(sqrt(2)*max(M(:)));
%         Npix_all = nan(rho_max,length(theta));
% %         R = nan(rho_max,length(theta));
%         Rnorm = nan(rho_max,length(theta));
    else
        % Maaike:
        % Fit surface F to values of I(x,y); Mf and Nf are grid and I is z-value but
        % in case of NaN in I, those are interpolated.
        % In case no NaN in I, F==I.
        
%         fprintf(['check datatype for scatteredInterpolant: ' , class(Mf), class(Nf), class(I) ])
        % convert to doubles (required for scatteredInterpolant)
%         Mf = double(Mf);
%         Nf = double(Nf);
        I = double(I);

        F = scatteredInterpolant(Mf(~Out),Nf(~Out),I(~Out), ...
            'linear', ...   % interpolation method
            'none');        % extrapolation method

        [M,N] = ndgrid(1:.5:m,1:.5:n); % expand grid resolution: 2x as many pixels.
        
        %%%%%% TO CHECK: the center value calculaged llike this is the
        %%%%%% corner(/max) value instead of center... ??? --> okay this
        %%%%%% doesnt actually matter
        center = [size(M,1)./2 size(N,2)./2]; % used to determine radius of circle

        
        % construct circle
%         mC = center(1).* cosd(theta) + center(1); % x = r*cos
%         nC = center(2).* sind(theta) + center(2); % y = r*sin
%         In = poly2mask(mC,nC,size(M,1),size(N,2));
%         theta = theta(1:180);
% 
%         Mc = M(In); Nc = N(In); % sample grid
%         center = [M(end)./2 N(end)./2];    
%         R = zeros(max((Nc.*2)-1),length(theta)); % matrix with each grid point as a row (therefore Nc*2; is eigenlijk Nc*Mc (?))
%         Npix_all = zeros(max((Nc.*2)-1),length(theta));

        % norm2: no circle
        rho_max = floor(sqrt(2)*max(M(:)));
        Npix_all = nan(rho_max,length(theta));
        R = nan(rho_max,length(theta));
        Rnorm = nan(rho_max,length(theta));
%         rho_max_test = 0;
        Iresampled = F(M,N); % original image on resampled grid.
        
        
        for i = theta
            % -- new rotated coordinate frame
%             Mnew = (sind(i).*(Mc-center(1)) - cosd(i).*(Nc-center(2))) +center(1);
%             Nnew = (-cosd(i).*(Mc-center(1)) + sind(i).*(Nc-center(2))) +center(2);

            % -- Maaike : correct (?) rotated frame
%             Mnew = (cosd(i).*(Mc-center(1)) + sind(i).*(Nc-center(2))) +center(1);
%             Nnew = (-sind(i).*(Mc-center(1)) + cosd(i).*(Nc-center(2))) +center(2);
            
            % -- rotated frame for square
            Mrot = (cosd(i).*(M-center(1)) + sind(i).*(N-center(2)));
            Nrot = (-sind(i).*(M-center(1)) + cosd(i).*(N-center(2)));
            Mnew3 = Mrot + abs( min(Mrot(:)) );
            Nnew3 = Nrot + abs( min(Nrot(:)) );
            % now try to extract the rigth samples for reprojected axis.
            
            
            
            % -- Mnew is new x-axis. So First re-order Mnew, then we can
            % integrate all Nnew values (new y-axis) for each step of Mnew
            % Also, bin Mnew to integer-values (such that for every theta,
            % the rotated coords are binned to the same values)
            % -- UPDATE: Nnew is new x-axis ha ha I switched them...
            [Nsort,idx_sort] = sort(Nnew3(:)); 
            Nrounded = round(Nsort); % round to nearest integer
%             Mrounded = round(Msort .* 4)/4; % roudn to nearest 0.25 
            [Nunique,~,i_Nsort] = unique(Nrounded); % round Msort to 4 decimals   
            
            
            % -- Determine the padding the projected bins (as every theta 
            % has different length of projection axis)
            count_nan = rho_max - length(Nunique); 
            if count_nan > 0 % padding required
                if mod(count_nan, 2) > 0 % unequal padding
                    pad_left = ceil(count_nan/2);
                    pad_right = floor(count_nan/2);
                else
                    pad_left = count_nan/2;
                    pad_right = pad_left;
                end
            else
                pad_left = 1;
                pad_rigth = 1;
            end
            
            
            % -- For circle method: Extract Inew values from fitted frame
%             Inew = F(Mnew,Nnew);     % image valuse interpolated (on new higher-res grid)
            
            % -- For square method: Extract Inew values from rotated frame
            % using idx_sort & idx_unique values
            Inew = Iresampled(idx_sort);
            
            if numel(Inew)==0
                R(:,i) = 0;
%                 keyboard;
            else

%                 try
%                     Nnew = interp2(Nf,Mf,double(Out),Nnew,Mnew, ...
%                         'nearest');             % interpolation NaN's
%                 catch
%                     keyboard;
%                 end
%                 Inew(isnan(Nnew)|Nnew==1) = NaN;
                
                % radon insensative to NaN's
                try
                    % splitapply performs funcion @func on all groups G
                    % defined in splitapply(@func,X,G). In Nc all values
                    % of 1:0,5:n (within circle) are included, with many
                    % duplicates. So this function averages the values back
                    % to 1 value per 1:0.5:n point.
                    
%                     R(:,i) = splitapply(@nanmean,Inew(:),(Nc(:).*2)-1); 
                    
                    % Maaike:
                    % Take @nansum over the bin, to actually get the
                    % integral (ipv @nanmean)
%                     to_group = (Nc(:).*2)-1;  % set all (half-point) coordinates to full integers; these must be grouped back to 1 value per coord
% %                     [~,groupID] = findgroups(to_group );
%                     Npix = groupcounts(  to_group ); % count element in each group (each x-point/coordinate)
%                     R(:,i) = splitapply(@nansum,Inew(:),to_group)./Npix; % take integral for each coordinate and divide by #pixels
%                     Npix_all(:,i) = Npix;

                    % Square radon_norm
                    Npix = groupcounts(i_Nsort); % number of elements that is integrated for each binned M-value
                    Npix_all(pad_left:pad_left+length(Npix)-1,i) = Npix; 
                    R( pad_left:pad_left+length(Npix)-1, i) = splitapply(@nansum,Inew,i_Nsort); 
                    Rnorm( pad_left:pad_left+length(Npix)-1, i) = splitapply(@nansum,Inew,i_Nsort)./Npix; % normalised value

                    % update use of nansum since matlab2022a
%                     R( pad_left:pad_left+length(Npix)-1, i) = splitapply(@sum,Inew,'omitnan',i_Nsort); 
%                     Rnorm( pad_left:pad_left+length(Npix)-1, i) = splitapply(@sum,Inew,'omitnan',i_Nsort)./Npix; % normalised value
                catch
                    keyboard;
                end
            end
        end
    end
%     toc;
% else
    % build in matlab funtion is probably faster
%     tic;
%     R = radon(I);
    
    % normalise radon output to number of pixels in integral
    % max width: [NOT WORKING CORRECTLY YET; so only use above optioon with normalisation]
%     xprime_maxwidth = size(Rmat,1);

    
%     toc;
% end

end

