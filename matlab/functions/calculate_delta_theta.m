function [delta_theta, R_resz] = calculate_delta_theta(alpha_c, Ivx , Ivy, Rv)
%Calculate the relative oriention of the detected crevasse angle w.r.t to
% the orientation of the local velocity field.
%   ALPHA_C is the crevasse angle w.r.t the image horizontal in range 
%   [-90 90]
%   IVX and IVY are the horizontal x and y velocity components in the same
%   area as ALPHA_C. The grid of IVX and IVY is fitted to the grid of
%   ALPHA_C.
%   RV is the MapCellsReference object of IVX and IVY.
%   The orientation of the principal strain (theta_p) is calculated as 
%   explained here: 
%   https://www.continuummechanics.org/principalstressesandstrains.html
%   DELTA_THETA is the difference between the two orientations, in range
%   [-180 180], indacting a parallel alignment for
%   DELTA_THETA = 0 and 180; or a perpendicular alignment for
%   DELTA_THETA = -90 or 90.
%
%   M. Izeboud, TU Delft, 2022

res = abs(Rv.CellExtentInWorldX); % velocity grid resolution
[~, theta_p] = calc_princip_strain_and_orientation(Ivx,Ivy,res,1);

if any( size(theta_p) ~= size(alpha_c) )

    scaleA = size(alpha_c,1)/size(theta_p,1); % should be <1, as we're downscaling VX and VY to radon-output-grid.

    [theta_p_resz,R_resz] = mapresize(theta_p,Rv,scaleA); %returns a raster B that is scale times the size of raster A. RA is a raster reference object that specifies the location and extent of data in A. mapresize returns the raster reference object RB that is associated with the returned raster B. By default, mapresize uses cubic interpolation.

    % if resize is not perfect, cut of 1 row/column
    theta_p_resz = fix_matrix_size(theta_p_resz,alpha_c); 

else 
    theta_p_resz = theta_p;
end

% -- angle of crevsse wrt velocity field
delta_theta = theta_p_resz - alpha_c; %[-180 to 180]
% meeaning: 0 and 180 = parallel
%         -90 and  90 = perpendicular


end


