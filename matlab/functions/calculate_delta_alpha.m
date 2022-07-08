function [delta_alpha, R_resz] = calculate_delta_alpha(alpha_c, Ivx , Ivy, Rvx)
%Calculate the relative oriention of the detected crevasse angle w.r.t to
% the orientation of the local velocity field.
%   ALPHA_C is the crevasse angle w.r.t the image horizontal in range 
%   [-90 90]
%   IVX and IVY are the horizontal x and y velocity components in the same
%   area as ALPHA_C. The grid of IVX and IVY is fitted to the grid of
%   ALPHA_C.
%   The orientation of the local velocity field is calculated using
%   pythagoras theorem, in range [-90 90].
%   DELTA_ALPHA is the difference between the two orientations, in range
%   [-180 180], indacting a parallel (along-flow) alignment for
%   DELTA_ALPHA= 0 and 180; or a perpendicular (acros-flow) alignment for
%   DELTA_ALPHA = -90 or 90.
%
%   M. Izeboud, TU Delft, 2022

if any( size(Ivx) ~= size(alpha_c) )

    scaleA = size(alpha_c,1)/size(Ivx,1); % should be <1, as we're downscaling VX and VY to radon-output-grid.

    [Ivx_resz,R_resz] = mapresize(Ivx,Rvx,scaleA); %returns a raster B that is scale times the size of raster A. RA is a raster reference object that specifies the location and extent of data in A. mapresize returns the raster reference object RB that is associated with the returned raster B. By default, mapresize uses cubic interpolation.
    [Ivy_resz,~] = mapresize(Ivy,Rvx,scaleA); 

    % if resize is not perfect, cut of 1 row/column
    Ivx_resz = fix_matrix_size(Ivx_resz,alpha_c); 
    Ivy_resz = fix_matrix_size(Ivy_resz,alpha_c); 

else 
    Ivx_resz = Ivx;
    Ivy_resz = Ivy;
end

% -- angle of crevsse wrt velocity field
alpha_v = atan(Ivy_resz./Ivx_resz)/(2*pi) * 360; % radians to degree, range [-90 90]; tan(alpha_v) = vy/vx

delta_alpha = alpha_v - alpha_c; %[-180 to 180]
% meeaning: 0 and 180 = along-flow. 
%         -90 and  90 = across-flow


end


