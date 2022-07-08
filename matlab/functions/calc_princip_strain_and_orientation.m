function [emax, theta_p] = calc_princip_strain_and_orientation(u,v,res,pix)
%calc_princip_strain_and_orientation returns the principle strain based on
%the velocity field components, u and v.
%   U and V are matrices of horizontal velocity components 
%   RES is the resolution of the grid
%   PIX is the size of the gradient, default to 1. So that du/dx is 
%   calculated for u(i) and u(i-PIX). 
% 
%   See also: https://www.continuummechanics.org/principalstressesandstrains.html
%   M. Izeboud, TU Delft, 2022

if nargin<3
pix=1;
end


dudx=(u-circshift(u,[0 pix]))/res;
dvdy=(v-circshift(v,[pix,0]))/res;
dudy=(u-circshift(u,[pix,0]))/res;
dvdx=(v-circshift(v,[0,pix]))/res;


exx = 0.5*(dudx+dudx);
eyy = 0.5*(dvdy+dvdy);
exy = 0.5*(dudy+dvdx);
emax = (exx+eyy)*0.5 + sqrt( (exx-eyy).^2 *0.25 + exy.^2 );
% emin = (exx+eyy)*0.5 - sqrt( (exx-eyy).^2 *0.25 + exy.^2 );

% calculate orientatioon of princip sstrain
theta_p = (atan(2*exy./(exx-eyy))/2)*360/(2*pi); % in degrees

end