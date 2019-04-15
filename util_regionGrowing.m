function [inxCollection,flag] = util_regionGrowing(I,ini_inx,t_dist,vol,bound_im,fg)

% Ali Abdollahzadeh 
% A.I. Virtanen Institute for Molecular Sciences
% University of Eastern Finland, Finland
% ali.abdollahzadeh at uef.fi
% -------------------------------------------------------------
% The method is described in 
% Abdollahzadeh, A., Belevich, I., Jokitalo, E., Tohka, J. & Sierra, A. 
% Automated 3D Axonal Morphometry of WhiteMatter.bioRxiv(2018). DOI 10.1101/239228.
% --------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software 
% for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies. The author and University of Eastern Finland make no representations
% about the suitability of this software for any purpose. It is 
% provided "as is" without express or implied warranty.
% -------------------------------------------------------------

[r,c,h] = size(I);
J = false(r,c,h); J(ini_inx) = 1;
flag = 0;

if fg == 1
    ngb = [r;-r;1;-1;(r*c);-(r*c)];
else
    ngb1 = [0; r-1; r; r+1; -r+1; -r; -r-1; 1; -1];
    ngb2 = ngb1 - (r*c);
    ngb3 = ngb1 + (r*c);
    ngb = [ngb1;ngb2;ngb3]; ngb = ngb(2:end);
end
nngb = length(ngb);
L_ini_inx = length(ini_inx);
rgn_mean = sum(I(ini_inx))/L_ini_inx;
inx = zeros(L_ini_inx,nngb);
L_rgn = L_ini_inx;
inxCollection = zeros(vol(2),1);
inxCollection(1:L_rgn) = ini_inx;
while L_ini_inx ~= 0
    rgn_mean_old = rgn_mean;
    for j = 1:nngb
        inx(:,j) = ini_inx + ngb(j);
    end
    inx = unique(inx);
    dist = abs(I(inx) - rgn_mean);
    ini_inx = inx((dist<t_dist) & J(inx)==0 & bound_im(inx)==0);
    J(ini_inx) = 1;
    L_ini_inx = length(ini_inx);
    old_intn = L_rgn*rgn_mean_old;
    inxCollection(L_rgn+1:L_rgn+L_ini_inx) =  ini_inx;
    L_rgn = L_rgn + L_ini_inx;
    if L_rgn>vol(2); flag = 2; J = false; break; end
    rgn_mean = (old_intn+sum(I(ini_inx))) / L_rgn;
    inx = zeros(L_ini_inx,nngb);
end
if L_rgn<vol(1); flag = 1; end
inxCollection = inxCollection(inxCollection~=0);
