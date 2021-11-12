function [dA,idx_good] = func_delta_areachange(Qriver,slr,s,r,w,bed_h,fr)

if nargin<7, fr=1; end
dA = (Qriver.*fr - (w.*(s-r).*0.5.*slr))./-bed_h;

idx_good = (~isnan(bed_h) & bed_h<0 & r<0); % & abs(dA./w)<1e3); %only appropriate deltas for calculation

end