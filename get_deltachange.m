function [dA,idx_good] = get_deltachange(Qriver,slr,s,r,w,bed_h,fr)

if nargin<7, fr=1; end
dA = (Qriver.*fr - (w.*(s-r).*0.5.*slr))./-bed_h;

idx_good = (~isnan(bed_h) & r<0 & abs(dA./w)<1e3); %only appropriate deltas for calculation

end