function [ds,dr] = get_deltachange(qs,qs_sus,slr,s,r,beta,alpha,bed_h,fr)

%do we need all of this? %assume profile is a reflection of the past 1000 years?

%how does qs compare to slr rate?
dr = -slr./beta; %m retreat of alluvial bedrock transition

qs_tot = qs_sus.*fr+qs;

%how much do you expect based on profile?
ds = (qs_tot - (s-r-dr).*slr)./slr;

idx = (slr>0 & ds<0);
ds(idx) = max(ds(idx),max(-slr(idx)./max(1e-4,alpha(idx)),-(slr(idx)-bed_h(idx))./beta(idx)));

idx = (slr>0 & ds>0);
ds(idx) = ((qs_tot(idx) - (s(idx)-r(idx)-dr(idx)).*slr(idx))./slr(idx)).*slr(idx)./(-bed_h(idx)+slr(idx));

idx = (slr<=0);
ds(idx) = -slr(idx)./alpha(idx) + qs_tot(idx)./(-bed_h(idx)+slr(idx));

