function get_retention
out = load('D:\Dropbox\WorldDeltas\scripts\GlobalDeltaData.mat','ee','MouthLon','MouthLat','delta_name','BasinID2','QRiver_prist');
load('D:\Dropbox\WorldDeltaLandArea\scripts\GlobalDeltaProfile.mat','s','r','qs_sus','w','alpha','r_h');

ds_obs = out.ee.net_aqua./w*1e6; %observed shoreline change (m/yr)

%inferred topset deposition from slr (m2/yr) and delta growth
qs_inv_topset_ds_f = @(x) (s(x)-r(x)).*ds_obs(x).*alpha(x);

%prof_inv_topset_f = @(x) (-r(x)./(s(x)-r(x)));
%prof_inv_topset_f = @(x) ((s(x)-r(x)).*w(x).*((s(x)-r(x)).*alpha(x))./6000);
%prof_inv_foreset_f = @(x) (s(x).*(1000.*psi(x))./6000);

%inferred from gross delta shape (like syvitski and saito) topset is
%created over 6000 years, so retained per year is volume/yr
prof_inv_topset_f = @(x) ((s(x)-r(x)).*r_h(x)./6000);

%fraction is mean of both estimates, then divided by m3/yr of fluvial input
fr = mean([(w.*qs_inv_topset_ds_f(true(size(w)))) w.*prof_inv_topset_f(true(size(w)))],2)./((365*24*3600*out.QRiver_prist./1600));
idx = ds_obs>0 & out.QRiver_prist>1 & ~isnan(fr) & r~=0;

out.fr = fr; out.qs_sus = qs_sus; out.delta_name = char(out.delta_name); out.idx = idx;

out = rmfield(out,'ee');

save('GlobalDeltaSedimentRetention.mat','-struct','out')


funits = {'dec deg','dec deg','','','ks/s','','m2/yr',''};
fmeta = {'delta river mouth longitude',...
    'delta river mouth latitude',...
    'delta name, where available',...
    'delta BasinID, HydroSheds basin ID',...
    'incoming suspended load sediment, pristine estimate',...
    'Retention efficiency: fraction of incoming suspended load sediment retained within delta topset and foreset',...
    'incoming suspended load sediment normalized by delta width',...
    'suitable for retention analysis',...
    };

create_netcdf('GlobalDeltaSedimentRetention.nc',out,funits,fmeta)
