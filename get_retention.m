function get_retention
out = load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','ee','MouthLon','MouthLat','delta_name','BasinID2','QRiver_prist');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaProfile.mat','s','r','qs_sus','w','alpha','r_h','bed_h');

ds_obs = out.ee.net_aqua./w*1e6; %observed shoreline change (m/yr)

%inferred topset deposition from delta growth and topset length
qs_inv_topset_ds_f = @(x) (s(x)-r(x)).*ds_obs(x).*alpha(x);
qs_inv_foreset_ds_f = @(x) (max(1,-bed_h(x)).*ds_obs(x));



%inferred from gross delta shape (like syvitski and saito) topset is
%created over 6000 years, so retained per year is volume/yr
prof_inv_topset_f = @(x) ((s(x)-r(x)).*r_h(x)*0.5./6000);
prof_inv_foreset_f = @(x) (s(x).*max(1,-bed_h(x)).*0.5./6000);

%fraction is mean of both estimates, then divided by m3/yr of fluvial input
out.fr = mean([(w.*qs_inv_topset_ds_f(true(size(w)))) w.*prof_inv_topset_f(true(size(w)))],2)./((365*24*3600*out.QRiver_prist./1600));
out.fr_f = mean([(w.*qs_inv_foreset_ds_f(true(size(w)))) w.*prof_inv_foreset_f(true(size(w)))],2)./((365*24*3600*out.QRiver_prist./1600));

out.qs_sus = qs_sus; out.delta_name = char(out.delta_name);

out.idx = (ds_obs>0 & out.QRiver_prist>1 & ~isnan(out.fr) & r~=0);

out = rmfield(out,'ee');

save('GlobalDeltaSedimentRetention.mat','-struct','out')

out.idx = double(out.idx);
funits = {'dec deg','dec deg','','','ks/s','','','m2/yr',''};
fmeta = {'delta river mouth longitude',...
    'delta river mouth latitude',...
    'delta name, where available',...
    'delta BasinID, HydroSheds basin ID',...
    'incoming suspended load sediment, pristine estimate',...
    'Retention efficiency: fraction of incoming suspended load sediment retained within delta topset',...
    'Retention efficiency: fraction of incoming suspended load sediment retained within delta foreset',...
    'incoming suspended load sediment normalized by delta width',...
    'suitable for retention analysis',...
    };

create_netcdf('GlobalDeltaSedimentRetention.nc',out,funits,fmeta)
