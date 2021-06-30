function get_sealevel
%% add slr to deltas, including projections and subsidence
out = load(['D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat'],'MouthLon','MouthLat','BasinID2','delta_name');
%out.MouthLon(out.MouthLon>180) = out.MouthLon(out.MouthLon>180)-360;
CoorImDelta = out.MouthLat + 1i*out.MouthLon;
ff = 'D:\OneDrive - Universiteit Utrecht\SeaLevelRise\';

%{ 
subsidence (%m/yr) from erkens
sub = ncread([ff 'subs_global2000-2014.nc'],'subsidence');
sub = (sub(:,:,end)-sub(:,:,1))./14; %m/yr
sub = sub([2161:end 1:2160],end:-1:1);
lat = flipud(ncread([ff 'subs_global2000-2014.nc'],'lat'));
lon = ncread([ff 'subs_global2000-2014.nc'],'lon');
lon = lon([2161:end 1:2160]);

subfilt = sub;

%add a few for which there is no subsidence(?!)
idx = [2824,5405]; %miss
subidx = [1e-2,4e-3];
idx = sub2ind(size(sub),max(1,round(out.MouthLon(idx)*12)),round((90+out.MouthLat(idx))*12));
subfilt(idx) = subidx;

subfilt(subfilt<1e-6) = nan;
subfilt = nanconv(subfilt,ones(5)./25,'same');
subfilt(isnan(subfilt)) = 0;

out.DeltaSub = subfilt(sub2ind(size(subfilt),max(1,round(out.MouthLon*12)),round((90+out.MouthLat)*12)));

out.DeltaSub([233,5390,4717,4789,6302,6473,2798,6598,5828,2681,3646,4397]) = [4e-3,1e-3,0.0028,0.001,0.032,1e-3,2e-3,1e-2,5e-3,1e-3,2e-3,5e-3];
%}

load([ff 'Fig3_data_GPS_Subsidence.mat'],'Fig3_data');
sub2 = max(-10,min(10,Fig3_data(:,5)))./1000;
%scatter(Fig3_data(:,1),Fig3_data(:,2),30,sub2,'filled')
CoorImSub = Fig3_data(:,2)+1i*(mod(Fig3_data(:,1)-1,360)+1);

%{
wr = abs(CoorImSub-rot90(CoorImDelta))<2;
out.DeltaSub = zeros(size(out.MouthLon));
for ii=1:length(CoorImDelta),
    if any(wr(:,ii)),
        out.DeltaSub(ii) = mean(-sub2(wr(:,ii)));
    else, 
        out.DeltaSub(ii) = 0;
    end
end
%}

[b,idx] = min(abs(CoorImSub-rot90(CoorImDelta)));

out.DeltaSub = -sub2(idx);
out.DeltaSub(b>1) = 0; %don't include subsidence data where distance to station >1


%{
%slr minus gia from past 50 years.
%gia effect from http://icdc.cen.uni-hamburg.de/las/
gia(:,:,1) = ncread('D:\GlobalDatasets\SeaLevelRise\9EAA985EE5D15A1A297A60B59C954B58_ferret_listing.nc','SELECTED_COMPONENTS'); %2100
gia(:,:,2) = ncread('D:\GlobalDatasets\SeaLevelRise\9EAA985EE5D15A1A297A60B59C954B58_ferret_listing_2007.nc','SELECTED_COMPONENTS'); %2007

gia = (gia(:,:,1)-gia(:,:,2))./(2100-2007); %gia contribution in m/yr

f = 'D:\GlobalDatasets\SeaLevelRise\recons_1950_2012_noib_seasrem.nc'; %http://www.cmar.csiro.au/sealevel/sl_hist_last_decades.html https://ws.data.csiro.au/collections/19291/data/4360561
blub = permute(double(ncread(f,'height')),[3 1 2]); %standard deviation = 0.3 mm/yr (see church et al 2004).


s = size(blub);
V = bsxfun(@power,(1:s(1))',0:1);
blub = V*pinv(V)*reshape(blub,s(1),[]);
blub = reshape(blub,s);
slr(:,25:155) = squeeze(blub(2,:,:)-blub(1,:,:))*12/1e3; %change to m/yr
slr(:,[1:24 156:180]) = nanmean(slr(:));
slr(isnan(slr)) = nanmean(slr(:));
slr = slr+gia;

slrsmooth = nanconv(slr,ones(5)./25,'same');

out.DeltaSLR = slrsmooth(sub2ind(size(slrsmooth),max(1,round(out.MouthLon)),round((90+out.MouthLat))));

[~,idx] = nanmin(abs(CoorImDelta-rot90(999i.*isnan(out.DeltaSLR)+CoorImDelta)),[],2);
out.DeltaSLR = out.DeltaSLR(idx);
%}

SLR = load([ff 'HR_Reconstruction_1900-2015.mat'],'REC','LRec','tTG');
blub = permute(SLR.REC(:,1020:end),[2 1]); %get 1985-2015
s = size(blub);
V = bsxfun(@power,(1:s(1))',0:1);
blub = V*pinv(V)*reshape(blub,s(1),[]);
blub = reshape(blub,s);
slr = ((blub(2,:)-blub(1,:))*12/1e3)'; %change to m/yr

%add hudson bay -10mm/yr and caspian sea = 2mm/yr (1985-2015)
%(https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL073958)
%slr(end+(1:2)) = [-0.01 0.002];
%SLR.LRec(end+(1:2),:) = [-85 60;51 41];

CoorImSLR = SLR.LRec(:,2)+1i*(mod(SLR.LRec(:,1)-1,360)+1);
[~,idx] = min(abs(CoorImSLR-rot90(CoorImDelta)));

out.DeltaSLR = slr(idx);

out.DeltaSLR_series(1:length(idx),2:(length(SLR.tTG)/12)) = diff(SLR.REC(idx,1:12:end),1,2)/1e3;
out.DeltaSLR_time = SLR.tTG(1:12:end)';


f = (['SROCC\rsl_ts_26.nc']); %http://icdc.cen.uni-hamburg.de/las/
[out.DeltaSLR_RCP26_2100,out.DeltaSLR_RCP26_tot,out.DeltaSLR_RCP26_low,out.DeltaSLR_RCP26_high,out.DeltaSLR_RCP26_series,out.DeltaSLR_RCP_time] = get_slr_from_source(ff,f,CoorImDelta);
f = (['SROCC\rsl_ts_45.nc']); %http://icdc.cen.uni-hamburg.de/las/
[out.DeltaSLR_RCP45_2100,out.DeltaSLR_RCP45_tot,out.DeltaSLR_RCP45_low,out.DeltaSLR_RCP45_high,out.DeltaSLR_RCP45_series] = get_slr_from_source(ff,f,CoorImDelta);
f = (['SROCC\rsl_ts_85.nc']); %http://icdc.cen.uni-hamburg.de/las/
[out.DeltaSLR_RCP85_2100,out.DeltaSLR_RCP85_tot,out.DeltaSLR_RCP85_low,out.DeltaSLR_RCP85_high,out.DeltaSLR_RCP85_series] = get_slr_from_source(ff,f,CoorImDelta);

%{
a = tight_subplot(1,3,0.05); 
edges = [-1e-2:5e-4:2e-2];
histogram(a(1),DeltaSub,edges), title(a(1),'Subsidence 2000-2014')
histogram(a(2),DeltaSLR,edges), title(a(2),'SLR Observations 1950-2012')
histogram(a(3),out.DeltaSLR_RCP26_2100,edges,'FaceColor','r'), hold on,
histogram(a(3),out.DeltaSLR_RCP45_2100,edges,'FaceColor','g'), 
histogram(a(3),out.DeltaSLR_RCP85_2100,edges,'FaceColor','b'), title(a(3),'SLR Projections 2090-2100')
legend('RCP 26','RCP 45','RCP 85')
arrayfun(@(x) (box(x,'on')),a)
arrayfun(@(x) (xlim(x,[-1.5e-2,+2e-2])),a)
arrayfun(@(x) (ylim(x,[0 6000])),a)
xlabel(a(1),'Rate (mm/yr)'),xlabel(a(2),'Rate (mm/yr)'),xlabel(a(3),'Rate (mm/yr)')
%}

%example timeseries

%plot(out.DeltaSLR_time,cumsum(-out.DeltaSLR_series(10000,:),2,'reverse'))
%hold on
%plot(out.DeltaSLR_RCP_time,cumsum(out.DeltaSLR_RCP26_series(10000,:)))

save('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData','-struct','out');

funits = {'dec deg','dec deg','','m/yr','m/yr','m/yr','yr','m/yr','m/yr','m/yr','m/yr','m/yr','yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr'};
fmeta = {'Delta location','Delta location','Basin ID2 from hydrosheds','Delta vertical land movement, doi:10.1038/s43017-020-00115-x',...
    'Delta sea-level rise 1985-2015, https://www.nature.com/articles/s41558-019-0531-8',...
    'Delta sea-level rise yearly 1900-2015, https://www.nature.com/articles/s41558-019-0531-8',...
    'Delta sea-level rise yearly time, https://www.nature.com/articles/s41558-019-0531-8',...
    'Delta sea-level rise RCP 2.6, 2081-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 2.6, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 2.6, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 2.6, upper bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 2.6, yearly 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP, yearly time, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, 2081-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, upper bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, yearly 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, 2081-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, upper bound, http://icdc.cen.uni-hamburg.de/las/'...
    'Delta sea-level rise RCP 8.5, yearly 2007-2100, http://icdc.cen.uni-hamburg.de/las/'};
out = rmfield(out,'delta_name');
create_netcdf('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.nc',out,funits,fmeta)

end
function [SLR_2100,SLR_tot,SLR_low,SLR_high,SLR_series,SLR_time] = get_slr_from_source(ff,f,CoorImDelta)

gia = ncread([ff 'SROCC\gia_mean.nc'],'rslrate'); %subtract GIA because it is in the VLM dataset

SLRraw = permute(double(ncread([ff f],'slr_md')),[3 1 2]); %m
SLRraw_high = permute(double(ncread([ff f],'slr_he')),[3 1 2]); %m slr
SLRraw_low = permute(double(ncread([ff f],'slr_le')),[3 1 2]); %m slr
SLR_time = ncread([ff f],'time')';

[Y,X] = meshgrid(ncread([ff f],'y'),ncread([ff f],'x'));

SLR_2100 = squeeze((SLRraw(end,:,:)-SLRraw(end-19,:,:))./19)-gia; %m/yr
SLR_tot = squeeze((SLRraw(end,:,:)-SLRraw(1,:,:))./length(SLR_time))-gia; %m/yr
SLR_high = squeeze((SLRraw_high(end,:,:)-SLRraw_high(1,:,:)))./length(SLR_time);
SLR_low = squeeze((SLRraw_low(end,:,:)-SLRraw_low(1,:,:)))./length(SLR_time);

ind_sea = find(~isnan(SLRraw(1,:,:)));
CoorImSLR = Y(ind_sea)+1i.*X(ind_sea);

[~,idx] = min(abs(CoorImSLR-rot90(CoorImDelta)));

SLR_2100 = SLR_2100(idx)';
SLR_tot = SLR_tot(idx)';
SLR_low = SLR_low(idx)';
SLR_high = SLR_high(idx)';

SLR_series = diff(SLRraw(:,idx),1,1)';
SLR_time(end) = [];

end
