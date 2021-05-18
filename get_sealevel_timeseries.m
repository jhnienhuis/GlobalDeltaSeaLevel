function get_sealevel_timeseries
%% add slr to deltas, including projections and subsidence
out = load('D:\Dropbox\WorldDeltas\scripts\GlobalDeltaData.mat','MouthLon','MouthLat','BasinID');
out.MouthLon(out.MouthLon>360) = out.MouthLon(out.MouthLon>360)-360;
CoorImDelta = out.MouthLat + 1i*out.MouthLon;
ff = 'D:\OneDrive - Universiteit Utrecht\SeaLevelRise\';

%subsidence (%m/yr)
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


%slr minus gia from past 50 years.
%gia effect from http://icdc.cen.uni-hamburg.de/las/
gia(:,:,1) = ncread([ff '9EAA985EE5D15A1A297A60B59C954B58_ferret_listing.nc'],'SELECTED_COMPONENTS'); %2100
gia(:,:,2) = ncread([ff '9EAA985EE5D15A1A297A60B59C954B58_ferret_listing_2007.nc'],'SELECTED_COMPONENTS'); %2007

gia = (gia(:,:,1)-gia(:,:,2))./(2100-2007); %gia contribution in m/yr

%http://www.cmar.csiro.au/sealevel/sl_hist_last_decades.html https://ws.data.csiro.au/collections/19291/data/4360561
blub = permute(double(ncread([ff 'recons_1950_2012_noib_seasrem.nc'],'height')),[3 1 2]); %standard deviation = 0.3 mm/yr (see church et al 2004).
t = ncread([ff 'recons_1950_2012_noib_seasrem.nc'],'year');
slr(:,:,25:155) = blub;
slr(:,:,[1:24 156:180]) = repmat(nanmean(slr,[2 3]),[1 360 49]); %add average slr to out of bounds

%convert to m/yr
slr2 = zeros(63,360,180);
for ii=1:63,
    slr2(ii,:,:) = mean(slr(t==(1949+ii),:,:),1);
end

slr2 = padarray(slr2,[2 2 2],'replicate');
slrsmooth = zeros(size(slr2));
for ii=1:size(slr2,1),
    slrsmooth(ii,:,:) = nanconv(squeeze(slr2(ii,:,:)),ones(5)./25,'same');
end
slrsmooth = convn(slrsmooth,ones(5,1,1)./5,'same');

slrsmooth = slrsmooth(3:end-2,3:end-2,3:end-2);


slrsmooth = nansum(cat(4,(diff(slrsmooth,1,1)./1000),repmat(permute(gia,[3 1 2]),[62 1 1])),4);
slrsmooth(:,isnan(gia)) = nan;
%imagesc(1:360,1:180,squeeze(mean(slrsmooth,1))'), axis xy, hold on
%scatter(out.MouthLon,out.MouthLat,10,nanmean(out.DeltaSLR,1),'filled')
%scatter(max(1,ceil(out.MouthLon)),ceil((90+out.MouthLat)))

x = ~isnan(slrsmooth(1,:,:));
[latx,lonx] = meshgrid(-89:90,1:360); coorx_i = latx+1i.*lonx;
slrsmooth = slrsmooth(:,x);
coorx_i = coorx_i(x);

slrsmooth(:,end+(1:2)) = repmat([-0.01 0.002],62,1);
coorx_i(end+(1:2)) = [60+275i;41+51i];

[~,xx] = min(abs(CoorImDelta-rot90(coorx_i)),[],2);
out.DeltaSLR = slrsmooth(:,xx)';

%out.DeltaSLR = slrsmooth(:,sub2ind([360 180],max(1,ceil(out.MouthLon)),ceil((90+out.MouthLat))))';
%[~,idx] = nanmin(abs(CoorImDelta-rot90(999i.*isnan(out.DeltaSLR(:,1))+CoorImDelta)),[],2);
%out.DeltaSLR = out.DeltaSLR(idx,:);


f = ([ff 'SROCC\rsl_ts_26.nc']); %http://icdc.cen.uni-hamburg.de/las/

out.DeltaSLR_RCP26 = get_slr_from_source(f,out.MouthLon,out.MouthLat,CoorImDelta);

f = ([ff 'SROCC\rsl_ts_45.nc']); %http://icdc.cen.uni-hamburg.de/las/

out.DeltaSLR_RCP45 = get_slr_from_source(f,out.MouthLon,out.MouthLat,CoorImDelta);

f = ([ff 'SROCC\rsl_ts_85.nc']); %http://icdc.cen.uni-hamburg.de/las/

out.DeltaSLR_RCP85 = get_slr_from_source(f,out.MouthLon,out.MouthLat,CoorImDelta);

out.t_historic = 1950:2011;
out.t_future = 2012:2100;

save GlobalDeltaSeaLevelTimeseries -struct out
%plot([out.t_historic out.t_future],cumsum([out.DeltaSLR(1000,:) out.DeltaSLR_RCP85(1000,:)]))
end
function [DeltaSLR_RCP] = get_slr_from_source(f,MouthLon,MouthLat,CoorImDelta);

slr = permute(double(ncread(f,'slr_md')),[3 1 2]); 
slr = diff(slr(5:end,:,:),1,1);%m/yr slr
DeltaSLR_RCP = slr(:,sub2ind([360 180],max(1,round(MouthLon)),round((90+MouthLat))))';
[~,idx] = nanmin(abs(CoorImDelta-rot90(999i.*isnan(DeltaSLR_RCP(:,1))+CoorImDelta)),[],2);
DeltaSLR_RCP = DeltaSLR_RCP(idx,:);


end
