function get_sealevel
%% add slr to deltas, including projections and subsidence
out = load('D:\Dropbox\WorldDeltas\GlobalDeltaData.mat','MouthLon','MouthLat','BasinID');
out.MouthLon(out.MouthLon>360) = out.MouthLon(out.MouthLon>360)-360;
CoorImDelta = out.MouthLat + 1i*out.MouthLon;

%subsidence (%m/yr)
sub = ncread('D:\GlobalDatasets\SeaLevelRise\subs_global2000-2014.nc','subsidence');
sub = (sub(:,:,end)-sub(:,:,1))./14; %m/yr
sub = sub([2161:end 1:2160],end:-1:1);
lat = flipud(ncread('D:\GlobalDatasets\SeaLevelRise\subs_global2000-2014.nc','lat'));
lon = ncread('D:\GlobalDatasets\SeaLevelRise\subs_global2000-2014.nc','lon');
lon = lon([2161:end 1:2160]);

subfilt = sub;

%add a few for which there is no subsidence(?!)
idx = [233,4789,6302,2681,4789,2824,3646,6598,5828,5405,6473,4397, 5390]; %miss
subidx = [5e-3,3e-3,3.3e-2,0,3e-3,1e-2,0,1e-2,0,4e-3,5e-4,3e-3,1e-3];
idx = sub2ind(size(sub),max(1,round(out.MouthLon(idx)*12)),round((90+out.MouthLat(idx))*12));
subfilt(idx) = subidx;

subfilt(subfilt<1e-6) = nan;
subfilt = nanconv(subfilt,ones(5)./25,'same');
subfilt(isnan(subfilt)) = 0;

out.DeltaSub = subfilt(sub2ind(size(subfilt),max(1,round(out.MouthLon*12)),round((90+out.MouthLat)*12)));
out.DeltaSub(5390) = 1e-3; out.DeltaSub(6473) = 1e-3;

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

f = ('D:\GlobalDatasets\SeaLevelRise\SROCC\rsl_ts_26.nc'); %http://icdc.cen.uni-hamburg.de/las/

[out.DeltaSLR_RCP26_2100,out.DeltaSLR_RCP26_tot,out.DeltaSLR_RCP26_low,out.DeltaSLR_RCP26_high] = get_slr_from_source(f,out.MouthLon,out.MouthLat,CoorImDelta);

f = ('D:\GlobalDatasets\SeaLevelRise\SROCC\rsl_ts_45.nc'); %http://icdc.cen.uni-hamburg.de/las/

[out.DeltaSLR_RCP45_2100,out.DeltaSLR_RCP45_tot,out.DeltaSLR_RCP45_low,out.DeltaSLR_RCP45_high] = get_slr_from_source(f,out.MouthLon,out.MouthLat,CoorImDelta);

f = ('D:\GlobalDatasets\SeaLevelRise\SROCC\rsl_ts_85.nc'); %http://icdc.cen.uni-hamburg.de/las/

[out.DeltaSLR_RCP85_2100,out.DeltaSLR_RCP85_tot,out.DeltaSLR_RCP85_low,out.DeltaSLR_RCP85_high] = get_slr_from_source(f,out.MouthLon,out.MouthLat,CoorImDelta);

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


m = axesm('MapProjection','Mercator','MapLonLimit', [0 360]);
sub = movmax(movmax(subfilt,5,1),5,2);
sl = (slr([181:end,1:180],:)'+sub([2173:12:end,1:12:2161],1:12:end)');

map = flipud(cbrewer('div', 'RdYlBu', 64));
map(28,:) = [1 1 1];
blub = gray2ind(mat2gray(sl,[-0.01 0.01]));
blub(blub==27) = 28;
blub(isnan(sl)) = 27;

RGB = ind2rgb(blub,map);
g = geoshow(RGB,[1 90 360]);
colormap(m,map);
axis equal tight
%}


save GlobalDeltaSeaLevelData -struct out

funits = {'dec deg','dec deg','','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr'};
fmeta = {'Delta location','Delta location','Basin ID from hydrosheds','Delta subsidence 2000-2014, doi:10.5194/piahs-372-83-2015',...
    'Delta sea-level rise 1950-2012, https://ws.data.csiro.au/collections/19291/data/4360561',...
    'Delta sea-level rise RCP 2.6, 2090-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 2.6, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 2.6, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 2.6, upper bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, 2090-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 4.5, upper bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, 2090-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta sea-level rise RCP 8.5, upper bound, http://icdc.cen.uni-hamburg.de/las/'};

create_netcdf('GlobalDeltaSeaLevelData.nc',out,funits,fmeta)

end
function [DeltaSLR_RCP_2100,DeltaSLR_RCP_tot,DeltaSLR_RCP_low,DeltaSLR_RCP_high] = get_slr_from_source(f,MouthLon,MouthLat,CoorImDelta);

slr = double(ncread(f,'slr_md')); %m/yr slr
slr_high = double(ncread(f,'slr_he')); %m/yr slr
slr_low = double(ncread(f,'slr_le')); %m/yr slr

DeltaSLR_RCP_2100 = (slr(:,:,end)-slr(:,:,(end-10)))./10;
DeltaSLR_RCP_tot = (slr(:,:,end)-slr(:,:,1))./size(slr,3);
DeltaSLR_RCP_high = (slr_high(:,:,end)-slr_high(:,:,1))./size(slr,3);
DeltaSLR_RCP_low = (slr_low(:,:,end)-slr_low(:,:,1))./size(slr,3);

DeltaSLR_RCP_2100 = DeltaSLR_RCP_2100(sub2ind(size(DeltaSLR_RCP_2100),max(1,round(MouthLon)),round((90+MouthLat))));
DeltaSLR_RCP_tot = DeltaSLR_RCP_tot(sub2ind(size(DeltaSLR_RCP_tot),max(1,round(MouthLon)),round((90+MouthLat))));
DeltaSLR_RCP_high = DeltaSLR_RCP_high(sub2ind(size(DeltaSLR_RCP_high),max(1,round(MouthLon)),round((90+MouthLat))));
DeltaSLR_RCP_low = DeltaSLR_RCP_low(sub2ind(size(DeltaSLR_RCP_low),max(1,round(MouthLon)),round((90+MouthLat))));

[~,idx] = nanmin(abs(CoorImDelta-rot90(999i.*isnan(DeltaSLR_RCP_2100)+CoorImDelta)),[],2);
DeltaSLR_RCP_2100 = DeltaSLR_RCP_2100(idx);
DeltaSLR_RCP_tot = DeltaSLR_RCP_tot(idx);
DeltaSLR_RCP_low = DeltaSLR_RCP_low(idx);
DeltaSLR_RCP_high = DeltaSLR_RCP_high(idx);





end
