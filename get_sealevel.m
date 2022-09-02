function get_sealevel
%% add slr to deltas, including projections and subsidence
out = load(['D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat'],'MouthLon','MouthLat','BasinID2','delta_name');
CoorImDelta = out.MouthLat + 1i*out.MouthLon;
ff = 'D:\OneDrive - Universiteit Utrecht\SeaLevelRise\';

%% Holocene SL

tt = [21 11 8 7 6];

for ii=1:length(tt),
    f = [ff 'anu_rsl\rsl.' num2str(tt(ii)) '.xyz'];
    gunzip([f '.gz'])
    fid = fopen(f); 
    SL = textscan(fid,'%f %f %f'); 
    fclose(fid);
    delete(f)
    
    F = scatteredInterpolant(SL{1},SL{2},SL{3});

    
    [xx,yy] = meshgrid(0:360,-90:90);
    RSL7 = F(xx,yy)./7;
    %sl_dif_8 = gray2ind(mat2gray(-RSL7,[-10 10]));

        
    SLt(:,ii) = -F(mod(out.MouthLon,360),out.MouthLat);
       
end

out.HoloceneSL_time = tt.*-1000;
out.HoloceneSL = SLt;

%% subsidence, incl GIA
load([ff 'Fig3_data_GPS_Subsidence.mat'],'Fig3_data');
sub2 = max(-10,min(10,Fig3_data(:,5)))./1000;
%scatter(Fig3_data(:,1),Fig3_data(:,2),30,sub2,'filled')
CoorImSub = Fig3_data(:,2)+1i*(mod(Fig3_data(:,1)-1,360)+1);
%find closest gps data for each delta
[b,idx] = min(abs(CoorImSub-rot90(CoorImDelta)));

out.DeltaSub = -sub2(idx);

%subtract GIA from GPS rates (GIA from sea-level reconstructions is likely to be better)
gia = ncread([ff 'SROCC\gia_mean.nc'],'rslrate'); %subtract GIA because it is in the VLM dataset
gia = gia(sub2ind(size(gia),floor(out.MouthLon),round(90+out.MouthLat)));
out.DeltaSub(b>1) = gia(b>1); %don't include subsidence data where distance to station >1

out.GIA = gia;

%% 20th century SL
SLR = load([ff 'HR_Reconstruction_1900-2015.mat'],'REC','LRec','tTG');
blub = permute(SLR.REC(:,1:end),[2 1]);
s = size(blub);
V = bsxfun(@power,(1:s(1))',0:1);
blub = V*pinv(V)*reshape(blub,s(1),[]);
blub = reshape(blub,s);
slr = ((blub(2,:)-blub(1,:))*12/1e3)'; %change to m/yr

CoorImSLR = SLR.LRec(:,2)+1i*(mod(SLR.LRec(:,1)-1,360)+1);
[~,idx] = min(abs(CoorImSLR-rot90(CoorImDelta)));


out.DeltaSLR = slr(idx);
out.DeltaSLR_series(1:length(idx),2:(length(SLR.tTG)/12)) = diff(SLR.REC(idx,1:12:end),1,2)/1e3;
out.DeltaSLR_time = SLR.tTG(1:12:end)';

%% future SLR

%SROCC data
f = (['SROCC\rsl_ts_26.nc']); %http://icdc.cen.uni-hamburg.de/las/
[out.DeltaSLR_RCP26_2100,out.DeltaSLR_RCP26_tot,out.DeltaSLR_RCP26_low,out.DeltaSLR_RCP26_high,out.DeltaSLR_RCP26_series,out.DeltaSLR_RCP_time] = get_slr_from_source(ff,f,CoorImDelta);
f = (['SROCC\rsl_ts_45.nc']); %http://icdc.cen.uni-hamburg.de/las/
[out.DeltaSLR_RCP45_2100,out.DeltaSLR_RCP45_tot,out.DeltaSLR_RCP45_low,out.DeltaSLR_RCP45_high,out.DeltaSLR_RCP45_series] = get_slr_from_source(ff,f,CoorImDelta);
f = (['SROCC\rsl_ts_85.nc']); %http://icdc.cen.uni-hamburg.de/las/
[out.DeltaSLR_RCP85_2100,out.DeltaSLR_RCP85_tot,out.DeltaSLR_RCP85_low,out.DeltaSLR_RCP85_high,out.DeltaSLR_RCP85_series] = get_slr_from_source(ff,f,CoorImDelta);

%add slangen timeseries up to 2300
[out.DeltaSLR_SSPt,out.DeltaSLR_SSP126,out.DeltaSLR_SSP245,out.DeltaSLR_SSP585] = get_ar6_from_source(out.MouthLat,out.MouthLon);


%plot(out.DeltaSLR_time,cumsum(-out.DeltaSLR_series(10000,:),2,'reverse'))
%hold on
%plot(out.DeltaSLR_RCP_time,cumsum(out.DeltaSLR_RCP26_series(10000,:)))

save('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData','-struct','out');

funits = {'dec deg','dec deg','','yr BP','m','m/yr','m/yr','m/yr','m/yr','yr','m/yr','m/yr','m/yr','m/yr','m/yr','yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','m/yr','yr','m/yr','m/yr','m/yr'};
fmeta = {'Delta location','Delta location','Basin ID2 from hydrosheds',...
    'Delta Holocene SLt, time in year BP',...
    'Delta Holocene SL, sea level relative to present from ANU ice model',...
    'Delta vertical land movement, incl GIA, doi:10.1038/s43017-020-00115-x',...
    'Delta modern GIA, http://icdc.cen.uni-hamburg.de/las/',...,
    'Delta SLR 1900-2015, https://www.nature.com/articles/s41558-019-0531-8',...
    'Delta SLR yearly 1900-2015, https://www.nature.com/articles/s41558-019-0531-8',...
    'Delta SLR yearly time, https://www.nature.com/articles/s41558-019-0531-8',...
    'Delta SLR RCP 2.6, 2081-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 2.6, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 2.6, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 2.6, upper bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 2.6, yearly 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP, yearly time, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 4.5, 2081-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 4.5, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 4.5, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 4.5, upper bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 4.5, yearly 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 8.5, 2081-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 8.5, 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 8.5, lower bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 8.5, upper bound, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR RCP 8.5, yearly 2007-2100, http://icdc.cen.uni-hamburg.de/las/',...
    'Delta SLR SSP times',...
    'Delta SLR SSP126, 2000-2050 2050-2100 2100-2200 2200-2300, IPCC AR6',...
    'Delta SLR SSP245, 2000-2050 2050-2100 2100-2200 2200-2300, IPCC AR6',...
    'Delta SLR SSP585, 2000-2050 2050-2100 2100-2200 2200-2300, IPCC AR6'
    };
out = rmfield(out,'delta_name');
create_netcdf('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.nc',out,funits,fmeta)

end

function [SLR_2100,SLR_tot,SLR_low,SLR_high,SLR_series,SLR_time] = get_slr_from_source(ff,f,CoorImDelta)

%gia = ncread([ff 'SROCC\gia_mean.nc'],'rslrate'); %subtract GIA because it is in the VLM dataset

SLRraw = permute(double(ncread([ff f],'slr_md')),[3 1 2]); %m
SLRraw_high = permute(double(ncread([ff f],'slr_he')),[3 1 2]); %m slr
SLRraw_low = permute(double(ncread([ff f],'slr_le')),[3 1 2]); %m slr
SLR_time = ncread([ff f],'time')';

[Y,X] = meshgrid(ncread([ff f],'y'),ncread([ff f],'x'));

SLR_2100 = squeeze((SLRraw(end,:,:)-SLRraw(end-19,:,:))./19); %m/yr
SLR_tot = squeeze((SLRraw(end,:,:)-SLRraw(1,:,:))./length(SLR_time)); %m/yr
SLR_high = squeeze((SLRraw_high(end,:,:)-SLRraw_high(1,:,:)))./length(SLR_time);
SLR_low = squeeze((SLRraw_low(end,:,:)-SLRraw_low(1,:,:)))./length(SLR_time);

ind_sea = find(~isnan(SLRraw(1,:,:)));
CoorImSLR = Y(ind_sea)+1i.*X(ind_sea);

[~,idx] = min(abs(CoorImSLR-rot90(CoorImDelta)));

SLR_2100 = SLR_2100(ind_sea(idx));
SLR_tot = SLR_tot(ind_sea(idx));
SLR_low = SLR_low(ind_sea(idx));
SLR_high = SLR_high(ind_sea(idx));
%gia = gia(ind_sea(idx));

SLR_series = diff(SLRraw(:,ind_sea(idx)),1,1)'; %-repmat(gia,1,length(SLR_time)-1)
SLR_time(end) = [];

end


function [t, DeltaSLR_SSP126,DeltaSLR_SSP245,DeltaSLR_SSP585] = get_ar6_from_source(lat,lon)

%lat = out.MouthLat;
lon = mod(lon,360);

load('D:\Drive\2022 DeltaSLR_AnnualRev\Data\gridded-SL-IPCC-AR6.mat')

lon_grid = mod(lon_grid,360);

lon_grid(361,:) = 360.*ones(size(lon_grid(1,:)));
lat_grid(361,:) = lat_grid(1,:);

k = dsearchn([lon_grid(:),lat_grid(:)],[lon lat]); 

t = [2050 2100 2200 2300];

DeltaSLR_SSP126 = permute(diff(cat(3,zeros(size(SL_126LC_2050_50)), SL_126LC_2050_50, SL_126LC_2100_50, SL_126LC_2200_50, SL_126LC_2300_50),1,3),[3 1 2]);
DeltaSLR_SSP245 = permute(diff(cat(3,zeros(size(SL_245LC_2050_50)), SL_245LC_2050_50, SL_245LC_2100_50, SL_245LC_2200_50, SL_245LC_2300_50),1,3),[3 1 2]);
DeltaSLR_SSP585 = permute(diff(cat(3,zeros(size(SL_585LC_2050_50)), SL_585LC_2050_50, SL_585LC_2100_50, SL_585LC_2200_50, SL_585LC_2300_50),1,3),[3 1 2]);

DeltaSLR_SSP126(:,361,:) = DeltaSLR_SSP126(:,1,:);
DeltaSLR_SSP245(:,361,:) = DeltaSLR_SSP245(:,1,:);
DeltaSLR_SSP585(:,361,:) = DeltaSLR_SSP585(:,1,:);

DeltaSLR_SSP126 = (DeltaSLR_SSP126(:,k)./[45.5; 50; 100; 100]/1000)';
DeltaSLR_SSP245 = (DeltaSLR_SSP245(:,k)./[45.5; 50; 100; 100]/1000)';
DeltaSLR_SSP585 = (DeltaSLR_SSP585(:,k)./[45.5; 50; 100; 100]/1000)';


end