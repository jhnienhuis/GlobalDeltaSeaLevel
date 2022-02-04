function get_deltaarea
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','Discharge_prist','QRiver_prist','channel_len','channel_len_lat','channel_len_lon','shelf_depth','BasinID2','MouthLon','MouthLat');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile','s')

%area proxy
delta_area = 1.07.*Discharge_prist.^1.1.*QRiver_prist.^0.45./max(50,-shelf_depth).*1e6; %delta area from syvitski2009

%width proxy
w = max(100,sqrt(delta_area./pi).*2);

%get data from edmonds et al., 2020
[ed_ID2,ed_area,ed_width,~,ed_mouth_latlon,ed_apex_latlon,ed_sho1_latlon,ed_sho2_latlon] = get_edmonds_data(BasinID2);

[~,ed_xx] = ismember(ed_ID2,BasinID2);

%then use these values also
w(ed_xx) = ed_width;
delta_area(ed_xx) = ed_area;

%extract delta apex/shoreline/mouth latitude/longitude
mouth_latlon = [MouthLat MouthLon];
apex_latlon_i = max(1,min(size(channel_len,2),size(channel_len,2)-sum(cummax(channel_len,2)>s,2)));
apex_latlon = sub2ind(size(channel_len),(1:length(channel_len))',apex_latlon_i);
apex_latlon = [channel_len_lat(apex_latlon) channel_len_lon(apex_latlon)];
a1 = atan2(MouthLat-apex_latlon(:,1),MouthLon-apex_latlon(:,2));
d1 = sqrt(sum((mouth_latlon-apex_latlon).^2,2));

%where s lat+lon could not be determined accurately, interpolate from first point
idx = apex_latlon_i==1 | (d1./km2deg(w./2000))>10;
apex_latlon(idx,:) = mouth_latlon(idx,:)+[-sin(a1(idx)) -cos(a1(idx))].*km2deg(min(10.*w(idx)./2000,s(idx)./1000));
a1 = atan2(MouthLat-apex_latlon(:,1),MouthLon-apex_latlon(:,2));
d1 = sqrt(sum((mouth_latlon-apex_latlon).^2,2));

d2 = sqrt(d1.^2+km2deg(w./2000).^2);
a2 = atan2(km2deg(w./2000),d1);

sho1_latlon = apex_latlon+[sin(a1+a2) cos(a1+a2)].*d2;
sho2_latlon = apex_latlon+[sin(a1-a2) cos(a1-a2)].*d2;


mouth_latlon(ed_xx,:) = ed_mouth_latlon;
apex_latlon(ed_xx,:) = ed_apex_latlon;
sho1_latlon(ed_xx,:) = ed_sho1_latlon;
sho2_latlon(ed_xx,:) = ed_sho2_latlon;

%ii=1:1000;
%plot([apex_latlon(ii,2) sho1_latlon(ii,2) mouth_latlon(ii,2) sho2_latlon(ii,2)]',[apex_latlon(ii,1) sho1_latlon(ii,1) mouth_latlon(ii,1) sho2_latlon(ii,1)]','o-')

save('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea','BasinID2','w','delta_area','mouth_latlon','apex_latlon','sho1_latlon','sho2_latlon');

%also save netcdf
out = struct('BasinID2',BasinID2,'w', w, 'delta_area',delta_area,'apex_latlon',apex_latlon,'sho1_latlon',sho1_latlon,'sho2_latlon',sho2_latlon,'mouth_latlon',mouth_latlon);

funits = {'','m','m2','dec deg','dec deg','dec deg','dec deg'};
fmeta = {'HydroSheds BasinID','Delta width','delta area','Latitude and longitude of delta apex',...
    'Latitude and longitude of left shoreline','Latitude and longitude of right shoreline',...
    'Latitude and longitude of the river mouth'};
create_netcdf('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.nc',out,funits,fmeta)


