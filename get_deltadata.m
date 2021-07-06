function get_deltadata
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','Discharge_prist','QRiver_prist','channel_len','channel_len_lat','channel_len_lon','shelf_len','shelf_lines','delta_name','shelf_depth','BasinID2','depth_upstream','MouthLon','MouthLat');
%get delta profile
beta = nan(size(channel_len,1),1);
alpha = beta;
bed_h = beta;
s = beta;
r = beta;
r_h=beta;
psi = beta;
t = beta;

fo_1 = fitoptions('poly2','Lower',[0,0,0],'Upper',[inf,inf,0]);
fo_2 = fittype('poly2');

for ii=1:length(s)
    if mod(ii,100)==1, ii, end

    [beta(ii),alpha(ii),s(ii),r(ii),bed_h(ii),psi(ii),r_h(ii),t(ii)] = get_deltaprofile(channel_len(ii,:),shelf_len(ii,:),shelf_lines,fo_1,fo_2,0);
end

%set sedimentary depth of delta wedge to be at least the channel depth
bed_h = min(bed_h,-depth_upstream);

%width
delta_area = 1.07.*Discharge_prist.^1.1.*QRiver_prist.^0.45./max(50,-shelf_depth).*1e6; %delta area from syvitski2009
w = max(100,sqrt(delta_area./pi).*2);

%three factors: beta, length, area

%get data from edmonds et al., 2020
[ed_ID2,ed_area,ed_width,ed_length,ed_mouth_latlon,ed_apex_latlon,ed_sho1_latlon,ed_sho2_latlon] = get_edmonds_data(BasinID2);

%beta, check against blum&tornqvist
bt_ID2 = [4346011,4165211,2098785,4301321,1248635,4267691];
[~,bt_xx] = ismember(bt_ID2,BasinID2);
bt_beta = [50,20,130,80,25,15].*1e-2./1e3;
%scatter(bt_beta,beta(bt_xx))

%length, check against blum&tornqvist and edmonds
bt_length = [40,90,100,90,150,350].*1e3;

[~,ed_xx] = ismember(ed_ID2,BasinID2);
%scatter(delta_area(ed_xx),ed_area)

%then use these values also
s(ed_xx) = ed_length;
w(ed_xx) = ed_width;
beta(bt_xx) = bt_beta;
delta_area(ed_xx) = ed_area;
%have calibrated values and they seem more or less fine. 


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

save('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile','w','r','s','t','psi','beta','alpha','bed_h','r_h','delta_area','mouth_latlon','apex_latlon','sho1_latlon','sho2_latlon');

%also save netcdf
out = struct('alpha', alpha,'beta', beta,'psi',psi,'r', r ,'s', s,'w', w,'bed_h', bed_h,'r_h',r_h,'delta_area',delta_area,'apex_latlon',apex_latlon,'sho1_latlon',sho1_latlon,'sho2_latlon',sho2_latlon,'mouth_latlon',mouth_latlon);

funits = {'','','','m','m','m','m','m','m2','dec deg','dec deg','dec deg','dec deg'};
fmeta = {'delta surface slope at river mouth', 'Delta basement slope','Delta foreset slope','Distance from delta center to alluvial-basement transition',...
    'Distance from delta center to the shoreline','Delta width','basement depth at shoreline','elevation of the alluvial-bedrock transition','delta area',...
    'Latitude and longitude of delta apex','Latitude and longitude of left shoreline','Latitude and longitude of right shoreline','Latitude and longitude of the river mouth'};
create_netcdf('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.nc',out,funits,fmeta)

