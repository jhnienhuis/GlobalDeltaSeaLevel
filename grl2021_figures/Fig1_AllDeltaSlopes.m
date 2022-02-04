%figure 1
out = load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLat','MouthLon');
%load('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaProfile.mat','bed_h')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaSeaLevelResponse.mat','idx')


axesm('Robinson','MapLonLimit', [0 360])
setm(gca,'Origin',[0 0])
land = shaperead('landareas.shp', 'UseGeoCoords', true);

for ii=1:length(land),
    [land(ii).Lat, land(ii).Lon] = reducem(land(ii).Lat', land(ii).Lon',1);
end

rivers = shaperead('worldrivers', 'UseGeoCoords', true);
hold on
geoshow(land, 'FaceColor', [1 0.92 0.8]);
axis tight

geoshow(rivers,'Color','blue');

scatterm(rem(out.MouthLat(idx)+180,360)-180,out.MouthLon(idx),'filled','r')
%plot(out.channel_len_lon',out.channel_len_lat','-b')
%plot(out.shelf_len_lon',out.shelf_len_lat','-g')
