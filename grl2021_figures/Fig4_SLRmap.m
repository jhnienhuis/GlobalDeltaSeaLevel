clr
out = load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat');

load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLon','MouthLat','BasinID2');
ee = load('D:\Drive\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat');



out.delta_change_obs = mean([ee.net_aqua ee.net_pekel],2).*1e6;
out.delta_change_obs(out.delta_change_obs==0) = 0.01;

m = {'obs','RCP26_2100','RCP45_2100','RCP85_2100'}; %1985_2015
a = tight_subplot(2,2);
land = shaperead('landareas.shp', 'UseGeoCoords', true);
[Lat, Lon] = reducem([land(:).Lat]', [land(:).Lon]',0.2);
land = geoshape(Lat,Lon,'Geometry','Polygon');

MouthLon(MouthLon>180) = MouthLon(MouthLon>180) - 360;

for ii=1:4,
    
axes(a(ii))

%rivers = shaperead('worldrivers', 'UseGeoCoords', true);
title(m{ii},'interpreter','none')
%setm(gca,'Origin',[0 180])
axesm('MapProjection','pcarree','MapLonLimit', [-180 180],'MapLatLimit',[-60 90])
axis tight
hold on
geoshow(land, 'FaceColor', [1 0.92 0.8]);
x = out.(['delta_change_' m{ii}]);
scatterm(MouthLat(out.idx),MouthLon(out.idx),5*abs(x(out.idx))./1e6,sign(x(out.idx)),'filled')
%colormap((cbrewer('div', 'RdYlGn', 64)))
colormap([0.5 0 0;0 0.8 0]);
%set(a(ii),'CLim',[-2 2])

if ii==3,
   scatterm([-20 -10 0],[-170 -170 -170],15*abs([1 10 50]),sign([1 10 50]),'filled','MarkerEdgeColor','k')
end
   
end


set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
%saveas(gcf,'Fig5_SLRmap.svg')