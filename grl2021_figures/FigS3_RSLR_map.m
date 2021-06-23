
ff = 'D:\OneDrive - Universiteit Utrecht\SeaLevelRise\';

%subsidence (%m/yr)
sub = ncread([ff 'subs_global2000-2014.nc'],'subsidence');
sub = (sub(:,:,end)-sub(:,:,1))./14; %m/yr
sub = sub([2161:end 1:2160],end:-1:1);
lat = flipud(ncread([ff 'subs_global2000-2014.nc'],'lat'));
lon = ncread([ff 'subs_global2000-2014.nc'],'lon');
lon = lon([2161:end 1:2160]);

subfilt = sub;
subfilt(subfilt<1e-6) = nan;
subfilt = nanconv(subfilt,ones(5)./25,'same');
subfilt(isnan(subfilt)) = 0;



SLR = load([ff 'HR_Reconstruction_1900-2015.mat'],'RECwithGIA','LRec');
blub = permute(SLR.RECwithGIA(:,1020:end),[2 1]); %get 1985-2015
s = size(blub);
V = bsxfun(@power,(1:s(1))',0:1);
blub = V*pinv(V)*reshape(blub,s(1),[]);
blub = reshape(blub,s);
slr = ((blub(2,:)-blub(1,:))*12/1e3)'; %change to m/yr
slr(end+(1:2)) = [-0.01 0.002];
SLR.LRec(end+(1:2),:) = [-85 60;51 41];
lat = SLR.LRec(:,2);
lon = SLR.LRec(:,1);

%scatter(lon,lat,10,slr,'filled')
%figure,

F = scatteredInterpolant(lon,lat,slr);
[xx,yy] = meshgrid(-180.5:179.5,-72.5:80.5);

sub = movmax(movmax(subfilt,5,1),5,2);
sl = F(xx,yy)+subfilt([2161:12:end,1:12:2161],210:12:2046)';


map = flipud(cbrewer('div', 'RdYlBu', 64));

map(28,:) = [1 1 1];
blub = gray2ind(mat2gray(sl,[-0.01 0.01]));
%blub(blub==27) = 28;
%blub(isnan(sl)) = 27;

RGB = ind2rgb(blub,map);

m = axesm('MapProjection','Mercator','MapLonLimit', [-180 180], 'MapLatLimit',[-72.5 80.5]);
g = geoshow(RGB,[1 80.5 -179.5]);

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(land)




setm(m,'MapLonLimit',[0 360])
colormap(m,map);
axis equal tight