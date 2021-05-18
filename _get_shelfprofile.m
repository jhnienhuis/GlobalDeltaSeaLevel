function get_shelfprofile
%get continental shelf depth test

%load world bathymetry
[Z, refvec] = etopo('D:\OneDrive - Universiteit Utrecht\GlobalDEM\',5,[-90 90]);

%hoy many km per cell latitude (simple!)
kmperlatcell = deg2km(1)./refvec(1);

%how many km per cell longitude (depends on latitude)
kmperloncell = kmperlatcell*cos(linspace(-0.5*pi,0.5*pi,size(Z,1)));

%find the coastline contour (read the contourc to see what the output looks like
C = contourc(Z,[0 0]);

%uncomment to plot 0 and 200m contour line


%find all instances of the coastline (all "islands")
% I = find(C(1,:)==0);
% 
% hold on
% for ii=1:length(I)-1,
%     plot(C(1,(1+I(ii)):(I(ii+1)-1)),C(2,(1+I(ii)):(I(ii+1)-1)))
% end
% 
% %plot 200 m contour
% C200 = contourc(Z,[-200 -200]);
% I = find(C200(1,:)==-200);
% hold on
% for ii=1:length(I)-1,
%     plot(C200(1,(1+I(ii)):(I(ii+1)-1)),C200(2,(1+I(ii)):(I(ii+1)-1)),'b')
% end

%loop through the entire coastline to find the nearest offshore contour
%lines
C2 = C(1,C(1,:)~=0)+1i*C(2,C(1,:)~=0);

%do only every 10 points along the coast.
C2 = C2(1:10:end);

%do profiles for these contour lines
c_lines = 0:-20:-400;

%pre-allocate matrices
d = zeros(length(C2),length(c_lines));
lat = zeros(length(C2),length(c_lines));
lon = zeros(length(C2),length(c_lines));



%get shoreline lat and lon
lat(:,1) = imag(C2);
lon(:,1) = real(C2);

%find minimum distance to the next contour line for every contour
for jj=2:length(c_lines),
    jj
    S = contourc(Z,[c_lines(jj) c_lines(jj)]);
    S2 = S(1,S(1,:)~=c_lines(jj))+1i*S(2,S(1,:)~=c_lines(jj));
    
    idx = zeros(length(C2),1);
    %for all the shoreline positions
    for ii=1:length(C2),
        
        %distance to contour in km (simple minimum of the complex magnitude
        %of the difference the two contour lines
        [d(ii,jj),idx(ii)] = min(abs((kmperloncell(uint16(imag(C2(ii))))+1i*kmperlatcell).*(S2-C2(ii))));
        
    end
    
    %use new contour locations as the next place to start finding minima.
    C2 = S2(idx);
    lat(:,jj) = imag(S2(idx));
    lon(:,jj) = real(S2(idx));
    
end

%turn matrix indices into latitude and longitude
lat = (lat./refvec(1))-90;
lon = lon./refvec(1);

d = cumsum(d,2);

%remove distances away from shoreline larger than 400km? (great lakes etc)
lat(d>500) = nan;
lon(d>500) = nan;
d(d>500) = nan;

%save everything into mat file
%save GlobalDeltaShelf c_lines d lat lon
%%

%load GlobalDeltaShelf
%see what all the transects are like!
%hold on
%for ii=1:length(lat),
%    plot(lon(ii,:),lat(ii,:));
%end

% plot a particular transect
% figure
% ii=200;
% plot(d(ii,:),c_lines,'-o')
% xlabel('Offshore Distance (km)')
% ylabel('Depth (m)')
% hold on

%plot pchip fitted line
% plot(0:0.1:d(ii,end),pchip(d(ii,:),c_lines,0:0.1:d(ii,end)),'r');

%what kind of metrics to establish from this? (needs improvement!)

%width in m of the shelf?
ShelfWidth = zeros(length(lat),1);

%depth of shelf break?
ShelfBreakDepth = zeros(length(lat),1);


for ii=1:length(lat)
    
    %don't do computation of shelfwidth and depth if all nans (>500km)
    l = find(~isnan(d(ii,:)),1,'last');
    if l<2, continue, end
    
    %fit a smooth line through the data
    SlopeSmooth = pchip(d(ii,1:l),c_lines(1:l),0:0.1:d(ii,l));
    
    %find where the second derivative is minimum
    [~,ShelfWidth(ii)] = min(diff(SlopeSmooth,2));
    
    %get the depth at which this derivative is at a minimum
    ShelfBreakDepth(ii) = SlopeSmooth(ShelfWidth(ii));
end
%convert back to m (smoothed line was at every 100m)
ShelfWidth = ShelfWidth.*100;

%simple way of estimating slope 
ShelfSlope = (c_lines(7)./d(:,7))./1000;
% figure
% histogram(log10(-ShelfSlope),100)

%uncomment these lines to get a complicated way to estimate the
%slope (derivative at the coast) (doesnt work as well)
%{
ShelfSlope = zeros(size(lat));
for ii=1:length(lat),
    pp = pchip(d(ii,:),c_lines);
    pp = polyder(pp.coefs(1,:));
    ShelfSlope(ii) = pp(end)./1000;
end
%}
%save GlobalDeltaShelf -append ShelfSlope ShelfWidth ShelfBreakDepth
%make plot of shelfwidth, shelf slope, and shelfbreakdepth
% figure
% a = tight_subplot(3,1,0.05,0.05);
% scatter(a(1),lon(:,1),lat(:,1),15,log10(-ShelfSlope),'filled')
% scatter(a(2),lon(:,1),lat(:,1),15,ShelfWidth,'filled'), set(gca,'CLim',[0 500])
% scatter(a(3),lon(:,1),lat(:,1),10,ShelfBreakDepth,'filled'), set(gca,'CLim',[-300 0])
% box(a(1),'on'),box(a(2),'on'),box(a(3),'on')
% set(a,'DataAspectRatio',[1 1 1],'XLim',[0 360],'YLim',[-90 90])
% colorbar(a(1)),colorbar(a(2)),colorbar(a(3))
% title(a(1),'ShelfSlope (log10)')
% title(a(2),'ShelfWidth (km)')
% title(a(3),'ShelfBreakDepth (m)')


%% add shelf to global delta dataset
out = load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat'],'MouthLon','MouthLat');
%load('GlobalDeltaShelf.mat',)

[~,idx] = min(abs(out.MouthLon'-lon(:,1)+1i*(out.MouthLat'-lat(:,1))));
shelf_len = d(idx,:);
shelf_len_lon = lon(idx,:);
shelf_len_lat = lat(idx,:);
shelf_depth = ShelfBreakDepth(idx);
shelf_width = ShelfWidth(idx);
shelf_slope = ShelfSlope(idx);

save GlobalDeltaShelfProfile shelf_len shelf_len_lon shelf_len_lat shelf_depth shelf_width shelf_slope

%s = geoshape(num2cell(out.shelf_len_lat,2),num2cell(out.shelf_len_lon,2),'Elevation',num2cell(out.shelf_depth,2));
%kmlwrite('ShelfTransects.kml',s,'name',num2cell(num2str(out.BasinID),2))

%get basin_depth_100k
y = 0:20:400;
for ii=1:length(shelf_len),
    x = shelf_len(ii,~isnan(shelf_len(ii,:)));
    if length(x)<2, basin_depth_50k(ii) = 10;
    else, basin_depth_50k(ii) = interp1(x,y(1:length(x)),50,'linear','extrap');
    end
end
basin_depth_50k = max(10,min(500,basin_depth_50k));

save GlobalDeltaShelfProfile basin_depth_50k -append
