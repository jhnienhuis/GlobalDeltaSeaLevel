function create_shapefile

%create shapefile with delta extent, and add a csv file with timeseries of
%sea-level/area change?

% load all deltas
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat'],'BasinID2','delta_name')

%load profile data
load([dropbox filesep 'github' filesep 'GlobalDeltaSeaLevel' filesep 'export_data' filesep 'GlobalDeltaProfile.mat'])

%load response data
slr = load([dropbox filesep 'github' filesep 'GlobalDeltaSeaLevel' filesep 'export_data' filesep 'GlobalDeltaSeaLevelData.mat']);


%load response data
res = load([dropbox filesep 'github' filesep 'GlobalDeltaSeaLevel' filesep 'export_data' filesep 'GlobalDeltaSeaLevelResponse.mat']);


delta_area_lat = cell(size(BasinID2));
delta_area_lon = cell(size(BasinID2));
% write to shapefiles
for ii=1:length(BasinID2),
    delta_area_lat{ii} = [apex_latlon(ii,1) sho1_latlon(ii,1) mouth_latlon(ii,1) sho2_latlon(ii,1)];
    delta_area_lon{ii} = [apex_latlon(ii,2) sho1_latlon(ii,2) mouth_latlon(ii,2) sho2_latlon(ii,2)];
end


%shapefile limited to 10 character attribute names..
p = geoshape(delta_area_lat(res.idx),delta_area_lon(res.idx),'BasinID2',double(BasinID2(res.idx)),'delta_name',char(delta_name(res.idx)));
p.Geometry = 'polygon';

p.dAslr_h = res.delta_change_1985_2015(res.idx);
p.dAslr_rcp26 = res.delta_change_RCP26_2100(res.idx);
p.dAslr_rcp45 = res.delta_change_RCP45_2100(res.idx);
p.dAslr_rcp85 = res.delta_change_RCP85_2100(res.idx);

fname = 'D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea';
shapewrite(p,fname)
zip([fname '_shp'],{[fname '.dbf'],[fname '.shx'],[fname '.shp']})
delete([fname '.dbf'],[fname '.shx'],[fname '.shp'])


%add sl historic timeseries as csv files
%adapt to show cumulative
slr.DeltaSLR_series = cumsum(-slr.DeltaSLR_series,2,'reverse');
slr.DeltaSLR_series = slr.DeltaSLR_series-slr.DeltaSLR_series(:,108); %put 2007 as reference year
slr.DeltaSLR_RCP26_series = cumsum(slr.DeltaSLR_RCP26_series,2);
slr.DeltaSLR_RCP45_series = cumsum(slr.DeltaSLR_RCP45_series,2);
slr.DeltaSLR_RCP85_series = cumsum(slr.DeltaSLR_RCP85_series,2);

T = array2table(slr.DeltaSLR_series(res.idx,:));
T.Properties.VariableNames = cellstr(num2str(slr.DeltaSLR_time'));
T = [table(BasinID2(res.idx)) T];
writetable(T,'GlobalDelta_DeltaSLR.csv')

T = array2table(slr.DeltaSLR_RCP26_series(res.idx,:));
T.Properties.VariableNames = cellstr(num2str(slr.DeltaSLR_RCP_time'));
T = [table(BasinID2(res.idx)) T];
writetable(T,'GlobalDelta_DeltaSLR_RCP26.csv')

T = array2table(slr.DeltaSLR_RCP45_series(res.idx,:));
T.Properties.VariableNames = cellstr(num2str(slr.DeltaSLR_RCP_time'));
T = [table(BasinID2(res.idx)) T];
writetable(T,'GlobalDelta_DeltaSLR_RCP45.csv')

T = array2table(slr.DeltaSLR_RCP85_series(res.idx,:));
T.Properties.VariableNames = cellstr(num2str(slr.DeltaSLR_RCP_time'));
T = [table(BasinID2(res.idx)) T];
writetable(T,'GlobalDelta_DeltaSLR_RCP85.csv')





%show 
%{
idx=5390;
plot(slr.DeltaSLR_time,slr.DeltaSLR_series(idx,:))
hold on
plot(slr.DeltaSLR_RCP_time,slr.DeltaSLR_RCP26_series(idx,:))
plot(slr.DeltaSLR_RCP_time,slr.DeltaSLR_RCP85_series(idx,:))
%}

