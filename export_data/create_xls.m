load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','delta_name','BasinID2','BasinArea','Discharge_prist');

load('D:\Dropbox\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat','net_pekel','net_aqua');

load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat')

T = table(delta_name(idx),int64(BasinID2(idx)),DeltaSLR(idx),DeltaSub(idx),DeltaSLR_RCP26_2100(idx),DeltaSLR_RCP45_2100(idx),DeltaSLR_RCP85_2100(idx),...
    delta_change_1985_2015(idx),delta_change_RCP26_2100(idx),delta_change_RCP45_2100(idx),delta_change_RCP85_2100(idx),...
'VariableNames',{'Delta Name','BasinID2','SLR_1985_2015 (m/yr)','VLM (m/yr, positive is downward)','SLR_2081_2100_RCP26 (m/yr)','SLR_2081_2100_RCP45 (m/yr)','SLR_2081_2100_RCP85 (m/yr)','DeltaChange_1985_2015_model (m2/yr, positive is land gain)','DeltaChange_2081_2100_RCP26 (m2/yr, positive is land gain)','DeltaChange_2081_2100_RCP45 (m2/yr, positive is land gain)','DeltaChange_2081_2100_RCP85 (m2/yr, positive is land gain)'});

writetable(T,'GlobalDeltaSeaLevel.xls');


idx = find(idx); 
[~,I] = maxk(Discharge_prist(idx),100); 
idx = idx(I);

T = table(delta_name(idx),int64(BasinID2(idx)),MouthLat(idx),MouthLon(idx),DeltaSLR(idx),DeltaSub(idx),DeltaSLR_RCP26_2100(idx),DeltaSLR_RCP45_2100(idx),DeltaSLR_RCP85_2100(idx),...
    mean([net_pekel(idx),net_aqua(idx)],2).*1e6, delta_change_1985_2015(idx),delta_change_RCP26_2100(idx),delta_change_RCP45_2100(idx),delta_change_RCP85_2100(idx),...
'VariableNames',{'Delta Name','BasinID2','Latitude','Longitude','SLR_1985_2015 (m/yr)','VLM (m/yr, positive is subsidence)','SLR_2081_2100_RCP26 (m/yr)','SLR_2081_2100_RCP45 (m/yr)','SLR_2081_2100_RCP85 (m/yr)','DeltaChange_1985_2015_obs (m2/yr, positive is land gain)','DeltaChange_1985_2015_model (m2/yr, positive is land gain)','DeltaChange_2081_2100_RCP26 (m2/yr, positive is land gain)','DeltaChange_2081_2100_RCP45 (m2/yr, positive is land gain)','DeltaChange_2081_2100_RCP85 (m2/yr, positive is land gain)'});

writetable(T,'GlobalDeltaSeaLevel_max100.xls');

load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','delta_name');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat')

T = table(delta_name,int64(BasinID2),delta_area,...
'VariableNames',{'Delta Name','BasinID2','Delta Area (m2)'});

writetable(T,'GlobalDeltaArea.xls');