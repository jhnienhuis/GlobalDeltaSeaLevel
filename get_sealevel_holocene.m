load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLon','MouthLat','BasinID2');


t = [0 2 4 6 9 13 17 21]';
SLt = zeros(length(MouthLon),length(t));
SDt = zeros(length(MouthLon),length(t));
f = ['D:\OneDrive - Universiteit Utrecht\SeaLevelRise\holoceneSL\']; %ice5g from Milne, Mitrovica


for ii=2:length(t),
    fid = fopen([f 'mean_rsl.' num2str(t(ii)) '.xyz']); 
    SL = textscan(fid,'%f %f %f'); 
    fclose(fid);
    
    fid = fopen([f 'sd.' num2str(t(ii)) '.xyz']); 
    SD = textscan(fid,'%f %f %f'); 
    fclose(fid);
    
    
    if ii==2,
        idx = dsearchn([SL{1} SL{2}],[MouthLon MouthLat]);
    end
    
    SLt(:,ii) = SL{3}(idx);
    
    SDt(:,ii) = SD{3}(idx);
    
       
end

save('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelHolocene','BasinID2','SLt','SDt','t');

scatter(MouthLon,MouthLat,30,SLt(:,8))
set(gca,'Clim',[-120 -100])