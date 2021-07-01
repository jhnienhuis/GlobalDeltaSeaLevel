%% sensitivity to delta width? medium
clr
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','QRiver_prist','delta_name');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat')
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel\')
w_r = (0.5:0.1:2).*w; land_pred = [];

for ii=1:size(w_r,2),
    
    land_pred(ii) = nansum(get_deltachange(365*24*3600*QRiver_dist./1600,DeltaSLR+DeltaSub,s,r,w_r(:,ii),bed_h).*idx);

end

figure,plot((0.5:0.1:2),land_pred./land_pred(6),'-o'), hold on, grid on

%% sensitivity to bed_depth? medium
bed_h_r = (0.5:0.1:2).*bed_h; land_pred = [];

for ii=1:size(w_r,2),
    
    land_pred(ii) = nansum(get_deltachange(365*24*3600*QRiver_dist./1600,DeltaSLR+DeltaSub,s,r,w,bed_h_r(:,ii)).*idx);

end
plot((0.5:0.1:2),land_pred./land_pred(6),'-o')
%% sensitivity to slr
slr_r = (0.5:0.1:2).*(DeltaSub+DeltaSLR); land_pred = [];

for ii=1:size(w_r,2),
    
    land_pred(ii) = nansum(get_deltachange(365*24*3600*QRiver_dist./1600,slr_r(:,ii),s,r,w,bed_h).*idx);

end

plot((0.5:0.1:2),land_pred./land_pred(6),'-o')
%% sensitivity to Qriver
qs_r = (0.5:0.1:2).*QRiver_dist; land_pred = [];

for ii=1:size(w_r,2),
    
    land_pred(ii) = nansum(get_deltachange(365*24*3600*qs_r(:,ii)./1600,DeltaSLR+DeltaSub,s,r,w,bed_h).*idx);

end

plot((0.5:0.1:2),land_pred./land_pred(6),'-o')

%% sensitivity to delta length
s_r = (0.5:0.1:2).*s; land_pred = [];

for ii=1:size(w_r,2),
    
    land_pred(ii) = nansum(get_deltachange(365*24*3600*QRiver_dist./1600,DeltaSLR+DeltaSub,s_r(:,ii),r,w,bed_h).*idx);

end

plot((0.5:0.1:2),land_pred./land_pred(6),'-o')

ylabel('Delta land change (compared to default estimate)')
xlabel('Parameter (compared to default estimate)')
legend('Delta Width','Delta Toe Depth','SLR','QRiver','Delta Length')

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'FigS2_deltachange_uncertainty.svg')