% bar plot linear decomposition of forcings (subsidence, dams, SLR)
clr
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','QRiver_prist','delta_name','MouthLon','MouthLat','BasinID2','shelf_depth','Discharge_prist');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat','idx')
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel\')

%idx mrd = 233 %rhine=3646 %mekong=5390 %volta=967 %nile 1168 %mahakam 6089, ganges = 4717


%simple back-of-the-envelope
%billion ton per year:
nansum(delta_area(idx).*1e-3).*1600./1e3./1e9

Qsc = {'QRiver_prist','min(QRiver_prist,QRiver_dist)','max(QRiver_prist,QRiver_dist)','QRiver_prist','QRiver_prist','QRiver_dist','QRiver_prist','QRiver_dist','QRiver_prist','QRiver_dist','QRiver_prist','QRiver_dist'};
Slrsc = {'zeros(size(w))','zeros(size(w))','zeros(size(w))','zeros(size(w))','DeltaSLR','DeltaSLR','DeltaSLR_RCP26_2100','DeltaSLR_RCP26_2100','DeltaSLR_RCP45_2100','DeltaSLR_RCP45_2100','DeltaSLR_RCP85_2100','DeltaSLR_RCP85_2100'};
Slrunsc = {'zeros(size(w))','zeros(size(w))','zeros(size(w))','zeros(size(w))','3e-4.*ones(size(w))','3e-4.*ones(size(w))','DeltaSLR_RCP26_high','DeltaSLR_RCP26_high','DeltaSLR_RCP45_high','DeltaSLR_RCP45_high','DeltaSLR_RCP85_high','DeltaSLR_RCP85_high'};
Subsc = {'zeros(size(w))','zeros(size(w))','zeros(size(w))','DeltaSub','zeros(size(w))','DeltaSub','zeros(size(w))','DeltaSub','zeros(size(w))','DeltaSub','zeros(size(w))','DeltaSub'};

for ii=1:12,
    SLRpred(ii) = nansum(get_deltachange(365*24*3600*eval(Qsc{ii})./1600,eval(Slrsc{ii})+eval(Subsc{ii}),s,r,w,bed_h,0.9).*idx);
    SLRpred_unc(ii) = nansum(std(get_deltachange_montecarlo(eval(Qsc{ii}),eval(Slrsc{ii}),eval(Subsc{ii}),eval(Slrunsc{ii}),s,r,w,bed_h),1,2).*idx)./sqrt(sum(idx));
    
end
% figure
SLRpred([2:5,7,9,11]) = SLRpred([2:5,7,9,11])-SLRpred(1);

figure
b = bar(SLRpred./1e6);
b.FaceColor = 'flat';
b.CData(7,:) = [0 1 0]; b.CData(8,:) = [0 1 0]; b.CData(9,:) = [1 1 0]; b.CData(10,:) = [1 1 0];  b.CData(11,:) = [1 0 0]; b.CData(12,:) = [1 0 0];
set(gca,'xticklabel',{'Pristine','Dams','Deforestation','Subsidence','SLR','Delta Change',...
    'RCP26 \newline 2090-2100','Change RCP26 \newline 2090-2100',...
    'RCP45 \newline 2090-2100','Change RCP45 \newline 2090-2100',...
    'RCP85 \newline 2090-2100','Change RCP85 \newline 2090-2100'})
hold on
errorbar(SLRpred./1e6,2.*SLRpred_unc./1e6)
ylabel('Delta Land Change (km2/yr)')

set(gcf,'PaperOrientation','Landscape')
title('Relative contribution of SLR and sediment on global delta land gain')

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'Fig4_Decomposition.svg')
