function Fig3_LandAreaChange
% continuous plot for SLR rates w/ pristine sediment, disturbed sediment, no sediment.
clr
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','QRiver_prist');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.mat','bed_h','w','r','s')
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel\')
addpath('D:\Dropbox\github\GlobalDeltaChange\')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat')


subplot(1,4,1)
binx = [-10:1:20].*1e-3;
hold on
histogram(DeltaSub(idx),binx)
histogram(DeltaSLR(idx),binx)
histogram(DeltaSLR_RCP26_2100(idx),binx)
histogram(DeltaSLR_RCP45_2100(idx),binx)
histogram(DeltaSLR_RCP85_2100(idx),binx)


subplot(1,4,2)
DeltaSLR_vec = linspace(-1e-2,2e-2,100);

for ii=1:100,
    [land_pred] = get_deltachange(365*24*3600*QRiver_prist./1600,DeltaSLR_vec(ii).*ones(size(w)),s,r,w,bed_h,0.9);
    SLRpred(ii) = nansum(land_pred(idx));
   
end

plot(DeltaSLR_vec.*1e3,SLRpred./1e6)
hold on


%DeltaSLR_vec = linspace(0,2e-2,100);
for ii=1:100,
    [land_pred] = get_deltachange(365*24*3600*QRiver_dist./1600,DeltaSLR_vec(ii).*ones(size(w)),s,r,w,bed_h,0.9);
    SLRpred(ii) = nansum(land_pred(idx));
end

plot(DeltaSLR_vec.*1e3,SLRpred./1e6)
hold on

%DeltaSLR_vec = linspace(0,2e-2,100);
for ii=1:100,
    [land_pred] = get_deltachange(zeros(size(w)),DeltaSLR_vec(ii).*ones(size(w)),s,r,w,bed_h,0.9);
    SLRpred(ii) = nansum(land_pred(idx));
   
end

plot(DeltaSLR_vec.*1e3,SLRpred./1e6)
hold on, grid on

legend('Pristine sediment supply','Modern sediment supply','No sediment supply')
xlabel('SLR Rate (mm/yr)'), ylabel('Delta Land Change (km2 yr-1)')

subplot(1,4,3)
SLRpred = nansum([delta_change_1985_2015 delta_change_RCP26_2100 delta_change_RCP45_2100 delta_change_RCP85_2100]);
SLRpred_unc = 2*nansum([delta_change_1985_2015_std delta_change_RCP26_2100_std delta_change_RCP45_2100_std delta_change_RCP85_2100_std])./sqrt(sum(idx));

bar(SLRpred./1e6)
hold on
errorbar(SLRpred./1e6,SLRpred_unc./1e6)
set(gca,'xticklabel',{'1985-2015','RCP26','RCP45','RCP85'})
ylabel('Delta Land Change (km2/yr)')


subplot(1,4,4)
SLRpred = nansum([delta_change_RCP26_tot delta_change_RCP45_tot delta_change_RCP85_tot]);
SLRpred_unc = 2*nansum([delta_change_RCP26_tot_std delta_change_RCP45_tot_std delta_change_RCP85_tot_std])./sqrt(sum(idx));

bar(SLRpred./1e6)
hold on
errorbar(SLRpred./1e6,SLRpred_unc./1e6)
set(gca,'xticklabel',{'RCP26','RCP45','RCP85'})
ylabel('Delta Land Change (2007-2100)')

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'Fig3_LandAreaChange.svg')

