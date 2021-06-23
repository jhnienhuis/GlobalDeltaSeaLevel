clr
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','QWave','QTide','QRiver_prist','delta_name','MouthLon','MouthLat','BasinID2','shelf_depth','Discharge_prist','depth_upstream');
ee = load('D:\Dropbox\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaSeaLevelData.mat','DeltaSLR','DeltaSub')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaProfile.mat','bed_h','s','r','w','delta_area')
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel\')
addpath('D:\Dropbox\github\GlobalDeltaChange\')



land_obs = mean([ee.net_pekel ee.net_aqua],2).*1e6;


fr=1;

[land_pred,idx] = get_deltachange(365*24*3600*QRiver_dist./1600,DeltaSLR+DeltaSub,s,r,w,bed_h,fr);  %has r2 of 0.38

%remove caspian sea
idx= idx & ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48);

%first, global delta land?
total_land = s.*0.5.*w; %edmonds 2020 847000km2
nansum(total_land(idx))./1e6

%how much global delta land loss
nansum(land_pred(idx))./1e6
sum(land_obs(idx))./1e6

%how good is this model?
corrcoef(land_obs(idx),land_pred(idx)).^2

%subplot(1,2,1)

scatter(land_pred(idx)./1e6,land_obs(idx)./1e6,(w(idx).*s(idx))./2e8,'filled','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k');

hold on, plot([-20 30],[-20 30],'--k')
%xlim([-10 20]),ylim([-10 20]), set(gca,'PlotBoxAspectRatio',[1 1 1])
box on, grid on, xlabel('Predicted delta area change (km^2/yr)')
ylabel('Observed delta area change (km^2/yr)')
RMSE = sqrt(mean((land_obs(idx)-land_pred(idx)).^2))./1e6

%idx = cellfun(@(x) (find(BasinID==x,1)),delta_name_id); table(delta_name,idx,DeltaSub(idx),ee.net_aqua(idx),land_pred(idx)./1e6,w(idx))
%{
subplot(1,2,2)
idx = find(idx);
for ii=1:100,
    [~,idxp] = maxk(QRiver_prist(idx),ii*10);
    c = corrcoef(land_obs(idx(idxp)),land_pred(idx(idxp))).^2;
    cx(ii) = c(1,2);
end
plot((1:100)*10/100,cx,'-o'), xlabel('Largest x% of deltas'),ylabel('r^2')
ylim([0 1]),set(gca,'PlotBoxAspectRatio',[1 1 1])
%}

%gname(delta_name(idx))

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
%saveas(gcf,'Fig2_Validation.svg')