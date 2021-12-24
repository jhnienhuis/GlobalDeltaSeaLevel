clr
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','QRiver_prist','QWave','QTide','delta_name','MouthLon','MouthLat','BasinID2','shelf_depth','Discharge_prist');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat')
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel\')
QTide(isnan(QTide)) = 1; QWave(isnan(QWave)| QWave==0) = 1;


subplot(1,2,1)
DeltaSLRlin = linspace(0,2e-2,100);
qs_tot = (365*24*3600*QRiver_dist./1600);
for ii=1:100,
    x = (365*24*3600*QRiver_dist./1600) - ((s-r).*w.*DeltaSLRlin(ii)./2);
    qs_left(ii) = sum(max(0,x(idx)));
    
    %x = qs_tot < ((s-r).*DeltaSLRlin(ii)./2);
    
    perc(ii) = sum(idx & x<0)./sum(idx);
   
end
xlim([0,0.02]),ylim([0 1])
plot(DeltaSLRlin,qs_left./qs_left(1)), hold on, plot(DeltaSLRlin,perc)


%remaining under RCP8.5
slr = DeltaSLR_RCP85_2100+DeltaSub;
perc = sum(qs_tot < (s-r).*slr./2)./numel(qs_tot);
qs_left = sum(max(0,(w.*(qs_tot - (s-r).*slr./2))))./sum(w.*qs_tot);

%no clear relation of land loss w/ delta size
%plot(accumarray(discretize(QRiver_prist(idx_good),logspace(log10(min(QRiver_prist(idx_good))),log10(max(QRiver_prist(idx_good))),100)),delta_change_RCP85_2100(idx_good)./delta_area(idx_good),[],@nanmean))
%accumarray(discretize(QRiver_prist,logspace(log10(min(QRiver_prist)),log10(max(QRiver_prist)),20)),delta_change_RCP85_2100./delta_area,[],@numel)

subplot(1,2,2)

%morphological distribution current SLR
QRiver_SLR = max(0,(QRiver_dist - (1600.*w.*(s-r).*(DeltaSub+DeltaSLR)./2/365/24/3600)));
[~,mor] = max([QRiver_dist,QWave,QTide],[],2);
histcounts(mor(idx))./sum(idx)

%morphological distribution future SLR RCP8.5
QRiver_SLR = max(0,(QRiver_dist - (1600.*w.*(s-r).*(DeltaSub+DeltaSLR_RCP85_2100)/365/24/3600)));
[~,mor] = max([QRiver_SLR,QWave,QTide],[],2);
histcounts(mor)

%plot effect in galloway triangle:

%current SLR
QRiver_SLR = max(0,(QRiver_dist - (1600.*w.*(s-r).*(DeltaSub+DeltaSLR)/365/24/3600./2)));

[QRiver_prist_log,QWave_prist_log,QTide_prist_log] = DeltaLogMaker(QRiver_SLR,QWave,QTide);

%future SLR
QRiver_SLR = max(0,(QRiver_dist - (1600.*w.*(s-r).*(DeltaSub+DeltaSLR_RCP45_2100)/365/24/3600./2)));
[QRiver_dist_log,QWave_dist_log,QTide_dist_log] = DeltaLogMaker(QRiver_SLR,QWave,QTide);

[~,x0,y0] =ternplot(QTide_prist_log,QRiver_prist_log,QWave_prist_log,'scatter','SizeData',5,'Marker','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4]);
[~,x1,y1] =ternplot(QTide_dist_log,QRiver_dist_log,QWave_dist_log,'scatter','SizeData',5,'Marker','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4]);
%close all
x_change = x1-x0;
y_change = y1-y0;
change = idx; % & sqrt(x_change.^2+y_change.^2)>0.05;
ii = [0 logspace(0,3,62) inf];
blub = flipud(cbrewer('div', 'RdYlBu', 64));

colormap(flipud(cbrewer('div', 'RdYlBu', 64)))
ternplot(QTide_prist_log(idx),QRiver_prist_log(idx),QWave_prist_log(idx),'scatter','filled')
hold on

for k= [2:64],
    iplot = QRiver_prist<ii(k) & QRiver_prist>ii(k-1);
    quiver(gca,x0(iplot & change),y0(iplot & change),x_change(iplot & change),y_change(iplot & change),0,'Color',blub(k,:),'LineWidth',max(1.5,k/30),'MaxHeadSize',0.1)
    %scatter(gca,x0(iplot & ~change),y0(iplot & ~change),k/10,blub(k,:),'filled')
end
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'Layer','top')
set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')
colormap(flipud(cbrewer('div', 'RdYlBu', 64)));
h = colorbar;
set(h,'YTick',linspace(0,1,4),'YTickLabel',cellstr(num2str((10.^(0:3))','%1.0f'))); ylabel(h,'QRiver (kgs^{-1})')
%axis off

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'Fig6_Galloway.svg')
