clr
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','QRiver_prist','QWave','QTide','delta_name','Discharge_prist','shelf_depth','BasinID2','Continent');
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat')
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel\')
addpath('D:\Dropbox\github\GlobalDeltaChange')

%numbers of deltas in study
sum(idx)
f = 365*24*3600./1600;
%representing x percentage of delta land
delta_area = s.*w*0.5;
nansum(delta_area(idx))./nansum(delta_area)

%yearly delta land loss for rcp scenarios
nansum(delta_change_RCP26_2100(idx))./1e6
2.*nansum(delta_change_RCP26_2100_std(idx))./1e6./sqrt(sum(idx))

nansum(delta_change_RCP85_2100(idx))./1e6
2.*nansum(delta_change_RCP85_2100_std(idx))./1e6./sqrt(sum(idx))

%cumulative land loss for rcp8.5
nansum(delta_change_RCP85_tot(idx))./1e6
2.*nansum(delta_change_RCP85_tot_std(idx))./1e6./sqrt(sum(idx))
%worst case delta land loss as a fraction of modern delta land
nansum(delta_change_RCP85_tot(idx))./nansum(delta_area(idx))

%percentage caused by slr compared to other effects

land_pred(1) = nansum(get_deltachange(365*24*3600*QRiver_prist./1600,DeltaSLR_RCP85_2100,s,r,w,bed_h,0.9).*idx)./1e6;
land_pred(2) = nansum(get_deltachange(365*24*3600*min(QRiver_prist,QRiver_dist)./1600,0,s,r,w,bed_h,0.9).*idx)./1e6;
land_pred(3) = nansum(get_deltachange(365*24*3600*QRiver_prist./1600,DeltaSub,s,r,w,bed_h,0.9).*idx)./1e6;
pristine = nansum(get_deltachange(365*24*3600*QRiver_prist./1600,0,s,r,w,bed_h,0.9).*idx)./1e6;
(pristine-land_pred(1))/sum(pristine-land_pred)

%number of deltas from doug vs proxy method:
[dd] = get_edmonds_data(BasinID2,0);
[dd] = ismember(BasinID2,dd);

sum(idx & dd)
sum(idx & ~dd)

nansum(delta_area(idx & dd))/nansum(delta_area)

%projections without sediment supply:
nansum(get_deltachange(zeros(size(s)),DeltaSub+DeltaSLR,s,r,w,bed_h,0.9).*idx)./1e6

%sensitivity large vs small deltas
fr_loss = delta_change_RCP85_tot./delta_area;
idx2 = (idx & fr_loss>-1 & fr_loss<1);
corrcoef(delta_area(idx2),fr_loss(idx2))

%arctic deltas
mean(DeltaSLR_RCP85_2100(idx & Continent==8))
mean(DeltaSLR_RCP85_2100(idx & Continent~=8))

nansum(delta_change_RCP85_tot(idx & Continent==8))./nansum(delta_area(idx & Continent==8))
nansum(delta_change_RCP85_tot(idx & Continent~=8))./nansum(delta_area(idx & Continent~=8))

nansum(QRiver_prist(idx & Continent==8))./nansum(delta_area(idx & Continent==8))
nansum(QRiver_prist(idx & Continent~=8))./nansum(delta_area(idx & Continent~=8))

%individual drivers, uninhibited delta growth
nansum(get_deltachange(365*24*3600*QRiver_prist./1600,zeros(size(s)),s,r,w,bed_h,0.9).*idx)./1e6

%rcp8.5 number of deltas in autoretreat
sum((f*QRiver_dist(idx))<(delta_area(idx).*DeltaSLR_RCP85_2100(idx)))./sum(idx)

%delta accomodation (m3/yr) compared to fluvial sed flux for 10mm/yr
sum(delta_area(idx).*1e-2)./sum(f.*QRiver_dist(idx))

%flux trapped:
1-(sum(max(0,f*QRiver_dist(idx)-(delta_area(idx).*DeltaSLR(idx))))./sum(f*QRiver_dist(idx)))
1-(sum(max(0,f*QRiver_dist(idx)-(delta_area(idx).*DeltaSLR_RCP85_2100(idx))))./sum(f*QRiver_dist(idx)))


%morphologic modern SLR;
[~,mor1] = max([QRiver_dist(idx)-(delta_area(idx).*DeltaSLR(idx)./f),QWave(idx),QTide(idx)],[],2);
[~,mor2] = max([QRiver_dist(idx)-(delta_area(idx).*DeltaSLR_RCP45_2100(idx)./f),QWave(idx),QTide(idx)],[],2);
(histcounts(mor2)-histcounts(mor1))./sum(mor1==1)

%fraction of change by delta size
[~,blub] = sort(delta_area(idx),'descend');
xx = abs(delta_change_RCP85_tot(idx));
nansum(xx(blub(1:round(length(blub)*0.01))))./nansum(xx)

