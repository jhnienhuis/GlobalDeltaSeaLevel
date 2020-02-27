function get_deltaresponse
load GlobalDeltaSeaLevelData
load GlobalDeltaProfile
load('D:\Dropbox\WorldDeltas\GlobalDeltaData.mat','QRiver_dist','Discharge_prist');


%yearly end-of-century SLR
slr_scenarios = {'DeltaSLR','DeltaSLR_RCP26_2100','DeltaSLR_RCP45_2100','DeltaSLR_RCP85_2100'};
slr_unc_scenarios =  {'3e-4.*ones(size(w))','DeltaSLR_RCP26_high','DeltaSLR_RCP45_high','DeltaSLR_RCP85_high'};

outcomes1 = {'delta_change_1985_2015','delta_change_RCP26_2100','delta_change_RCP45_2100','delta_change_RCP85_2100'};
outcomes2 = {'delta_change_1985_2015_std','delta_change_RCP26_2100_std','delta_change_RCP45_2100_std','delta_change_RCP85_2100_std'};

SLRpred = zeros(size(slr_scenarios));
SLRpred_unc = SLRpred;

for jj=1:length(slr_scenarios),
    slr = eval(slr_scenarios{jj})+DeltaSub;
    [ds,dr] = get_deltachange(qs,365*24*3600*QRiver_dist./1600./w,slr,s,r,beta,alpha,bed_h,1);
    SLRpred(jj) = nansum(w.*ds./1e6);
    
    land_area = get_deltachange_montecarlo(eval(slr_scenarios{jj}),eval(slr_unc_scenarios{jj}),w,DeltaSub,bed_h,QRiver_dist,Discharge_prist,rab,beta,s,r,alpha);
    SLRpred_unc(jj) = std(sum(land_area,1)./1e6);
    
    out.(outcomes1{jj}) = w.*ds./1e6;
    out.(outcomes2{jj}) = std(land_area,1,2)./1e6;
end


%cumulative SLR
slr_scenarios = {'DeltaSLR_RCP26_tot','DeltaSLR_RCP45_tot','DeltaSLR_RCP85_tot'};
slr_unc_scenarios =  {'3e-4.*ones(size(w))','DeltaSLR_RCP26_high','DeltaSLR_RCP45_high','DeltaSLR_RCP85_high'};

outcomes1 = {'delta_change_RCP26_tot','delta_change_RCP45_tot','delta_change_RCP85_tot'};
outcomes2 = {'delta_change_RCP26_tot_std','delta_change_RCP45_tot_std','delta_change_RCP85_tot_std'};

SLRpred = zeros(size(slr_scenarios));
SLRpred_unc = SLRpred;
lt = length(ncread('D:\GlobalDatasets\SeaLevelRise\SROCC\rsl_ts_26.nc','time'));

for jj=1:length(slr_scenarios),
    slr = eval(slr_scenarios{jj})+DeltaSub;
    [ds,dr] = get_deltachange(qs,365*24*3600*QRiver_dist./1600./w,slr,s,r,beta,alpha,bed_h,1);
    SLRpred(jj) = nansum(w.*ds./1e6).*lt;
    
    land_area = get_deltachange_montecarlo(eval(slr_scenarios{jj}),eval(slr_unc_scenarios{jj}),w,DeltaSub,bed_h,QRiver_dist,Discharge_prist,rab,beta,s,r,alpha).*lt;
    SLRpred_unc(jj) = std(sum(land_area,1)./1e6);
    
    out.(outcomes1{jj}) = w.*ds./1e6.*lt;
    out.(outcomes2{jj}) = std(land_area,1,2)./1e6;
    
end

save GlobalDeltaSeaLevelResponse -struct out

%also save netcdf
funits = repmat({'km2/yr'},1,14);
fmeta = repmat({'Land Area Change prediction per delta','Land area change standard deviation (from Monte Carlo-assessment)'},1,7);
create_netcdf('GlobalDeltaSeaLevelResponse.nc',out,funits,fmeta)