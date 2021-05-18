function get_deltaresponse
load('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaSeaLevelData')
load('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaProfile','s','r','w','bed_h')
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist');

%yearly end-of-century SLR
slr_scenarios = {'DeltaSLR','DeltaSLR_RCP26_2100','DeltaSLR_RCP45_2100','DeltaSLR_RCP85_2100'};
slr_unc_scenarios =  {'3e-4.*ones(size(w))','DeltaSLR_RCP26_high','DeltaSLR_RCP45_high','DeltaSLR_RCP85_high'};

outcomes1 = {'delta_change_1985_2015','delta_change_RCP26_2100','delta_change_RCP45_2100','delta_change_RCP85_2100'};
outcomes2 = {'delta_change_1985_2015_std','delta_change_RCP26_2100_std','delta_change_RCP45_2100_std','delta_change_RCP85_2100_std'};

ff = 365*24*3600/1600;
fr = 0.9;

%remove caspian sea and other non deltas
[~,idx] = get_deltachange(ff.*QRiver_dist,DeltaSLR+DeltaSub,s,r,w,bed_h,fr);
out.idx= idx & ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48);


for jj=1:length(slr_scenarios),
    slr = eval(slr_scenarios{jj})+DeltaSub;
    dA = get_deltachange(ff.*QRiver_dist,slr,s,r,w,bed_h,fr);
    dA(~out.idx) = nan;   
    dA_unc = get_deltachange_montecarlo(QRiver_dist,eval(slr_scenarios{jj}),DeltaSub,eval(slr_unc_scenarios{jj}),s,r,w,bed_h);
    dA_unc(~out.idx,:) = nan;
    
    out.(outcomes1{jj}) = dA;
    out.(outcomes2{jj}) = double(std(dA_unc,1,2));
    nansum(dA)
end




%cumulative SLR
slr_scenarios = {'DeltaSLR_RCP26_tot','DeltaSLR_RCP45_tot','DeltaSLR_RCP85_tot'};
slr_unc_scenarios =  {'DeltaSLR_RCP26_high','DeltaSLR_RCP45_high','DeltaSLR_RCP85_high'};

outcomes1 = {'delta_change_RCP26_tot','delta_change_RCP45_tot','delta_change_RCP85_tot'};
outcomes2 = {'delta_change_RCP26_tot_std','delta_change_RCP45_tot_std','delta_change_RCP85_tot_std'};

lt = length(ncread('D:\OneDrive - Universiteit Utrecht\SeaLevelRise\SROCC\rsl_ts_26.nc','time'));

for jj=1:length(slr_scenarios),
    slr = eval(slr_scenarios{jj})+DeltaSub;
    dA = get_deltachange(ff.*QRiver_dist,slr,s,r,w,bed_h,fr);
    dA(~out.idx) = nan;   
    dA_unc = get_deltachange_montecarlo(QRiver_dist,eval(slr_scenarios{jj}),DeltaSub,eval(slr_unc_scenarios{jj}),s,r,w,bed_h);
    dA_unc(~out.idx,:) = nan;
    
    out.(outcomes1{jj}) = dA.*lt;
    out.(outcomes2{jj}) = double(std(dA_unc,1,2).*lt);
    
end



save('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaSeaLevelResponse','-struct','out')

out.idx = double(out.idx);

%also save netcdf
funits = [{' '},repmat({'m2/yr'},1,14)];
fmeta = [{'included deltas'},repmat({'Land Area Change prediction per delta','Land area change standard deviation (from Monte Carlo-assessment)'},1,7)];
create_netcdf('GlobalDeltaSeaLevelResponse.nc',out,funits,fmeta)