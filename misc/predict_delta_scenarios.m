% load stuff
load('GlobalDeltaSeaLevelTimeseries.mat')
load('D:\Dropbox\WorldDeltaLandArea\scripts\GlobalDeltaProfile.mat','qs','bed_h','s','r','psi','beta','w','alpha');
load Qriver_scenarios 
load('D:\Dropbox\WorldDeltas\scripts\GlobalDeltaData.mat','delta_name')
qs_bar = cat(3,qs_bar(:).QRiver);

%qs_scenario (1 through 12)
qs_scenario=12
qs_sus = 365*24*3600*qs_bar(:,:,qs_scenario)./1600./w; %m2/yr

%time period
t = 1980:2099; %yr

%slr scenario (26,45,85)
slr = [DeltaSLR(:,t_historic>=t(1)), DeltaSLR_RCP26(:,t_future<=t(end))]; %m/yr




delta_change = zeros(size(slr));
ds=0; dr=0;

for ii=1:length(t),
    %update delta size after each iteration
    s = s+ds; r = r+dr;
    %delta shoreline and upstream change in m/yr
    [ds,dr] = get_deltachange(qs,qs_sus(:,ii),slr(:,ii)+DeltaSub,s,r,beta,alpha,psi,bed_h,1);
    %delta shoreline change in m2/yr
    delta_change(:,ii) = w.*ds;
    
end

%all deltas
%plot(t,cumsum(sum(delta_change,1)./1e6))
%xlabel('time (yr)'), ylabel('delta area, cumulative change since 1980 (km^2)'), hold on

%specific delta
plot(t,cumsum(delta_change(delta_name=='Ebro',:)./1e6))
xlabel('time (yr)'), ylabel('delta area, cumulative change since 1980 (km^2)'), hold on

%tune bedload down