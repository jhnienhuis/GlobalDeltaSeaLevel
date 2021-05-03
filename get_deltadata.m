function get_deltadata
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','Discharge_prist','QRiver_prist','channel_len','shelf_len','shelf_lines','delta_name','shelf_depth','MouthLon','MouthLat','BasinID2');
%get delta profile
beta = nan(size(channel_len,1),1);
alpha = beta;
bed_h = beta;
s = beta;
r = beta;
r_h=beta;
psi = beta;
t = beta;

fo_1 = fitoptions('poly2','Lower',[0,0,0],'Upper',[inf,inf,0]);
fo_2 = fittype('poly2');

for ii=1:length(s)
    if mod(ii,100)==1, ii, end
    
    [beta(ii),alpha(ii),s(ii),r(ii),bed_h(ii),psi(ii),r_h(ii),t(ii)] = get_deltaprofile(channel_len(ii,:),shelf_len(ii,:),shelf_lines,fo_1,fo_2,0);
end

%width
delta_area = 1.07.*Discharge_prist.^1.1.*QRiver_prist.^0.45./max(50,-shelf_depth).*1e6; %delta area from syvitski2009
w = sqrt(delta_area./pi).*2;

%three factors: beta, length, area

%get data from edmonds et al., 2020
[ed_ID2,ed_area,ed_width,ed_length] = get_edmonds_data(BasinID2);

%beta, check against blum&tornqvist
bt_ID2 = [4346011,4165211,2098785,4301321,1248635,4267691];
[~,bt_xx] = ismember(bt_ID2,BasinID2);
bt_beta = [50,20,130,80,25,15].*1e-2./1e3;
%scatter(bt_beta,beta(bt_xx))

%length, check against blum&tornqvist and edmonds
bt_length = [40,90,100,90,150,350].*1e3;

[~,ed_xx] = ismember(ed_ID2,BasinID2);


%then use these values also
s(ed_xx) = ed_length;
w(ed_xx) = ed_width;
beta(bt_xx) = bt_beta;
delta_area(ed_xx) = ed_area;

%have calibrated values and they seem more or less fine. 

save('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaProfile','w','r','s','t','psi','beta','alpha','bed_h','r_h','delta_area');

%also save netcdf
out = struct('alpha', alpha,'beta', beta,'psi',psi,'r', r ,'s', s,'w', w,'bed_h', bed_h,'r_h',r_h,'delta_area',delta_area);

funits = {'','','','m','m','m','m','m','m2'};
fmeta = {'delta surface slope at river mouth', 'Delta basement slope','Delta foreset slope','Distance from delta center to alluvial-basement transition',...
    'Distance from delta center to the shoreline','Delta width','basement depth at shoreline','elevation of the alluvial-bedrock transition','delta area'};
create_netcdf('D:\Dropbox\github\GlobalDeltaSeaLevel\GlobalDeltaProfile.nc',out,funits,fmeta)

