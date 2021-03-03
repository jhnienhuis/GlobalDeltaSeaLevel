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
fo_1 = fitoptions('poly2','Lower',[0,0,0],'Upper',[inf,inf,0]);
fo_2 = fittype('poly2');
for idx=1:length(s)
    if mod(idx,100)==1, 
        idx, 
    end
    
    [beta(idx),alpha(idx),s(idx),r(idx),bed_h(idx),psi(idx),r_h(idx)] = get_deltaprofile(channel_len(idx,:),shelf_len(idx,:),shelf_lines,fo_1,fo_2,0);
end

delta_area_proxy = 1.07.*Discharge_prist.^1.1.*QRiver_prist.^0.45./max(100,-shelf_depth).*1e6; %delta area from syvitski2009

w = sqrt(delta_area_proxy./pi)*2;

%get data from edmonds et al., 2020
[delta_area,delta_width,delta_length] = get_delta_area(Discharge_prist,QRiver_prist,delta_name,shelf_depth,MouthLon,MouthLat,BasinID2);

corrcoef(w(~isnan(delta_width)),delta_width(~isnan(delta_width)))
corrcoef(s(~isnan(delta_length)),delta_length(~isnan(delta_length)))

w(~isnan(delta_width)) = delta_width(~isnan(delta_width));
s(~isnan(delta_length)) = delta_length(~isnan(delta_length));



save GlobalDeltaProfile psi w r s beta alpha bed_h r_h delta_area

%also save netcdf
out = struct('alpha', alpha,'beta', beta,'psi',psi,'r', r ,'s', s,'w', w,'bed_h', bed_h,'r_h',r_h,'delta_area',delta_area);

funits = {'','','','m','m','m','m','m','m2'};
fmeta = {'delta surface slope at river mouth', 'Delta basement slope','Delta foreset slope','Distance from delta center to alluvial-basement transition',...
    'Distance from delta center to the shoreline','Delta width','basement depth at shoreline','elevation of the alluvial-bedrock transition','delta area'};
create_netcdf('GlobalDeltaProfile.nc',out,funits,fmeta)

