function [ed_BasinID2,ed_area,ed_width,ed_length,ed_mouth_latlon,ed_apex_latlon,ed_sho1_latlon,ed_sho2_latlon,ed_apex_ele] = get_edmonds_data(BasinID2,get_elevation)

data = xlsread('Edmondsetal2020_NatCom_suppdata.xlsx','B2:K2177');
ed_area = xlsread('Edmondsetal2020_NatCom_suppdata.xlsx','Q2:Q2177').*1e6; %area in m2
ed_BasinID2 = xlsread('Edmondsetal2020_NatCom_suppdata.xlsx','AF2:AF2177');

ed_BasinID2(isnan(ed_area)) = 0;

[~,idx] = ismember(ed_BasinID2,BasinID2);

idx = find(idx);
ed_BasinID2 = ed_BasinID2(idx);

if nargin==1,
    get_elevation=0;
end

%ed_area = zeros(size(idx));
%for ii=1:length(idx)
%ed_area(ii) = areaint([data(idx(ii),1) data(idx(ii),5) data(idx(ii),3) data(idx(ii),7)] ,[data(idx(ii),2) data(idx(ii),6) data(idx(ii),4) data(idx(ii),8)],referenceSphere('earth'));
%end


ed_area = ed_area(idx);

ed_width = zeros(size(idx));
for ii=1:length(idx)
ed_width(ii) = distance(data(idx(ii),5),data(idx(ii),6),data(idx(ii),7),data(idx(ii),8),referenceSphere('earth'));
end
ed_length = zeros(size(idx));
for ii=1:length(idx)
ed_length(ii) = distance(data(idx(ii),1),data(idx(ii),2),data(idx(ii),3),data(idx(ii),4),referenceSphere('earth'));
end

ed_mouth_latlon = data(idx,1:2); 
ed_apex_latlon = data(idx,3:4);
ed_sho1_latlon = data(idx,5:6);
ed_sho2_latlon = data(idx,7:8);

%ed_mouth_latlon(:,2) = mod(ed_mouth_latlon(:,2)-1,360)+1;
%ed_apex_latlon(:,2) = mod(ed_apex_latlon(:,2)-1,360)+1;
%ed_sho1_latlon(:,2) = mod(ed_sho1_latlon(:,2)-1,360)+1;
%ed_sho2_latlon(:,2) = mod(ed_sho2_latlon(:,2)-1,360)+1;


if get_elevation,
    load('D:\OneDrive - Universiteit Utrecht\GlobalDEM\SRTM15plus_int8.mat','cz');
    cz = cz(:,[43201:end 1:43200]);
    lat_grid = int64(linspace(-90*240+0.5,(90*240)-0.5,180*240));
    lon_grid = int64(linspace((-180*240)+0.5,(180*240)-0.5,360*240));
    [~,idx_lat] = ismember(round(data(idx,3)*240),lat_grid);
    [~,idx_lon] = ismember(round(data(idx,4)*240),lon_grid);
    
    idx_lin = sub2ind(size(cz),idx_lat,idx_lon);
    
    ed_apex_ele = cz(idx_lin);
else,
    ed_apex_ele = 0;
end


%{
old method
delta_area_proxy = 1.07.*Discharge_prist.^1.1.*QRiver_prist.^0.45./max(100,-shelf_depth).*1e6; %delta area from syvitski2009
%xx = sort(delta_area,'descend'); yy = sort(delta_area_proxy,'descend');
%plot(xx(1:2000),yy(1:2000)), grid on, axis equal
%set(gca,'XScale','log','YScale','log')


%do the largest delta first, then work backwards
[~,delta_idx] = sort(doug_delta_area,'descend');
idx2 = zeros(size(doug_delta_area));


%define matching parameter
delta_coor = MouthLon+1i.*MouthLat;
dougdelta_coor = data(:,2)+1i.*data(:,1);
%do separate function for large and small delta
blub = abs(rot90(dougdelta_coor)-delta_coor);
%blub(:,delta_idx(1:200)) = abs(rot90(dougdelta_coor(delta_idx(1:200)))-delta_coor) + abs(log10(doug_delta_area(delta_idx(1:200)))-log10(delta_area_proxy)) + (abs(rot90(dougdelta_coor(delta_idx(1:200)))-delta_coor)>4).*10;
%blub(:,delta_idx(201:end)) = abs(rot90(dougdelta_coor(delta_idx(201:end)))-delta_coor) + abs(log10(doug_delta_area(delta_idx(201:end)))-log10(delta_area_proxy)) + (abs(rot90(dougdelta_coor(delta_idx(201:end)))-delta_coor)>1).*10;

%histogram(min(blub,[],2))
%


for ii=delta_idx
    [d,idx2(ii)] = min(blub(:,ii));
    %d
    %don't match if matching parameter exceeds 50
    if d>5, idx2(ii) = 0; continue, end
    
    %set blub to infinite to make sure a POC flux is not matched twice
    blub(idx2(ii),:) = inf;
    doug_BasinID2(ii) = BasinID2(idx2(ii));
    
    
end


plot(delta_coor,'o'), hold on, plot(dougdelta_coor,'or')
hold on
plot([data(idx2>0,2) MouthLon(idx2(idx2>0))]',[data(idx2>0,1) MouthLat(idx2(idx2>0))]')

sum(QRiver_prist(idx2(idx2>0)))./sum(QRiver_prist)

histogram(doug_delta_area(idx2>0)-delta_area_proxy(idx2(idx2>0)))

scatter(doug_delta_area(idx2>0),delta_area_proxy(idx2(idx2>0)))
set(gca,'XScale','log','YScale','log')


delta_area = nan(size(delta_name));
delta_width = nan(size(delta_name));
delta_length= nan(size(delta_name));

delta_area(idx2(idx2>0)) = doug_delta_area(idx2>0);
delta_width(idx2(idx2>0)) = doug_delta_width(idx2>0);
delta_length(idx2(idx2>0)) = doug_delta_length(idx2>0);

kmlwrite('jaap_matching',MouthLat,MouthLon,'name',string(BasinID2),'description',string(QRiver_prist))


blubs = idx2(idx2>0);
blub2 = find(idx2>0);
s = [];
for ii=1:length(blubs),
    slat{ii} = [MouthLat(blubs(ii)),data(blub2(ii),1)];
    slon{ii} = [MouthLon(blubs(ii)),data(blub2(ii),2)];
end
s = geoshape(slat,slon);    
s.Name = string(doug_id(idx2>0));
s.Geometry = 'line';
kmlwrite('jaap_matching2',s)
%}





