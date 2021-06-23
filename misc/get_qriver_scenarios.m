% get Qriver scenarios for Frances
p = shaperead([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.shp']);
p_BasinID2 = [p(:).BasinID2];
load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat'],'BasinID2','BasinArea');
%resolution is 240 per degree (15s)
res = 240;
%resolution WBMSed is 10 per degree (6 min)
res_wbm = 10;
remfun = @(lon) (rem(res_wbm*360-1+lon,res_wbm*360)+1);

for jj=1:12,
    qs_bar(jj).m = matfile(['D:\OneDrive - Universiteit Utrecht\WBMSed\Qs_bar' num2str(jj) '_Yield.mat']);
end

%a = load('D:\OneDrive - Universiteit Utrecht\WBMSed\Prist_Yield.mat');
%b = load('D:\OneDrive - Universiteit Utrecht\WBMSed\Dist_Yield.mat');
%WBMsed = permute(qs_bar(1).m.WBMsed,[3 1 2]);
%area per cell of WBMSed in km2
areapercell = repmat(6371.^2.*2*pi/360/res_wbm*(sin(deg2rad(-(90-(1/res_wbm)):(1/res_wbm):90))-sin(deg2rad(-90:(1/res_wbm):(90-(1/res_wbm))))),360*res_wbm,1);


cellareabasins = zeros(length(BasinID2),1);
QRiver_bar1 = zeros(length(BasinID2),120);
QRiver = zeros(length(BasinID2),1);
for jj=1:12,
    WBMsed = qs_bar(jj).m.WBMsed;
for ii=1:length(BasinID2),
    
    if mod(ii,100)==1, ii, end
    idx = find(p_BasinID2==BasinID2(ii)); idx = idx(1);
    %get yield
    a = false(res_wbm*360,res_wbm*180);
    a(sub2ind(size(a),floor(remfun(res_wbm*p(idx).X(~isnan(p(idx).Y)))),floor(res_wbm*(90+p(idx).Y(~isnan(p(idx).Y)))))) = 1;
    
    if any(a(1,:)),
        abasin = circshift(imfill(circshift(a,360*res_wbm/2),8,'holes'),-(360*res_wbm/2));
    else,
        abasin = imfill(a,8,'holes');
    end
    cellareabasins(ii) = sum(areapercell(abasin(:)));
    
    %[xx,yy] = find(abasin);
    %idx = sub2ind([3600,1800,120],xx,yy,1:120);
    
    
    
        %qs_bar(jj).QRiver(ii,1:120) = sum(sum(WBMsed(:,abasin))).*out.BasinArea(ii)./cellareabasins(ii);
        qs_bar(jj).QRiver(ii,1:120) = sum(WBMsed(:,abasin),2).*BasinArea(ii)./cellareabasins(ii);
    
    
    %QRiver(ii) = sum(sum(b.WBMsed(abasin)));
end
end


save Qriver_scenarios qs_bar

%%
clr
load Qriver_scenarios
QRiver = cat(3,qs_bar(:).QRiver);
out = load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat']);
plot(QRiver(233,:,1))