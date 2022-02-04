function get_deltaprofile
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','channel_len','shelf_len','shelf_lines','BasinID2','Hs','depth_mouth');

load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelHolocene','BasinID2','SLt');

SL_LGM = SLt(:,end);

%get delta profile
beta = nan(size(channel_len,1),1);
alpha = beta;
s = beta;
r = beta;
r_h = beta;

fo_1 = fitoptions('poly2','Lower',[0,0,0],'Upper',[inf,inf,0]);
fo_2 = fittype('poly2');

for ii=1:length(s)
    if mod(ii,100)==1, ii, end

    [beta(ii),alpha(ii),s(ii),r(ii), r_h(ii)] = profile_func(channel_len(ii,:),shelf_len(ii,:),shelf_lines,fo_1,fo_2,SL_LGM(ii));
end

%set sedimentary depth of delta wedge to be at least the channel depth
bed_h = double(-max(min(50,depth_mouth),10.*Hs));

%three factors: beta, length, area

%get data from edmonds et al., 2020
[ed_ID2,~,~,ed_length] = get_edmonds_data(BasinID2);

%beta, check against blum&tornqvist
bt_ID2 = [4346011,4165211,2098785,4301321,1248635,4267691];
[~,bt_xx] = ismember(bt_ID2,BasinID2);
bt_beta = [50,20,130,80,25,15].*1e-2./1e3;
%scatter(bt_beta,beta(bt_xx)),axis equal,ylim([0 inf]),xlim([0 inf])

%length, check against blum&tornqvist and edmonds
bt_length = [40,90,100,90,150,350].*1e3;

[~,ed_xx] = ismember(ed_ID2,BasinID2);
%scatter(ed_length,s(ed_xx)),axis equal,ylim([0 inf]),xlim([0 inf]),
%set(gca,'xscale','log'),set(gca,'yscale','log')

%then use these values also
s(ed_xx) = ed_length;
beta(bt_xx) = bt_beta;
%have calibrated values and they seem more or less fine. 

save('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile','r','s','beta','alpha','bed_h','r_h');

%also save netcdf
out = struct('BasinID2',BasinID2,'alpha', alpha,'beta', beta,'r', r ,'s', s,'bed_h', bed_h,'r_h',r_h);

funits = {'','m','m','m','m','m','m'};
fmeta = {'BasinID2','delta surface slope at river mouth', 'Delta basement slope','Distance from delta center to alluvial-basement transition',...
    'Distance from delta center to the shoreline','basement depth at shoreline','elevation at alluvial-basement transition'};
create_netcdf('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.nc',out,funits,fmeta)



end


function [beta,alpha,s,r,r_h] = profile_func(l_x,shelf_x,shelf_h,fo_1,fo_2,LGM)

if any(isnan(shelf_x)),
    shelf_x(isnan(shelf_x)) = interp1(shelf_h(~isnan(shelf_x)),shelf_x(~isnan(shelf_x)),shelf_h(isnan(shelf_x)),'linear','extrap');
end

[~,sh1] = min(abs(shelf_h(2:end)-LGM));
shelf_h = shelf_h(sh1+1);
shelf_x = shelf_x(sh1+1);


%if all(isnan(shelf_x)),
%    shelf_x = 120/1e-3;
%    shelf_h = -120;
%if any(isnan(shelf_x)),
%    sh1 = find(~isnan(shelf_x(1:find(shelf_h>=-120,1,'last'))),1,'last');
%    shelf_x = shelf_x(sh1);
%    shelf_h = shelf_h(sh1);
%end
%}
%get rid of first zeros (mouth height)
x_nan = find(isnan(l_x),1); if isempty(x_nan), x_nan=length(l_x); end
shore = find(l_x>0,1)-1; if shore==0, l_x = [0 l_x]; shore=1; end

if max(l_x)==0 | x_nan<4 | (x_nan-shore)<5 | LGM>0,
    beta = nan; alpha = nan; s = nan; r = nan; r_h = nan;
    return,
end

l_x = l_x(shore(1):x_nan-1)';
l_h = (0:length(l_x)-1)';



sl = nan(2,length(l_h)-1);
for ii=2:length(l_h)-1,
    sl(:,ii) = [ones(3,1), [fliplr(shelf_x)'; -l_x(ii+(0:1))]]\[fliplr(shelf_h)'; l_h(ii+(0:1))];
    
    

end

[bed_h,ii] = min(sl(1,:));

if bed_h>=0, %no sedimentary wedge
    beta = -sl(2,ii);
    alpha = beta;
    s=nan;
    r=nan;
    bed_h=nan; 
    r_h = 0;
else, %sedimentary wedge; get river delta values
    
    %fit 2nd degree polynomial
    l_delta = coeffvalues(fit(l_x(1:(ii+1)),l_h(1:(ii+1)),fo_2,fo_1));
    beta = -sl(2,ii);   
    alpha = l_delta(2);
    r_h = l_h(ii+1);
    s = max(1,min(l_x(ii),-bed_h/beta));
    r = min(0,-l_x(ii)+s);
    if r==0 && x_nan==30, r=-s; end %didn't find end of delta, assume big delta (yellow fix)
end

%{
    figure
    hold on
    grid on
    plot(-l_x,l_h,'ok')
    plot(shelf_x,shelf_h,'ob')
    plot([-s 0],[0 0]);
    plot([-s+r -s],[0 0]);
    plot(0:-1:-l_x(ii),polyval(l_delta,0:l_x(ii)))
    xlabel('Along-stream (m)')
    ylabel('Elevation (m)')
%}
end


