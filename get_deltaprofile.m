function [beta,alpha,s,r,bed_h,psi,r_h] = get_deltaprofile(l_x,shelf_x,shelf_h,fo_1,fo_2,pl)

%old bathy
sh1 = 1; sh2 = find(shelf_h<=-40,1);

%get rid of first zeros (mouth height)
x_nan = find(isnan(l_x),1); if isempty(x_nan), x_nan=length(l_x); end
shore = find(l_x>0,1)-1; if shore==0, l_x = [0 l_x]; shore=1; end

%if shelf_len doesn't exist, put in 1e-3
s_nan = min(sh2,find(isnan(shelf_x),1)-1); 
if isempty(s_nan), s_nan=sh2; elseif s_nan<sh2, s_nan=sh2; shelf_x = -shelf_h./1e-3; end




%estimate r, s and rab based on profile
if max(l_x)==0 | x_nan<4 | (x_nan-shore)<5,
    beta = nan; alpha = nan; s = nan; r = nan; bed_h = nan; psi=nan; r_h=nan;
    return,
end

l_x = l_x(shore(1):x_nan-1)';
l_h = (0:length(l_x)-1)';

shelf_x = shelf_x(sh1:s_nan);
shelf_h = shelf_h(sh1:s_nan);

%l_basement = coeffvalues(fit([fliplr(shelf_x)'; -l_x(end+(-4:0))],[fliplr(shelf_h)'; l_h(end+(-4:0))],'poly1'));
sl = nan(2,length(l_h)-1);
for ii=2:length(l_h)-1,
    sl(:,ii) = [ones(length(shelf_x)+1,1), [fliplr(shelf_x(2:end))'; -l_x(ii+(0:1))]]\[fliplr(shelf_h(2:end))'; l_h(ii+(0:1))];
end

[bed_h,ii] = min(sl(1,:));

%no sedimentary wedge
if bed_h>=0,
    beta = -sl(2,ii);
    alpha = beta;
    s=1;
    r=0;
    bed_h=0;   
    psi = -(fliplr(shelf_x)'\fliplr(shelf_h)');
    r_h = 0;
else,
    %sedimentary wedge; get river delta values
    
    %fit 2nd degree polynomial
    l_delta = coeffvalues(fit(l_x(1:(ii+1)),l_h(1:(ii+1)),fo_2,fo_1));
    
    %get basement slope from mean of local fit and global fit
    beta = mean([-sl(2,ii),l_delta(1)*2*l_x(ii)+l_delta(2)]);
    
    %get thickness of alluvial deposit under river mouth from mean of
    %global and local fit:
    bed_h = mean([sl(1,ii),(l_h(ii)-(l_x(ii)*(l_delta(1)*2*l_x(ii)+l_delta(2))))]);
    alpha = l_delta(2);
    s = max(1,min(l_x(ii),-bed_h/beta));
    r = min(0,-l_x(ii)+s);
        
    psi = -(fliplr(shelf_x)'\fliplr(shelf_h)');
    %height of alluvial bedrock transition
    r_h = l_h(ii+1);
    
end

if pl == 1;
    disp(['r = ' num2str(r)])
    disp(['s = ' num2str(s)])
    disp(['bedh = ' num2str(bed_h)])
    disp(['beta = ' num2str(beta)])
    disp(['alpha = ' num2str(alpha)])
    disp(['rab1 = ' num2str(rab1)])
    disp(['rab2 = ' num2str(rab2)])
    
    figure
    hold on
    grid on
    plot(-l_x,l_h,'ok')
    plot(shelf_x,shelf_h,'ob')
    plot([fliplr(shelf_x(2:end))'; -l_x(ii+(0:1))],[fliplr(shelf_x(2:end))'; -l_x(ii+(0:1))].*sl(2,ii)+sl(1,ii));
    plot(0:-1:-l_x(ii),polyval(l_delta,0:l_x(ii)))
    xlabel('Along-stream (m)')
    ylabel('Elevation (m)')
end
%}