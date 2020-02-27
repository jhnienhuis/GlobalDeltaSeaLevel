function [rab1,rab2,beta,alpha,s,r,bed_h] = get_deltaprofile(x,shelf_len,fo_1,fo_2,pl)

%get rid of first zeros (mouth height)
l_nan = find(isnan(x),1);
shore = find(x>0,1)-1;

%estimate r, s and rab based on profile
if max(x)==0 | l_nan<4 | (l_nan-shore)<5,
    rab1 = nan; rab2 = nan; beta = nan; alpha = nan; s = nan; r = nan; bed_h = nan;
    return,
end

x = x(shore(1):l_nan-1)';
l_h = (0:length(x)-1)';
%fit 2nd degree polynomial
l_fit = coeffvalues(fit(x,l_h,fo_2,fo_1));
%get basement slope
beta = l_fit(1)*2*x(end)+l_fit(2);
%get depth of basement transition
bed_h = -(beta*x(end))+polyval(l_fit,x(end));

%fit another polynomial ont he first 5 values, if poor fit use polynomial
%of all values
[l_fit2,gof] = fit(x(1:5),l_h(1:5),fo_2,fo_1);
if gof.adjrsquare<0.25 | (l_fit2.p1*x(end)^2+l_fit2.p2*x(end))>50, l_fit2 = l_fit;
else, l_fit2 = coeffvalues(l_fit2);
end


l_fit2(2) = min(l_fit2(2),beta);
l_fit_y= polyval(l_fit2,x);
s = -bed_h/beta;
r = min(0,s-(beta-l_fit2(2))/2/l_fit2(1)); %maximum of same slope R, and intersect R
in = -x(find(l_fit_y<(beta*x+bed_h),1))+s;
if in>r, r = in; end


rab1 = (l_fit2(1)*(s-r)+0.5*l_fit2(2))*2/beta; %always unity

if r>-1 || -r>s , rab1 = 0; r = max(r,-s); end

% offshore
%
%sea_fit = coeffvalues(fit(fliplr(-shelf_len(idx,1:3)*1000)',sea_h(4:6)','poly1',fo_2));
%sea_fit = polyfit(fliplr(-shelf_len(idx,1:3)*1000),sea_h(4:6),1);

%get rab based on delta gross shape
l_toe = bed_h/(beta-(20/shelf_len(2)));
if l_toe<0, l_toe = max(s+shelf_len(2),(bed_h+20)/beta); end %intersect w/ bathy

fr = (-r./(s-r));
rab2 = fr./(fr+(1-fr)*(l_toe./(l_toe-r)));

alpha = (l_fit_y(2)-l_fit_y(1))./x(2);

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
plot(-x,l_h,'ok')
plot(-x,polyval(l_fit,x))
plot(-x,l_fit_y,'g')
plot(-x,beta*x+bed_h)
plot(fliplr(shelf_len(1:6)),-100:20:0,'-ob')
xlabel('Along-stream (m)')
ylabel('Elevation (m)')
end
    %}