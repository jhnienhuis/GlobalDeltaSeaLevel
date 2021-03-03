function land_area = get_deltachange_montecarlo(QRiver_dist,DeltaSLR,DeltaSub,SLR_unc,s,r,w,bed_h)

%f_lin = @(x) (0.25+(1.5*rand(x,1)));
f_mul = @(x) (10.^(-0.5+1*rand(x,1)));

land_area = nan(length(bed_h),2000,'single');
len = numel(w);


for ii=1:2000,
    DeltaSLR_r = DeltaSLR+(SLR_unc.*randn(numel(w),1)); %church et al 2004
    DeltaSub_r = DeltaSub.*f_mul(len);
    w_r = w.*f_mul(len);
    bed_r = bed_h.*f_mul(len);
    QRiver_dist_r = QRiver_dist.*(0.62+(2*0.38.*rand(numel(w),1)));  %+/- 38% from cohen et al 2013
    s_r = s.*f_mul(len);
    r_r = r.*f_mul(len);
           
    land_area(:,ii) = get_deltachange(365*24*3600*QRiver_dist_r./1600,DeltaSub_r+DeltaSLR_r,s_r,r_r,w_r,bed_r);
        
end

%[~,idx] = get_deltachange(365*24*3600*QRiver_dist./1600,DeltaSub+DeltaSLR,s,r,w,bed_h);
%land_area(~idx,:) = nan;

end
