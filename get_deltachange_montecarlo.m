function land_area = get_deltachange_montecarlo(QRiver_dist,Discharge_prist,DeltaSLR,DeltaSub,SLR_unc,s,r,beta,alpha,psi,bed_h,rab,w)

land_area = zeros(length(beta),5000,'single');
for ii=1:5000,
    DeltaSLR_r = DeltaSLR+(SLR_unc.*randn(numel(w),1)); %church et al 2004
    DeltaSub_r = DeltaSub.*(0.5+rand(numel(w),1));
    w_r = w.*(0.5+rand(numel(w),1));
    bed_r = bed_h.*(0.5+rand(numel(w),1));
    QRiver_dist_r = QRiver_dist.*(0.62+(2*0.38.*rand(numel(w),1)));  %+/- 38% from cohen et al 2013
    nu = 0.1 * 365*24*3600*Discharge_prist./w_r; %paola et al, %m2/yr %get nu from field data, way to undercontrained
    qs = rab.*nu.*beta;
    qs_sus = 365*24*3600*QRiver_dist_r./1600./w_r; %suspended flux per m coastline
    [ds,dr] = get_deltachange(qs,qs_sus,DeltaSub_r+DeltaSLR_r,s,r,beta,alpha,psi,bed_r,1);
    land_area(:,ii) = w_r.*ds;
end
end
