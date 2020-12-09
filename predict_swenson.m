%model variables
tau_star = 1; %shield stress
R = 0.6; %density of water vs sediment
Cf = 1e-2; %drag coefficient fluvial
Cfs = 1e-2; %drag coefficient marine

beta = 1e-2; %fraction of basin covered by channel (do channel width / delta width)

q_wf = 1.*3600*24*365; %m2/yr (water discharge from WBMsed)
If = 1e-2; %frequency (get from WBMsed) (or, get discharge at this frequency)

Is = 1e-5; %storm frequency (get from wavewatch)
H_b = 6; %wave breaking height (m)

D = 50;


h_b = H_b/gamma_b;
e_ss = 1e-2;
gamma_b = 0.6;
W_s = 1e-2.*3600*24*365; %m/yr %keep in seconds very fine sand..
g = 10*3600*24*365.*3600*24*365; %m/yr^2 keep in seconds

nu_0 = 0.3.*3600*24*365;
q_wo = 1e-1;
upsilon = 1e6;
q_so = 500; %sediment supply at m2/yr
upsilon = 1/20.*R.*If.*q_wf.*tau_star./sqrt(Cf)

L = 1e5;
L2 = D*upsilon./q_so
tau_b = 1e4;
tau_b2 = (D*upsilon./q_so).^2./(upsilon)

%add tidal currents?
%v0 is a tidal current... seems appropriate

%dimensionless parameters



hb_star = h_b./D
k_star = (2./3./pi.*Cfs*e_ss.*(gamma_b.^3).*sqrt(g.*h_b)./W_s).*(gamma_b.^2.*sqrt(Cf)./tau_star).*(Is.*g.*h_b.^2./(If.*q_wf.*W_s))

nu_star = (2./3./pi.*Cfs*e_ss.*(gamma_b.^3).*sqrt(g.*h_b)./W_s).*R.*(Is.*nu_0.*h_b./q_so)

%% get wave breaking/storm frequency
gamma_b = 0.6;
g = 9.81;

hb = @(tp,hs) ((hs./gamma_b.*(sqrt(g./4./pi.*tp))./(g.^(1/4))).^(4/5));
v0 = out.Discharge_tide./(out.width_mouth.*out.depth_mouth);
q_wf = out.Discharge_prist./out.width_upstream;

hb_star = hb(out.Tp,out.Hs)./(psi.*10000);