N = 100;
dx = 1./N;
dt = 0.001;
T = 10000;
Sbase = tan(deg2rad(12.5));
Sa = 1e-1;
Sfi = 1e-2; %initial height of foreset
slr = 0; %slr rate
qw = 4.36; %cm2/s

qsf = 0.904; %cm2/s

%initial conditions
eta_si = 0; %initial height of shoreline

xi_i = 0; %m
ss = 10; %m
su = 0; %m

eta_bi = (Sbase-Sfi).*ss./(1-(Sbase./Sa)); %initial height of foreset

sb=ss+eta_bi./Sa;

x = linspace(0,ss,N+1);

eta = Sfi.*(ss-x); %eq. 13a

plot(x,eta,'-o');

for ii=1:T,
    
    
    
    xb = (x-su)./(ss-su);
    S(1) = (1./(ss-su)).*(eta(1)-eta(2))./dx;
    S(2:N) = (1./(ss-su)).*(eta(1:(N-1))-eta(3:(N+1)))./(2*dx);
    
    S_r(2:(N+1)) = (1./(ss-su)).*(eta(1:N)-eta(2:(N+1)))./dx;
    qs(2:(N+1)) = (qw.*12.3).*S_r(2:(N+1)).^2.24;
    qs(1) = qsf;
    %qs(N+1) = %interp??
    
    dsu = -1./(Sbase-S(1)).*(qsf-qs(2))./dx.*dt./(ss-su);
    
    dss = dt./(Sa+S(N)).*((qs(N+1)./(sb-ss))+(qs(N+1)-qs(N))./dx./(ss-su));
    
    dsb = 1./(Sa-Sbase).*((dss.*(Sa-S(N))-dt.*(qs(N+1)-qs(N))./dx./(ss-su)));
    
    eta(1) = eta(1) + Sbase./(Sbase-S(1)).*(qsf-qs(2))./dx.*dt./(ss-su);
    eta(N) = xi_i+(slr.*ii);
    eta(1:N) = eta(1:N) - (xb(1:N).*dss+(1-xb(1:N)).*dsu).*S.*dt + (qs(1:N)-qs(2:(N+1)))./dx.*dt./(ss-su);
    
    su = su+dsu;
    ss = ss+dss;
    sb = sb+dsb;
    
    plot(x*(ss-su)+su,eta,'-o'), drawnow, %pause(0.001)
    
end