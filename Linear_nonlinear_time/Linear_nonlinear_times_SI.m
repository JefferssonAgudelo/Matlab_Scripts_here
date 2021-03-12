%------------------------------------------------------------------
% This is how to define the linear and non-linear times
%------------------------------------------------------------------
% The linear time is tl = 1/w where omega is the frequency of the
% waves that are solutions to the disspersion relation. There is a
% different definition for each type of wave. Alfven or fast waves.

% Plasma Constants in SI to make it consistent with the simulation units
%-------------------------------------------------------------------
mu0=1;
eps0=1;
kb=1;
c=1;
mi=1;
mime=100;
me= mi/mime;
qi=1;
qe=-qi;
di=1;
 
% Plasma parameters
%-------------------------------------------------------------------
Va= 0.1;
B0 = Va;

%Ions 
ni=1;
Ti=1;
Ti_per=1;
Ti_par_per = 1; 
Ti_par = Ti_per / Ti_par_per; 
Vth_i_per = sqrt(2 * Ti_per / mi);
Vth_i_par = sqrt(2 * Ti_par / mi);
beta_i= (2*mu0*ni*kb*Ti) / (B0*B0);
beta_i_per= (2*mu0*ni*kb*Ti_per) / (B0*B0);
beta_i_par= (2*mu0*ni*kb*Ti_par) / (B0*B0);
Omega_i = (qi * B0) / (mi);
omega_pi = sqrt( (ni*qi*qi) /( mi*eps0 ));
lambda_Di = sqrt( (eps0 * kb * Ti) / (ni*qi*qi) ) ;
rho_i = Vth_i_per / Omega_i ;

%Va = B0/sqrt(4*pi*ni*mi);

% Electrons
ne=1;
Te=1;
Te_per=1;
Te_par_per = 1; 
Te_par = Te_per / Te_par_per; 
Vth_e_per = sqrt(2 * Te_per / me);
Vth_e_par = sqrt(2 * Te_par / me);
neni = ne/ni;
memi = 1/mime;
beta_e= (2*mu0*ne*kb*Te) / (B0*B0);
beta_e_per= (2*mu0*ne*kb*Te_per) / (B0*B0);
beta_e_par= (2*mu0*ne*kb*Te_par) / (B0*B0);
Omega_e = (qe * B0) / (me) ;
omega_pe = sqrt( (ne*qe*qe) /( me*eps0 ));
lambda_De = sqrt( (eps0*kb * Te) / (ne*qe*qe) ) ;
rho_e = Vth_e_per / Omega_e ;

d_i=rho_i ./ sqrt(beta_i);

%Ion acustic scaling
Vs = sqrt(Te / mi); 
Rs = Vs / Omega_i;
VsVa = Vs./Va; 
%-------------------------------------------------------------------


%These are the values of the normalization parameters in the simulation
%-------------------------------------------------------------------
d_i=1; ni=1; ne=1; Va= 0.1; B0 = Va;
beta_i=1;
beta_e=1;
rho_i = sqrt(beta_i) * d_i;
Ti_0 = (beta_i * B0*B0)/(2*ni);
Te_0 = (beta_e * B0*B0)/(2*ne);
Ti=Ti_0;
Te=Te_0;

%-------------------------------------------------------------------

%Here is where I am going to include the input for k_per and k_par
%-------------------------------------------------------------------
load ('/Volumes/PSC_DiRAC_DATA/nonlineartime_images/Vi_cell_fourier_spectrum.mat')
spectravi=spectravis{10};
P2D=spectravi.P2D;
k_per=spectravi.kper;
k_par=spectravi.kpar;
%k_per = 1;


%k_par = 1;

k = sqrt( k_per.*k_per .* ones(size(P2D)) + k_par'.*k_par' .* ones(size(P2D)));
%-------------------------------------------------------------------
%-------------------------------------------------------------------
% For Alfven waves in the inertial range. (eq. 163)

w_aw = k_par' .* Va ;

v_aw_phase = (w_aw.* ones(size(P2D))) ./ k_par' ;

tl_aw = 1 ./ (k_par' .* v_aw_phase);   
%-------------------------------------------------------------------

%-------------------------------------------------------------------
% For Kinetic Alfven waves (eq. 10.179) Baumjohann Treumann (There is a problem with the dimentions)
TeTi = Te./Ti;

%w_kaw_p = +( k_par' .* Va  ).* sqrt( 1 +  ((k_per ).^2 .* rho_i.^2) .*( (3/4) + ( TeTi ) ));
%w_kaw_m = -( k_par' .* Va  ).* sqrt( 1 +  ((k_per ).^2 .* rho_i.^2) .*( (3/4) + ( TeTi ) ));

w_kaw_p = +( k_par' .* Va .* ones(size(P2D)) ).* sqrt( 1 +  ((k_per .* ones(size(P2D)) ).^2 .* rho_i.^2) .*( (3/4) + ( TeTi ) ));
w_kaw_m = -( k_par' .* Va .* ones(size(P2D)) ).* sqrt( 1 +  ((k_per .* ones(size(P2D)) ).^2 .* rho_i.^2) .*( (3/4) + ( TeTi ) ));

v_kaw_phase_p = (w_kaw_p) ./ k_par';
v_kaw_phase_m = (w_kaw_m) ./ k_par';

tl_kaw_p = 1 ./ (k_par' .* v_kaw_phase_p);
tl_kaw_m = 1 ./ (k_par' .* v_kaw_phase_m);


%-------------------------------------------------------------------

%-------------------------------------------------------------------
% For Alfve/Ion-Cyclotron waves (eq. 168)

%w_IC_p = (k_par*d_i).^2 * Omega_i * (sqrt( 1 + (2/(k_par * d_i)).^2 ) - 1); 
%w_IC_m = (k_par*d_i).^2 * Omega_i * (-sqrt( 1 + (2/(k_par * d_i)).^2 ) - 1); 

w_IC_p = (k_par' .* ones(size(P2D)) .* d_i).^2 .* Omega_i .* (sqrt( 1 + (2 ./ (k_par' .*ones(size(P2D)) .* d_i)).^2 ) - 1); 
w_IC_m = (k_par' .* ones(size(P2D)) .* d_i).^2 .* Omega_i .* (-sqrt( 1 + (2 ./ (k_par' .*ones(size(P2D)) .* d_i)).^2 ) - 1); 


v_IC_phase_p = w_IC_p ./ k_par';
v_IC_phase_m = w_IC_m ./ k_par';

tl_IC_p = 1 ./ (k .* v_IC_phase_p);
tl_IC_m = 1 ./ (k .* v_IC_phase_m);
%-------------------------------------------------------------------


%-------------------------------------------------------------------
% For Ion-acustic waves (eq. 165)
%-------------------------------------------------------------------

% Now calculating the non linear time
%-------------------------------------------------------------------

u_k_per_1=(sqrt((P2D) .* k_per ) .* k_per) .* ones(size(P2D));
u_k_per_2=(sqrt((P2D./rms(P2D)) .* k_per ) .* k_per) .* ones(size(P2D));
tnl_v_1 = 1 ./u_k_per_1; 
tnl_v_2 = 1 ./u_k_per_2;




% Making plots
%--------------------------------------------------------------------

    cd /Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here
    f1=figure(1);
    ax(1)=subplot(3,3,1);
    [C,h]=contourf(k_par,k_per,log10(tl_aw'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(ax(1),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{l,aw}$','Interpreter','latex')
    hold on
    x=log10([0.1 50]);
    y=(3/2)*x;
    %loglog(10.^x,10.^y,'--k')
    ylim([0.28 52])
    txt = '\bf{3/2 \rightarrow}';
    %text(0.1,0.5,txt)
    ylim([0.28 52])
    hold off
    
    ax(2)=subplot(3,3,2);
    [C,h]=contourf(k_par,k_per,log10(tl_kaw_p'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(ax(2),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{l,kaw,p}$','Interpreter','latex')
    hold on
    x3=log10([0.1 70]);
    y3=(3)*x3;
    %loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt2 = '\bf{\leftarrow 3}';
    %text(5,28,txt2)
    ylim([0.28 52])
    hold off
    
    ax(3)=subplot(3,3,3);
    [C,h]=contourf(k_par,k_per,log10(abs(tl_kaw_m)'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(ax(3),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}|\tau_{l,kaw,m}|$','Interpreter','latex')
    hold on
    x3=log10([0.1 70]);
    y3=(3)*x3;
    %loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt2 = '\bf{\leftarrow 3}';
    %text(5,28,txt2)
    ylim([0.28 52])
    hold off
    
     ax(4)=subplot(3,3,4);
    [C,h]=contourf(k_par,k_per,log10(tl_IC_p'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(ax(4),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{l,IC,p}$','Interpreter','latex')
    hold on
    x3=log10([0.1 70]);
    y3=(3)*x3;
    %loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt2 = '\bf{\leftarrow 3}';
    %text(5,28,txt2)
    ylim([0.28 52])
    hold off
    
    ax(5)=subplot(3,3,5);
    [C,h]=contourf(k_par,k_per,log10(abs(tl_IC_m')),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(ax(5),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}|\tau_{l,IC,m}|$','Interpreter','latex')
    hold on
    x3=log10([0.1 70]);
    y3=(3)*x3;
    %loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt2 = '\bf{\leftarrow 3}';
    %text(5,28,txt2)
    ylim([0.28 52])
    hold off
    
    ax(6)=subplot(3,3,6);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(tnl_v_1'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(ax(6),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{nl_1}$','Interpreter','latex')
    hold on
    %x3=log10([0.1 70]);
    y3=(3)*x3;
    loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    %txt2 = '\bf{\leftarrow 3}';
    text(5,28,txt2)
    ylim([0.28 52])
    hold off

    ax(7)=subplot(3,3,7);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(tnl_v_2'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(ax(7),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{nl_2} \ P2D/rms(P2D)$','Interpreter','latex')
    hold on
    x3=log10([0.1 70]);
    y3=(3)*x3;
    %loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt2 = '\bf{\leftarrow 3}';
    %text(5,28,txt2)
    ylim([0.28 52])
    hold off
