% This is to plot and analyse the NHDS output data.


%Change the directory
%-----------------------------------------------------------------------
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/NHDS/NHDS-master';
path = '/Users/jeffersson_agudelo/Documents/CB104_local_data/NHDS/NHDS-master';

% Removing the first line of the file
%--------------------------------------------------------------------------
%# kk, theta, omega, gamma, Re(Ey/Ex), Im(Ey/Ex), Re(Ez/Ex), Im(Ez/Ex), energy, quality
%--------------------------------------------------------------------------
%wave_file = load('output_kaw.in_wave.dat');

filename='output_kaw.in_wave.dat';
nn=4000*77 + 1;
wave_file = readmatrix(filename,'Range',[1 1 nn 4]);



%cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

k_nhds = wave_file(:,1);
theta_nhds = wave_file(:,2);
w_nhds = wave_file(:,3);
gamma_nhds = wave_file(:,4);

theta_rad_nhds = wave_file(:,2) * (pi/180);
k_par_nhds = k_nhds.*cos(theta_rad_nhds);
k_per_nhds = k_nhds.*sin(theta_rad_nhds);

f1=figure(1);
pointsize = 5;
scatter(k_nhds,w_nhds,pointsize,theta_nhds,'f')
hold on
set(gca,'XScale','log','YScale','log','FontSize',18)
cb = colormap(jet);
%caxis([0.0001 2]);
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '\theta \circ';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
xlabel('$k / d_i$','Interpreter','latex')
ylabel('$\omega / \Omega_{i}$','Interpreter','latex')
title('$\omega$ vs $k$','Interpreter','latex')
hold off

f2=figure(2);
pointsize = 5;
scatter(k_nhds,abs(gamma_nhds),pointsize,theta_nhds,'f')
hold on
set(gca,'XScale','log','YScale','log','FontSize',18)
cb = colormap(jet);
%caxis([0.0001 2]);
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '\theta \circ';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
xlabel('$k / d_i$','Interpreter','latex')
ylabel('$| \gamma | / \Omega_{i}$','Interpreter','latex')
title('$| \gamma |$ vs $k$','Interpreter','latex')
hold off

saveas(f1,'omega_k_all.png');
saveas(f2,'gamma_k_all.png');
%------------------------------------------------------------------------  



%------------------------------------------------------------------------
[kpar2D_s, kper2D_s]=meshgrid( kpar,kper) ;
k2D_s=sqrt(kpar2D_s.*kpar2D_s + kper2D_s.*kper2D_s);   
theta_s = atan(kper2D_s./kpar2D_s) .* (180/pi);


f3=figure(3);
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'normal';

nexttile(tl)
[C,h]=contourf(kpar2D_s,kper2D_s,k2D_s);
set(gca,'XScale','log','YScale','log','FontSize',18)
cb = colormap(jet);
%caxis([1 27]);
hcb=colorbar;
xlabel('$k / d_i$','Interpreter','latex')
ylabel('$k_\| d_i$','Interpreter','latex')
title('$|k|$','Interpreter','latex')

nexttile(tl)
[C,h]=contourf(kpar2D_s,kper2D_s,theta_s);
set(gca,'XScale','log','YScale','log','FontSize',18)
cb = colormap(jet);
%caxis([1 27]);
hcb=colorbar;
xlabel('$k / d_i$','Interpreter','latex')
ylabel('$k_\| d_i$','Interpreter','latex')
title('$\theta \circ$','Interpreter','latex')

nexttile(tl)
[C,h]=contourf(kpar2D_s,kper2D_s,kper2D_s./kpar2D_s);
set(gca,'XScale','log','YScale','log','FontSize',18)
cb = colormap(jet);
%caxis([1 27]);
hcb=colorbar;
xlabel('$k / d_i$','Interpreter','latex')
ylabel('$k_\| d_i$','Interpreter','latex')
title('$k_\perp / k_\|$','Interpreter','latex')


% Left to right, top to bottom
n_stp=100;
DeltaN=200;
f6=figure(6);
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
j1=1;
j2=n_stp;
for i = 1:4
    nexttile(tl)
    [kpar2D, kper2D]=meshgrid( k_par_nhds(j1:j2),k_per_nhds(j1:j2)) ;
    contourf(kpar2D,kper2D,kper2D./kpar2D);
    %imagesc(magic(2*i))
    set(gca,'XScale','log','YScale','log','FontSize',10)
    j1 = DeltaN*i;
    j2 = n_stp + DeltaN*i;
end
title(tl,'$k_\perp / k_\| $','Interpreter','latex','FontSize',18);
xlabel(tl,'$k_\| d_i$','Interpreter','latex','FontSize',18);
ylabel(tl,'$k_\perp d_i$','Interpreter','latex','FontSize',18);
%------------------------------------------------------------------------



%Things to compare and understand
%------------------------------------------------------------------------
theta=linspace(0,27,90); %degrees
theta_rad=pi*theta/180;
cos_theta=cos(theta_rad);
sin_theta=sin(theta_rad);

f1=figure(1);
plot(theta,1.9*cos_theta,'b');
hold on
plot(theta,2.*cos_theta,'r');
hold off

f2=figure(2);
plot(theta,1.9*sin_theta,'b');
hold on
plot(theta,2.*sin_theta,'r');
hold off

f3=figure(3)
plot(2*cos_theta, 2*sin_theta)
hold on
plot(1.9*cos_theta, 1.9*sin_theta)
%------------------------------------------------------------------------



% These are things that were useful at some point
%-------------------------------------------------------

%{
Usage:
h=subaxis(rows,cols,cellno[,settings])
h=subaxis(rows,cols,cellx,celly[,settings])
h=subaxis(rows,cols,cellx,celly,spanx,spany[,settings])

SETTINGS: Spacing,SpacingHoriz,SpacingVert
Padding,PaddingRight,PaddingLeft,PaddingTop,PaddingBottom
Margin,MarginRight,MarginLeft,MarginTop,MarginBottom
Holdaxis
%}

f5=figure(5);
title('$\omega$ vs $k$','Interpreter','latex')
for i=1:9
h1=subaxis(3,3,i,'SV',0.04,'SH',0.04,'MR',0.05,'ML',0.05);
imagesc(magic(2*i))
%set(gca,'XScale','log','YScale','log','FontSize',18)
end
%-----------------------------------------------------------
