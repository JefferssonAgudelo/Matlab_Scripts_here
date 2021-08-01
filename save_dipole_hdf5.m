%--------------------------------------------------------------------------
% This is to create the simple theoretical case in the .h5 format
%--------------------------------------------------------------------------
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
% Lets define the grid
%--------------------------------------------------------------------------
Lx=40; Ly=40; Lz=40;
Nx = 161; Ny = 161; Nz = 161;
nx = Nx-1; ny = Ny-1; nz = Nz-1;
dx = Lx/(nx); dy = Ly/(ny); dz = Lz/(nz);  

%Lx=8; Ly=9; Lz=10;
%nx = 80; ny = 90; nz = 100;

%xl=linspace(-Lx,Lx,nx);
%yl=linspace(-Ly,Ly,ny);
%zl=linspace(-Lz,Lz,nz);

xl=linspace(0,Lx,nx);
yl=linspace(0,Ly,ny);
zl=linspace(0,Lz,nz);


%rl=sqrt(xl.^2 + yl.^2 + zl.^2); This doesn't work
[XL, YL, ZL]=meshgrid(xl,yl,zl);

% For the streamlines
fr=1;
NN=200;  %Do this just once

ax=xl(1); bx=xl(end);
ay=yl(1); by=yl(end);
az=zl(1); bz=zl(end);
xstart = ax + (bx-ax)*rand(NN,1);
ystart = ay + (by-ay)*rand(NN,1);
zstart = az + (bz-az)*rand(NN,1);

fs=14; % Font size
sli=1; % slide 
%--------------------------------------------------------------------------


% One Dipole magnetic field cartesian
%--------------------------------------------------------------------------
xp = Lx/2; yp = Ly/2; zp = Lz/2; 
M_i=1;
M=M_i*ones(size(XL));
%--------------------------------------------------------------------------
%RL=sqrt(XL.^2 + YL.^2 + ZL.^2);
%Bx1d=3.*M.*XL.*ZL./(RL.^5);
%By1d=3.*M.*YL.*ZL./(RL.^5);
%Bz1d=M.*(3.* (ZL.^2) - RL.^2)./(RL.^5);

RL=sqrt((XL-xp).^2 + (YL-yp).^2 + (ZL-zp).^2);
Bx1d=3.*M.*(XL-xp).*(ZL-zp)./(RL.^5);
By1d=3.*M.*(YL-yp).*(ZL-zp)./(RL.^5);
Bz1d=M.*(3.* ((ZL-zp).^2) - RL.^2)./(RL.^5);
%--------------------------------------------------------------------------
XYZ1=stream3(XL,YL,ZL,Bx1d,By1d,Bz1d,xstart,ystart,zstart);
f31=figure(31);
hline=streamline(XYZ1);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('Dipole field lines','Interpreter','latex','FontSize',fs)

%--------------------------------------------------------------------------
f21=figure(21);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bx1d(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(By1d(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bx1d(sli,:,:).^2 + By1d(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bx1d(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(By1d(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bx1d(:,sli,:).^2 + By1d(:,sli,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
sgtitle('One Dipole, $B_{x} = 3Mxz/r^{5}; B_{y} = 3Myz/r^{5}; B_{z} = M(3z^{2} -r^(2))/r^{5}; M =$' + string(M_i),'Interpreter','latex','FontSize',fs)


%divB = divergence(XL,YL,ZL,Bx1d,By1d,Bz1d); 
%f341=figure(341);
%histogram(divB)

%ax=xl(1)+7.5; bx=xl(end)-7.5;
%ay=yl(1)+7.5; by=yl(end)-7.5;
%az=zl(1)+7.5; bz=zl(end)-7.5;
%xstart = ax + (bx-ax)*rand(NN,1);
%ystart = ay + (by-ay)*rand(NN,1);
%zstart = az + (bz-az)*rand(NN,1);

%dx=xl(2)-xl(1); dy=yl(2)-yl(1); dz=zl(2)-zl(1); 
%[sx,sy,sz] = meshgrid(xl(1):8*dx:xl(end),yl(1):8*dy:yl(end),zl(83));
%XYZ2=stream3(XL,YL,ZL,Bx,By,Bz,sx,sy,sz);
%--------------------------------------------------------------------------


% This is for the Harris current sheets
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%Role of electron physics in the development of turbulent magnetic reconnection 
% in collisionless plasmas. Daughton, 2011
%--------------------------------------------------------------------------
B0=1*ones(size(XL));
lambda_i=4;
lambda=lambda_i*ones(size(XL));
%--------------------------------------------------------------------------
xp = Lx/2; yp = Ly/2; zp = Lz/2; 

Bx = B0 .* tanh((ZL-zp)./lambda); 
By = B0;
Bz = zeros(size(XL));


%--------------------------------------------------------------------------
XYZ2=stream3(XL,YL,ZL,Bx,By,Bz,xstart,ystart,zstart);
f32=figure(32);
hline=streamline(XYZ2);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('One Harris field lines','Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
f22=figure(22);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bx(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(By(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bx(sli,:,:).^2 + By(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bx(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(By(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bx(:,sli,:).^2 + By(:,sli,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
sgtitle('Single Harris, $B_{x} = B_{0}\tanh(z/\lambda); B_{y}=B_{0}; B_{z}=0$ ; $\lambda =$'...
    + string(lambda_i) + '$d_{i}$' ,'Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------

% Force free Harris current sheet
% Energy release and transfer in guide field reconnection. Birn & Hesse,
% 2010
%--------------------------------------------------------------------------
Bxf = B0 .* tanh((ZL-zp)./lambda); 
Byf = B0 ./ cosh((ZL-zp)./lambda);
Bzf = zeros(size(XL));
%--------------------------------------------------------------------------
XYZ3=stream3(XL,YL,ZL,Bxf,Byf,Bzf,xstart,ystart,zstart);
f33=figure(33);
hline=streamline(XYZ3);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('One Harris force free field lines','Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------
f23=figure(23);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bxf(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(Byf(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bxf(sli,:,:).^2 + Byf(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bxf(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(Byf(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bxf(:,sli,:).^2 + Byf(:,sli,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
sgtitle('Single Harris force free, $B_{x} = B_{0}\tanh \left( z/\lambda \right); B_{y}=B_{0} / \cosh(z/\lambda); B_{z}=0$ ; $\lambda =$'...
    + string(lambda_i) + '$d_{i}$' ,'Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------


% Multiple harris current sheets force free
% Energy dissipation and entropy in collisionless plasma. Senbei Du, 2020
%--------------------------------------------------------------------------
d_ac_i= 6;
d_ac = d_ac_i*ones(size(XL));
Bg = 1.*B0;

Bxfm = B0 .* tanh( (d_ac ./(pi.*lambda)) .* sin( (pi.*(ZL-zp))./d_ac ) ); 
Byfm = B0 .* sqrt( 1 + (Bg./B0).^2 - ( tanh( (d_ac ./(pi.*lambda)) .* sin( (pi.*(ZL-zp))./d_ac ) ).^2 ) );
Bzfm = zeros(size(XL));
%--------------------------------------------------------------------------
XYZ4=stream3(XL,YL,ZL,Bxfm,Byfm,Bzfm,xstart,ystart,zstart);
f34=figure(34);
hline=streamline(XYZ4);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('Multiple Harris force free field lines','Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
f24=figure(24);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bxfm(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(Byfm(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bxfm(sli,:,:).^2 + Byfm(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
h=colorbar;
yt=get(h,'XTick');
set(h,'XTickLabel',sprintf('%2.2f',yt));
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bxfm(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(Byfm(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bxfm(:,sli,:).^2 + Byfm(:,sli,:).^2) )');
%hc=pcolor(yl,zl,squeeze(sqrt(Bxfm(:,sli,:).^2 ))');
set(hc,'edgecolor','none')
colormap(BWR);
cb=colorbar;
%set(cb, 'Ticks', [-1, 0 1], 'TickLabels', {'-1', '0', '1'})
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
sgtitle('Multple Harris force free, $$B_{x} = B_{0}\tanh \left( \frac{d}{\pi \lambda} \sin\frac{\pi z}{d} \right); B_{y}=B_{0}\sqrt{1 + \left( \frac{B_{g}}{B_{0}} \right)^2 - \tanh^{2}\left( \frac{d}{\pi \lambda} \sin\frac{\pi z}{d} \right) }; B_{z}=0 ; \lambda =$$'...
    + string(lambda_i) + '$$d_{i}; d = $$' + string(d_ac_i) + '$$ d_{i}; B_{g} = B_{0}$$' ,'Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------
%}


% This is for the flux ropes
% THREE-DIMENSIONAL SIMULATION STUDY OF FLUX ROPE DYNAMICS IN THE SOLAR CORONA. S. Inoue, 2006
%--------------------------------------------------------------------------
% Elongated dipole
%--------------------------------------------------------------------------
M=0.5*ones(size(XL));
L0_i = (Ly/2);
L0 = L0_i.*ones(size(XL));
d0_i =(Lz/10);
d0 =d0_i.*ones(size(XL));

R_d = sqrt((YL-L0).^2 + (ZL + d0).^2);
Bx_d = zeros(size(XL)); 
By_d = (M./(2*pi)) .* ( ( (YL - L0).^2 - (ZL + d0).^2 )./ (R_d).^4 );
Bz_d = (M./(2*pi)) .*  ( (2.*(YL - L0).*(ZL + d0))./ (R_d).^4 );

%xp = Lx/2; yp = Ly/2; zp = Lz/2; 
%R_d = sqrt((YL-L0).^2 + (ZL + d0).^2);
%Bx_d = zeros(size(XL)); 
%By_d = (M./(2*pi)) .* ( ( (YL - L0).^2 - (ZL + d0).^2 )./ (R_d).^4 );
%Bz_d = (M./(2*pi)) .*  ( (2.*(YL - L0).*(ZL + d0))./ (R_d).^4 );

%--------------------------------------------------------------------------
XYZ5=stream3(XL,YL,ZL,Bx_d,By_d,Bz_d,xstart,ystart,zstart);
f35=figure(35);
hline=streamline(XYZ5);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('Elongated dipole field lines','Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------
f25=figure(25);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bz_d(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{z}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(By_d(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bz_d(sli,:,:).^2 + By_d(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
h=colorbar;
yt=get(h,'XTick');
set(h,'XTickLabel',sprintf('%2.2f',yt));
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{yz}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bz_d(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{z}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(By_d(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bz_d(:,sli,:).^2 + By_d(:,sli,:).^2) )');
%hc=pcolor(yl,zl,squeeze(sqrt(Bxfm(:,sli,:).^2 ))');
set(hc,'edgecolor','none')
colormap(BWR);
cb=colorbar;
%set(cb, 'Ticks', [-1, 0 1], 'TickLabels', {'-1', '0', '1'})
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{yz}$','Interpreter','latex','FontSize',fs)
sgtitle('Elongated dipole, $$B_{x} = 0; B_{y}=\frac{M}{2\pi} \frac{(y-L_{0})^{2} - (z+d_{0})^{2} }{r_{d}^{4}}; B_{z}=\frac{M}{2\pi} \frac{(y-L_{0}) (z+d_{0}) }{r_{d}^{4}}; r_{d} =\sqrt{(y-L_{0})^{2} + (z+d_{0})^{2}}$$' ,'Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------

% Two elongated dipoles
%--------------------------------------------------------------------------
M=0.5*ones(size(XL));
L0_i1 = (Ly/10); L01 = L0_i1.*ones(size(XL));
d0_i1 =(-4*Lz/10); d01 =d0_i1.*ones(size(XL));
R_d1 = sqrt((YL-L01).^2 + (ZL + d01).^2);
Bx_d1 = zeros(size(XL)); 
By_d1 = (M./(2*pi)) .* ( ( (YL - L01).^2 - (ZL + d01).^2 )./ (R_d1).^4 );
Bz_d1 = (M./(2*pi)) .*  ( (2.*(YL - L01).*(ZL + d01))./ (R_d1).^4 );

L0_i2 = (7*Ly/10); L02 = L0_i2.*ones(size(XL));
d0_i2 = (-6*Lz/10); d02 =d0_i2.*ones(size(XL));
R_d2 = sqrt((YL-L02).^2 + (ZL + d02).^2);
Bx_d2 = zeros(size(XL)); 
By_d2 = (M./(2*pi)) .* ( ( (YL - L02).^2 - (ZL + d02).^2 )./ (R_d2).^4 );
Bz_d2 = (M./(2*pi)) .*  ( (2.*(YL - L02).*(ZL + d02))./ (R_d2).^4 );

Bx_d12 = Bx_d1 + Bx_d2;
By_d12 = By_d1 + By_d2;
Bz_d12 = Bz_d1 + Bz_d2;
%--------------------------------------------------------------------------
XYZ6=stream3(XL,YL,ZL,Bx_d12,By_d12,Bz_d12,xstart,ystart,zstart);
f36=figure(36);
hline=streamline(XYZ6);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('Two elongated dipoles field lines','Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------
f26=figure(26);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bz_d12(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{z}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(By_d12(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bz_d12(sli,:,:).^2 + By_d12(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
h=colorbar;
yt=get(h,'XTick');
set(h,'XTickLabel',sprintf('%2.2f',yt));
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{yz}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bz_d12(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{z}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(By_d12(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bz_d12(:,sli,:).^2 + By_d12(:,sli,:).^2) )');
%hc=pcolor(yl,zl,squeeze(sqrt(Bxfm(:,sli,:).^2 ))');
set(hc,'edgecolor','none')
colormap(BWR);
cb=colorbar;
%set(cb, 'Ticks', [-1, 0 1], 'TickLabels', {'-1', '0', '1'})
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{yz}$','Interpreter','latex','FontSize',fs)
sgtitle('Two elongated dipoles, $$B_{x} = 0; B_{y}=\frac{M}{2\pi} \frac{(y-L_{0})^{2} - (z+d_{0})^{2} }{r_{d}^{4}}; B_{z}=\frac{M}{2\pi} \frac{(y-L_{0}) (z+d_{0}) }{r_{d}^{4}}; r_{d} =\sqrt{(y-L_{0})^{2} + (z+d_{0})^{2}}$$' ,'Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Real One flux rope
%--------------------------------------------------------------------------
M=0.5*ones(size(XL));
L0_i = (Ly/2); L0 = L0_i.*ones(size(XL));
d0_i =(Lz/4); d0 =d0_i.*ones(size(XL));

x1 = 2.5045;
r0_i=4;%0.3;
r0=r0_i.*ones(size(XL));
Bc = 1.0;

I0 = 2.*pi.*r0.*Bc.*besselj(1,x1);

h0_i =(Lz/4); h0 =h0_i.*ones(size(XL)); %???
M2 = M./(I0.*d0);
h0p = d0.*( (M2 -1) + sqrt(M2.^2 - 2.*M2) );
h0m = d0.*( (M2 -1) - sqrt(M2.^2 - 2.*M2) );
%h0 = h0m;

R_d = sqrt((YL-L0).^2 + (ZL + d0).^2);
rc = sqrt((YL-L0).^2 + (ZL - h0).^2);
alpha = 24.045;

alpharc = alpha.*rc;
J0 = besselj(0,alpharc);
J1 = besselj(1,alpharc);

% Define the conditions
%--------------------------------------------------------------------------
cnd1 = ((0 < rc) & (rc < r0_i)); 
cnd2 = (rc >= r0_i);
%--------------------------------------------------------------------------

% Calculate the magnetic field in both cases
%--------------------------------------------------------------------------
Bx_fr_1 = Bc.*J0; 
By_fr_1 = (-Bc.*J1 .*(ZL - h0)) ./ sqrt( (YL - L0).^2 + (ZL - h0).^2 ) + ...
        (I0./(2.*pi)) .* ( (ZL + h0) ./ ( (YL - L0).^2 + (ZL + h0).^2 ) ) ;
Bz_fr_1 = (Bc.*J1 .*(YL - L0)) ./ sqrt( (YL - L0).^2 + (ZL - h0).^2 ) - ...
        (I0./(2.*pi)) .* ( (YL - L0) ./ ( (YL - L0).^2 + (ZL + h0).^2 ) ) ;

Bx_fr_2 = Bc.*J0; 
By_fr_2 = (-I0./(2.*pi)) .* ( (ZL - h0) ./ ( (YL - L0).^2 + (ZL - h0).^2 ) ) + ...
        (I0./(2.*pi)) .* ( (ZL + h0) ./ ( (YL - L0).^2 + (ZL + h0).^2 ) ) ;
Bz_fr_2 = (I0./(2.*pi)) .* ( (YL - L0) ./ ( (YL - L0).^2 + (ZL - h0).^2 ) ) - ...
        (I0./(2.*pi)) .* ( (YL - L0) ./ ( (YL - L0).^2 + (ZL + h0).^2 ) ) ;
%--------------------------------------------------------------------------

% Get magnetic field satisfying the conditions 
%--------------------------------------------------------------------------
Bx_fr_12 = Bx_fr_1.*cnd1 + Bx_fr_2.*cnd2; 
By_fr_12 = By_fr_1.*cnd1 + By_fr_2.*cnd2; 
Bz_fr_12 = Bz_fr_1.*cnd1 + Bz_fr_2.*cnd2; 

Bx_fr_12_2 = Bx_fr_2; 
By_fr_12_2 = By_fr_2; 
Bz_fr_12_2 = Bz_fr_2; 
%--------------------------------------------------------------------------
% The condition is fullfiled where is it?
%{
sz = size(Bx_fr_1);
[I1_1,I2_1,I3_1] = ind2sub(sz,cnd1);

for i=1:length(yl)
f347 = figure(347);
hc=pcolor(xl,zl,squeeze(Bz_fr_12(i,:,:))');
%dum=rand(size(Bx_fr)).*cnd2;
%dum=rc;
%hc=pcolor(xl,zl,squeeze(dum(i,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$both RC$','Interpreter','latex','FontSize',fs)
end
%}
%subplot(2,3,2)
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
XYZ7=stream3(XL,YL,ZL,Bx_fr_12_2,By_fr_12_2,Bz_fr_12_2,xstart,ystart,zstart);
f37=figure(37);
hline=streamline(XYZ7);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('Flux rope field lines','Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
f27=figure(27);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bz_fr_12(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{z}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(By_fr_12(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bz_fr_12(sli,:,:).^2 + By_fr_12(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
h=colorbar;
yt=get(h,'XTick');
set(h,'XTickLabel',sprintf('%2.2f',yt));
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{yz}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bz_fr_12(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{z}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(By_fr_12(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bz_fr_12(:,sli,:).^2 + By_fr_12(:,sli,:).^2) )');
%hc=pcolor(yl,zl,squeeze(sqrt(Bxfm(:,sli,:).^2 ))');
set(hc,'edgecolor','none')
colormap(BWR);
cb=colorbar;
%set(cb, 'Ticks', [-1, 0 1], 'TickLabels', {'-1', '0', '1'})
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{yz}$','Interpreter','latex','FontSize',fs)
sgtitle('Flux rope, conditions $$(r_{c} > r_{0}); (r_{c} < r_{0})$$, Inoue and Kusano, 2006, Eq 12 to 17' ,'Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------

%
%--------------------------------------------------------------------------
% Multiple Dipoles magnetic field cartesian
%--------------------------------------------------------------------------
xpl = randi([2 (Lx-2)],1,8); 
ypl = randi([2 (Ly-2)],1,8); 
zpl = randi([2 (Lz-2)],1,8); 
M_i=1.0;
M=M_i*ones(size(XL));
%--------------------------------------------------------------------------

Bx1d_m_i=zeros(size(Bx1d)); By1d_m_i=zeros(size(Bx1d)); Bz1d_m_i=zeros(size(Bx1d));
dumx = Bx1d_m_i; dumy = By1d_m_i; dumz = Bz1d_m_i;
%B1d_mc=cell(3,8);
for i=1:8    
    xp =xpl(i); yp =ypl(i); zp = zpl(i);
    RL=sqrt((XL-xp).^2 + (YL-yp).^2 + (ZL-zp).^2);
    Bx1d_m_i=3.*M.*(XL-xp).*(ZL-zp)./(RL.^5);
    By1d_m_i=3.*M.*(YL-yp).*(ZL-zp)./(RL.^5);
    Bz1d_m_i=M.*(3.* ((ZL-zp).^2) - RL.^2)./(RL.^5);    
    dumx = dumx + Bx1d_m_i;
    dumy = dumy + By1d_m_i;
    dumz = dumz + Bz1d_m_i;
   % B1d_mc{1,i} = Bx1d_m_i; B1d_mc{2,i} = By1d_m_i; B1d_mc{3,i} = Bz1d_m_i;  
end

Bx1dm=dumx;
By1dm=dumy;
Bz1dm=dumz;
%--------------------------------------------------------------------------
XYZ1=stream3(XL,YL,ZL,Bx1dm,By1dm,Bz1dm,xstart,ystart,zstart);
f39=figure(39);
hline=streamline(XYZ1);
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs); xlim([xl(1),xl(end)])
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',fs); ylim([yl(1),yl(end)])
zlabel('$z/d_{i}$','Interpreter','latex','FontSize',fs); zlim([zl(1),zl(end)])
view(3);
title('Multiple dipoles field lines','Interpreter','latex','FontSize',fs)

%--------------------------------------------------------------------------
f29=figure(29);
subplot(2,3,1)
hc=pcolor(xl,zl,squeeze(Bx1dm(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,2)
hc=pcolor(xl,zl,squeeze(By1dm(sli,:,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,3)
hc=pcolor(xl,zl,squeeze(sqrt(Bx1dm(sli,:,:).^2 + By1dm(sli,:,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
subplot(2,3,4)
hc=pcolor(yl,zl,squeeze(Bx1dm(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{x}$','Interpreter','latex','FontSize',fs)
subplot(2,3,5)
hc=pcolor(yl,zl,squeeze(By1dm(:,sli,:))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{y}$','Interpreter','latex','FontSize',fs)
subplot(2,3,6)
hc=pcolor(yl,zl,squeeze(sqrt(Bx1dm(:,sli,:).^2 + By1dm(:,sli,:).^2))');
set(hc,'edgecolor','none')
colormap(BWR);
colorbar
ylabel('$z/d_{i}$','Interpreter','latex','FontSize',fs)
xlabel('$y/d_{i}$','Interpreter','latex','FontSize',fs)
title('$B_{xy}$','Interpreter','latex','FontSize',fs)
sgtitle('One Dipole, $B_{x} = 3Mxz/r^{5}; B_{y} = 3Myz/r^{5}; B_{z} = M(3z^{2} -r^(2))/r^{5}; M =$' + string(M_i),'Interpreter','latex','FontSize',fs)


%--------------------------------------------------------------------------

% Save the plots
%--------------------------------------------------------------------------
%{
cd '/Volumes/PSC_DiRAC_DATA/MagneTORE/Synthetic_field_configurations'
saveas(f21,strcat('dipole_components.png'));
saveas(f31,strcat('dipole_field_lines.png'));
saveas(f22,strcat('single_harris_components.png'));
saveas(f32,strcat('single_harris_field_lines.png'));
saveas(f23,strcat('single_harris_FF_components.png'));
saveas(f33,strcat('single_harris_FF_field_lines.png'));
saveas(f24,strcat('Multiple_harris_FF_components.png'));
saveas(f34,strcat('Multiple_harris_FF_field_lines.png'));
saveas(f25,strcat('Elongated_dipole_components.png'));
saveas(f35,strcat('Elongated_dipole_field_lines.png'));
saveas(f26,strcat('Two_Elongated_dipole_components.png'));
saveas(f36,strcat('Two_Elongated_dipole_field_lines.png'));
saveas(f27,strcat('Flux_rope_components.png'));
saveas(f37,strcat('Flux_rope_field_lines.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%}
%--------------------------------------------------------------------------

% Check the divergence
%--------------------------------------------------------------------------
% Dipole Bx1d, By1d ,Bz1d
% single Harris Bx ,By, Bz
% single Harris FF Bxf ,Byf, Bzf
% Multiple Harris FF Bxfm ,Byfm, Bzfm
% Elongated dipole Bx_d ,By_d, Bz_d
% Two Elongated dipoles Bx_d12 ,By_d12, Bz_d12
% Flux rope Bx_fr_12 ,By_fr_12, Bz_fr_12

%divB_h = divergence(XL,YL,ZL,Bx,By,Bz);

DIVB=cell(1,9);
DIVB{1,1} = divergence(XL,YL,ZL,Bx1d,By1d,Bz1d);
DIVB{1,2} = divergence(XL,YL,ZL,Bx ,By, Bz);
DIVB{1,3} = divergence(XL,YL,ZL,Bxf ,Byf, Bzf);
DIVB{1,4} = divergence(XL,YL,ZL,Bxfm ,Byfm, Bzfm);
DIVB{1,5} = divergence(XL,YL,ZL,Bx_d ,By_d, Bz_d);
DIVB{1,6} = divergence(XL,YL,ZL,Bx_d12 ,By_d12, Bz_d12);
DIVB{1,7} = divergence(XL,YL,ZL,Bx_fr_12, By_fr_12, Bz_fr_12);
DIVB{1,8} = divergence(XL,YL,ZL,Bx_fr_12_2, By_fr_12_2, Bz_fr_12_2);
DIVB{1,9} = divergence(XL,YL,ZL,Bx1dm,By1dm,Bz1dm);

f3411=figure(3411);
subplot(3,3,1)
dump=DIVB{1,1}; histogram(dump,100); title('Dipole')
subplot(3,3,2)
dump=DIVB{1,2}; histogram(dump,100); title('Harris')
subplot(3,3,3)
dump=DIVB{1,3}; histogram(dump,100); title('Harris FF')
subplot(3,3,4)
dump=DIVB{1,4}; histogram(dump,100); title('M Harris FF')
subplot(3,3,5)
dump=DIVB{1,5}; histogram(dump,100); title('E dipole')
subplot(3,3,6)
dump=DIVB{1,6}; histogram(dump,100); title('Two E dipoles')
subplot(3,3,7)
dump=DIVB{1,7}; histogram(dump,100); title('Flux rope')
subplot(3,3,8)
dump=DIVB{1,8}; histogram(dump,100); title('E dipole 2')
subplot(3,3,9)
dump=DIVB{1,9}; histogram(dump,100); title('Multiple dipoles')
sgtitle('$$\nabla \cdot B$$' ,'Interpreter','latex','FontSize',fs)

%--------------------------------------------------------------------------

%


% This works
%
%--------------------------------------------------------------------------
% Create the data file
cd '/Volumes/PSC_DiRAC_DATA/MagneTORE/Synthetic_field_configurations'

for i=1:9
    if (i==1)
        filename='dipole';
        filenameh5 = 'dipole.h5'; 
        Bxh5 = permute(Bx1d,[2,1,3]);
        Byh5 = permute(By1d,[2,1,3]);
        Bzh5 = permute(Bz1d,[2,1,3]);
    elseif (i==2)
        filename='harris';
        filenameh5 = 'harris.h5'; 
        Bxh5 = permute(Bx,[2,1,3]);
        Byh5 = permute(By,[2,1,3]);
        Bzh5 = permute(Bz,[2,1,3]);
    elseif (i==3)
        filename='harris_FF';
        filenameh5 = 'harris_FF.h5'; 
        Bxh5 = permute(Bxf,[2,1,3]);
        Byh5 = permute(Byf,[2,1,3]);
        Bzh5 = permute(Bzf,[2,1,3]);
    elseif (i==4)
        filename='Multiple_harris_FF';
        filenameh5 = 'Multiple_harris_FF.h5'; 
        Bxh5 = permute(Bxfm,[2,1,3]);
        Byh5 = permute(Byfm,[2,1,3]);
        Bzh5 = permute(Bzfm,[2,1,3]);
    elseif (i==5)
        filename='Elongated_dipole';
        filenameh5 = 'Elongated_dipole.h5'; 
        Bxh5 = permute(Bx_d,[2,1,3]);
        Byh5 = permute(By_d,[2,1,3]);
        Bzh5 = permute(Bz_d,[2,1,3]);    
    elseif (i==6)
        filename='Two_Elongated_dipoles';
        filenameh5 = 'Two_Elongated_dipoles.h5'; 
        Bxh5 = permute(Bx_d12,[2,1,3]);
        Byh5 = permute(By_d12,[2,1,3]);
        Bzh5 = permute(Bz_d12,[2,1,3]);
    elseif (i==7)
        filename='Flux_rope';
        filenameh5 = 'Flux_rope.h5'; 
        Bxh5 = permute(Bx_fr_12,[2,1,3]);
        Byh5 = permute(By_fr_12,[2,1,3]);
        Bzh5 = permute(Bz_fr_12,[2,1,3]); 
    elseif (i==8)
        filename='Elongated_dipole_2';
        filenameh5 = 'Elongated_dipole_2.h5'; 
        Bxh5 = permute(Bx_fr_12_2,[2,1,3]);
        Byh5 = permute(By_fr_12_2,[2,1,3]);
        Bzh5 = permute(Bz_fr_12_2,[2,1,3]); 
    elseif (i==9)
        filename='Multiple_dipoles';
        filenameh5 = 'Multiple_dipoles.h5'; 
        Bxh5 = permute(Bx1dm,[2,1,3]);
        Byh5 = permute(By1dm,[2,1,3]);
        Bzh5 = permute(Bz1dm,[2,1,3]); 
    end

%if isfile('dipole.h5'), delete dipole.h5; end

ssA=size(permute(Bx1d,[2,1,3]));
xh5=permute(XL,[2,1,3]); yh5=permute(YL,[2,1,3]); zh5=permute(ZL,[2,1,3]);
% y,x,z

%{
%h5create(filenameh5,'/g1/XL',ssA)
%h5create(filenameh5,'/g1/YL',ssA)
%h5create(filenameh5,'/g1/ZL',ssA)
h5create(filenameh5,'/g2/Bx',ssA)
h5create(filenameh5,'/g2/By',ssA)
h5create(filenameh5,'/g2/Bz',ssA)
% write the data file
%h5write(filenameh5, '/g1/XL', xh5)
%h5write(filenameh5, '/g1/YL', yh5)
%h5write(filenameh5, '/g1/ZL', zh5)
h5write(filenameh5, '/g2/Bx', Bxh5)
h5write(filenameh5, '/g2/By', Byh5)
h5write(filenameh5, '/g2/Bz', Bzh5)
%}

%if isfile(filenameh5), delete string(filenameh5); end %this doesn't work
%yet
% Create the h5 data file
h5create(filenameh5,'/Bx',ssA)
h5create(filenameh5,'/By',ssA)
h5create(filenameh5,'/Bz',ssA)
% Write the h5 data file
h5write(filenameh5, '/Bx', Bxh5)
h5write(filenameh5, '/By', Byh5)
h5write(filenameh5, '/Bz', Bzh5)

% This only works saving the scalar components
% Also I have to change by hand the resolution in the xdmfwrite.m file % <------
xdmfwrite(filenameh5,ssA,'double');
end



%--------------------------------------------------------------------------

% This is to create the xmdf file. This try did not work
%--------------------------------------------------------------------------
%{
fileID = fopen( strcat(filename,'.xmf'),'w');
fprintf(fileID,'<?xml version="1.0" ?>\n');
fprintf(fileID,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fileID,'<Xdmf Version="2.0">\n');
fprintf(fileID,'<Domain>\n');
fprintf(fileID,'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n');
fprintf(fileID,'<Grid Name="XYZ" GridType="Uniform">\n');
fprintf(fileID,'<Time Type="Single" Value="0.00000" />\n');
fprintf(fileID,'<Topology TopologyType="3DCoRectMesh" NumberOfElements="%d %d %d"/>\n',Nz,Ny,Nx);
fprintf(fileID,'<Geometry GeometryType="Origin_DxDyDz">\n');
fprintf(fileID,'<DataItem Name="Origin" DataType="Float" Dimensions="3" Format="XML"> "%d %d %d" </DataItem>\n',zl(1),yl(1),xl(1));
fprintf(fileID,'<DataItem Name="DxDyDz" DataType="Float" Dimensions="3" Format="XML"> "%d %d %d" </DataItem>\n',dz,dy,dx);
fprintf(fileID,'</Geometry>\n');
fprintf(fileID,'<Attribute Name="B" AttributeType="Vector" Center="Cell">\n');
fprintf(fileID,'<DataItem Dimensions="%d %d %d 3" NumberType="Float" Precision="8" Format="HDF">\n',nz,ny,nx);
fprintf(fileID,'%s:/g2\n',filenameh5);
%:/B
fprintf(fileID,'</DataItem>\n');
fprintf(fileID,'</Attribute>\n');
fprintf(fileID,'</Grid>\n');
fprintf(fileID,'</Grid>\n');
fprintf(fileID,'</Domain>\n');
fprintf(fileID,'</Xdmf>\n');
fclose(fileID);
%}
%--------------------------------------------------------------------------

% This is the other example which didn't work either. the volume problem
%--------------------------------------------------------------------------
%{
% Define dimension of the 3D-array
Lx = 1.0;  nx = 121; 
Ly = 1.0;  ny = 111;
Lz = 1.0;  nz = 101;
% Build axis of the mesh
dx = Lx/(nx-1);  x = 0:dx:Lx;
dy = Ly/(ny-1);  y = 0:dy:Ly;
dz = Lz/(nz-1);  z = 0:dz:Lz;
% Build mesh rectangular grid
[X,Y,Z] = meshgrid(y,x,z);
% Initial condition of pressure field
pressure = 1.0*exp(-( (X-0.5).^2+(Y-0.5).^2+(Z-0.5).^2 ) ./ (0.31)^2 );
% Create geometry file
if isfile('geometry.h5'), delete geometry.h5; end
% x-axis
h5create('geometry.h5','/x_nodes',nx);
h5write('geometry.h5','/x_nodes',x);
h5disp('geometry.h5');
% y-axis
h5create('geometry.h5','/y_nodes',ny);
h5write('geometry.h5','/y_nodes',y);
h5disp('geometry.h5');
% z-axis
h5create('geometry.h5','/z_nodes',nz);
h5write('geometry.h5','/z_nodes',z);
h5disp('geometry.h5');
% Create data file
if isfile('fields_it0.h5'), delete fields_it0.h5; end
% Pressure field
h5create('fields_it0.h5','/p',[nx,ny,nz]);
h5write('fields_it0.h5','/p',pressure);
%--------------------------------------------------------------------------
h5write('fields_it0.h5','/p',pressure);
%--------------------------------------------------------------------------
h5disp('fields_it0.h5');
% Numbers precision
precision = 8; % bit for real double 
time = 0.0; % assuming t=0
it = 0; % assuming this is the initial condition

% This doesnt produce a volume option and the variables are saved in point
% arrays
%--------------------------------------------------------------------------
% write associated XDMF file
fileID = fopen('openWithParaView.xmf','w');
fprintf(fileID,'<?xml version="1.0" ?>\n');
fprintf(fileID,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fileID,'<Xdmf Version="2.0">\n');
fprintf(fileID,'<Domain>\n');
fprintf(fileID,'<Grid Name="mesh" GridType="Uniform">\n');
fprintf(fileID,'<Topology TopologyType="3DRectMesh" NumberOfElements="%d %d %d"/>\n',nz,ny,nx);
fprintf(fileID,'<Geometry GeometryType="VXVYVZ">\n');
fprintf(fileID,'<DataItem Name="coordx" Dimensions="%g" NumberType="Float" Precision="%g" Format="HDF">\n',nx,precision);
fprintf(fileID,'geometry.h5:/x_nodes\n');
fprintf(fileID,'</DataItem>\n');
fprintf(fileID,'<DataItem Name="coordy" Dimensions="%g" NumberType="Float" Precision="%g" Format="HDF">\n',ny,precision);
fprintf(fileID,'geometry.h5:/y_nodes\n');
fprintf(fileID,'</DataItem>\n');
fprintf(fileID,'<DataItem Name="coordz" Dimensions="%g" NumberType="Float" Precision="%g" Format="HDF">\n',nz,precision);
fprintf(fileID,'geometry.h5:/z_nodes\n');
fprintf(fileID,'</DataItem>\n');
fprintf(fileID,'</Geometry>\n');
fprintf(fileID,'<Time TimeType="Single" Value="%g"/>\n',time);
%fprintf(fileID,'<Attribute Name="p" AttributeType="Scalar" Center="Node">\n');
fprintf(fileID,'<Attribute Name="p" AttributeType="Scalar" Center="Cell">\n');
fprintf(fileID,'<DataItem Dimensions="%d %d %d" NumberType="Float" Precision="%d" Format="HDF">\n',nz,ny,nx,precision);
fprintf(fileID,'fields_it%g.h5:/p\n',it);
fprintf(fileID,'</DataItem>\n');
fprintf(fileID,'</Attribute>\n');
fprintf(fileID,'</Grid>\n');
fprintf(fileID,'</Domain>\n');
fprintf(fileID,'</Xdmf>\n');
fclose(fileID);
%}
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
%
% Trying conditions
%{
%--------------------------------------------------------------------------
A=[1,2,3,4;5,6,7,8;9,10,11,12]; % 3x4 matrix array

aa = find(A>6); % this returns the list of indices nxm. This goes rows firts
B1 = A(A>6) %this returns a 1s array of the values in those indices
B = A.*(A>6) %this returns the matrix with the elements fullfiling the condition
aaa=((A>6) & (A<=11))
B = A.*aaa
C = A.*(A<=6)
AA=B+C

sz = size(A);
[I1,I2] = ind2sub(sz,aa);
%--------------------------------------------------------------------------


f34=figure(34)
%histogram(squeeze(sqrt(Bxfm(:,sli,:).^2 + Byfm(:,sli,:).^2) )',50);
histogram(squeeze(Byfm(:,sli,:)  )',50);

divB_h = divergence(XL,YL,ZL,Bx,By,Bz); 
f3411=figure(3411);
histogram(divB_h)

% For the streamlines
fr=1;
NN=200;  %Do this just once


%dx=xl(2)-xl(1); dy=yl(2)-yl(1); dz=zl(2)-zl(1); 
%[sx,sy,sz] = meshgrid(xl(1):8*dx:xl(end),yl(1):8*dy:yl(end),zl(83));
%XYZ2=stream3(XL,YL,ZL,Bx,By,Bz,sx,sy,sz);


XYZ2=stream3(XL,YL,ZL,Bx,By,Bz,xstart,ystart,zstart);

f3411=figure(3411);
histogram(divB_h)

f3421=figure(3421);
hline=streamline(XYZ2);
view(3);

f3422=figure(3422);
pcolor(Bx(:,1,:));
view(3);
%--------------------------------------------------------------------------

f343=figure(343);
scatter3(xstart,ystart,zstart); 
view(3);


%--------------------------------------------------------------------------

mydata = rand(100,100,400);


f88=figure(88)
pcolor(Bz(:,:,4))


pt=[50,50,100];
rl=sqrt((xl-pt(1)).*(xl-pt(1)) + (yl-pt(2)).*(yl-pt(2)) + (zl-pt(3)).*(zl-pt(3)));
Rhat=[XL./RL, YL./RL, ZL./RL];
ui=[1,0,0]; uj=[0,1,0]; uk=[0,0,1];

%}
