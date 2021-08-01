% separatrix.m
clear all; close all;
hold on
for phimax=pi/5:pi/5:pi
phidot=@(phi)sqrt(2*(cos(phi)-cos(phimax))); % this is a way to define the function before defining the arrays!
phi=-phimax:phimax/500:phimax;
om=phidot(phi);
if phimax==pi
    plot(phi,om,'k',phi,-om,'k-','LineWidth',2);
else
plot(phi,om,'k-.',phi,-om,'k-.','LineWidth',2);
end
end
xlim([-pi,pi]); ylim([-2.2,2.2]);
xlabel('\phi'); ylabel('(d\phi/dt)/\Omega_s')
text(1.35,1.7,'Separatrix','FontSize',16)
set(gca,'xtick',[-pi,-pi/2,0,pi/2,pi],'fontsize',16, ...
'xticklabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})


%-------------------------------------------------------------------------
%{
% The point where the reconnection is happening is around 
xp_mr = xll(99); % this is in the new RF
yp_mr = yll(127);
%--------------------------------------------------------------------------

% This are the projections for the magnetic field
%ubx = slice_Gp{1,1};   vby = slice_Gp{1,2}; wbz = slice_Gp{1,3};
ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};

%ubx = B_Gp{1,2};   vby = B_Gp{1,1}; wbz = B_Gp{1,3}; %Is it this way? it doesn't look that way

bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

% This is for the electron velocities vectors
%Vev_x = slice_Gp{2,1}; Vev_y = slice_Gp{2,2}; Vev_z = slice_Gp{2,3}; 
Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,2}; 

%Vev_x = ve_Gp{1,2}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,3}; 

Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

Vivm=sqrt(Viv_x.*Viv_x + Viv_y.*Viv_y + Viv_z.*Viv_z);
Viv_x=Viv_x.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));
Viv_y=Viv_y.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));

%NN=50; xstart = (max(xll)/1.2)*rand(NN,1); ystart = (max(yll)/1.2)*rand(NN,1); %Do this just once

NN=40; 
xstart1 = (resx)*ones(1,NN); ystart1 = linspace(0,yll(end),NN); 
xstart2 = linspace(0,xll(end),NN); ystart2 = (resy)*ones(1,NN); 
xstart3 = (xll(end)-resx)*ones(1,NN); ystart3 = linspace(0,yll(end),NN); 
xstart4 = linspace(0,xll(end),NN); ystart4 = (yll(end)-resy)*ones(1,NN); 
xstart = horzcat(xstart1,xstart2,xstart3,xstart4);
ystart = horzcat(ystart1,ystart2,ystart3,ystart4);

%}
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
[X,Y] = meshgrid(phi,phi);
f = sin(Y).*cos((cos(X)));
%--------------------------------------------------------------------------

f32=figure(32)
histogram(sin(phi),50)

% Get the critical points (ubx = 0;   vby = 0; wbz =0;)
% Define the zero
%--------------------------------------------------------------------------
zero_d=3e-2;
f_cr_nan=isnan(f);
f_cr_i=((abs(f) < zero_d));
fXY_cr_i=((abs(X) < zero_d) & (abs(Y) < zero_d));
%--------------------------------------------------------------------------

% To get the Jacobian 
%--------------------------------------------------------------------------
[f_x,f_y] = gradient(f);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
f3=figure(3);
h1=subplot(2,2,1);
hc=pcolor(f);
colorbar;colormap(h1,BWR);
set(hc,'edgecolor','none')
%surf(f)
h2=subplot(2,2,2);
%contour(f)
surf(f)
h1=subplot(2,2,3);
hc=pcolor(f.*f_cr_i);
colorbar;colormap(h1,BWR);
set(hc,'edgecolor','none')
%surf(f)
h1=subplot(2,2,4);
hc=pcolor(f.*fXY_cr_i);
colorbar;colormap(h1,BWR);
set(hc,'edgecolor','none')
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% large_oscillation_period.m
clear all; close all
hold off
phihat=0:0.02:pi;
m=sin(phihat/2).^2;
Tp=(2/pi)*ellipke(m);
plot(phihat,Tp,'k','LineWidth',2)
xlim([0,pi]); ylim([0,5.7])
xlabel({'$\hat{\phi}$'},'Interpreter','latex')
ylabel('T_p/T_s')
set(gca,'xtick',[0,pi/2,pi],'fontsize',16, ...
'xticklabels',{'0','\pi/2','\pi'})


% pendulumtracker.m
function xout=pendulumtracker(x,omega,dt)
k2=(0.5*x(2)/omega)^2+sin(0.5*x(1))^2; k=sqrt(k2);
if (x(1)>pi) x(1)=x(1)-2*pi; end
if (x(1)<-pi) x(1)=x(1)+2*pi; end
s=1; if (x(1)<0) s=-s; x(1)=-x(1); end
s1=1; if (x(2)<0) s1=-s1; end
if (k>1) % outside separatrix; check this part next
% kelf=ellipticK(1/k2); % slow, uses symbolic
kelf=ellipke(1/k2); % much faster
trev=2*kelf/(k*omega);
t0=mod(dt,trev);
% tmp=s1*k*omega*t0+s*ellipticF(0.5*x(1),1/k2); % slow, uses symbolic
tmp=s1*k*omega*t0+s*elliptic12(0.5*x(1),1/k2); % much faster
[sn,cn,dn]=ellipj(tmp,1/k2);
if (abs(tmp) > kelf) sn=-sn; end
xout(1)=2*asin(sn);
xout(2)=2*s1*omega*k*dn;
else % inside separatrix
% trev=4*ellipticK(k2)/omega; % slow, uses symbolic
trev=4*ellipke(k2)/omega; % much faster
t0=mod(dt,trev);
z0=asin(min(1,sin(0.5*x(1))/k));
% tmp=s1*omega*t0+s*ellipticF(z0,k2); % slow, uses symbolic
tmp=s1*omega*t0+s*elliptic12(z0,k2); % much faster
[sn,cn,dn]=ellipj(tmp,k2);
xout(1)=2*asin(k*sn);
xout(2)=2*s1*omega*k*cn;
end


% This is another set of data
%--------------------------------------------------------------------------
N = 64 ;
xl=-2:4/(N-1) :2 ;
yl = -2:4/(N-1 ) : 2;
[ x , y]=meshgrid(xl, yl) ;
z = x .* exp(-x .^ 2 - y .^ 2 ) ;
[ px , py ] = gradient ( z , 4 / (N-1) ,4/(N-1) ) ;
quiver (px , py ) ;

ubx = px; vby= py;
% Now lets get the Jacobian
[ubx_x,ubx_y] = gradient(ubx);
[vby_x,vby_y] = gradient(vby);

zero_d=0.05;
ubvb_cr_i=((abs(ubx) < zero_d) & (abs(vby) < zero_d));

%ubx_x1 = ubx_x.*ubvb_cr_i;ubx_y1 = ubx_y.*ubvb_cr_i;  % this in matrix representation
%vby_x1 = vby_x.*ubvb_cr_i;vby_y1 = vby_y.*ubvb_cr_i;
ubx_x2 = ubx_x(ubvb_cr_i); ubx_y2 = ubx_y(ubvb_cr_i);  % this in array representation
vby_x2 = vby_x(ubvb_cr_i); vby_y2 = vby_y(ubvb_cr_i);

l2e=length(ubx_x2);
EIG=cell(2,l2e);
for i=1:l2e
    Jac = [ubx_x2(i), ubx_y2(i);vby_x2(i), vby_y2(i)];
    
    [Vast,Dast] = eig(Jac);
    EIG{1,i} = Vast;
    EIG{2,i} = Dast;
end    

vec=EIG{1,2};

vec1=vec(:,1) % vector 1
vec2=vec(:,2) % vector 2

vecx=vec(1,:) % x components
vecy=vec(2,:) % y components

% I need to have the coordinates for the vectors 
% bbb5=stream2(x,y,ubx,vby,xx,yy);

%
%--------------------------------------------------------------------------
NN=20; 
xx = linspace(xl(1),xl(end),NN);
yy = linspace(yl(1),yl(end),NN);
[xx yy] = meshgrid(xx,yy);
aaa5=stream2(x,y,ubx,vby,xx,yy);

xstart1 = (xl(1))*ones(1,NN); ystart1 = linspace(yl(1),yl(end),NN); 
xstart2 = linspace(xl(1),xl(end),NN); ystart2 = (yl(1))*ones(1,NN); 
xstart3 = (xl(end))*ones(1,NN); ystart3 = linspace(yl(1),yl(end),NN); 
xstart4 = linspace(xl(1),xl(end),NN); ystart4 = (yl(end))*ones(1,NN); 

xstart = horzcat(xstart1,xstart2,xstart3,xstart4);
ystart = horzcat(ystart1,ystart2,ystart3,ystart4);

aaa1=stream2(x,y,ubx,vby,xstart,ystart); %This is are the streamlines of the field 
aaa2=stream2(x,y,-ubx,-vby,xstart,ystart); %This is are the streamlines of the field 
aaa3=stream2(x,y,ubx,-vby,xstart,ystart);
aaa4=stream2(x,y,-ubx,vby,xstart,ystart);


%aaa2=stream2(x,y,-ubx,-vby,xl,yl); %This is are the streamlines of the field 
%--------------------------------------------------------------------------
figure213=figure(213);
h1=subplot(2,2,1);
dum_p=z.*ubvb_cr_i;
hc=pcolor(xl,yl,dum_p); set(hc,'edgecolor','none')
lim_yp=max(max(max(abs(dum_p))));  caxis([-lim_yp lim_yp]);
colorbar
h3=subplot(2,2,2);
h2=quiver (px , py ) ;
set(h2,'AutoScale','on', 'AutoScaleFactor', 2, 'Color', 'k');
h3=subplot(2,2,3);
%hold on
hlines1=streamline(aaa5);
set(hlines1,'LineWidth',1,'Color', 'k');
hold on
%
hlines2=streamline(aaa1);set(hlines2,'LineWidth',1.3,'Color', 'r');
hlines3=streamline(aaa2);set(hlines3,'LineWidth',1.3,'Color', 'm');
hlines4=streamline(aaa3);set(hlines4,'LineWidth',1.3,'Color', 'c');
hlines5=streamline(aaa4);set(hlines5,'LineWidth',1,'Color', 'k');
%scatter(xstart,ystart)
%}
%xlim([xll(1) xll(end)]);ylim([yll(1) yll(end)]);
hold off

fid = fopen( ' field.vtk ' , 'w' ) ;
fprintf ( fid , '# vtk DataFile Version 4.0\ n ' ) ;
fprintf ( fid , 'Our dataset \n ' ) ;
fprintf ( fid , 'ASCII\n ' ) ;
fprintf ( fid , 'DATASET STRUCTURED POINTS\n ' ) ;
fprintf ( fid , 'DIMENSIONS %d %d 1\n ' ,N,N) ;
fprintf ( fid , 'ORIGIN 0.000 0.000 0.00\ n ' ) ;
fprintf ( fid , 'SPACING 1 1 1\n ' ) ;
fprintf ( fid , 'POINT DATA %d\n ' , NN) ;
fprintf ( fid , 'VECTORS vectorsfloat \n ' ) ;
for i =1:N*N
fprintf ( fid , '%f %f 0.000\ n ' , px( i ) , py( i ) ) ;
end
fclose ( fid ) ;

%--------------------------------------------------------------------------
