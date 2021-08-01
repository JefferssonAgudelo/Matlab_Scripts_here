%--------------------------------------------------------------------------
%This is to do only plots no calculations and it goes with the
%Computing_energy_reconnection4. This scripts is made only for the part in
%the original reference frame
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% These are the plots including the streamlines
%--------------------------------------------------------------------------
X = B_Gp{1,3};
%--------------------------------------------------------------------------
% built the grid to do the right plots
% this migth ruin the plots as it was used for the untranformed quantities
size_X=size(X);
xll=linspace(1,size_X(2),size_X(2))*0.06; %Jeff flag it works with (1,size_x(2))
yll=linspace(1,size_X(1),size_X(1))*0.06;
[XLL, YLL] = meshgrid(xll,yll);

resx = abs(xll(2)-xll(1)); 
resy = abs(yll(2)-yll(1)); 
%--------------------------------------------------------------------------

% The point where the reconnection is happening is around 
xp_mr = xll(99); % this is in the new RF
yp_mr = yll(127);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% This are the projections for the magnetic field
ubx = B_Gp{1,2};   vby = B_Gp{1,1}; wbz = B_Gp{1,3};

%--------------------------------------------------------------------------

% This is to check a dummy field. The plot needs to be run after the ccc
% streamlines have been calculated.
%--------------------------------------------------------------------------
%{
N = 100 ;
xl=linspace(1,6,N);%-2:4/(N-1) :2 ;
yl=linspace(1,6,N);% -2:4/(N-1 ) : 2;
[ x , y]=meshgrid(xl, yl) ;
z = (x -3 ).* exp(-(x -3).^ 2 - (y-3) .^ 2 ) ;
[ px , py ] = gradient ( z , 4 / (N-1) ,4/(N-1) ) ;
ubx = px; vby= py; wbz =px./px;
XLL=x ; YLL=y; xll=xl; yll=yl; 
resx = abs(xll(2)-xll(1)); 
resy = abs(yll(2)-yll(1)); 
NN=20; 
xx = linspace(xl(1),xl(end),NN);
yy = linspace(yl(1),yl(end),NN);
[xx yy] = meshgrid(xx,yy);
aaa51=stream2(x,y,ubx,vby,xx,yy);

%-------------------------------------------------------------
f780=figure(780)
h1=subplot(2,2,1);
h2=quiver (px , py ) ;
set(h2,'AutoScale','on', 'AutoScaleFactor', 2, 'Color', 'k');
h3=subplot(2,2,2);
dum_p = abs(real(EIGV1).*saddle_cond12./real(EIGV1));
hc=pcolor(dum_p');
%hc=pcolor(ubx);
colorbar
set(hc,'edgecolor','none'); colormap(jet(2))
h3=subplot(2,2,3);
hlines1=streamline(aaa51);
set(hlines1,'LineWidth',1,'Color', 'k');
h3=subplot(2,2,4);
hlines1=streamline(ccc12);
set(hlines1,'LineWidth',1.5,'Color', 'm');
%-------------------------------------------------------------
%}
%--------------------------------------------------------------------------

bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
%ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
%vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

% This is for the electron velocities vectors
Vev_x = ve_Gp{1,2}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,3}; 
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
%Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
%Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

Viv_x = vi_Gp{1,3}; Viv_y = vi_Gp{1,1}; Viv_z = vi_Gp{1,2}; 
Vivm=sqrt(Viv_x.*Viv_x + Viv_y.*Viv_y + Viv_z.*Viv_z);
%Viv_x=Viv_x.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));
%Viv_y=Viv_y.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));


% In this par reference frame there is no projection
%ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};
%Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,2};
%Viv_x = vi_Gp{1,3}; Viv_y = vi_Gp{1,1}; Viv_z = vi_Gp{1,2}; 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


%NN=50; xstart = (max(xll)/1.2)*rand(NN,1); ystart = (max(yll)/1.2)*rand(NN,1); %Do this just once

NN=40; 
xstart1 = (resx)*ones(1,NN); ystart1 = linspace(0,yll(end),NN); 
xstart2 = linspace(0,xll(end),NN); ystart2 = (resy)*ones(1,NN); 
xstart3 = (xll(end)-resx)*ones(1,NN); ystart3 = linspace(0,yll(end),NN); 
xstart4 = linspace(0,xll(end),NN); ystart4 = (yll(end)-resy)*ones(1,NN); 
xstart = horzcat(xstart1,xstart2,xstart3,xstart4);
ystart = horzcat(ystart1,ystart2,ystart3,ystart4);
aaa=stream2(XLL,YLL,ubx,vby,xstart,ystart);

NN4=10; 
xx = linspace(xll(10),xll(end-10),NN4);
yy = linspace(yll(10),yll(end-10),NN4);
[xx, yy] = meshgrid(xx,yy);
aaa5=stream2(XLL,YLL,ubx,vby,xx,yy);

NN5=10; 
xx5 = linspace(xll(10),xll(end-10),NN5);
yy5 = linspace(yll(10),yll(end-10),NN5);
[xx5, yy5] = meshgrid(xx5,yy5);

aaa=stream2(XLL,YLL,ubx,vby,xstart,ystart);
aaam=stream2(XLL,YLL,-ubx,-vby,xstart,ystart);

%veve=stream2(XLL,YLL,Vev_x,Vev_y,xstart,ystart);
%vivi=stream2(XLL,YLL,Viv_x,Viv_y,xstart,ystart);
veve5=stream2(XLL,YLL,Vev_x,Vev_y,xx5,yy5);
vivi5=stream2(XLL,YLL,Viv_x,Viv_y,xx5,yy5);
bbbb5=stream2(XLL,YLL,ubx,vby,xx5,yy5);

veve5m=stream2(XLL,YLL,-Vev_x,-Vev_y,xx5,yy5);
vivi5m=stream2(XLL,YLL,-Viv_x,-Viv_y,xx5,yy5);
bbbb5m=stream2(XLL,YLL,-ubx,-vby,xx5,yy5);

%--------------------------------------------------------------------------
%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
dCp=25;
colorFVBl=[0 0 0 0.5];
sv=3;
asf=3;
%--------------------------------------------------------------------------
% This is the plot for the energy terms in the kinetic equation
%--------------------------------------------------------------------------

% 2D plots for the kinetic energy rate terms
%
%--------------------------------------------------------------------------
GPijf{1,1}=dtKe_in_Gp_2{1,1}; GPijf{1,2}=ve_grad_Ke_Gp{1,1};
GPijf{1,3}=Ke_Div_ve_Gp{1,1}; GPijf{1,4}=ve_Div_Pet_Gp{1,1}; 
GPijf{1,5}=-qe_ne_ve_E_Gp{1,1};  
f4=figure(4);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
h1=subaxis(1,5,1,'SV',0,'SH',0.0,'MR',0,'ML',0,'PL',0.025,'PR',0);
dum_p=GPijf{1,1}; %dum_p2 = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h1,BWR); set(hc,'edgecolor','none'); 
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -2;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs);
ylh=ylabel('$x / d_{i}$','Interpreter','latex','FontSize',fs);
xlh.Position(2) = xlh.Position(2) + 0.4; ylh.Position(1) = ylh.Position(1) + 0.4; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%xlh.Position(2) = xlh.Position(2) + abs(xlh.Position(2) * 0.1);
title('$\partial \varepsilon^{k}_{e} / \partial t $','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h2=subaxis(1,5,2,'SV',0,'SH',0.000,'MR',0,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,2}; 
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-1 5e-1]); %caxis([-lim_yp lim_yp]); 
colormap(h2,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex'; ax.YTick = [];
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$u_{e} \cdot \nabla \varepsilon^{k}_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h3=subaxis(1,5,3,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,3}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-1 5e-1]); %caxis([-lim_yp lim_yp]); 
colormap(h3,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$\varepsilon^{k}_{e}\nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h4=subaxis(1,5,4,'SV',0,'SH',0.00,'MR',0,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,4}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-2e3 2e3]); %caxis([-lim_yp lim_yp]); 
colormap(h4,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$u_{e} \cdot \nabla \cdot \overline{P}_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h5=subaxis(1,5,5,'SV',0,'SH',0.00,'MR',0.001,'ML',0,'PL',0.02,'PR',0.005);
dum_p=GPijf{1,5}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-4e-0 4e-0]); %caxis([-lim_yp lim_yp]); 
colormap(h5,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$-q_{e}n_{e}E\cdot u_{e}$','Interpreter','latex','FontSize',fs);
hold off
clearvars GPijf
%--------------------------------------------------------------------------


GPijf{1,1}=Zernitani_e_Gp{1,1}; GPijf{1,2}=p_theta_e_Gp{1,1};
GPijf{1,3}=PiD_e_Gp{1,1};   
f41=figure(41);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
h1=subaxis(1,3,1,'SV',0,'SH',0.0,'MR',0,'ML',0,'PL',0.025,'PR',0);
dum_p=GPijf{1,1}; %dum_p2 = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
%caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h1,BWR); set(hc,'edgecolor','none'); 
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -2;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs);
ylh=ylabel('$x / d_{i}$','Interpreter','latex','FontSize',fs);
xlh.Position(2) = xlh.Position(2) + 0.4; ylh.Position(1) = ylh.Position(1) + 0.4; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%xlh.Position(2) = xlh.Position(2) + abs(xlh.Position(2) * 0.1);
title('$D_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h2=subaxis(1,3,2,'SV',0,'SH',0.000,'MR',0,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,2}; 
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-1e3 1e3]); %caxis([-lim_yp lim_yp]); 
colormap(h2,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex'; ax.YTick = [];
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$p\theta_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h3=subaxis(1,3,3,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,3}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-1e3 1e3]); %caxis([-lim_yp lim_yp]); 
colormap(h3,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$y / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$PiD_{e}$','Interpreter','latex','FontSize',fs);
hold off

clearvars GPijf
%--------------------------------------------------------------------------
