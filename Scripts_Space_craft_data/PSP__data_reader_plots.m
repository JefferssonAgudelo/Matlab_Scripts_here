% Read CDF data from PSP

% PSP
%-------------------------------------------------------------------------
addpath '/Applications/matlab_cdf380_patch-64'

cd '/Volumes/PSC_DiRAC_DATA/DATA_SolO_PSP/PSP/DATA/20191012';
path_psp = '/Volumes/PSC_DiRAC_DATA/DATA_SolO_PSP/PSP/DATA/20191012';

S_psp = dir(fullfile(path_psp,'*.cdf'));
N_psp=numel(S_psp); % number of files to use
H_psp=zeros(N_psp);
i=1;

fileID_MAG_psp =  S_psp(1).name;
PSP_MAG=spdfcdfread(S_psp(1).name);
Info_MAG_psp=spdfcdfinfo(fileID_MAG_psp);

psp_epoch=PSP_MAG{1};
psp_B=PSP_MAG{2};

% Normalization units to SI
nT2T=1e-9;
psp_B = psp_B*nT2T;
psp_Bx=psp_B(:,1);
psp_By=psp_B(:,2);
psp_Bz=psp_B(:,3);
psp_Bmag=sqrt(psp_Bx.*psp_Bx + psp_By.*psp_By + psp_Bz.*psp_Bz);


fileID_SPC_psp =  S_psp(2).name;
PSP_SPC=spdfcdfread(S_psp(2).name);
Info_SPC_psp=spdfcdfinfo(fileID_SPC_psp);

psp_spc_epoch=PSP_SPC{1};
psp_vi=PSP_SPC{25};
% Normalization units to SI
kmtom=1000;
psp_vi=psp_vi*kmtom;
psp_vix=psp_vi(:,1);
psp_viy=psp_vi(:,2);
psp_viz=psp_vi(:,3);
psp_vimag=sqrt(psp_vix.*psp_vix + psp_viy.*psp_viy + psp_viz.*psp_viz);

psp_np = PSP_SPC{19};
% Normalization units to SI
cm_32m_3=(0.01)^3;
psp_np=psp_np*cm_32m_3;


%-------------------------------------------------------------------------



%-------------------------------------------------------------------------
% Normalization
%B are normalized to the background magnetic field 
%B0=1;
%n is normalise to the proton density n0=1
%n0=1;
%Vi normalised to the proton alfven speed V_a=0.1c The alfven speed of
%electrons is Vae=10Vai
%Vai=1; Vae=1;
%E = vB, as B0=0.1 and v in units of c, E0= B0*Vai/c
E0=1;
%T=(beta/2)B^2, Thus, as beta=1 
T0=1;
%The current density J0=1?
%J0=Vai;

mu0=4*pi*1e-7;
mi=1.67e-27;
%-------------------------------------------------------------------------
fr=120;
n=1;
d=psp_epoch(fr*n:fr*(n+1)); %d(isinf(d)|isnan(d)) = 0; 
hx = psp_Bx(fr*n:fr*(n+1)); %hx(isinf(hx)|isnan(hx)) = 0;
hy = psp_By(fr*n:fr*(n+1)); %hy(isinf(hy)|isnan(hy)) = 0;
hz = psp_Bz(fr*n:fr*(n+1)); %hz(isinf(hz)|isnan(hz)) = 0;
hm=sqrt(hx.*hx + hy.*hy + hz.*hz);  %psp_Bmag(fr*n:fr*(n+1)); %hm(isinf(hm)|isnan(hm)) = 0;
B0=mean(hm);

%B0=1;
%dateFormat='HH:MM:SS';
dateFormat='ss';
xlow = d(1);
xup = d(fr);
%-------------------------------------------------------------------------
f2=figure(2);
%-------------------------------------------------------------------------
subplot(4,1,1),
plot(d, hm/B0, '-k','LineWidth',2);set(gca,'FontSize',16)
set(gca,'ycolor','k') 
ylabel('$$\frac{|B|}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
xlabel('$$t(s) (\ r \sim 60d_{i})$$','Interpreter','latex','FontSize',18)
xlim([xlow xup]);
%yt1 = 1*max(hm/B0)/3; yt2 = 2*max(hm/B0)/3; 
yticks([1]);
datetick('x',dateFormat,'keepticks')
%datetick('x',dateFormat)
%rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.4])
%-------------------------------------------------------------------------
subplot(4,1,2);
plot(d,hz/B0, '-k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k')
yticks([-0.75 -0.72]);
%yline(0,'--'); 
ylabel('$$\frac{B_{z}}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%-------------------------------------------------------------------------
subplot(4,1,3);
plot(d,hy/B0, '-k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k') 
yticks([-0.48 -0.42]);
%yline(0,'--'); 
ylabel('$$\frac{B_{y}}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%-------------------------------------------------------------------------
subplot(4,1,4);
plot(d,hx/B0, '-k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k') 
%yline(0,'--'); 
%ylabel('$$\frac{B_{x}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
ylabel('$$\frac{B_{x}}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
yticks([0.48 0.52]);
title('$$20191012;\ \langle v_{sw} \rangle =450 km/s;\ \langle |B| \rangle =6.4 nT;\ R \approx 0.78 AU $$','Interpreter','latex','FontSize',18)
%-------------------------------------------------------------------------
ha=get(gcf,'children');
set(ha(1),'position',[.12 .74 .8 .21]); set(ha(2),'position',[.12 .53 .8 .21]);
set(ha(3),'position',[.12 .32 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
cd '/Volumes/PSC_DiRAC_DATA/';
saveas(f2,'Bxyzm_20191012.png');
%--------------------------------------------------------------------------


% video settings
frames = 500;
vfrate = 25; % frame rate of video
vname = 'videoname';
% figure settings
fhand = figure('visible','off');
set(gcf,'color','w');
set(gcf,'Renderer','zbuffer');
axis vis3d;
rotate3d on;
% video settings
writerObj = VideoWriter(vname,'MPEG-4');
writerObj.FrameRate = vfrate;
writerObj.Quality = 100;
open(writerObj);
% video frames
x = d; %rand(100,1); % some random example data to plot
y = hx/B0; %rand(100,1); % some random example data to plot
for i = 1:200
    % do something to change your figure (i.e. plot new data)
    cla; % clear the axis between frames to avoid plotting on top of the old data
    plot(x(1:i),y(1:i),'k'); % this is just an example
    axis([0 1 0 1]); % this just keeps the axis consistent between frames, you will have to personalise it to your data
    
    % get ready for snapshot
    drawnow;
    cframe = getframe(gcf); % take shot
    writeVideo(writerObj,cframe); % add it to video
end
close(writerObj); % close the video file so it's no longer 'in use' by Matlab

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------

%axis([0,100,-1,1])
%f3=figure(3);
x = double(d);
y = double(hm/B0);
h = animatedline; %(x,y,'Color','k','LineWidth',2);
%datetick('x',dateFormat,'keepticks')
for k = 1:60%length(x)
    addpoints(h,x(k),y(k));
    drawnow
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%A = A( ~any( isnan( A ) | isinf( A ), 2 ),: )

frac=10;
d2=psp_spc_epoch(1:frac); %d2(isinf(d2)|isnan(d2)) = 0;
d3=d2;
vix = psp_vix(1:frac); x=vix; mask = (isnan(x) | isinf(x)); x(mask) = []; vix = x;
viy = psp_viy(1:frac); x=viy; mask = (isnan(x) | isinf(x)); x(mask) = []; viy = x;
viz = psp_viz(1:frac); x=viz; mask = (isnan(x) | isinf(x)); x(mask) = []; viz = x;
vim = psp_vimag(1:frac); x=vim; mask = (isnan(x) | isinf(x)); x(mask) = []; vim = x; 
d3(mask) = []; dvim = d3;
 
%Direccion
%kra 29c # 1a-05, Santa Isabel, 
dt = datestr(d(100080),'mmmm dd, yyyy HH:MM:SS.FFF AM');

ni=psp_np(1:frac); %ni(isinf(ni)|isnan(ni)|ni<0) = 0; 
d3=psp_spc_epoch(1:frac);
x=ni; mask = (isnan(x) | isinf(x) | x<0);
x(mask) = []; ni = x; d3(mask) = []; dni = d3;

n0=mean(ni);
Vai=B0./sqrt(mu0*n0*mi)'; Vai=Vai';
%Vai=Vai(1:frac); Vai(isinf(Vai)|isnan(Vai)) = 0; 

%hx = vix; hy = viy; hz = viz; hm = vim;
%vix = hx; viy = hy; viz = hz; vim = hm;
%vim=psp_np(1:frac);

xpos=31.5;%25.5;
xlpos=30;
dateFormat=15; %21;

f1=figure(1);
%-------------------------------------------------------------------------
subplot(4,1,1),
yyaxis left
line1 = plot(dvim,vim/Vai,'*r','LineWidth',1); set(gca,'FontSize',16)
%xlim([0 xlpos]);
maxi=abs(max(vim/Vai)); mini=abs(min(vim/Vai));
%ylim1=max(maxi,mini); ylim([0 ylim1]);
%yticks([0.2 0.4 0.6]);
set(gca,'ycolor','k') 
ylabel('$$\frac{|v_{i}|}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
yyaxis right
line2 = plot(d, hm/B0, '*k','LineWidth',1);set(gca,'FontSize',16)
%xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
%xlim([0 xlpos]);
%ylim([0 max(hm/B0)]) 
%%yticks([1 2]);
set(gca,'ycolor','k') 
ylabel('$$\frac{|B|}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
xlabel('$r/d_{i}$','Position',[median(d)+0.2 -0.1],'Interpreter','latex','FontSize',18)
datetick('x',dateFormat,'keepticks')
%datetick('x',dateFormat)
%rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.4])
%-------------------------------------------------------------------------
subplot(4,1,2);
yyaxis left
plot(d2,viz/Vai,'*r','LineWidth',1),set(gca,'xtick',[],'FontSize',16)
%xline(jmax_pos,'--');  xline(hmin_pos,'--');
yline(0,'--'); 
%xlim([0 xlpos]);
ylim1=max(abs(max(viz/Vai)), abs(min(viz/Vai)));
%ylim([-ylim1 ylim1]);
set(gca,'ycolor','k') 
ylabel('$$\frac{v_{iz}}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
yyaxis right
plot(d,hz/B0, '*k','LineWidth',1),set(gca,'xtick',[],'FontSize',16)
%xline(jmax_pos,'--');  xline(hmin_pos,'--');
yline(0,'--'); 
%xlim([0 xlpos]);  
ylimh=max(abs(max(hz/B0)), abs(min(hz/B0)));
%ylim([-ylimh ylimh]);
%%yticks([0 1]);
set(gca,'ycolor','k') 
ylabel('$$\frac{B_{z}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
%-------------------------------------------------------------------------
subplot(4,1,3);
yyaxis left
plot(d2,viy/Vai,'*r','LineWidth',1),set(gca,'xtick',[],'FontSize',16)
%xline(jmax_pos,'--');  xline(hmin_pos,'--');
yline(0,'--');
%xlim([0 xlpos]); 
ylimi=max(abs(max(viy/Vai)), abs(min(viy/Vai)));
%ylim([-ylimi ylimi]);
set(gca,'ycolor','k') 
ylabel('$$\frac{v_{iy}}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
yyaxis right
plot(d,hy/B0, '*k','LineWidth',1),set(gca,'xtick',[],'FontSize',16)
%xline(jmax_pos,'--');  xline(hmin_pos,'--');
yline(0,'--'); 
%xlim([0 xlpos]); 
ylimh=max(abs(max(hy/B0)), abs(min(hy/B0)));
%ylim([-ylimh ylimh]);
%%yticks([-0.5 0.5]);%6 8 10]);
set(gca,'ycolor','k') 
ylabel('$$\frac{B_{y}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
%-------------------------------------------------------------------------
subplot(4,1,4);
yyaxis left
plot(d2,vix/Vai,'*r','LineWidth',1),set(gca,'xtick',[],'FontSize',16)
%xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
yline(0,'--'); 
%xline(12,'--'); 
%xlim([0 xlpos]); 
ylimi=max(abs(max(vix/Vai)), abs(min(vix/Vai)));
%ylim([-ylimi ylimi]);
%%yticks([-0.1 0.1]);
set(gca,'ycolor','k') 
ylabel('$$\frac{v_{ix}}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
yyaxis right
plot(d,hx/B0, '*k','LineWidth',1),set(gca,'xtick',[],'FontSize',16)
%xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
yline(0,'--'); 
%xlim([0 xlpos]); 
%%yticks([-0.5 0.5]); 
ylimh=max(abs(max(hx/B0)), abs(min(hx/B0)));
%ylim([-ylimh ylimh]);
set(gca,'ycolor','k') 
ylabel('$$\frac{B_{x}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
%-------------------------------------------------------------------------
ha=get(gcf,'children');
set(ha(1),'position',[.12 .74 .8 .21]); set(ha(2),'position',[.12 .53 .8 .21]);
set(ha(3),'position',[.12 .32 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
%legend({'$v_{i}$', '$B$' },'Location','northwest','Interpreter','latex')
hL = legend([line1,line2],{'$v_{i}$', '$B$'},'Interpreter','latex');
set(hL);

 
