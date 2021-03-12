% Read CDF data from Solar Orbiter and PSP

% Solar orbiter
%-------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/DATA_SolO_PSP/SolO/20200728';
path_so = '/Volumes/PSC_DiRAC_DATA/DATA_SolO_PSP/SolO/20200728';

S_so = dir(fullfile(path_so,'*.cdf'));
N_so=numel(S_so); % number of files to use
H_so=zeros(N_so);
i=15;
disp(strcat('Computing step ...',S_so(i).name))
% Read files  
%{
1. solo_L1_swa-eas1-NM3D
2. solo_L1_swa-eas1-strahlc
3. solo_L1_swa-eas-2Dbusrtc
4. solo_L1_swa-eas-OnbPartMoms
5. solo_L1_swa-pas-3d
6. solo_L1_swa-pas-cal
7. solo_L1_swa-pas-hsk
8. solo_L1_swa-pas-mom
9. solo_L2_swa-pas-eflux
10. solo_L2_swa-pas-grnd-mom
11. solo_L2_swa-pas-vdf
12. solo_L2a_swa-eas-2DBurstc
13. solo_L2a_swa-eas1-NM3D-psd
14. solo_L2a_swa-eas1-strahlc
15. solo_LL02_mag
%}

fileID_MAG_so =  S_so(1).name; %change the 1 per i
%Solo_eas1_NM3D = cdfread(S_so(1).name); %This works with 1,2,3,4,12,13,14 
addpath '/Applications/matlab_cdf380_patch-64'

Solo_MAG=spdfcdfread(S_so(1).name);
Info_MAG_so=spdfcdfinfo(fileID_MAG_so);
so_mag_epoch=Solo_MAG{1};
so_B=Solo_MAG{2};
% Normalization units to SI
nT2T=1e-9;
so_B = so_B*nT2T;
so_Bx=so_B(:,1);
so_By=so_B(:,2);
so_Bz=so_B(:,3);
so_Bmag=sqrt(so_Bx.*so_Bx + so_By.*so_By + so_Bz.*so_Bz);
%-------------------------------------------------------------------------


%Solo SWA PAS
%--------------------------------------------------------------------------
fileID_swa_pas_gm =  S_so(6).name;
swa_pas_gm=spdfcdfread(S_so(6).name);
Info_swa_pas_gm=spdfcdfinfo(fileID_swa_pas_gm);
swa_pas_epo=swa_pas_gm{1};
swa_pas_gm_v_nrt=swa_pas_gm{8}; %NRT frame
swa_pas_gm_p_nrt=swa_pas_gm{10}; % Energy Density
% Normalization units to SI
km2m=1000;
so_vnrt = swa_pas_gm_v_nrt*km2m;
so_vix=so_vnrt(:,1);
so_viy=so_vnrt(:,2);
so_viz=so_vnrt(:,3);
so_vimag=sqrt(so_vix.*so_vix + so_viy.*so_viy + so_viz.*so_viz);

swa_pas_gm_n=swa_pas_gm{6};


%Cheking the epoch to see where they match and to define an interval
%-------------------------------------------------------------------------
tmag = datetime(so_mag_epoch,'ConvertFrom','datenum');
tswa = datetime(swa_pas_epo,'ConvertFrom','datenum');

% the resolution of the moments calculated is too small.
%4 points
%8045 {'28-Jul-2020 06:34:05'}
%8047 {'28-Jul-2020 06:34:13'}
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
fr=200;
n=1;
dm=so_mag_epoch(fr*n:fr*(n+1)); %d(isinf(d)|isnan(d)) = 0; 
hx = so_Bx(fr*n:fr*(n+1)); %hx(isinf(hx)|isnan(hx)) = 0;
hy = so_By(fr*n:fr*(n+1)); %hy(isinf(hy)|isnan(hy)) = 0;
hz = so_Bz(fr*n:fr*(n+1)); %hz(isinf(hz)|isnan(hz)) = 0;
hm=sqrt(hx.*hx + hy.*hy + hz.*hz);  %psp_Bmag(fr*n:fr*(n+1)); %hm(isinf(hm)|isnan(hm)) = 0;
B0=mean(hm);

%----------------------------------------------------------
frs=1; ns=8045; dns=4; % this is chosen by hand to match mag data with pas data

ds=swa_pas_epo(frs*ns:frs*(ns+dns)); %d(isinf(d)|isnan(d)) = 0; 
vix = so_vix(frs*ns:frs*(ns+dns)); %hx(isinf(hx)|isnan(hx)) = 0;
viy = so_viy(frs*ns:frs*(ns+dns)); %hy(isinf(hy)|isnan(hy)) = 0;
viz = so_viz(frs*ns:frs*(ns+dns)); %hz(isinf(hz)|isnan(hz)) = 0;
vim=sqrt(vix.*vix + viy.*viy + viz.*viz);  %psp_Bmag(fr*n:fr*(n+1)); %hm(isinf(hm)|isnan(hm)) = 0;
vimean=mean(vim);
vsw=vimean;

cm_32m_3=(0.01)^3;
ni=cm_32m_3*swa_pas_gm_n(frs*ns:frs*(ns+dns));
nim=mean(ni);

R=0.7;

mu0=4*pi*1e-7;
ep0=8.854e-12;
q=1.602e-19;
mi=1.67e-27;
c=299792458; %m/s
wpi=sqrt( (nim*q*q) / (mi*ep0) );
di=c/wpi;

Time10di = (10 * di) / (vsw); % in seconds
Time10di_d = Time10di/(60*60*24); 

Via=B0/sqrt(mu0 * mi * nim);

deltat=seconds(tswa(ns+dns)-tswa(ns));
frac_di = (deltat*vsw) / di;

%-----------------------------------------------------------------

% Plots magnetic field
%B0=1;
dateFormat='HH:MM:SS';
%dateFormat='ss';
xlow = dm(1);
xup = dm(fr);
f21=figure(21); %string(round(frac_di,3))
%-------------------------------------------------------------------------
subplot(4,1,1),
plot(dm, hm/B0, '-k','LineWidth',2);set(gca,'FontSize',16)
set(gca,'ycolor','k') 
ylabel('$$\frac{|B|}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
xlabel('$$t(s) (\ r \sim' + string(frac_di) + 'd_{i})$$','Interpreter','latex','FontSize',18)
xlim([xlow xup]);
%yt1 = 1*max(hm/B0)/3; yt2 = 2*max(hm/B0)/3; 
%yticks([1]);
datetick('x',dateFormat,'keepticks')
%datetick('x',dateFormat)
%rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.4])
%-------------------------------------------------------------------------
subplot(4,1,2);
plot(dm,hz/B0, '-k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k')
%yticks([-0.75 -0.72]);
%yline(0,'--'); 
ylabel('$$\frac{B_{z}}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%-------------------------------------------------------------------------
subplot(4,1,3);
plot(dm,hy/B0, '-k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k') 
%yticks([-0.48 -0.42]);
%yline(0,'--'); 
ylabel('$$\frac{B_{y}}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%-------------------------------------------------------------------------
subplot(4,1,4);
plot(dm,hx/B0, '-k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k') 
%yline(0,'--'); 
%ylabel('$$\frac{B_{x}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
ylabel('$$\frac{B_{x}}{\langle |B| \rangle}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%yticks([0.48 0.52]);
title('$$20200728;\ \langle v_{sw} \rangle =' + string(round(vsw/1000,2))+'km/s;\ \langle |B| \rangle ='+ string(round(B0*10^(9),2)) + ' nT;\ R \approx' +string(R) +'AU $$','Interpreter','latex','FontSize',18)
%-------------------------------------------------------------------------
ha=get(gcf,'children');
set(ha(1),'position',[.12 .74 .8 .21]); set(ha(2),'position',[.12 .53 .8 .21]);
set(ha(3),'position',[.12 .32 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);

%cd '/Volumes/PSC_DiRAC_DATA/';
%saveas(f21,'Bxyzm_20200924_SO.png');
%--------------------------------------------------------------------------


%Plots moments velocity
dateFormat='HH:MM:SS';
%dateFormat='ss';
xlow = ds(1);
xup = ds(dns);
f22=figure(22);
%-------------------------------------------------------------------------
subplot(4,1,1),
plot(ds, vim/Via, '-r','LineWidth',2);set(gca,'FontSize',16)
set(gca,'ycolor','k') 
ylabel('$$\frac{|v_{i}|}{v_{Ai}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
xlabel('$$t(s) (\ r \sim 60d_{i})$$','Interpreter','latex','FontSize',18)
xlim([xlow xup]);
%yt1 = 1*max(hm/B0)/3; yt2 = 2*max(hm/B0)/3; 
%yticks([1]);
datetick('x',dateFormat,'keepticks')
%datetick('x',dateFormat)
%rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.4])
%-------------------------------------------------------------------------
subplot(4,1,2);
plot(ds,viz/Via, '-r','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k')
%yticks([-0.75 -0.72]);
%yline(0,'--'); 
ylabel('$$\frac{v_{iz}}{v_{Ai}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%-------------------------------------------------------------------------
subplot(4,1,3);
plot(ds,viy/Via, '-r','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k') 
%yticks([-0.48 -0.42]);
%yline(0,'--'); 
ylabel('$$\frac{v_{iy}}{v_{Ai}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%-------------------------------------------------------------------------
subplot(4,1,4);
plot(ds,vix/Via, '-r','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
set(gca,'ycolor','k') 
%yline(0,'--'); 
%ylabel('$$\frac{B_{x}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
ylabel('$$\frac{v_{ix}}{v_{Ai}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
%rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
xlim([xlow xup]);
%yticks([0.48 0.52]);
title('$$20200728;\ \langle v_{sw} \rangle =450 km/s;\ \langle |B| \rangle =8.3 nT;\ R \approx 1.0 AU $$','Interpreter','latex','FontSize',18)
%-------------------------------------------------------------------------
ha=get(gcf,'children');
set(ha(1),'position',[.12 .74 .8 .21]); set(ha(2),'position',[.12 .53 .8 .21]);
set(ha(3),'position',[.12 .32 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);

%-------------------------------------------------------------------------





%-------------------------------------------------------------------------
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


