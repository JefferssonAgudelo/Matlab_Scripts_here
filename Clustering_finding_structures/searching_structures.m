%#-----------------------------------------------------------------------
%This is an script that is not meant to be run but to go trough it running
%specific part. The ram memory is not enough to run it at once. Besides,
%the amount of output will be unmanagible.
%#----------------------------------------------------------------------

cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data'; %This is important because the xdmf files are in that directory
path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
i=19;
%for i=5:6
    disp(strcat('Computing step ...',S(i).name))
    
    % Read files 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
  
    B_1st = info.Groups(18).Name;
    V_1st = info.Groups(25).Name;
    n_1st = info.Groups(23).Name;
    J_1st = info.Groups(19).Name;
    
    Bx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    By=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    Bz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));
    
    %Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    %By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    %Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    %{
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'));
    ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    
    Bx=(vix.*ni + vex.*ne)./(ni+ne); 
    By=(viy.*ni + vey.*ne)./(ni+ne);
    Bz=(viz.*ni + vez.*ne)./(ni+ne);
    clearvars vix viy viz vex vey vez ni ne
%}
    %Bxy = sqrt(Bx.^2 + By.^2 + (Bz-0.1).^2) ;
    Bxy = sqrt(Bx.^2 + By.^2 + Bz.^2) ; % this is the one for the velocity
        
    clearvars Bx By Bz 
      
    Bxy_mean=mean(Bxy,'all');
    Bxy2_mean=mean(Bxy.*Bxy, 'all');
    Bxy_std=sqrt(Bxy2_mean - Bxy_mean^2);
    Bxyth_m1=Bxy_mean+1*Bxy_std;
    Bxyth_m2=Bxy_mean+2*Bxy_std;
    Bxyth_m3=Bxy_mean+3*Bxy_std;
    Bxyth_m4=Bxy_mean+4*Bxy_std;
    
    Bxy_above = Bxy > Bxyth_m4;
    Bxy_stats = regionprops3(Bxy_above,'all');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First to calculate and save the individual statistics of the variables along
    % the box    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Current, velocities and Temperatures
    
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    filename_1 = 'Jxyz_geometry_m4';
    save(filename_1,'Bxy_stats');
        
    %{
    % To make the statistics about the geometrical properties
    %Lstats = regionprops3(B_zero,'all'); L=bwlabel(B_zero(:,:,1)); CC = bwconncomp(B_zero); Cstats = regionprops3(CC);    
    %%% Load the file 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Contents of workspace before loading file:')
    whos
    
    disp('Bxy_geometry.mat:')
    whos('-file','Bxy.mat')
    load('Bxy_geometry.mat')
    disp('Contents of workspace after loading file:')
    whos    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    
    %----------------------------------------------------------------------
    % This is how to filter and delete rows in Table and matrix
    %yourtable(yourtable.column2 == 0, :) = [];
    %yourmatrix(yourmatrix(:, 2) == 0, :) = [];
    
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    %Bxy_stats=load('Bxyz_geometry_m2.mat');    
    Bxy_stats=load('Bxyz_geometry_m2_initial.mat');
    Bxy_stats=Bxy_stats.Bxy_stats;
 
    %Bxy_stats=J2Vie2Tie_stats;
    %Bxy_stats=JmViemTie4_stats;
    %Bxy_stats2=Bxy_stats;
    
    Bxy_stats(Bxy_stats.Volume < 16*16*16, :) = [];   
    Bxy_pal=Bxy_stats.PrincipalAxisLength;
    Bxy_angle=Bxy_stats.Orientation;
    Bxy_positions=0.06*Bxy_stats.Centroid;
    
    Bxy_image=Bxy_stats.Image;
    Bxy_image1=cell2mat(Bxy_stats.Image(1));
    Bxy_image2=cell2mat(Bxy_stats.Image(2));
    
 %{
    f11=figure(11);
    img = permute(Bxy_image1, [2 1 3]);
    img = flip(img, 2); % former x-axis
    img = flip(img, 3); 
    volshow(img)
    
    f12=figure(12);
    volshow(Bxy_image)
  %}  
   
    %Epar_positions=Bxy_positions;
    %Epar_axisdi=0.06*Bxy_pal;
    
    Asp_ratio=Bxy_pal(:,1)./(sqrt(Bxy_pal(:,2).^2 + Bxy_pal(:,3).^2));
    rpar=Bxy_pal(:,1)*0.06;
    Diameter=sqrt(Bxy_pal(:,2).^2 + Bxy_pal(:,3).^2)*0.06;
    
    Asp_mean=mean(Asp_ratio); Asp_std=std(Asp_ratio);
    Asp_median = median(Asp_ratio);
    rpar_mean=mean(rpar); rpar_std=std(rpar);
    rpar_median = median(rpar);
    Diameter_mean=mean(Diameter); Diameter_std=std(Diameter);
    Diameter_median = median(Diameter);
    
    Asp_ratio_1=Asp_ratio;
    rpar_1 = rpar;
    Diameter_1 =  Diameter;
    Asp_mean_1 = Asp_mean; Asp_std_1 = Asp_std; Asp_median_1 = Asp_median;
    rpar_mean_1 = rpar_mean; rpar_std_1 = rpar_std; rpar_median_1 = rpar_median;
    Diameter_mean_1 = Diameter_mean; Diameter_std_1 = Diameter_std; Diameter_median_1 =Diameter_median;
    
    %Asp_ratio_2=Asp_ratio;
    %rpar_2 = rpar;
    %Diameter_2 =  Diameter;
    %Asp_mean_2 = Asp_mean; Asp_std_2 = Asp_std; Asp_median_2 = Asp_median;
    %rpar_mean_2 = rpar_mean; rpar_std_2 = rpar_std; rpar_median_2 = rpar_median;
    %Diameter_mean_2 = Diameter_mean; Diameter_std_2 = Diameter_std; Diameter_median_2 =Diameter_median;
    
    % 16*16*16
    Asp_ratio_3=Asp_ratio;
    rpar_3 = rpar;
    Diameter_3 =  Diameter;
    Asp_mean_3 = Asp_mean; Asp_std_3 = Asp_std; Asp_median_3 = Asp_median;
    rpar_mean_3 = rpar_mean; rpar_std_3 = rpar_std; rpar_median_3 = rpar_median;
    Diameter_mean_3 = Diameter_mean; Diameter_std_3 = Diameter_std; Diameter_median_3 =Diameter_median;
    
    
    
    %Doing the plots
    %------------------------------------------------------------------
    %{
    f1=figure(1);
    %loglog(Diameter*0.06,'*k')
    scatter3(J2Vie2Tie_positions(:,1),J2Vie2Tie_positions(:,2),J2Vie2Tie_positions(:,3))
    set(gca,'FontSize',18)
    hold on
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('3D Positions')
    hold off
    
    figure1=figure(1);
    loglog(Bxy_pal(:,1)*0.06,'*k')
    set(gca,'FontSize',18)
    hold on
    loglog(Bxy_pal(:,2)*0.06,'*r')
    loglog(Bxy_pal(:,3)*0.06,'*b')
    title('principal axis lengths')
    hold off
    
    figure2=figure(2);
    plot(Bxy_angle(:,1),'*k')
    set(gca,'FontSize',18)
    hold on
    plot(Bxy_angle(:,2),'*r')
    plot(Bxy_angle(:,3),'*b')
    title('Angles')
    hold off
    %}
    
    %{
    edges1 = linspace(0,40,41);
    f3=figure(3);
    subplot(3,2,1)
    h1=histogram(Bxy_pal(:,1)*0.06, edges1);
    h1.FaceColor='k';
    set(gca,'YScale','log','FontSize',16)
    title('$Axis_{1} \ \omega_{pi}t=120$','Interpreter','latex','FontSize',16)
    xlabel('$r_{1}/d_i$','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex')
    subplot(3,2,2)
    h2=histogram(Bxy_pal(:,2)*0.06, edges1);
    h2.FaceColor='r';
    set(gca,'YScale','log','FontSize',16)
    title('$Axis_{2}$','Interpreter','latex','FontSize',16)
    xlabel('$r_{2}/d_i$','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex')
    subplot(3,2,3)
    h3=histogram(Bxy_pal(:,3)*0.06, edges1);
    h3.FaceColor='b';
    set(gca,'YScale','log','FontSize',16)
    title('$Axis_{3}$','Interpreter','latex','FontSize',16)
    xlabel('$r_{3}/d_i$','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex')
    subplot(3,2,4)
    h4=histogram(Diameter*0.06, edges1);
    h4.FaceColor='g';
    set(gca,'YScale','log','FontSize',16)
    title('$Diameter$','Interpreter','latex','FontSize',16)
    xlabel('$r_{\parallel}/r_{\perp}$','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex')
    subplot(3,2,5)
    h5=histogram(Asp_ratio, edges1);
    h5.FaceColor='m';
    set(gca,'YScale','log','FontSize',16)
    title('$Aspec \ ratio$','Interpreter','latex','FontSize',16)
    xlabel('$r_{\parallel}/r_{\perp}$','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex')
    %}
    
    % This is the plot for the first paper about LD, Lpar, ratio
    % distribution
    %----------------------------------------------------------------------
    Asp_ratio_i=Asp_ratio_3;
    Diameter_i=Diameter_3;
    rpar_i=rpar_3;
    
    edges1 = linspace(0,40,41);
    edges2 = linspace(0,20,21);
    f31=figure(31);
    subplot(3,1,1)
    %h1=histogram(Bxy_pal(:,1)*0.06, edges1);
    h1=histogram(rpar_i, edges1);
    h1.Normalization='pdf';
    h1.FaceColor='r';
    set(gca,'YScale','log','FontSize',18)
    %title('$Axis_{1} \ \omega_{pi}t=0$','Interpreter','latex','FontSize',18)
    xlabel('$L_{\parallel}/d_{i}$','Interpreter','latex')
    ylabel('$PDF$','Interpreter','latex')
    ylim([0.001 1])
    subplot(3,1,2)
    %h4=histogram(Diameter*0.06, edges2);
    h4=histogram(Diameter_i, edges2);
    h4.Normalization='pdf';
    h4.FaceColor='b';
    set(gca,'YScale','log','FontSize',18)
    %title('$Diameter$','Interpreter','latex','FontSize',16)
    xlabel('$L_{D}/d_{i}$','Interpreter','latex')
    ylabel('$PDF$','Interpreter','latex')
    ylim([0.001 1])
    xlim([1 10])
    subplot(3,1,3)
    h5=histogram(Asp_ratio_i, edges1);
    h5.Normalization='pdf';
    h5.FaceColor='k';
    set(gca,'YScale','log','FontSize',18)
    %title('$Aspec \ ratio$','Interpreter','latex','FontSize',16)
    xlabel('$L_{\parallel}/L_{D}$','Interpreter','latex')
    ylabel('$PDF$','Interpreter','latex')
    ylim([0.001 1])
    xlim([1 20])
    
    cd '/Volumes/PSC_DiRAC_DATA';
    saveas(f31,'LpLd_asp_Bxyz_2000_log_no_vol.png');
    %----------------------------------------------------------------------
    
    %Comparing using the z-test
    %----------------------------------------------------------------------   
    Z_asp = (Asp_mean_1 - Asp_mean_2) / sqrt(Asp_std_2^2 + Asp_std_1^2) ;
    Z_rpar = (rpar_mean_1 - rpar_mean_2) / sqrt(rpar_std_2^2 + rpar_std_1^2);
    Z_Diameter = (Diameter_mean_1 - Diameter_mean_2) / sqrt(Diameter_std_2^2 + Diameter_std_1^2) ;

    %Z_rpar = 1.2623;
    %rpar_mean_1 = 16.3333; rpar_mean_2 = 2.1603; 
    %rpar_median_1 = 14.4516; rpar_median_2 = 0.6429;
    %rpar_std_1 = 8.3213; rpar_std_2 = 5.0804; 
    
    %Z_Diameter = 0.6421;
    %Diameter_mean_1 = 1.5515; Diameter_mean_2 = 0.6242; 
    %Diameter_median_1 = 1.3580; Diameter_median_2 = 0.3664;
    %Diameter_std_1 = 0.9497; Diameter_std_2 = 0.7157;  
    
    %Z_asp = 1.1193;
    %Asp_mean_1 = 11.0065; Asp_mean_2 = 2.5451; 
    %Asp_median_1 = 10.6259; Asp_median_2 = 1.9756;
    %Asp_std_1 = 7.0649; Asp_std_2 = 1.9441; 
    

    Z_asp_31 = (Asp_mean_1 - Asp_mean_3) / sqrt(Asp_std_3^2 + Asp_std_1^2) 
    Z_rpar_31 = (rpar_mean_1 - rpar_mean_3) / sqrt(rpar_std_3^2 + rpar_std_1^2)
    Z_Diameter_31 = (Diameter_mean_1 - Diameter_mean_3) / sqrt(Diameter_std_3^2 + Diameter_std_1^2) 
    
    %Z_rpar_31 = -0.0058;
    %rpar_mean_3 = 14.9729; 
    %rpar_median_3 = 12.4015;
    %rpar_std_3 = 9.0059;  
    
    %Z_Diameter_31 = -0.6485;
    %Diameter_mean_3 = 3.1356; 
    %Diameter_median_3 = 2.4501; 
    %Diameter_std_3 = 2.2504; 
    
    %Z_asp_31 = 0.7410
    %Asp_mean_3 = 5.4593; 
    %Asp_median_3 = 4.9743; 
    %Asp_std_3 = 2.4759; 
    
    
    
 
    %----------------------------------------------------------------------
    
    f4=figure(4);
    subplot(2,2,1)
    h1=histogram(Bxy_angle(:,1),80);
    h1.FaceColor='k';
    set(gca,'YScale','log','FontSize',16)
    title('$Angle_{1} \ \omega_{pi}t =120$','Interpreter','latex','FontSize',16)
    xlabel('$\theta_{1}^{\circ} $','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex')
    subplot(2,2,2)
    h2=histogram(Bxy_angle(:,2),80);
    h2.FaceColor='r';
    set(gca,'YScale','log','FontSize',16)
    title('$Angle_{2}$','Interpreter','latex','FontSize',16)
    xlabel('$\theta_{2}^{\circ} $','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex','FontSize',16)
    subplot(2,2,3)
    h3=histogram(Bxy_angle(:,3),80);
    h3.FaceColor='b';
    set(gca,'YScale','log','FontSize',16)
    title('$Angle_{3}$','Interpreter','latex','FontSize',16)
    xlabel('$\theta_{3}^{\circ} $','Interpreter','latex','FontSize',16)
    ylabel('$Counts \ J_{xyz}$','Interpreter','latex','FontSize',16)
 
    
    
    %---------------------------------------------------------
    cd '/Volumes/PSC_DiRAC_DATA';
    %-------------------------------------------------------------------------
    saveas(f3,'axis_Jxyz4_initial.png');
    saveas(f31,'LpLd_asp_Jxyz4_120.png');
    saveas(f4,'angles_Jxyz4_initial.png');
    
   
 
    
 
    
  

 
    %----------------------------------------------------------------------
    %cd /disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics/2D_1D_PSD
    cd '/Volumes/PSC_DiRAC_DATA/tmp_images';
    % Save the plots
    %-------------------------------------------------------------------------
    saveas(f1,strcat(S(i).name,'2Dspectrum_Bpd_CB04.png'));


     
    %cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'; 
    cd '/Volumes/PSC_DiRAC_DATA/DATACB104';
    
end




% Sample data
[latGrid,lonGrid] = meshgrid(25:45,125:145);
T = table(latGrid(:),lonGrid(:),randi([0,35],size(latGrid(:))),...
  'VariableNames',{'lat','lon','temp'});
% Discretize with 5 degree
[latGroupID, latEdge] = discretize(T.lat, 25:5:45);
[lonGroupID, lonEdge] = discretize(T.lon, 125:5:145);
% Add 5 degree mesh ID to the table T
[C,~,ic] = unique([latGroupID,lonGroupID],'rows');
T.meshID = ic;
% Calculate mean temperature in 5 degree mesh
Output = table(latEdge(C(:,1))',lonEdge(C(:,2))',...
  splitapply(@mean,T.temp,T.meshID),...
  'VariableNames',{'latEdge','lonEdge','aveTemp'});
% Visualize the output
heatmap(Output,'latEdge','lonEdge','ColorVariable','aveTemp');


[x,y,z] = meshgrid(1:50,1:50,1:50);

bw1xy= sqrt((x-10).^2 / 16 + (y-15).^2 / 9 ) < 2;
bw2xy = sqrt((x-10).^2 + (y-30).^2) < 5;
theta3=-pi/4;
bw3xy = sqrt( ((x-40)*cos(theta3) + (y-40)*sin(theta3)).^2 / 9 + ...
    ((x-40)*sin(theta3) - (y-40)*cos(theta3)).^2 / 20) < 2;
theta4=pi/4;
bw4xy = sqrt( ((x-40)*cos(theta4) + (y-10)*sin(theta4)).^2 / 9 + ...
    ((x-40)*sin(theta4) - (y-10)*cos(theta4)).^2 / 20) < 1;
bwxy = bw1xy | bw2xy | bw3xy | bw4xy;
f67=figure(67);
pcolor(bwxy(:,:,2))

s = regionprops3(bwxy,"Centroid","PrincipalAxisLength","orientation");
centers = s.Centroid;
s.PrincipalAxisLength
s.Orientation


% 3D ellipsoids
%angles of inclination
a1=0*pi/4; b1=pi/2;
a2=0*3*pi/4; b2=pi/4;
Sa1=sin(a1); Sa2=sin(a2); Ca1=cos(a1) ; Ca2=cos(a2);
Sb1=sin(b1); Sb2=sin(b2); Cb1=cos(b1) ; Cb2=cos(b2);
%centres
h1=10; j1=20; k1=30; 
h2=35; j2=35; k2=30;
%variables
s1=(x-h1)*Ca1*Cb1 + (y-j1)*Sa1*Cb1 + (z-k1)*Sb1;
t1=(y-j1)*Ca1 - (x-h1)*Sa1;
u1=(z-k1)*Cb1 - (y-j1)*Sa1*Sb1 - (x-h1)*Ca1*Sb1;

s2=(x-h2)*Ca2*Cb2 + (y-j2)*Sa2*Cb2 + (z-k2)*Sb2;
t2=(y-j2)*Ca2 - (x-h2)*Sa2;
u2=(z-k2)*Cb2 - (y-j2)*Sa2*Sb2 - (x-h2)*Ca2*Sb2;

% Equations of the ellipsoids
bw1 = sqrt((s1).^2 / 9  + (t1).^2 / 9 + (u1).^2 / 49) < 1;
bw2 = sqrt((s2).^2 / 9  + (t2).^2 / 9 + (u2).^2 /49) < 1;
bw = bw1 | bw2; % this seems to joing the two binary files

f68=figure(68);
pcolor(squeeze(bw(1:45,:,30))) % This plot x in the vertical axis


s3d = regionprops3(bw,"Centroid","PrincipalAxisLength","orientation",'EigenVectors','EigenValues');
centers3d = s3d.Centroid;
s3d.PrincipalAxisLength
s3d.Orientation
evectors=s3d.EigenVectors;
evalues=s3d.EigenValues;

%--------------------------------------------------------------------------

%This is to try to stimate the scalling law.
%--------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    
Bxy_stats2=load('Bxyz_geometry_m2.mat');  
Bxy_stats3=load('Bxyz_geometry_m2.mat');  

%Bxy_stats2=load('Bxyz_geometry_m2_initial.mat');   
Bxy_stats2=Bxy_stats2.Bxy_stats;

Bxy_stats3=Bxy_stats3.Bxy_stats;

%Bxy_stats2(Bxy_stats2.Volume < 16*16*16, :) = []; 
Bxy_stats2(Bxy_stats2.Volume < 8, :) = [];
%Bxy_stats2(Bxy_stats2.Volume > 16*16*16, :) = []; %to filter one more time
%removing the V>1di
Bxy_pal2=Bxy_stats2.PrincipalAxisLength;
  
rpar2=Bxy_pal2(:,1)*0.06;
rper1=Bxy_pal2(:,2)*0.06;
rper2=Bxy_pal2(:,3)*0.06;
rper12=sqrt(rper1.*rper1 + rper2.*rper2);
ratio2=rpar2./rper2;

%Bxy_stats2(Bxy_stats2.Volume < 16*16*16, :) = []; 
Bxy_stats3(Bxy_stats3.Volume < 16*16*16, :) = [];
Bxy_pal3=Bxy_stats3.PrincipalAxisLength;
rpar11=Bxy_pal3(:,1)*0.06;
rper11_1=Bxy_pal3(:,2)*0.06;
rper11_2=Bxy_pal3(:,3)*0.06;
rper112=sqrt(rper11_1.*rper11_1 + rper11_2.*rper11_2);



%{
logpar=log10(rpar2);logper=log10(rper12);
logpar1=logpar((0.08<logper) & (logper<0.6));
logper1=logper((0.08<logper) & (logper<0.6));
logpar2=logpar((0.8<logper) & (logper<1));
logper2=logper((0.8<logper) & (logper<1));
logpar3=logpar((0.8 < logper) & (logper < 2));
logper3=logper((0.8 < logper) & (logper < 2));
figure666=figure(666);
scatter(rper1,rpar2)
figure668=figure(668);
scatter(logper2,logpar2)
hold on
%scatter(logper2,logpar2)
scatter(logper3,logpar3)
hold off
%}

%dstep1=0.1; dstep2=0.46; dstep3=max(rper12)*1.2;
dstep1=0.1; dstep2=0.79; dstep3=max(rper12)*1.2;

%time 0
%rpar2_1=rpar2((dstep1<=rper12) & (rper12<=dstep3)); %2
%rper2_1=rper12((dstep1<=rper12) & (rper12<=dstep3));

%time 120
rpar2_1=rpar2((dstep1<=rper12) & (rper12<=dstep2)); %2
rper2_1=rper12((dstep1<=rper12) & (rper12<=dstep2));

rpar2_111=rpar11((dstep1<=rper12) & (rper12<=dstep2)); %2 This is ot working
rper2_111=rper11((dstep1<=rper12) & (rper12<=dstep2));

rpar2_2=rpar2((dstep2<rper112) & (rper112<=dstep3));
rper2_2=rper12((dstep2<rper112) & (rper112<=dstep3));

logpar2=log(rpar2_2);
logper2=log(rper2_2);

f123=figure(123);
scatter(log10(rper2_1), log10(rpar2_1),'b')
%scatter(log10(rper12), log10(rpar2),'b')
hold on 
%scatter(log10(rper112), log10(rpar11),'+')
%scatter(log10(rper12), log10(rpar2),'k')
%scatter(log10(rper2_2), log10(rpar2_2),'r')
%scatter(log10(rper12), log10(rpar2),'+')
hold off


f1=figure(1);
%scatter(rper12,rpar2,15,'+') %time 0 
%scatter(rper2_1,rpar2_1,10,'f') %time 120
h1=scatter(rper12,rpar2,15,'f') ;%time 120
hold on
%scatter(rper2_2,rpar2_2,10,'f')
h2=scatter(rper112,rpar11,15,'f');
set(gca,'XScale','log','YScale','log','FontSize',20)
xlabel('$L_{D}/d_{i}$','Interpreter','latex')
ylabel('$L_{\parallel}/d_{i}$','Interpreter','latex')
x=log10([dstep1 dstep2]); y=(0.67)*x + 0.3;%-0.8;
x2=log10([dstep2 dstep3]); y2=(1.2)*x2 + 0.15 ;
x3=log10([dstep1 dstep3]); y3=(0.67)*x3 + 0.81; %1.1 ;
x4=log10([dstep1 dstep3]); y4=x4 +0.28;
%x4=log10([dstep1 dstep3]); y4=x4 +1; %initial
%x4=log10([dstep1 dstep3]); y4=x4 +0.02;
%loglog(10.^x,10.^y,'--k','LineWidth',2)
%loglog(10.^x2,10.^y2,'--k','LineWidth',2)
h3=loglog(10.^x3,10.^y3,'--k','LineWidth',2);
h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
h4=loglog(10.^x4,10.^y4,'--r','LineWidth',2);
h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend([h1 h2],{'$V \leq d_{i}^3$' '$V > d_{i}^3$'},'Interpreter','latex')
%title('$\log_{10}{P^{\mathbf{B}}_{2D}}$','Interpreter','latex')
xlim([dstep1 dstep3])
%xlim([0.1 100])
ylim([0.1 100])
txt1 = '$L_{\parallel} \sim L_{D}$';
txt2 = '$L_{\parallel} \sim L_{D}^{0.7}$';
txt3 = '$L_{\parallel} \sim L_{D}^{0.66}$';
text(0.3,0.15,txt1,'Interpreter','latex','FontSize',20)
text(1,50,txt2,'Interpreter','latex','FontSize',20)
%text(1,50,txt3,'Interpreter','latex','FontSize',20)
hold off

cd '/Volumes/PSC_DiRAC_DATA';
%-------------------------------------------------------------------------
saveas(f1,'L_par_per_scaling_120_B255_2.png');


f22=figure(22);
%scatter(rper2_1,rpar2_1,'f') %time 0 
scatter(1./rpar2_1,1./rper2_1,'f') %time 120
hold on
scatter(1./rpar2_2,1./rper2_2,'f') 
set(gca,'XScale','log','YScale','log','FontSize',20)
xlabel('$d_{i}/L_{\parallel}$','Interpreter','latex')
ylabel('$d_{i}/L_{D}$','Interpreter','latex')
x=log10([dstep1 dstep2]); y=(0.67)*x ;%-0.8;
x2=log10([dstep2 dstep3]); y2=(1.5)*x2 + 0.3 ;
x3=log10([dstep1 dstep3]); y3=(0.67)*x3 + 1 ;
x4=log10([dstep1 dstep3]); y4=(1.5)*x4 + 0.7 ;
loglog(10.^x,10.^y,'--k','LineWidth',2)
loglog(10.^x2,10.^y2,'--k','LineWidth',2)
loglog(10.^x4,10.^y4,'--r','LineWidth',2)
%loglog(10.^x3,10.^y3,'--k','LineWidth',2)
%title('$\log_{10}{P^{\mathbf{B}}_{2D}}$','Interpreter','latex')
xlim([dstep1 dstep3])
ylim([0.1 100])
txt1 = '$L_{\parallel} \sim L_{D}^{0.67}$';
txt2 = '$L_{\parallel} \sim L_{D}^{1.5}$';
txt3 = '$L_{\parallel} \sim L_{D}^{0.67}$';
text(0.3,0.15,txt1,'Interpreter','latex','FontSize',20)
text(2,2,txt2,'Interpreter','latex','FontSize',20)
%text(2,4,txt3,'Interpreter','latex','FontSize',20)
hold off