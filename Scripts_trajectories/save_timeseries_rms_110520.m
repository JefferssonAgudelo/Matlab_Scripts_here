%#-----------------------------------------------------------------------
%To calculate the rms values and store in a single file
%#----------------------------------------------------------------------
%diary myDiaryFile
%cd '/disk/plasma2/jaa/CB8WAVES/CB8waves_04'; %Set the directory of the files 
%path = '/disk/plasma2/jaa/CB8WAVES/CB8waves_04';

cd '/Volumes/PSC_DiRAC_DATA/DATACB104'; %This is important because the xdmf files are in that directory
path = '/Volumes/PSC_DiRAC_DATA/DATACB104';

S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
mag_rms=zeros(N,6);
i=1;
for i=16:33
    disp(strcat('Computing step ...',S(i).name)) %6 is 2000
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    %T_1st = info.Groups(1).Name;
    %E_1st = info.Groups(17).Name;
    B_1st = info.Groups(18).Name;
    %J_1st = info.Groups(19).Name;
    V_1st = info.Groups(25).Name;
    %n_1st = info.Groups(23).Name;
    %diary off
    %s2='/hx/p0/3d';
    %sx = strcat(B_1st,s2);
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    %Ex=h5read(fileID,strcat(E_1st,'/ex/p0/3d'));
    %Ey=h5read(fileID,strcat(E_1st,'/ey/p0/3d'));
    %Ez=h5read(fileID,strcat(E_1st,'/ez/p0/3d'));
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    %vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    %vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    %vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'));
    %ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    %ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    %Jx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    %Jy=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    %Jz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));
    
    %------------------------------------------------------------------
    % Calculate the magnitudes and rms
    Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2); Bmag2=(Bx.^2 + By.^2 + Bz.^2);
    %Emag = sqrt(Ex.^2 + Ey.^2 + Ez.^2); Emag2 =(Ex.^2 + Ey.^2 + Ez.^2);
    %Jmag = sqrt(Jx.^2 + Jy.^2 + Jz.^2); Jmag2 = (Jx.^2 + Jy.^2 + Jz.^2);
    vimag = sqrt(vix.^2 + viy.^2 + viz.^2); vimag2 = (vix.^2 + viy.^2 + viz.^2);
    clearvars Bx By Bz 
    %Ex Ey Ez 
    %clearvars Jx Jy Jz
    clearvars vix viy viz
    
    Bmag = Bmag(~isnan(Bmag)); Bmag2 = Bmag2(~isnan(Bmag2));
    %Emag = Emag(~isnan(Emag)); Emag2 = Emag2(~isnan(Emag2));
    %Jmag = Jmag(~isnan(Jmag)); Jmag2 = Jmag2(~isnan(Jmag2));
    vimag = vimag(~isnan(vimag)); vimag2 = vimag2(~isnan(vimag2));
    
    Bmag_mean=mean(Bmag,'all'); Bmag_mean2=mean(Bmag2,'all'); 
    %Emag_mean=mean(Emag,'all'); Emag_mean2=mean(Emag2,'all'); 
    %Jmag_mean=mean(Jmag,'all'); Jmag_mean2=mean(Jmag2,'all'); 
    vimag_mean=mean(vimag,'all'); vimag_mean2=mean(vimag2,'all'); 
    
    clearvars Bmag Bmag2 
    %Emag Emag2 Jmag Jmag2  
    clearvars vimag vimag2
    
    B_rms= sqrt(Bmag_mean2 - (Bmag_mean)^2);
    %E_rms= sqrt(Emag_mean2 - (Emag_mean)^2);
    %J_rms= sqrt(Jmag_mean2 - (Jmag_mean)^2); 
    vi_rms= sqrt(vimag_mean2 - (vimag_mean)^2); 
    
    % Bmag, B_rms, Emag, E_rms, Jmag, J_rms,
    mag_rms(i,1)=Bmag_mean; mag_rms(i,2)=B_rms;
    %mag_rms(i,3)=Emag_mean; mag_rms(i,4)=E_rms;
    %mag_rms(i,5)=Jmag_mean; mag_rms(i,6)=J_rms;
    mag_rms(i,3)=vimag_mean; mag_rms(i,4)=vi_rms;
    
    %mag_rmsi=mag_rms;
    %mag_rms=[mag_rms mag_rmsi];
    %-------------------------------------------------------------------
    filename='mag_rms_2.mat';
    save(filename, 'mag_rms');
%    save(filename,'Bmag_mean','B_rms','Emag_mean','E_rms','Jmag_mean','J_rms','-append','-nocompression')
end


whos('-file','mag_rms_2.mat')

Bm1=magBEJ_rms(:,1);
Bm1_rms=magBEJ_rms(:,2);
Em1=magBEJ_rms(:,3);
Em1_rms=magBEJ_rms(:,4);
Jm1=magBEJ_rms(:,5);
Jm1_rms=magBEJ_rms(:,6);
Vim1=vim_rms(:,3);
Vim1_rms=vim_rms(:,4);


Bm1=BviJE2(:,1);
Bm1_rms=BviJE2(:,2);
Vim1=BviJE2(:,3);
Vim1_rms=BviJE2(:,4);
Jm1=mag_rms(:,5);
Jm1_rms=mag_rms(:,6);


BviJE2=[Bm1 Bm1_rms Vim1 Vim1_rms Jm1 Jm1_rms Em1 Em1_rms];
BviJE3=[Bm1 Bm1_rms Vim1 Vim1_rms Jm1 Jm1_rms];

filename='BviJE3.mat';
    save(filename, 'BviJE3');

    mi2me=100;
    UE=0.5*Emag.^2;
    UM=0.5*Bmag.^2;
    Ki = 0.5*vimag.^2;
    Ke = (1/mi2me)*0.5*vemag.^2;
    Kith = 0.5*vithmag.^2;
    Keth = (1/mi2me)*0.5*vethmag.^2;
    Utotal = UE + UM +Ki + Ke + Kith+ Keth;
    Bmag2 = Bx.^2 +By.^2 +Bz.^2;
    Bmag2_2 = Bmag2(:);
    Brms_1 = sqrt(mean(Bmag2_2));
    
  %To plot rms values and also make large the text
    %{ 
    f111=figure(111);
    plot(Rmsvalues.t_pi,Rmsvalues.B_rms);
    hold on;
    plot(Rmsvalues.t_pi,Rmsvalues.E_rms);
    plot(Rmsvalues.t_pi,Rmsvalues.J_rms);
    xlabel('$t \omega_{pi}$','Interpreter','latex')
    %ylabel('$P1D_\|$','Interpreter','latex')
    legend({'$B_{rms}$','$E_{rms}$','$J_{rms}$'},'Interpreter','latex')
    set(gca,'FontSize',15)
    hold off;    
    saveas(f111,'RMS_values.png');
    %}
    
  
    %aa=load('mag_rms.txt', '-ascii');

    
    %----------------------------------------------------------------------
    %cd /disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics/2D_1D_PSD
    %cd '/Volumes/PSC_DiRAC_DATA/tmp_images';
    % Save the plots
    %-------------------------------------------------------------------------
    %saveas(f1,strcat(S(i).name,'2Dspectrum_Bpd_CB04_not3.png'));
    %saveas(f5,strcat(S(i).name,'1_perDspectrum_CB04_Bvipd_2.png'));
    %saveas(f6,strcat(S(i).name,'1_parDspectrum_CB04_Bvipd_2.png'));
    %saveas(f7,strcat(S(i).name,'1_parDspectrum_CB04_Bparper.png'));    
    %cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'; 
    %cd '/Volumes/PSC_DiRAC_DATA/DATACB104';
