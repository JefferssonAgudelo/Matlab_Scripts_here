%#-----------------------------------------------------------------------
% This file is to make the plots related with the calculation of the
% spectrum
% This is the file that I use to checke the different runs changing
%parameters
%#----------------------------------------------------------------------

% Save the workspace in this directory 
cd /Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis;

%#----------------------------------------------------------------------
% This is to generate the spectra. One this has been done. There is no need
% to do it again
%#----------------------------------------------------------------------
%predir = 'DATACB103_1'; % DATACB103_1 
predir ='DATACB103_thin_8';
%predir = 'DATACB103_2L';
prepath = '/Volumes/PSC_DiRAC_DATA/'; %changing_ppc/, changing_mime/ 

path = strcat(prepath,predir);
cd(path) 

S = dir(fullfile(path,'*_p000000.h5'));
N=numel(S); % number of files to use
H=zeros(N);
%i=18;
i=1;
%for i=1:3
    disp(strcat('Computing step ...',S(i).name)) %6 is 2000   
    fileID =  S(i).name; %change the 1 per i
    h5disp(fileID);
    info = h5info(fileID);
%    B_1st = info.Groups(18).Name;
    B_1st = info.Groups(16).Name;
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
     % important to check the resolution 
    res=17.77/length(Bx(:,1,1)); %ppc res=0.06, mime res=0.04, VA2c_008 res=0.05  
    %res=24/length(Bx(:,1,1));
    res_x=res; res_y=res; res_z=res;
    threshold=pi/res; 
    spectraB=CalcSpectra_CB(Bx,By,Bz,res,threshold);
    assignin('base',strcat('spectraB_',string(i)), spectraB);
    clearvars Bx By Bz; 
    cd(path) 
    %{
    cd(path) 
    V_1st = info.Groups(25).Name;
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
    spectravi=CalcSpectra_CB(vix,viy,viz,res,threshold);
    assignin('base',strcat('spectravi_',predir), spectravi);
    clearvars vix viy viz; 
    cd(path)
    %}
    %------------------------------------------------------------------
    %{
    spectraB_1 = spectraB_DATACB103_1;
    spectravi_1 = spectravi_DATACB103_1;
    spectraB_2 = spectraB_DATACB104;
    spectravi_2 = spectravi_DATACB104;
    spectraB_1 = spectraB_DATACB103_2L;
    spectravi_1 = spectravi_DATACB103_2L;
    spectraB_2 = spectraB2L_DATACB103_2L;
    spectravi_2 = spectravi2L_DATACB103_2L;
    %}
    
% Making plots ppc
%#----------------------------------------------------------------------    
    
    % Plots Magnetic field spectra
    %#----------------------------------------------------------------------    
    %axes('NextPlot','replacechildren', 'ColorOrder',C);
    
    f7=figure(7);
    %C = linspecer(10);
    %----------------------------------------------------------------------------------
    %
    lg=loglog(spectraB_1.kper,spectraB_1.P1Dper/rms(spectraB_1.P1Dper,'all'),...
        spectraB_2.kper,spectraB_2.P1Dper/rms(spectraB_2.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\perp}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\perp}d_{e}=1$$'; 
    txt5 = '$$ k_{\perp}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1428,'--k'); 
    %xline(23.5702,'--k');
    legend({'$10^{-3}$','$10^{-4}$'},'Interpreter','latex')
    xlim([0.26 pi/0.06]);  
    %ylim([10^(-5) 100])
    %}
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_\perp$','Interpreter','latex','FontSize',20)
    hold off

   
    f8=figure(8);
    %----------------------------------------------------------------------------------
    %
    lg=loglog(spectraB_1.kpar,spectraB_1.P1Dpar/rms(spectraB_1.P1Dpar,'all'),...
        spectraB_2.kpar,spectraB_2.P1Dpar/rms(spectraB_2.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; 
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\|}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\|}d_{e}=1$$'; 
    txt5 = '$$ k_{\|}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.28 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k');
    legend({'$10^{-3}$','$10^{-4}$'},'Interpreter','latex')
    xlim([0.05 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_{\|}$','Interpreter','latex','FontSize',20)
    hold off

    %#----------------------------------------------------------------------
   
    % Plots ion velocity spectra
    %#----------------------------------------------------------------------     
    
    f71=figure(71);
   %C = linspecer(10);
    %
    lg=loglog(spectravi_1.kper,spectravi_1.P1Dper/rms(spectravi_1.P1Dper,'all'),...
        spectravi_2.kper,spectravi_2.P1Dper/rms(spectravi_2.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\perp}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\perp}d_{e}=1$$'; 
    txt5 = '$$ k_{\perp}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1428,'--k'); 
    %xline(23.5702,'--k');
    legend({'$10^{-3}$','$10^{-4}$'},'Interpreter','latex')
    xlim([0.26 pi/0.06]);  
    %ylim([10^(-5) 100])
    %}
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{v}_{i}}1D_\perp$','Interpreter','latex','FontSize',20)
    hold off

    f81=figure(81);
    %----------------------------------------------------------------------------------
    %
    lg=loglog(spectravi_1.kpar,spectravi_1.P1Dpar/rms(spectravi_1.P1Dpar,'all'),...
        spectravi_2.kpar,spectravi_2.P1Dpar/rms(spectravi_2.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; 
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\|}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\|}d_{e}=1$$'; 
    txt5 = '$$ k_{\|}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.28 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k');
    legend({'$10^{-3}$','$10^{-4}$'},'Interpreter','latex')
    xlim([0.05 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{vi}}1D_{\|}$','Interpreter','latex','FontSize',20)
    hold off
    %--------------------------------------------------------------------------------

    
    % Save the plots
    %-------------------------------------------------------------------------
    cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis/Images';    
    saveas(f7,strcat(S(i).name,'1_perDspectrum_B_0304.png'));
    saveas(f8,strcat(S(i).name,'1_parDspectrum_B_0304.png'));
    saveas(f71,strcat(S(i).name,'1_perDspectrum_vi_0304.png'));
    saveas(f81,strcat(S(i).name,'1_parDspectrum_vi_0304.png'));
    %cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/ppc_50';
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Making plots ppc
%#----------------------------------------------------------------------    
    
    % Plots Magnetic field spectra
    %#----------------------------------------------------------------------    
    %axes('NextPlot','replacechildren', 'ColorOrder',C);
    
    f7=figure(7);
    %C = linspecer(10);
    %----------------------------------------------------------------------------------
    %
    lg=loglog(spectraB_1.kper,spectraB_1.P1Dper/rms(spectraB_1.P1Dper,'all'),...
        spectraB_2.kper,spectraB_2.P1Dper/rms(spectraB_2.P1Dper,'all'),...
        spectraB_3.kper,spectraB_3.P1Dper/rms(spectraB_3.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\perp}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\perp}d_{e}=1$$'; 
    txt5 = '$$ k_{\perp}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1428,'--k'); 
    %xline(23.5702,'--k');
    legend({'$10^{-3}$','$10^{-4}$','$10^{-3},2L$'},'Interpreter','latex')
    xlim([0.26 pi/0.06]);  
    %ylim([10^(-5) 100])
    %}
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_\perp$','Interpreter','latex','FontSize',20)
    hold off

   
    f8=figure(8);
    %----------------------------------------------------------------------------------
    %
    lg=loglog(spectraB_1.kpar,spectraB_1.P1Dpar/rms(spectraB_1.P1Dpar,'all'),...
        spectraB_2.kpar,spectraB_2.P1Dpar/rms(spectraB_2.P1Dpar,'all'),...
        spectraB_3.kpar,spectraB_3.P1Dpar/rms(spectraB_3.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\|}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\|}d_{e}=1$$'; 
    txt5 = '$$ k_{\|}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.28 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k');
    legend({'$10^{-3}$','$10^{-4}$','$10^{-3},2L$'},'Interpreter','latex')
    xlim([0.05 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_{\|}$','Interpreter','latex','FontSize',20)
    hold off

    %#----------------------------------------------------------------------
   
    % Plots ion velocity spectra
    %#----------------------------------------------------------------------     
    
    f71=figure(71);
   %C = linspecer(10);
    %----------------------------------------------------------------------------------
    lg=loglog(spectravi_1.kper,spectravi_1.P1Dper/rms(spectravi_1.P1Dper,'all'),...
        spectravi_2.kper,spectravi_2.P1Dper/rms(spectravi_2.P1Dper,'all'),...
        spectravi_3.kper,spectravi_3.P1Dper/rms(spectravi_3.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\perp}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\perp}d_{e}=1$$'; 
    txt5 = '$$ k_{\perp}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1428,'--k'); 
    %xline(23.5702,'--k');
    legend({'$10^{-3}$','$10^{-4}$','$10^{-3},2L$'},'Interpreter','latex')
    xlim([0.26 pi/0.06]);  
    %ylim([10^(-5) 100])
    %}
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{v}_{i}}1D_\perp$','Interpreter','latex','FontSize',20)
    hold off

    f81=figure(81);
    %----------------------------------------------------------------------------------
    %
    lg=loglog(spectravi_1.kpar,spectravi_1.P1Dpar/rms(spectravi_1.P1Dpar,'all'),...
        spectravi_2.kpar,spectravi_2.P1Dpar/rms(spectravi_2.P1Dpar,'all'),...
        spectravi_3.kpar,spectravi_3.P1Dpar/rms(spectravi_3.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\|}\lambda_{D}=1$$'; 
    txt3 = '$$ k_{\|}d_{e}=1$$'; 
    txt5 = '$$ k_{\|}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.28 , 10,txt1,'Interpreter','latex','FontSize',20); 
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k');
    legend({'$10^{-3}$','$10^{-4}$','$10^{-3},2L$'},'Interpreter','latex')
    xlim([0.05 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{vi}}1D_{\|}$','Interpreter','latex','FontSize',20)
    hold off
    %--------------------------------------------------------------------------------

    
    % Save the plots
    %-------------------------------------------------------------------------
    cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis/Images';    
    saveas(f7,strcat(S(i).name,'1_perDspectrum_B_0304_2L.png'));
    saveas(f8,strcat(S(i).name,'1_parDspectrum_B_0304_2L.png'));
    saveas(f71,strcat(S(i).name,'1_perDspectrum_vi_0304_2L.png'));
    saveas(f81,strcat(S(i).name,'1_parDspectrum_vi_0304_2L.png'));
    %cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/ppc_50';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------------------
    mime=100;
    spectraB_i = spectraB_4;
    f1=figure(1);
    [C,h]=contourf(spectraB_i.kpar,spectraB_i.kper,log10(spectraB_i.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([0 13]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',20)
    xlabel('$k_\| d_i$','Interpreter','latex','FontSize',20)
    ylabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    title('$P_{\mathbf{B}}2D_{-3,2L}$','Interpreter','latex','FontSize',20)
    hold on
    xlim([0.05 52.35]);    ylim([0.35 52.35]);
    %xlim([0.1 52.35]);    ylim([0.35 52.35]);
    yline(sqrt(mime),'--k');
    xline(sqrt(mime),'--k');
    hold off
    cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis/Images'; 
    saveas(f1,'1_perDspectrum_B_03_2L.png');
    %#----------------------------------------------------------------------

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f5=figure(5);
    loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper,'all'), 'k',spectravi.kper,spectravi.P1Dper/rms(spectravi.P1Dper,'all'), 'r')
    set(gca,'FontSize',18)
    x2=log10([1 10]); y2=x2.*(-3)+0.4; x3=log10([0.099 1]); y3=x3.*(-1.7);
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x2,10.^y2,'--k')
    %loglog(10.^x3,10.^y3,'--b')
    loglog([x33 x33],y1,'--')
    loglog([x44 x44],y1,'--')
    loglog([x55 x55],y1,'--')
    % Text in the plot
    txt3 = '$$ k_{\perp}d_e=1$$';
    txt4 = '$$ k_{\perp}d_i=1$$';
    txt5 = '$$ k_{\perp}^{-3}$$';
    xline(0.52,'--k')
    %txt5 = '$$ k_{\perp}^{-2.33}$$'; this is the teorethical value for the
    %magnetic field
    txt6 = '$$ k_{\perp}^{-1.7}$$';
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    text(4.5 , 30,txt3,'Interpreter','latex','FontSize',15); 
    text(1.1 , 30,txt4,'Interpreter','latex','FontSize',15); 
    text(1.1 , 0.2,txt5,'Interpreter','latex','FontSize',15);
    %text(0.35 , 16,txt6,'Interpreter','latex','FontSize',15); 
    text(16 , 0.3,txt8,'Interpreter','latex','FontSize',15);
    legend({'$P_{\tilde{B}}1D_\perp$','$P_{\tilde{v}_i}1D_\perp$'},'Interpreter','latex')
    xlabel('$k_\perp d_i$','Interpreter','latex')
    ylabel('$P_{\psi}1D_\perp$','Interpreter','latex')
    xlim([0.26 52.36]); ylim([10^(-5) 100])
    hold off
    
    
    f6=figure(6);
    loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), '-*k',spectravi.kpar,spectravi.P1Dpar/rms(spectravi.P1Dpar), '-*r')
    %loglog(spectraB.kpar,spectraB.P1Dpar,'k', spectravi.kpar,spectravi.P1Dpar,'r')
    set(gca,'FontSize',18)
    x=log10([0.099 1]);
    x2=log10([1 10]);
    x3=log10([5 30]);
    %y=x.*(-2.3)+12.2;
    y=x.*(-2)-1  ;%+12.2;
    y2=x2.*(-5);%+12.2;
    y3=x3.*(-0.3);%+10.6;
    x33=1;    x44=10;   x55=14.2857;
    xlim([0.05 52.36]); ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x,10.^y,'--b')
    loglog([x33 x33],y1,'--')
    loglog([x44 x44],y1,'--')
    loglog([x55 x55],y1,'--')
    % Text in the plot
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$P_{\psi}1D_\|$','Interpreter','latex')
    txt3 = '$$ k_{||}d_e=1$$'; txt4 = '$$ k_{||}d_i=1$$';    txt5 = '$$ k_{||}^{-5}$$';
    txt6 = '$$ k_{||}^{-0.3}$$';    txt7 = '$$ k_{||}^{-2}$$';    txt8 = '$$ k_{||}\lambda_D=1$$';
    text(3 , 10,txt3,'Interpreter','latex','FontSize',15); text(1.1 , 30,txt4,'Interpreter','latex','FontSize',15)
    %text(2 , 0.4,txt5,'Interpreter','latex'); %text(18 , 1,txt6,'Interpreter','latex')
    text(0.5 , 7,txt7,'Interpreter','latex','FontSize',15); text(16 , 0.3,txt8,'Interpreter','latex','FontSize',15)
    legend({'$P_{\tilde{B}}1D_\|$','$P_{\tilde{v}_i}1D_\|$'},'Interpreter','latex')
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    
    f9=figure(9);
    %loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'r')
    scatter(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'r')
    f91=figure(91);
    loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k')
    %scatter(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'r')
    %----------------------------------------------------------------------
    cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/Images';
    % Save the plots
    %-------------------------------------------------------------------------
    saveas(f1,strcat(S(i).name,'2Dspectrum_Bpd_Isotropic_ppc50.png'));
    saveas(f5,strcat(S(i).name,'1_perDspectrum_Isotropic_ppc50.png'));
    saveas(f6,strcat(S(i).name,'1_parDspectrum_Isotropic_ppc50.png'));
    
    saveas(f7,strcat(S(i).name,'1_perDspectrum_B.png'));
    saveas(f8,strcat(S(i).name,'1_parDspectrum_B.png'));
    saveas(f71,strcat(S(i).name,'1_perDspectrum_vi.png'));
    saveas(f81,strcat(S(i).name,'1_parDspectrum_vi.png'));
    %cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/ppc_50';
    
%end
