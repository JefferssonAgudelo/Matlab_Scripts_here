%#-----------------------------------------------------------------------
% This file is to make the plots related with the calculation of the
% spectrum
% This is the file that I use to checke the different runs changing
%parameters
%#----------------------------------------------------------------------


%#----------------------------------------------------------------------
% This is to generate the spectra. One this has been done. There is no need
% to do it again
%#----------------------------------------------------------------------
predir = 'VA2c_008'; % ppc_100, mime_500  
prepath = '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_VA2c/'; %changing_ppc/, changing_mime/ 
path = strcat(prepath,predir);
cd(path) 

S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
%i=3; %ppc_50
i=1;
%for i=1:3
    disp(strcat('Computing step ...',S(i).name)) %6 is 2000   
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    B_1st = info.Groups(18).Name;
    %V_1st = info.Groups(25).Name;
    %n_1st = info.Groups(23).Name;
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    %{
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    Ex=h5read(fileID,strcat(B_1st,'/ex/p0/3d'));
    Ey=h5read(fileID,strcat(B_1st,'/ey/p0/3d'));
    Ez=h5read(fileID,strcat(B_1st,'/ez/p0/3d'));
    %}
    %------------------------------------------------------------------
    
    % important to check the resolution 
    res=12/length(Bx); %ppc res=0.06, mime res=0.04, VA2c_008 res=0.05  
    res_x=res; res_y=res; res_z=res;
    threshold=pi/res; 
   
    cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
    % kpar (1,99), kper (1,142), P2D (99,142) P1par (99,1) P1per (1,142)
    spectraB=CalcSpectra_CB(Bx,By,Bz,res,threshold);
    %spectravi=CalcSpectra_CB(vix,viy,viz,res,threshold);
    %spectrani=CalcSpectra_CB_n(ni,res,threshold); %here this is a different file
    %spectrane=CalcSpectra_CB_n(ne,res,threshold);
    
    assignin('base',strcat('spectraB_',predir), spectraB)
    %assignin('base',strcat('spectravi_',predir), spectravi)
    %assignin('base',strcat('spectrani_',predir), spectrani)
    %assignin('base',strcat('spectrane_',predir), spectrane)
%#----------------------------------------------------------------------    
    

%#----------------------------------------------------------------------    

%{
spectraB_1=spectraB_ppc_50;
spectraB_2=spectraB_ppc_100;
spectraB_3=spectraB_ppc_150;
spectraB_4=spectraB_ppc_200;
spectraB_5=spectraB_ppc_250;
spectraB_6=spectraB_ppc_300;
spectravi_1=spectravi_ppc_50;
spectravi_2=spectravi_ppc_100;
spectravi_3=spectravi_ppc_150;
spectravi_4=spectravi_ppc_200;
spectravi_5=spectravi_ppc_250;
spectravi_6=spectravi_ppc_300;


spectraB_1=spectraB_mime_1;
spectraB_2=spectraB_mime_10;
spectraB_3=spectraB_mime_50;
spectraB_4=spectraB_mime_100;
spectraB_5=spectraB_mime_150;
spectraB_6=spectraB_mime_200;
spectraB_7=spectraB_mime_250;
spectraB_8=spectraB_mime_300;
spectraB_9=spectraB_mime_400;
spectraB_10=spectraB_mime_500;
spectravi_1=spectravi_mime_1;
spectravi_2=spectravi_mime_10;
spectravi_3=spectravi_mime_50;
spectravi_4=spectravi_mime_100;
spectravi_5=spectravi_mime_150;
spectravi_6=spectravi_mime_200;
spectravi_7=spectravi_mime_250;
spectravi_8=spectravi_mime_300;
spectravi_9=spectravi_mime_400;
spectravi_10=spectravi_mime_500;


spectraB_1=spectraB_mime_100; %0.06
spectraB_2=spectraB_VA2c_008; %0.08
spectraB_3=spectraB_ppc_100; %0.1

spectravi_1=spectravi_mime_100; %0.06
spectravi_2=spectravi_VA2c_008; %0.08
spectravi_3=spectravi_ppc_100; %0.1
%}
%#----------------------------------------------------------------------    


% Making plots ppc
%#----------------------------------------------------------------------    
    
    % Plots Magnetic field spectra
    %#----------------------------------------------------------------------    
    %axes('NextPlot','replacechildren', 'ColorOrder',C);
    
    f7=figure(7);
    %C = linspecer(10);
    %----------------------------------------------------------------------------------
    %VA2cs    %{
    lg=loglog(spectraB_1.kper,spectraB_1.P1Dper/rms(spectraB_1.P1Dper,'all'),...
        spectraB_2.kper,spectraB_2.P1Dper/rms(spectraB_2.P1Dper,'all'),...
        spectraB_3.kper,spectraB_3.P1Dper/rms(spectraB_3.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\perp}\lambda_{D0.1}=1$$'; txt2 = '$$ k_{\perp}\lambda_{D0.08}=1$$';
    txt3 = '$$ k_{\perp}d_{e}=1$$'; txt4 = '$$ k_{\perp}\lambda_{D0.06}=1$$';
    txt5 = '$$ k_{\perp}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    text(17.7 , 30,txt2,'Interpreter','latex','FontSize',20);
    text(23.6 , 10,txt4,'Interpreter','latex','FontSize',20);
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k'); xline(17.6777,'--k'); xline(23.5702,'--k');
    legend({'$VA2c_{0.06}$','$VA2c_{0.08}$','$VA2c_{0.1}$'},'Interpreter','latex')
    xlim([0.52 pi/0.04]);  
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %ppc
    %{
    lg=loglog(spectraB_1.kper,spectraB_1.P1Dper/rms(spectraB_1.P1Dper,'all'),...
        spectraB_2.kper,spectraB_2.P1Dper/rms(spectraB_2.P1Dper,'all'),...
        spectraB_3.kper,spectraB_3.P1Dper/rms(spectraB_3.P1Dper,'all'),...
        spectraB_4.kper,spectraB_4.P1Dper/rms(spectraB_4.P1Dper,'all'),...
        spectraB_5.kper,spectraB_5.P1Dper/rms(spectraB_5.P1Dper,'all'),...
        spectraB_6.kper,spectraB_6.P1Dper/rms(spectraB_6.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    xline(1,'--k'); xline(23.6,'--k'); 
    xline(sqrt(100),'--k');  
    txt1 = '$$ k_{\perp}d_{i}=1$$'; txt2 = '$$ k_{\perp}d_{e}=1$$';
    txt3 = '$$ k_{\perp}\lambda_{D}=1$$'; 
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(23.6 , 30,txt3,'Interpreter','latex','FontSize',20);
    legend({'$ppc_{50}$','$ppc_{100}$','$ppc_{150}$','$ppc_{200}$',...
       '$ppc_{250}$','$ppc_{300}$'},'Interpreter','latex')
    xlim([0.52 pi/0.06]); 
    ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %mime
    %{
    lg=loglog(spectraB_1.kper,spectraB_1.P1Dper/rms(spectraB_1.P1Dper,'all'),...
        spectraB_2.kper,spectraB_2.P1Dper/rms(spectraB_2.P1Dper,'all'),...
        spectraB_3.kper,spectraB_3.P1Dper/rms(spectraB_3.P1Dper,'all'),...
        spectraB_4.kper,spectraB_4.P1Dper/rms(spectraB_4.P1Dper,'all'),...
        spectraB_5.kper,spectraB_5.P1Dper/rms(spectraB_5.P1Dper,'all'),...
        spectraB_6.kper,spectraB_6.P1Dper/rms(spectraB_6.P1Dper,'all'),...
        spectraB_7.kper,spectraB_7.P1Dper/rms(spectraB_7.P1Dper,'all'),'-*',...
        spectraB_8.kper,spectraB_8.P1Dper/rms(spectraB_8.P1Dper,'all'),'-*',...
        spectraB_9.kper,spectraB_9.P1Dper/rms(spectraB_9.P1Dper,'all'),'-*',...
        spectraB_10.kper,spectraB_10.P1Dper/rms(spectraB_10.P1Dper,'all'),'-*');
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; lg(9).LineWidth = 2;
    lg(10).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    xline(1,'--k'); xline(sqrt(10),'--k'); 
    xline(sqrt(100),'--k');  xline(sqrt(300),'--k'); xline(sqrt(500),'--k');
    txt1 = '$$ k_{\perp}d_{e1}=1$$'; txt2 = '$$ k_{\perp}d_{e10}=1$$';
    txt3 = '$$ k_{\perp}d_{e100}=1$$'; txt4 = '$$ k_{\perp}d_{e300}=1$$';
    txt5 = '$$ k_{\perp}d_{e500}=1$$';
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(10) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(sqrt(300) , 30,txt4,'Interpreter','latex','FontSize',20);
    text(sqrt(500) , 10,txt5,'Interpreter','latex','FontSize',20);
    legend({'$mime_{1}$','$mime_{10}$','$mime_{50}$','$mime_{100}$',...
            '$mime_{150}$','$mime_{200}$','$mime_{250}$','$mime_{300}$',...
            '$mime_{400}$', '$mime_{500}$'},'Interpreter','latex')
    xlim([0.52 pi/0.05]); 
    %ylim([10^(-5) 100])
      %}
    %----------------------------------------------------------------------------------
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_\perp$','Interpreter','latex','FontSize',20)
    hold off

   
    f8=figure(8);
    %----------------------------------------------------------------------------------
    %VA2c
    %{
    lg=loglog(spectraB_1.kpar,spectraB_1.P1Dpar/rms(spectraB_1.P1Dpar,'all'),...
        spectraB_2.kpar,spectraB_2.P1Dpar/rms(spectraB_2.P1Dpar,'all'),...
        spectraB_3.kpar,spectraB_3.P1Dpar/rms(spectraB_3.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    txt1 = '$$ k_{\|}\lambda_{D0.1}=1$$'; txt2 = '$$ k_{\|}\lambda_{D0.08}=1$$';
    txt3 = '$$ k_{\|}d_{e}=1$$'; txt4 = '$$ k_{\|}\lambda_{D0.06}=1$$';
    txt5 = '$$ k_{\|}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    text(17.7 , 30,txt2,'Interpreter','latex','FontSize',20);
    text(23.6 , 10,txt4,'Interpreter','latex','FontSize',20);
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k'); xline(17.6777,'--k'); xline(23.5702,'--k');
    legend({'$VA2c_{0.06}$','$VA2c_{0.08}$','$VA2c_{0.1}$'},'Interpreter','latex')
    xlim([0.52 pi/0.04]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %ppc
    %{
    lg=loglog(spectraB_1.kpar,spectraB_1.P1Dpar/rms(spectraB_1.P1Dpar,'all'),...
        spectraB_2.kpar,spectraB_2.P1Dpar/rms(spectraB_2.P1Dpar,'all'),...
        spectraB_3.kpar,spectraB_3.P1Dpar/rms(spectraB_3.P1Dpar,'all'),...
        spectraB_4.kpar,spectraB_4.P1Dpar/rms(spectraB_4.P1Dpar,'all'),...
        spectraB_5.kpar,spectraB_5.P1Dpar/rms(spectraB_5.P1Dpar,'all'),...
        spectraB_6.kpar,spectraB_6.P1Dpar/rms(spectraB_6.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on 
    xline(1,'--k'); xline(23.6,'--k'); 
    xline(sqrt(100),'--k'); 
    txt1 = '$$ k_{\|}d_{i}=1$$'; txt2 = '$$ k_{\|}d_{e}=1$$';
    txt3 = '$$ k_{\|}\lambda_{D}=1$$'; 
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(23.6 , 30,txt3,'Interpreter','latex','FontSize',20);
    legend({'$ppc_{50}$','$ppc_{100}$','$ppc_{150}$','$ppc_{200}$',...
       '$ppc_{250}$','$ppc_{300}$'},'Interpreter','latex')
    xlim([0.52 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %mime
    %{
    lg=loglog(spectraB_1.kpar,spectraB_1.P1Dpar/rms(spectraB_1.P1Dpar,'all'),...
        spectraB_2.kpar,spectraB_2.P1Dpar/rms(spectraB_2.P1Dpar,'all'),...
        spectraB_3.kpar,spectraB_3.P1Dpar/rms(spectraB_3.P1Dpar,'all'),...
        spectraB_4.kpar,spectraB_4.P1Dpar/rms(spectraB_4.P1Dpar,'all'),...
        spectraB_5.kpar,spectraB_5.P1Dpar/rms(spectraB_5.P1Dpar,'all'),...
        spectraB_6.kpar,spectraB_6.P1Dpar/rms(spectraB_6.P1Dpar,'all'),...
        spectraB_7.kpar,spectraB_7.P1Dpar/rms(spectraB_7.P1Dpar,'all'),'-*',...
        spectraB_8.kpar,spectraB_8.P1Dpar/rms(spectraB_8.P1Dpar,'all'),'-*',...
        spectraB_9.kpar,spectraB_9.P1Dpar/rms(spectraB_9.P1Dpar,'all'),'-*',...
        spectraB_10.kpar,spectraB_10.P1Dpar/rms(spectraB_10.P1Dpar,'all'),'-*');
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; lg(9).LineWidth = 2;
    lg(10).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    xline(1,'--k'); xline(sqrt(10),'--k'); 
    xline(sqrt(100),'--k');  xline(sqrt(300),'--k'); xline(sqrt(500),'--k');
    txt1 = '$$ k_{\|}d_{e1}=1$$'; txt2 = '$$ k_{\|}d_{e10}=1$$';
    txt3 = '$$ k_{\|}d_{e100}=1$$'; txt4 = '$$ k_{\|}d_{e300}=1$$';
    txt5 = '$$ k_{\|}d_{e500}=1$$';
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(10) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(sqrt(300) , 30,txt4,'Interpreter','latex','FontSize',20);
    text(sqrt(500) , 10,txt5,'Interpreter','latex','FontSize',20);
    legend({'$mime_{1}$','$mime_{10}$','$mime_{50}$','$mime_{100}$',...
            '$mime_{150}$','$mime_{200}$','$mime_{250}$','$mime_{300}$',...
            '$mime_{400}$', '$mime_{500}$'},'Interpreter','latex')
    xlim([0.52 pi/0.05]); 
    %ylim([10^(-5) 100])
     %}
    %----------------------------------------------------------------------------------
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_{\|}$','Interpreter','latex','FontSize',20)
    xlim([0.52 pi/0.04]); ylim([10^(-5) 100])
    hold off

    %#----------------------------------------------------------------------
   
    % Plots ion velocity spectra
    %#----------------------------------------------------------------------     
    
    f71=figure(71);
    %----------------------------------------------------------------------------------
    %VA2c
    %{
    lg=loglog(spectravi_1.kper,spectravi_1.P1Dper/rms(spectravi_1.P1Dper,'all'),...
        spectravi_2.kper,spectravi_2.P1Dper/rms(spectravi_2.P1Dper,'all'),...
        spectravi_3.kper,spectravi_3.P1Dper/rms(spectravi_3.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on 
    txt1 = '$$ k_{\perp}\lambda_{D0.1}=1$$'; txt2 = '$$ k_{\perp}\lambda_{D0.08}=1$$';
    txt3 = '$$ k_{\perp}d_{e}=1$$'; txt4 = '$$ k_{\perp}\lambda_{D0.06}=1$$';
    txt5 = '$$ k_{\perp}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    text(17.7 , 30,txt2,'Interpreter','latex','FontSize',20);
    text(23.6 , 10,txt4,'Interpreter','latex','FontSize',20);
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k'); xline(17.6777,'--k'); xline(23.5702,'--k');
    legend({'$VA2c_{0.06}$','$VA2c_{0.08}$','$VA2c_{0.1}$'},'Interpreter','latex')
    xlim([0.52 pi/0.04]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %ppc
    %{
    lg=loglog(spectravi_1.kper,spectravi_1.P1Dper/rms(spectravi_1.P1Dper,'all'),...
        spectravi_2.kper,spectravi_2.P1Dper/rms(spectravi_2.P1Dper,'all'),...
        spectravi_3.kper,spectravi_3.P1Dper/rms(spectravi_3.P1Dper,'all'),...
        spectravi_4.kper,spectravi_4.P1Dper/rms(spectravi_4.P1Dper,'all'),...
        spectravi_5.kper,spectravi_5.P1Dper/rms(spectravi_5.P1Dper,'all'),...
        spectravi_6.kper,spectravi_6.P1Dper/rms(spectravi_6.P1Dper,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    xline(1,'--k'); xline(23.6,'--k'); 
    xline(sqrt(100),'--k'); 
    txt1 = '$$ k_{\perp}d_{i}=1$$'; txt2 = '$$ k_{\perp}d_{e}=1$$';
    txt3 = '$$ k_{\perp}\lambda_{D}=1$$'; 
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(23.6 , 30,txt3,'Interpreter','latex','FontSize',20);
    legend({'$ppc_{50}$','$ppc_{100}$','$ppc_{150}$','$ppc_{200}$',...
       '$ppc_{250}$','$ppc_{300}$'},'Interpreter','latex')
    xlim([0.52 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %mime
    %{
    lg=loglog(spectravi_1.kper,spectravi_1.P1Dper/rms(spectravi_1.P1Dper,'all'),...
        spectravi_2.kper,spectravi_2.P1Dper/rms(spectravi_2.P1Dper,'all'),...
        spectravi_3.kper,spectravi_3.P1Dper/rms(spectravi_3.P1Dper,'all'),...
        spectravi_4.kper,spectravi_4.P1Dper/rms(spectravi_4.P1Dper,'all'),...
        spectravi_5.kper,spectravi_5.P1Dper/rms(spectravi_5.P1Dper,'all'),...
        spectravi_6.kper,spectravi_6.P1Dper/rms(spectravi_6.P1Dper,'all'),...
        spectravi_7.kper,spectravi_7.P1Dper/rms(spectravi_7.P1Dper,'all'),'-*',...
        spectravi_8.kper,spectravi_8.P1Dper/rms(spectravi_8.P1Dper,'all'),'-*',...
        spectravi_9.kper,spectravi_9.P1Dper/rms(spectravi_9.P1Dper,'all'),'-*',...
        spectravi_10.kper,spectravi_10.P1Dper/rms(spectravi_10.P1Dper,'all'),'-*');
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; lg(9).LineWidth = 2;
    lg(10).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    xline(1,'--k'); xline(sqrt(10),'--k'); 
    xline(sqrt(100),'--k');  xline(sqrt(300),'--k'); xline(sqrt(500),'--k');
    txt1 = '$$ k_{\perp}d_{e1}=1$$'; txt2 = '$$ k_{\perp}d_{e10}=1$$';
    txt3 = '$$ k_{\perp}d_{e100}=1$$'; txt4 = '$$ k_{\perp}d_{e300}=1$$';
    txt5 = '$$ k_{\perp}d_{e500}=1$$';
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(10) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(sqrt(300) , 30,txt4,'Interpreter','latex','FontSize',20);
    text(sqrt(500) , 10,txt5,'Interpreter','latex','FontSize',20);
    legend({'$mime_{1}$','$mime_{10}$','$mime_{50}$','$mime_{100}$',...
            '$mime_{150}$','$mime_{200}$','$mime_{250}$','$mime_{300}$',...
            '$mime_{400}$', '$mime_{500}$'},'Interpreter','latex')
    xlim([0.52 pi/0.05]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{v}_{i}}1D_\perp$','Interpreter','latex','FontSize',20)
    xlim([0.52 pi/0.04]); ylim([10^(-5) 100])
    hold off

    f81=figure(81);
    %----------------------------------------------------------------------------------
    %VA2c
    %{
    lg=loglog(spectravi_1.kpar,spectravi_1.P1Dpar/rms(spectravi_1.P1Dpar,'all'),...
        spectravi_2.kpar,spectravi_2.P1Dpar/rms(spectravi_2.P1Dpar,'all'),...
        spectravi_3.kpar,spectravi_3.P1Dpar/rms(spectravi_3.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on 
    txt1 = '$$ k_{\|}\lambda_{D0.1}=1$$'; txt2 = '$$ k_{\|}\lambda_{D0.08}=1$$';
    txt3 = '$$ k_{\|}d_{e}=1$$'; txt4 = '$$ k_{\|}\lambda_{D0.06}=1$$';
    txt5 = '$$ k_{\|}d_{i}=1$$'; 
    text(1 , 30,txt5,'Interpreter','latex','FontSize',20);
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(14.15 , 10,txt1,'Interpreter','latex','FontSize',20); 
    text(17.7 , 30,txt2,'Interpreter','latex','FontSize',20);
    text(23.6 , 10,txt4,'Interpreter','latex','FontSize',20);
    xline(1,'--k');  
    xline(sqrt(100),'--k'); 
    xline(14.1421,'--k'); xline(17.6777,'--k'); xline(23.5702,'--k');
    legend({'$VA2c_{0.06}$','$VA2c_{0.08}$','$VA2c_{0.1}$'},'Interpreter','latex')
    xlim([0.52 pi/0.04]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %ppc
    %{
    lg=loglog(spectravi_1.kpar,spectravi_1.P1Dpar/rms(spectravi_1.P1Dpar,'all'),...
        spectravi_2.kpar,spectravi_2.P1Dpar/rms(spectravi_2.P1Dpar,'all'),...
        spectravi_3.kpar,spectravi_3.P1Dpar/rms(spectravi_3.P1Dpar,'all'),...
        spectravi_4.kpar,spectravi_4.P1Dpar/rms(spectravi_4.P1Dpar,'all'),...
        spectravi_5.kpar,spectravi_5.P1Dpar/rms(spectravi_5.P1Dpar,'all'),...
        spectravi_6.kpar,spectravi_6.P1Dpar/rms(spectravi_6.P1Dpar,'all'));
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    xline(1,'--k'); xline(23.6,'--k'); 
    xline(sqrt(100),'--k'); 
    txt1 = '$$ k_{\|}d_{i}=1$$'; txt2 = '$$ k_{\|}d_{e}=1$$';
    txt3 = '$$ k_{\|}\lambda_{D}=1$$'; 
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(23.6 , 30,txt3,'Interpreter','latex','FontSize',20);
    legend({'$ppc_{50}$','$ppc_{100}$','$ppc_{150}$','$ppc_{200}$',...
       '$ppc_{250}$','$ppc_{300}$'},'Interpreter','latex')
    xlim([0.52 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    %mime
    %{
    lg=loglog(spectravi_1.kpar,spectravi_1.P1Dpar/rms(spectravi_1.P1Dpar,'all'),...
        spectravi_2.kpar,spectravi_2.P1Dpar/rms(spectravi_2.P1Dpar,'all'),...
        spectravi_3.kpar,spectravi_3.P1Dpar/rms(spectravi_3.P1Dpar,'all'),...
        spectravi_4.kpar,spectravi_4.P1Dpar/rms(spectravi_4.P1Dpar,'all'),...
        spectravi_5.kpar,spectravi_5.P1Dpar/rms(spectravi_5.P1Dpar,'all'),...
        spectravi_6.kpar,spectravi_6.P1Dpar/rms(spectravi_6.P1Dpar,'all'),...
        spectravi_7.kpar,spectravi_7.P1Dpar/rms(spectravi_7.P1Dpar,'all'),'-*',...
        spectravi_8.kpar,spectravi_8.P1Dpar/rms(spectravi_8.P1Dpar,'all'),'-*',...
        spectravi_9.kpar,spectravi_9.P1Dpar/rms(spectravi_9.P1Dpar,'all'),'-*',...
        spectravi_10.kpar,spectravi_10.P1Dpar/rms(spectravi_10.P1Dpar,'all'),'-*');
    lg(1).LineWidth = 2; lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; lg(9).LineWidth = 2;
    lg(10).LineWidth = 2;
    set(gca,'FontSize',20) 
    hold on
    xline(1,'--k'); xline(sqrt(10),'--k'); 
    xline(sqrt(100),'--k');  xline(sqrt(300),'--k'); xline(sqrt(500),'--k');
    txt1 = '$$ k_{\|}d_{e1}=1$$'; txt2 = '$$ k_{\|}d_{e10}=1$$';
    txt3 = '$$ k_{\|}d_{e100}=1$$'; txt4 = '$$ k_{\|}d_{e300}=1$$';
    txt5 = '$$ k_{\|}d_{e500}=1$$';
    text(1 , 30,txt1,'Interpreter','latex','FontSize',20); 
    text(sqrt(10) , 30,txt2,'Interpreter','latex','FontSize',20); 
    text(sqrt(100) , 30,txt3,'Interpreter','latex','FontSize',20);
    text(sqrt(300) , 30,txt4,'Interpreter','latex','FontSize',20);
    text(sqrt(500) , 10,txt5,'Interpreter','latex','FontSize',20);
    legend({'$mime_{1}$','$mime_{10}$','$mime_{50}$','$mime_{100}$',...
            '$mime_{150}$','$mime_{200}$','$mime_{250}$','$mime_{300}$',...
            '$mime_{400}$', '$mime_{500}$'},'Interpreter','latex')
    xlim([0.52 pi/0.05]); 
    %ylim([10^(-5) 100])
    %}
    %---------------------------------------------------------------------------------- 
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{v}_{i}}1D_{\|}$','Interpreter','latex','FontSize',20)
    xlim([0.52 pi/0.04]); ylim([10^(-5) 100])
    hold off
%--------------------------------------------------------------------------------

    
    % Save the plots
    %-------------------------------------------------------------------------
    cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis/Images';    
    saveas(f7,strcat(S(i).name,'1_perDspectrum_B_ppc.png'));
    saveas(f8,strcat(S(i).name,'1_parDspectrum_B_mime.png'));
    saveas(f71,strcat(S(i).name,'1_perDspectrum_vi_mime.png'));
    saveas(f81,strcat(S(i).name,'1_parDspectrum_vi_mime.png'));
    %cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/ppc_50';
    
    
    %----------------------------------------------------------------------
    mime=100;
    spectraB_i = spectraB_3;
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
    title('$P_{\mathbf{B}}2D_{VA2c0.1}$','Interpreter','latex','FontSize',20)
    hold on
    ylim([0.52 62.8])
    yline(sqrt(mime),'--k');
    xline(sqrt(mime),'--k');
    hold off
    cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis/Images'; 
    saveas(f1,'1_perDspectrum_B_VA2c01.png');
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
