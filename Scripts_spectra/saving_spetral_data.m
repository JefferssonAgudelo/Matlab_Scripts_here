%#-----------------------------------------------------------------------
% This file is to make the plots related with the calculation of the
% spectrum
% This is the file that I use to checke the different runs changing
%parameters
%#----------------------------------------------------------------------

% Save the workspace in this directory 
%cd /Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis;

%#----------------------------------------------------------------------
% This is to generate the spectra. One this has been done. There is no need
% to do it again
%#----------------------------------------------------------------------
%predir = 'DATACB103_1'; % DATACB103_1 
%predir ='DATACB103_thin_8_dia';
%predir ='103_thinner';
predir = 'DATACB103_1';
prepath = '/Volumes/PSC_DiRAC_DATA/'; %changing_ppc/, changing_mime/ 

path = strcat(prepath,predir);
cd(path) 

S = dir(fullfile(path,'*_p000000.h5'));
N=numel(S); % number of files to use
H=zeros(N);
%i=18;
i=1;
for i=1:9
    disp(strcat('Computing step ...',S(i).name)) %6 is 2000   
    fileID =  S(i).name; %change the 1 per i
    h5disp(fileID);
    info = h5info(fileID);
    %B_1st = info.Groups(18).Name;
    B_1st = info.Groups(16).Name;
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
     % important to check the resolution 
    res=17.77/length(Bx(:,1,1)); %ppc res=0.06, mime res=0.04, VA2c_008 res=0.05  
    %res=24/length(Bx(:,1,1));
    %res=12/length(Bx(:,1,1)); %beta01
    res_x=res; res_y=res; res_z=res;
    threshold=pi/res; 
    spectraB=CalcSpectra_CB(Bx,By,Bz,res,threshold);
    assignin('base',strcat('spectraB_double_',string(i)), spectraB);
    clearvars Bx By Bz; 
    cd(path)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------------------
    mime=100;
    %spectraB_i = spectraB_4;
    f1=figure(1);
    %[C,h]=contourf(spectraB.kpar,spectraB.kper,log10(spectraB.P2D'),'LevelStep',0.5);
    [C,h]=contourf(spectraB_thi_2.kpar,spectraB_thi_2.kper,log10(spectraB_thi_2.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([0 13]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',20)
    xlabel('$k_\| d_i$','Interpreter','latex','FontSize',20)
    ylabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    title('$P_{\mathbf{B}}2D_{16}$','Interpreter','latex','FontSize',20)
    hold on
    %xlim([0.05 52.35]);    ylim([0.35 52.35]);
    xlim([0.06 52.35]);    ylim([0.35 52.35]);
    yline(sqrt(mime),'--k');
    xline(sqrt(mime),'--k');
    hold off
    %cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis/Images'; 
    saveas(f1,strcat('1_perDspectrum_B_double_',string(i),'.png'));
    %#----------------------------------------------------------------------
    
end
    


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
    
    
    %----------------------------------------------------------------------------------
    %
    f7=figure(7);
    %C = linspecer(10);
    %{
    lg=loglog(spectraB_thi_2.kper,spectraB_thi_2.P1Dper/rms(spectraB_thi_2.P1Dper,'all'),...
        spectraB_thi_3.kper,spectraB_thi_3.P1Dper/rms(spectraB_thi_3.P1Dper,'all'),...
        spectraB_thi_4.kper,spectraB_thi_4.P1Dper/rms(spectraB_thi_4.P1Dper,'all'),...
        spectraB_thi_5.kper,spectraB_thi_5.P1Dper/rms(spectraB_thi_5.P1Dper,'all'),...
        spectraB_thi_6.kper,spectraB_thi_6.P1Dper/rms(spectraB_thi_6.P1Dper,'all'));
    %}
    %{    
    lg=loglog(spectraB_2.kper,spectraB_2.P1Dper/rms(spectraB_2.P1Dper,'all'),...
        spectraB_3.kper,spectraB_3.P1Dper/rms(spectraB_3.P1Dper,'all'),...
        spectraB_4.kper,spectraB_4.P1Dper/rms(spectraB_4.P1Dper,'all'),...
        spectraB_5.kper,spectraB_5.P1Dper/rms(spectraB_5.P1Dper,'all'),...
        spectraB_6.kper,spectraB_6.P1Dper/rms(spectraB_6.P1Dper,'all'),...
        spectraB_7.kper,spectraB_7.P1Dper/rms(spectraB_7.P1Dper,'all'),...
        spectraB_8.kper,spectraB_8.P1Dper/rms(spectraB_8.P1Dper,'all'),...
        spectraB_9.kper,spectraB_9.P1Dper/rms(spectraB_9.P1Dper,'all'),...
        spectraB_10.kper,spectraB_10.P1Dper/rms(spectraB_10.P1Dper,'all'),...
        spectraB_11.kper,spectraB_11.P1Dper/rms(spectraB_11.P1Dper,'all'),...
        spectraB_12.kper,spectraB_12.P1Dper/rms(spectraB_12.P1Dper,'all'),...
        spectraB_13.kper,spectraB_13.P1Dper/rms(spectraB_13.P1Dper,'all'));
    lg(1).LineWidth = 2; 
    lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; lg(9).LineWidth = 2;
    lg(10).LineWidth = 2; lg(11).LineWidth = 2; lg(12).LineWidth = 2;
    lg(13).LineWidth = 2;
    %}
    
    lg=loglog(spectraB_double_2.kper,spectraB_double_2.P1Dper/rms(spectraB_double_2.P1Dper,'all'),...
        spectraB_double_3.kper,spectraB_double_3.P1Dper/rms(spectraB_double_3.P1Dper,'all'),...
        spectraB_double_4.kper,spectraB_double_4.P1Dper/rms(spectraB_double_4.P1Dper,'all'),...
        spectraB_double_5.kper,spectraB_double_5.P1Dper/rms(spectraB_double_5.P1Dper,'all'),...
        spectraB_double_6.kper,spectraB_double_6.P1Dper/rms(spectraB_double_6.P1Dper,'all'),...
        spectraB_double_7.kper,spectraB_double_7.P1Dper/rms(spectraB_double_7.P1Dper,'all'),...
        spectraB_double_8.kper,spectraB_double_8.P1Dper/rms(spectraB_double_8.P1Dper,'all'),...
        spectraB_double_9.kper,spectraB_double_9.P1Dper/rms(spectraB_double_9.P1Dper,'all'));
    lg(1).LineWidth = 2; 
    lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; 
    %}
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
    %{
legend({'$t_{1}$','$t_{2}$','$t_{3}$','$t_{4}$','$t_{5}$',...
        '$t_{6}$','$t_{7}$','$t_{8}$','$t_{9}$','$t_{10}$',...
        '$t_{11}$','$t_{12}$'},'Interpreter','latex')
    %}
    %legend({'$t_{1}$','$t_{2}$','$t_{3}$','$t_{4}$','$t_{5}$'},'Interpreter','latex')
    legend({'$t_{1}$','$t_{2}$','$t_{3}$','$t_{4}$','$t_{5}$',...
        '$t_{6}$','$t_{7}$','$t_{8}$'},'Interpreter','latex')
    xlim([0.35 pi/0.06]);  
    ylim([10^(-7) 100])
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_\perp$','Interpreter','latex','FontSize',20)
    hold off
    %}

   
    
    %----------------------------------------------------------------------------------
    %
    f8=figure(8);
    %{
    lg=loglog(spectraB_thi_2.kpar,spectraB_thi_2.P1Dpar/rms(spectraB_thi_2.P1Dpar,'all'),...
        spectraB_thi_3.kpar,spectraB_thi_3.P1Dpar/rms(spectraB_thi_3.P1Dpar,'all'),...
        spectraB_thi_4.kpar,spectraB_thi_4.P1Dpar/rms(spectraB_thi_4.P1Dpar,'all'),...
        spectraB_thi_5.kpar,spectraB_thi_5.P1Dpar/rms(spectraB_thi_5.P1Dpar,'all'),...
        spectraB_thi_6.kpar,spectraB_thi_6.P1Dpar/rms(spectraB_thi_6.P1Dpar,'all'));
   %}
    %{
        lg=loglog(spectraB_2.kpar,spectraB_2.P1Dpar/rms(spectraB_2.P1Dpar,'all'),...
        spectraB_3.kpar,spectraB_3.P1Dpar/rms(spectraB_3.P1Dpar,'all'),...
        spectraB_4.kpar,spectraB_4.P1Dpar/rms(spectraB_4.P1Dpar,'all'),...
        spectraB_5.kpar,spectraB_5.P1Dpar/rms(spectraB_5.P1Dpar,'all'),...
        spectraB_6.kpar,spectraB_6.P1Dpar/rms(spectraB_6.P1Dpar,'all'),...
        spectraB_7.kpar,spectraB_7.P1Dpar/rms(spectraB_7.P1Dpar,'all'),...
        spectraB_8.kpar,spectraB_8.P1Dpar/rms(spectraB_8.P1Dpar,'all'),...
        spectraB_9.kpar,spectraB_9.P1Dpar/rms(spectraB_9.P1Dpar,'all'),...
        spectraB_10.kpar,spectraB_10.P1Dpar/rms(spectraB_10.P1Dpar,'all'),...
        spectraB_11.kpar,spectraB_11.P1Dpar/rms(spectraB_11.P1Dpar,'all'),...
        spectraB_12.kpar,spectraB_12.P1Dpar/rms(spectraB_12.P1Dpar,'all'),...
        spectraB_13.kpar,spectraB_13.P1Dpar/rms(spectraB_13.P1Dpar,'all'));
    lg(1).LineWidth = 2; 
    lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; lg(9).LineWidth = 2;
    lg(10).LineWidth = 2; lg(11).LineWidth = 2; lg(12).LineWidth = 2;
        %}
     %
       lg=loglog(spectraB_double_2.kpar,spectraB_double_2.P1Dpar/rms(spectraB_double_2.P1Dpar,'all'),...
       spectraB_double_3.kpar,spectraB_double_3.P1Dpar/rms(spectraB_double_3.P1Dpar,'all'),...
       spectraB_double_4.kpar,spectraB_double_4.P1Dpar/rms(spectraB_double_4.P1Dpar,'all'),...
       spectraB_double_5.kpar,spectraB_double_5.P1Dpar/rms(spectraB_double_5.P1Dpar,'all'),...
       spectraB_double_6.kpar,spectraB_double_6.P1Dpar/rms(spectraB_double_6.P1Dpar,'all'),...
       spectraB_double_7.kpar,spectraB_double_7.P1Dpar/rms(spectraB_double_7.P1Dpar,'all'),...
       spectraB_double_8.kpar,spectraB_double_8.P1Dpar/rms(spectraB_double_8.P1Dpar,'all'),...
       spectraB_double_9.kpar,spectraB_double_9.P1Dpar/rms(spectraB_double_9.P1Dpar,'all'));
    lg(1).LineWidth = 2; 
    lg(2).LineWidth = 2; lg(3).LineWidth = 2;
    lg(4).LineWidth = 2; lg(5).LineWidth = 2; lg(6).LineWidth = 2;
    lg(7).LineWidth = 2; lg(8).LineWidth = 2; 
        %}
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
    %legend({'$t_{1}$','$t_{2}$','$t_{3}$','$t_{4}$','$t_{5}$',...
     %   '$t_{6}$','$t_{7}$','$t_{8}$','$t_{9}$','$t_{10}$',...
     %   '$t_{11}$','$t_{12}$'},'Interpreter','latex')
    %legend({'$t_{1}$','$t_{2}$','$t_{3}$','$t_{4}$','$t_{5}$'},'Interpreter','latex')
    legend({'$t_{1}$','$t_{2}$','$t_{3}$','$t_{4}$','$t_{5}$',...
        '$t_{6}$','$t_{7}$','$t_{8}$'},'Interpreter','latex')
    xlim([0.06 pi/0.06]); 
    %ylim([10^(-5) 100])
    %----------------------------------------------------------------------------------
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_{\|}$','Interpreter','latex','FontSize',20)
    hold off
    %}

    %#----------------------------------------------------------------------
   
    
    % Plots Magnetic field spectra
    %#----------------------------------------------------------------------    
    %axes('NextPlot','replacechildren', 'ColorOrder',C);
    
    f72=figure(72);
    %C = linspecer(10);
    %----------------------------------------------------------------------------------
    lg=loglog(spectraB_double_9.kper,spectraB_double_9.P1Dper/rms(spectraB_double_9.P1Dper,'all'),...
        spectraB_thi_6.kper,spectraB_thi_6.P1Dper/rms(spectraB_thi_6.P1Dper,'all'));
    %lg=loglog(spectraB_DATACB103_1.kper,spectraB_DATACB103_1.P1Dper/rms(spectraB_DATACB103_1.P1Dper,'all'),...
        %spectraB_beta01_1.kper,spectraB_beta01_1.P1Dper/rms(spectraB_beta01_1.P1Dper,'all')); 
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
    legend({'$16$','$8$'},'Interpreter','latex')
    %legend({'$\beta = 1$','$\beta = 0.1$'},'Interpreter','latex')
    xlim([0.35 pi/0.06]);  
    %ylim([10^(-5) 100])
    %}
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_\perp$','Interpreter','latex','FontSize',20)
    hold off

   
    f82=figure(82);
    %----------------------------------------------------------------------------------
    %
    lg=loglog(spectraB_double_9.kpar,spectraB_double_9.P1Dpar/rms(spectraB_double_9.P1Dpar,'all'),...
        spectraB_thi_6.kpar,spectraB_thi_6.P1Dpar/rms(spectraB_thi_6.P1Dpar,'all'));
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
    legend({'$16$','$8$'},'Interpreter','latex')
    xlim([0.06 pi/0.06]); 
    %ylim([10^(-5) 100])
    %}
    %----------------------------------------------------------------------------------
    xlabel('$k_{\|} d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D_{\|}$','Interpreter','latex','FontSize',20)
    hold off

    %#----------------------------------------------------------------------
    
    
    
    
    % Save the plots
    %-------------------------------------------------------------------------
    %cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/Analysis/Images';    
    saveas(f7,strcat(S(i).name,'1_perDspectrum_B_double.png'));
    saveas(f8,strcat(S(i).name,'1_parDspectrum_B_double.png'));
    saveas(f72,strcat(S(i).name,'1_perDspectrum_B_16_8.png'));
    saveas(f82,strcat(S(i).name,'1_parDspectrum_B_16_8.png'));
  
    
    
    
    
    
    

