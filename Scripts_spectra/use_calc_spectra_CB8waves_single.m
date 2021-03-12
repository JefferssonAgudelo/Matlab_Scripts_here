%#-----------------------------------------------------------------------
%C= 10^{-3}; beta =1
%#----------------------------------------------------------------------
%diary myDiaryFile
%cd '/disk/plasma2/jaa/CB8WAVES/CB8waves_04'; %Set the directory of the files 
%path = '/disk/plasma2/jaa/CB8WAVES/CB8waves_04';
%cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/ppc_50'; %This is important because the xdmf files are in that directory
%path = '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/ppc_50';

cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

%cd '/Volumes/PSC_DiRAC_DATA/RUNS2020/REAL_8_small_08_cor';
%path = '/Volumes/PSC_DiRAC_DATA/RUNS2020/REAL_8_small_08_cor';

S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
i=30; %3, 22, 19

    disp(strcat('Computing step ...',S(i).name)) %6 is 2000
    
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    %T_1st = info.Groups(1).Name;
    E_1st = info.Groups(17).Name;
    B_1st = info.Groups(18).Name;
    %J_1st = info.Groups(19).Name;
    V_1st = info.Groups(25).Name;
    n_1st = info.Groups(23).Name;
    %diary off
    %s2='/hx/p0/3d';
    %sx = strcat(B_1st,s2);
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    %
    Ex=h5read(fileID,strcat(E_1st,'/ex/p0/3d'));
    Ey=h5read(fileID,strcat(E_1st,'/ey/p0/3d'));
    Ez=h5read(fileID,strcat(E_1st,'/ez/p0/3d'));
    
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    
    vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'));
    
    ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    
    ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    %}
    %jx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    %jy=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    %jz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));
    
    %------------------------------------------------------------------
    %resx=18/296; %0.06081
    %resy=18/296;
    %resz=63/1024; %0.06152
    
    res=0.06;
    threshold=52;% 52;
    
    %res=13/400;
    %threshold=96;% 52; %100; % Check if by changing this the extension in the kper will be reduced
    
    %for 103
    %res=0.07;
    %threshold=44;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalization
    %B are normalized to the background magnetic field 
    B0=0.1;%0.07;%0.1;
    %n is normalise to the proton density n0=1
    n0=1;
    %Vi normalised to the proton alfven speed V_a=0.1c The alfven speed of
    %electrons is Vae=10Vai
    Vai=0.1; Vae=1;
    %E = vB, as B0=0.1 and v in units of c, E0= B0*Vai/c
    E0=0.01;
    %T=(beta/2)B^2, Thus, as beta=1 
    T0=0.005;
    %The current density J0=1?
    J0=Vai;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
    %cd '/disk/plasma2/jaa/CB8WAVES/Processing_data_Scripts/Matlab_Scripts/';
    cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here/Scripts_spectra';
    
    %spectraB=CalcSpectra_CB(Bx,By,Bz,res,threshold);
    %Btotal=sqrt(Bx.*Bx + By.*By + Bz.*Bz);
    %Bmean=mean(Btotal,'all');
    %clearvars Btotal 
    %spectraB_240=CalcSpectra_CB_local(Bx/Bmean,By/Bmean,Bz/Bmean,res,threshold);
    spectraB_120=CalcSpectra_CB_local(Bx/B0,By/B0,Bz/B0,res,threshold);
    clearvars Bx By Bz 
    %spectraB=CalcSpectra_CB_local(Bx,By,Bz,res,threshold);
    
    Etotal=sqrt(Ex.*Ex + Ey.*Ey + Ez.*Ez);
    Emean=mean(Etotal,'all');
    clearvars Etotal 
%    spectraE=CalcSpectra_CB_local(Ex/Emean,Ey/Emean,Ez/Emean,res,threshold);
    spectraE=CalcSpectra_CB_local(Ex/E0,Ey/E0,Ez/E0,res,threshold);
    clearvars Ex Ey Ez 
    
    Vitotal=sqrt(vix.*vix + viy.*viy + viz.*viz);
    Vimean=mean(Vitotal,'all');
    clearvars Vitotal
%    spectravi=CalcSpectra_CB_local(vix/Vimean,viy/Vimean,viz/Vimean,res,threshold);
    spectravi=CalcSpectra_CB_local(vix/Vai,viy/Vai,viz/Vai,res,threshold);
    clearvars vix viy viz
    %spectravi=CalcSpectra_CB_local(vix,viy,viz,res,threshold);
    
    Vetotal=sqrt(vex.*vex + vey.*vey + vez.*vez);
    Vemean=mean(Vetotal,'all');
    clearvars Vetotal
%    spectrave=CalcSpectra_CB_local(vex/Vemean,vey/Vemean,vez/Vemean,res,threshold);
    spectrave=CalcSpectra_CB_local(vex/Vai,vey/Vai,vez/Vai,res,threshold);
    clearvars vex vey vez
    
    nimean=mean(ni,'all');
    %spectrani=CalcSpectra_CB_n(ni/nimean,res,threshold);
    spectrani_0=CalcSpectra_CB_n(ni/n0,res,threshold);
    clearvars ni
    
    nemean=mean(ne,'all');
    %spectrane=CalcSpectra_CB_n(ne/nemean,res,threshold);
    spectrane=CalcSpectra_CB_n(ne/n0,res,threshold);
    clearvars ne
   
    %spectraJ=CalcSpectra_CB(Jx,Jy,Jz,res,threshold);
    %cd '/disk/plasma2/jaa/CB8WAVES/CB8waves_03' %I must bear this in mind
    %cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/';
    %Bz2=Bz-0.1;
    %spectraB2=CalcSpectra_CB(Bx,By,Bz2,res,threshold);
    
    % Making plots
    %----------------------------------------------------------------------
    [kp2J, kpJ]=meshgrid(spectraB_120.kpar,spectraB_120.kper);
    %P2D=spectraB_120.P2D./(kpJ');
    %spectraB=spectraB_nkper;
    
    spectraB240=spectraB_120;
    
    spectraB=spectraB_120;
    
    P2D=spectraB.P2D/max(max(spectraB.P2D));
    
    f91=figure(91);
    [C,h]=contourf(spectraB.kpar,spectraB.kper,log10(P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([-12 0]);
    %caxis([-25 -15]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}{P^{\mathbf{B}}_{2D}}$','Interpreter','latex')
    %title('$\log_{10}{P^{\mathbf{B}_{t}}_{2D}}$','Interpreter','latex')
    %x=log10([0.06 52]); y=(2/3)*x + 0.2;
    x2=log10([0.21 1]); y2=x2 + 1.05;
    x3=log10([0.06 0.2]); y3=(3/2)*x3 + 1.4;
    hold on
    %scatter(1./rpar2_1,1./rper2_1,'k') %time 120
    %hold on
    %scatter(1./rpar2_2,1./rper2_2,'k')
    %loglog(10.^x,10.^y,'--k','LineWidth',2)
    %loglog(10.^x2,10.^y2,'--k','LineWidth',2)
    loglog(10.^x3,10.^y3,'--k','LineWidth',2)
    ylim([0.28 52])
    xlim([0.06 52])
    %ylim([0.48 96])
    %txt = '\bf{3/2 \rightarrow}';
    txt = '$k_{\perp} \sim k_{\parallel}^{3/2}$';
    txt2 = '$k_{\perp} \sim k_{\parallel}$';
    %text(0.2,0.5,txt,'FontSize',15)
    text(0.1,0.4,txt,'Interpreter','latex','FontSize',20)
    %text(0.15,7,txt2,'Interpreter','latex','FontSize',20)
    yline(10,'--k','LineWidth',1);
    %yline(2.17,'--r','LineWidth',1);
    xline(1,'--k','LineWidth',1);
    %xline(12.5,'--k','LineWidth',1);
    hold off
    
    
    cd '/Volumes/PSC_DiRAC_DATA';
    %-------------------------------------------------------------------------
    saveas(f91,'CB104_log10_P2D_12_slope_3.png');
    
    f33=figure(33);
    [C,h]=contourf(spectraB.kpar,spectraB.kper,log10(P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([-12 0]);
    %caxis([-25 -15]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}{P^{\mathbf{B}}_{2D}}$','Interpreter','latex')
    %title('$\log_{10}{P^{\mathbf{B}_{t}}_{2D}}$','Interpreter','latex')
    x=log10([0.06 52]); y=(2/3)*x +0.2;
    x2=log10([0.06 52]); y2=x2 +0.2;
    x3=log10([0.06 52]); y3=(3/2)*x3 +0.2;
    hold on
    scatter(1./rpar2_1,1./rper2_1,6,'f','k') %time 120
    %hold on
    scatter(1./rpar2_2,1./rper2_2,6,'f','k')
    loglog(10.^x,10.^y,'--k','LineWidth',2)
    loglog(10.^x2,10.^y2,'--b','LineWidth',2)
    loglog(10.^x3,10.^y3,'--r','LineWidth',2)
    ylim([0.28 52])
    xlim([0.06 52])
    %ylim([0.48 96])
    %txt = '\bf{2/3 \rightarrow}';
    txt = '$k_{\perp} \sim k_{\parallel}^{3/2}$';
    text(1,30,txt,'Interpreter','latex','FontSize',20)
    tx2 = '$k_{\perp} \sim k_{\parallel}$';
    text(8,20,tx2,'Interpreter','latex','FontSize',20)
    tx2 = '$k_{\perp} \sim k_{\parallel}^{2/3}$';
    text(5.0,3,tx2,'Interpreter','latex','FontSize',20)
    %txt2 = '\bf{\leftarrow 3}';
    %text(0.2,0.5,txt,'FontSize',15)
    %text(5,28,txt2)
    %yline(13.65,'--k','LineWidth',1);
    %yline(2.17,'--r','LineWidth',1);
    %xline(2,'--r','LineWidth',1);
    %xline(12.5,'--k','LineWidth',1);
    hold off
    
    cd '/Volumes/PSC_DiRAC_DATA';
    %-------------------------------------------------------------------------
    %saveas(f1,'CB103_log10_P2D_1500.png');
    saveas(f33,'CB104_log10_P2D_2000_points.png');
    
    
    
    %{
    f12=figure(12);
    [C,h]=contourf(spectraB.kper,spectraB.kpar,log10(spectraE.P2D),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    ylabel('$k_\| d_i$','Interpreter','latex')
    xlabel('$k_\perp d_i$','Interpreter','latex')
    title('B','Interpreter','latex')
    x=log10([0.26 1]);
    x3=log10([1 70]);
    y=(2/3)*(x-0.314)+(1/5);
    y3=(1/3)*(x3-1)+(1/3);
    %y1=get(gca,'ylim'); %x1=10; %x2=get(gca,'xlim'); y2=10;
    %grid on
    hold on
    loglog(10.^(x),10.^y,'--k')
    loglog(10.^(x3),10.^y3,'--b')
    %loglog([x1 x1],y1,'--k') %loglog(x2,y2*ones(size(x)),'--k')
    %pbaspect([1 1 1])
    %%xlim([0.1 50])
    ylim([0.28 50])
    %legend('2/3','1/2')
    txt = '\bf{\leftarrow 2/3}';
    txt2 = '\bf{\downarrow 1/3}';
    text(0.8,0.7,txt)
    text(20,4,txt2)
    hold off
       
    
    f2=figure(2);
    [C,h]=contourf(spectraE.kpar,spectraE.kper,log10(spectraE.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('E','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    %{
    f31=figure(31);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(spectravi.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colorbar
    colormap(hot)
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$v_i$','Interpreter','latex')
    ylim([1 70])
    pbaspect([1 1 1])
    %}
    f3=figure(3);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(spectravi.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$v_i$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    %ylim([1 70])
    %pbaspect([1 1 1])
    hold off
    
    
    f4=figure(4);
    [C,h]=contourf(spectrave.kpar,spectrave.kper,log10(spectrave.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$v_e$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    
    f41=figure(41);
    [C,h]=contourf(spectrani.kpar,spectrani.kper,log10(spectrani.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    colormap(hot)
    colorbar
    set(gca,'XScale','log','YScale','log')
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$n_i$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 50])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    %}
    
    f51=figure(51);
    %loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper), 'k', spectraE.kper,spectraE.P1Dper/rms(spectraE.P1Dper), 'b',spectravi.kper,spectravi.P1Dper/rms(spectravi.P1Dper), 'r', spectrave.kper,spectrave.P1Dper/rms(spectrave.P1Dper), 'g', spectrani.kper,spectrani.P1Dper/rms(spectrani.P1Dper), 'c',spectrane.kper,spectrane.P1Dper/rms(spectrane.P1Dper), 'm')
    loglog(spectraB.kper,(spectraB.P1Dper./spectraB.kper)/rms(spectraB.P1Dper), 'k',spectravi.kper,(spectravi.P1Dper./spectravi.kper)/rms(spectravi.P1Dper), 'r')
    set(gca,'FontSize',18)
    %loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper), 'k',spectraJ.kper,spectraJ.P1Dper/rms(spectraJ.P1Dper), 'b')
    %x2=log10([1 16]); y2=x2.*(-3.2)+14.3; x3=log10([16 50]); y3=x3.*(-0.3)+10.8;
    x2=log10([1 10]); y2=x2.*(-2.33); x3=log10([0.099 1]); y3=x3.*(-1.7);
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x2,10.^y2,'--k')
    loglog(10.^x3,10.^y3,'--b')
    loglog([x33 x33],y1,'--')
    loglog([x44 x44],y1,'--')
    loglog([x55 x55],y1,'--')
    % Text in the plot
    txt3 = '$$ k_{\perp}d_e=1$$';
    txt4 = '$$ k_{\perp}d_i=1$$';
    txt5 = '$$ k_{\perp}^{-2.33}$$';
    txt6 = '$$ k_{\perp}^{-1.7}$$';
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    text(5 , 30,txt3,'Interpreter','latex','FontSize',15); text(1.1 , 30,txt4,'Interpreter','latex','FontSize',15); 
    text(1.1 , 0.1,txt5,'Interpreter','latex','FontSize',15);
    text(0.35 , 16,txt6,'Interpreter','latex','FontSize',15); text(16 , 0.3,txt8,'Interpreter','latex','FontSize',15);
    legend({'$P_{\tilde{B}}1D_\perp$','$P_{\tilde{v}_i}1D_\perp$'},'Interpreter','latex')
    xlabel('$k_\perp d_i$','Interpreter','latex')
    ylabel('$P_{\psi}1D_\perp$','Interpreter','latex')
    %xlim([0.3 52.36]); ylim([10^(-5) 100])
    hold off
    
    
    f61=figure(61);
    %loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k', spectraE.kpar,spectraE.P1Dpar/rms(spectraE.P1Dpar), 'b',spectravi.kpar,spectravi.P1Dpar/rms(spectravi.P1Dpar), 'r', spectrave.kpar,spectrave.P1Dpar/rms(spectrave.P1Dpar), 'g', spectrani.kpar,spectrani.P1Dpar/rms(spectrani.P1Dpar), 'c',spectrane.kpar,spectrane.P1Dpar/rms(spectrane.P1Dpar), 'm')
    loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k',spectravi.kpar,spectravi.P1Dpar/rms(spectravi.P1Dpar), 'r')
    %loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k',spectraJ.kpar,spectraJ.P1Dpar/rms(spectraJ.P1Dpar), 'b')
    set(gca,'FontSize',18)
    x=log10([0.099 1]);
    x2=log10([1 10]);
    x3=log10([5 30]);
    %y=x.*(-2.3)+12.2;
    y=x.*(-2)-1  ;%+12.2;
    y2=x2.*(-5);%+12.2;
    y3=x3.*(-0.3);%+10.6;
    x33=1;    x44=10;   x55=14.2857;
    xlim([0.3 52.36]); ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x,10.^y,'--b')
    %loglog(10.^x2,10.^y2,'--k')
    %loglog(10.^x3,10.^y3,'--')
    loglog([x33 x33],y1,'--')
    loglog([x44 x44],y1,'--')
    loglog([x55 x55],y1,'--')
    % Text in the plot
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$P_{\psi}1D_\|$','Interpreter','latex')
    txt3 = '$$ k_{||}d_e=1$$'; txt4 = '$$ k_{||}d_i=1$$';    txt5 = '$$ k_{||}^{-5}$$';
    txt6 = '$$ k_{||}^{-0.3}$$';    txt7 = '$$ k_{||}^{-2}$$';    txt8 = '$$ k_{||}\lambda_D=1$$';
    text(5 , 30,txt3,'Interpreter','latex','FontSize',15); text(1.1 , 30,txt4,'Interpreter','latex','FontSize',15)
    %text(2 , 0.4,txt5,'Interpreter','latex'); %text(18 , 1,txt6,'Interpreter','latex')
    text(0.5 , 7,txt7,'Interpreter','latex','FontSize',15); text(16 , 0.3,txt8,'Interpreter','latex','FontSize',15)
    %title('$P_{\psi}1D_\|$','Interpreter','latex')
    %legend({'$B$','$E$','$v_i$','$v_e$','$n_i$'},'Interpreter','latex') 
    %legend({'$\tilde{B}$','$\tilde{v}_i$'},'Interpreter','latex')
    %legend({'$\tilde{B}$','$\tilde{J}$'},'Interpreter','latex')
    legend({'$P_{\tilde{B}}1D_\|$','$P_{\tilde{v}_i}1D_\|$'},'Interpreter','latex')
    hold off
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f5=figure(5);
    %loglog(spectraB.kper,spectraB.P1Dper/Bper_norm, 'k',spectraB.kper,spectraB.P1Dper/Bper_norm, 'r')
    %loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper,'all'),'-k','LineWidth',2)
    loglog(spectraB.kper,spectraB.P1Dper,'-k','LineWidth',2)
    hold on
    %loglog(spectravi.kper,spectravi.P1Dper/rms(spectravi.P1Dper,'all'),'-r','LineWidth',2)
    loglog(spectravi.kper,spectravi.P1Dper,'-r','LineWidth',2)
    %loglog(spectraB.kper,spectraB.P1Dper/rms(spectraB.P1Dper,'all'), 'k',spectravi.kper,spectravi.P1Dper/rms(spectravi.P1Dper,'all'), 'r')
    set(gca,'FontSize',20)
    x2=log10([1 10]); y2=x2.*(-3)+0.4; x3=log10([0.099 1]); y3=x3.*(-1.7);
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    loglog(10.^x2,10.^y2,'--b')
    %loglog(10.^x3,10.^y3,'--b')
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    % Text in the plot
    txt3 = '$$ k_{\perp}d_e=1$$';
    txt4 = '$$ k_{\perp}d_i=1$$';
    txt5 = '$$ k_{\perp}^{-3}$$';
    %txt5 = '$$ k_{\perp}^{-2.33}$$'; this is the teorethical value for the
    %magnetic field
    txt6 = '$$ k_{\perp}^{-1.7}$$';
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    text(4.5 , 30,txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 30,txt4,'Interpreter','latex','FontSize',20); 
    text(1.1 , 0.2,txt5,'Interpreter','latex','FontSize',20);
    %text(0.35 , 16,txt6,'Interpreter','latex','FontSize',15); 
    text(16 , 0.3,txt8,'Interpreter','latex','FontSize',20);
    legend({'$P_{\tilde{B}}1D_\perp$','$P_{\tilde{v}_i}1D_\perp$'},'Interpreter','latex','FontSize',20)
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{\psi}}1D_\perp$','Interpreter','latex','FontSize',20)
    %xlim([0.26 52.36]); ylim([10^(-5) 100])
    hold off
    
    
    f6=figure(6);
    %loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar), 'k')
    loglog(spectraB.kpar,spectraB.P1Dpar/rms(spectraB.P1Dpar,'all'),'-k','LineWidth',2)
    hold on
    loglog(spectravi.kpar,spectravi.P1Dpar/rms(spectravi.P1Dpar,'all'),'-r','LineWidth',2)
    set(gca,'FontSize',20)
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
    loglog(10.^x,10.^y,'--b')
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    % Text in the plot
    txt3 = '$$ k_{||}d_e=1$$'; txt4 = '$$ k_{||}d_i=1$$';    txt5 = '$$ k_{||}^{-5}$$';
    txt6 = '$$ k_{||}^{-0.3}$$';    txt7 = '$$ k_{||}^{-2}$$';    txt8 = '$$ k_{||}\lambda_D=1$$';
    text(3 , 2,txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 30,txt4,'Interpreter','latex','FontSize',20)
    %text(2 , 0.4,txt5,'Interpreter','latex'); %text(18 , 1,txt6,'Interpreter','latex')
    text(0.3 , 2,txt7,'Interpreter','latex','FontSize',20); 
    text(10.5 , 0.3,txt8,'Interpreter','latex','FontSize',20)
    xlabel('$k_\| d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{\psi}}1D_\|$','Interpreter','latex','FontSize',20)
    legend({'$P_{\tilde{B}}1D_\|$','$P_{\tilde{v}_i}1D_\|$'},'Interpreter','latex','FontSize',20)
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Including the paricle density
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kper=spectraB.kper;%+spectraB.kper(2); % to correct the zero
    kpar=spectraB.kpar;%+spectraB.kpar(2);
    kper0=kper(2);
    
    N1=400; N2 = 400 ; N3 =2016;
    P1D_perB=(spectraB.P1Dper)./(max(spectraB.P1Dper));
    P1D_pervi=(spectravi.P1Dper)./(max(spectravi.P1Dper));
    P1D_perni=(spectrani.P1Dper)./(max(spectrani.P1Dper));
    %P1D_perni_0=(spectrani_0.P1Dper)./(max(spectrani_0.P1Dper));
    
    f55=figure(55);
    loglog(kper,P1D_perB,'-k','LineWidth',2)
    hold on
    loglog(kper,P1D_pervi, '-r','LineWidth',2) 
    loglog(kper,P1D_perni, '-b','LineWidth',2)
    %loglog(kper,P1D_perni_0, '--','LineWidth',1)
    set(gca,'FontSize',20)
    x2=log10([1.8 7]);     y2=x2.*(-3)-0.24; 
    x3=log10([7 20]);     y3=x3.*(-4)+0.5;
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    loglog(10.^x2,10.^y2,'--b','LineWidth',2)
    loglog(10.^x3,10.^y3,'--b','LineWidth',2)
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    % Text in the plot
    txt3 = '$$ k_{\perp}d_e=1$$';
    txt4 = '$$ k_{\perp}d_i=1$$';
    txt5 = '$$ k_{\perp}^{-3}$$';
    txt7 = '$$ k_{\perp}^{-4}$$'; %this is the teorethical value for the
    %magnetic field
    txt6 = '$$ k_{\perp}^{-1.7}$$';
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    text(4.5 , 10e-5,txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 10e-5,txt4,'Interpreter','latex','FontSize',20); 
    text(1.1 , 10e-5,txt5,'Interpreter','latex','FontSize',20);
    text(1.1 , 10e-5,txt7,'Interpreter','latex','FontSize',20);
    %text(0.35 , 16,txt6,'Interpreter','latex','FontSize',15); 
    text(16 , 10e-5,txt8,'Interpreter','latex','FontSize',20);
    %legend({'$P^{\tilde{B}}_{1D_\perp}$','$P^{\tilde{v}_i}_{1D_\perp}$',...
     %   '$P^{\tilde{n}_i}_{1D_\perp}$'},'Interpreter','latex','FontSize',20) 
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P^{\psi}_{1D_\perp}$','Interpreter','latex','FontSize',20)
    xlim([0.26 52.36]); ylim([10^(-5) 10^0])
    legend({'$\mathbf{B}$','$\mathbf{v}_{i}$','$n_{i}$'},'Interpreter','latex','FontSize',20)
    hold off
    
    P1D_parB=(spectraB.P1Dpar)/(max(spectraB.P1Dpar));
    P1D_parvi=(spectravi.P1Dpar)/(max(spectravi.P1Dpar));
    P1D_parni=(spectrani.P1Dpar)/(max(spectrani.P1Dpar));
    %P1D_parni_0=(spectrani_0.P1Dpar)/(max(spectrani_0.P1Dpar));
    
    f66=figure(66);
    loglog(kpar,P1D_parB,'-k','LineWidth',2)
    hold on
    loglog(kpar,P1D_parvi, '-r','LineWidth',2) 
    loglog(kpar,P1D_parni, '-b','LineWidth',2)
    %loglog(kpar,P1D_parni_0, '--','LineWidth',1)
    set(gca,'FontSize',20)
    %x=log10([0.05 1]);    y=x.*(-2) + -3;%+12.2;     %y=x.*(-2.3)+12.2;
    x=log10([0.1 0.3]);    y=x.*(-2) + -2.8;
    x2=log10([0.4 2]);    y2=x2.*(-2.5)-3;%+12.2;
    x3=log10([2 4]);    y3=x3.*(-4)-2.5;
    x33=1;    x44=10;   x55=14.2857;
    y1=get(gca,'ylim');
    loglog(10.^x,10.^y,'--b', 'LineWidth',2)
    loglog(10.^x2,10.^y2,'--b', 'LineWidth',2)
    loglog(10.^x3,10.^y3,'--b', 'LineWidth',2)
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    % Text in the plot
    txt3 = '$$ k_{||}d_e=1$$'; txt4 = '$$ k_{||}d_i=1$$';   
    txt5 = '$$ k_{||}^{-2}$$';
    txt6 = '$$ k_{||}^{-2.5}$$';    
    txt7 = '$$ k_{||}^{-4}$$';   
    txt8 = '$$ k_{||}\lambda_D=1$$';
    text(3 , 10^(-5),txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 10^(-5),txt4,'Interpreter','latex','FontSize',20)
    text(2 , 10^(-5),txt5,'Interpreter','latex','FontSize',20);
    text(18 , 10^(-5),txt6,'Interpreter','latex','FontSize',20)
    text(0.3 , 10^(-5),txt7,'Interpreter','latex','FontSize',20); 
    text(10.5 , 10^(-5),txt8,'Interpreter','latex','FontSize',20)
    xlabel('$k_\| d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P^{\psi}_{1D_\|}$','Interpreter','latex','FontSize',20)
%    legend({'$P_{\tilde{B}}1D_\|$','$P_{\tilde{v}_i}1D_\|$','$P_{\tilde{n}_i}1D_\|$'},'Interpreter','latex','FontSize',20)
    legend({'$\mathbf{B}$','$\mathbf{v}_{i}$','$n_{i}$'},'Interpreter','latex','FontSize',20)
    xlim([0.05 52.36]); ylim([10^(-6) 10^(0)]);
    hold off
        
    f555=figure(555);
    loglog(spectrani.kper,(spectrani.P1Dper/rms(spectrani.P1Dper,'all'))./(spectraB.P1Dper/rms(spectraB.P1Dper,'all')), 'k', ...
           spectrani.kper,(spectrani.P1Dper/rms(spectrani.P1Dper,'all'))./(spectravi.P1Dper/rms(spectravi.P1Dper,'all')), 'r',...
           spectravi.kper,(spectravi.P1Dper/rms(spectravi.P1Dper,'all'))./(spectraB.P1Dper/rms(spectraB.P1Dper,'all')), 'm')
    set(gca,'FontSize',20)
    x2=log10([1 10]); y2=x2.*(-3)+0.4; x3=log10([0.099 1]); y3=x3.*(-1.7);
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    hold on
    %loglog(10.^x2,10.^y2,'--k')
    %loglog(10.^x3,10.^y3,'--b')
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    % Text in the plot
    txt3 = '$$ k_{\perp}d_e=1$$';
    txt4 = '$$ k_{\perp}d_i=1$$';
    txt5 = '$$ k_{\perp}^{-3}$$';
    %txt5 = '$$ k_{\perp}^{-2.33}$$'; this is the teorethical value for the
    %magnetic field
    txt6 = '$$ k_{\perp}^{-1.7}$$';
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    text(4.5 , 30,txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 30,txt4,'Interpreter','latex','FontSize',20); 
    text(16 , 0.3,txt8,'Interpreter','latex','FontSize',20);
    legend({'$\frac{P_{\tilde{n}_{i}}1D_\perp}{P_{\tilde{B}}1D_\perp}$',...
            '$\frac{P_{\tilde{n}_{i}}1D_\perp}{P_{\tilde{v}_{i}}1D_\perp}$',...
            '$\frac{P_{\tilde{v}_{i}}1D_\perp}{P_{\tilde{B}}1D_\perp}$'},'Interpreter','latex','FontSize',20)
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$\frac{P_{\tilde{\psi}_{1}}1D_\perp}{P_{\tilde{\psi}_{2}}1D_\perp}$','Interpreter','latex','FontSize',20)
    xlim([0.26 52.36]); %ylim([10^(-5) 100])
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Including the paricle density
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f5555=figure(5555);
    loglog(spectraB.kper,(spectraB.P1Dper/max(spectraB.P1Dper)),'LineWidth',2)
    hold on
    loglog(spectraE.kper,(spectraE.P1Dper/max(spectraE.P1Dper)),'LineWidth',2)
    loglog(spectravi.kper,(spectravi.P1Dper/max(spectravi.P1Dper)),'LineWidth',2)
    loglog(spectrave.kper,(spectrave.P1Dper/max(spectrave.P1Dper)),'LineWidth',2) 
    loglog(spectrani.kper,(spectrani.P1Dper/max(spectrani.P1Dper)),'LineWidth',2)
    loglog(spectrane.kper,(spectrane.P1Dper/max(spectrane.P1Dper)),'LineWidth',2)
    set(gca,'FontSize',20)
    x2=log10([1 10]); 
    y2=x2.*(-3)-0.2; 
    x3=log10([0.099 1]); 
    y3=x3.*(-1.7);
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x2,10.^y2,'--b')
    %loglog(10.^x3,10.^y3,'--b')
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    legend({'$\mathbf{B}$','$\mathbf{E}$','$\mathbf{v}_{i}$','$\mathbf{v}_{e}$',...
        '$n_{i}$','$n_{e}$'},'Interpreter','latex','FontSize',20)
    % Text in the plot
    txt3 = '$$ k_{\perp}d_e=1$$';
    txt4 = '$$ k_{\perp}d_i=1$$';
    txt5 = '$$ k_{\perp}^{-3.3}$$';
    %txt5 = '$$ k_{\perp}^{-2.33}$$'; this is the teorethical value for the
    %magnetic field
    txt6 = '$$ k_{\perp}^{-1.7}$$';
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    text(4.5 , 10^(-5),txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 10^(-5),txt4,'Interpreter','latex','FontSize',20); 
    text(1.1 , 10^(-5),txt5,'Interpreter','latex','FontSize',20);
    %text(0.35 , 16,txt6,'Interpreter','latex','FontSize',15); 
    text(16 , 10^(-5),txt8,'Interpreter','latex','FontSize',20);
    %legend({'$P^{\tilde{B}}_{1D_\perp}$','$P^{\tilde{v}_i}_{1D_\perp}$',...
     %   '$P^{\tilde{n}_i}_{1D_\perp}$'},'Interpreter','latex','FontSize',20) 
    xlabel('$k_\perp d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P^{\psi}_{1D_\perp}$','Interpreter','latex','FontSize',20)
    xlim([0.26 52.36]); ylim([10^(-5) 1])
    hold off
    
    f6666=figure(6666);
    loglog(spectraB.kpar,(spectraB.P1Dpar/max(spectraB.P1Dpar)),'LineWidth',2)
    hold on
    loglog(spectraE.kpar,(spectraE.P1Dpar/max(spectraE.P1Dpar)),'LineWidth',2)
    loglog(spectravi.kpar,(spectravi.P1Dpar/max(spectravi.P1Dpar)),'LineWidth',2) 
    loglog(spectrave.kpar,(spectrave.P1Dpar/max(spectrave.P1Dpar)),'LineWidth',2) 
    loglog(spectrani.kpar,(spectrani.P1Dpar/max(spectrani.P1Dpar)),'LineWidth',2)
    loglog(spectrane.kpar,(spectrane.P1Dpar/max(spectrane.P1Dpar)),'LineWidth',2)
    set(gca,'FontSize',20)
    x=log10([0.099 1]);
    x2=log10([1 10]);
    x3=log10([5 30]);
    %y=x.*(-2.3)+12.2;
    y=x.*(-2)-3  ;%+12.2;
    y2=x2.*(-5);%+12.2;
    y3=x3.*(-0.3);%+10.6;
    x33=1;    x44=10;   x55=14.2857;
    xlim([0.05 52.36]); ylim([10^(-5) 1]);
    y1=get(gca,'ylim');
    loglog(10.^x,10.^y,'--b')
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    % Text in the plot
    txt3 = '$$ k_{||}d_e=1$$'; txt4 = '$$ k_{||}d_i=1$$';    txt5 = '$$ k_{||}^{-5}$$';
    txt6 = '$$ k_{||}^{-0.3}$$';    txt7 = '$$ k_{||}^{-2}$$';    txt8 = '$$ k_{||}\lambda_D=1$$';
    text(3 , 10^(-5),txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 10^(-5),txt4,'Interpreter','latex','FontSize',20)
    %text(2 , 0.4,txt5,'Interpreter','latex'); %text(18 , 1,txt6,'Interpreter','latex')
    text(0.3 , 10^(-5),txt7,'Interpreter','latex','FontSize',20); 
    text(10.5 , 10^(-5),txt8,'Interpreter','latex','FontSize',20)
    xlabel('$k_\| d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P^{\psi}_{1D_\|}$','Interpreter','latex','FontSize',20)
    legend({'$\mathbf{B}$','$\mathbf{E}$','$\mathbf{v}_{i}$','$\mathbf{v}_{e}$',...
        '$n_{i}$','$n_{e}$'},'Interpreter','latex','FontSize',20)
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Bper_norm=sqrt(sum(spectraB.P1Dper.*spectraB.P1Dper, 'all'));
    Bpar_norm=sqrt(sum(spectraB.P1Dpar.*spectraB.P1Dpar, 'all'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Interpolation to compare the spectrums parallel and perpendicular in
    % the same plot. So far it ends up as a plot where the x label is of
    % the kpar lenght but the points kdi corresponds to kperpdi and not
    % kparalleldi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tlinkper=spectraB.kper;
    xlinkper=spectraB.P1Dper/rms(spectraB.P1Dper,'all');
    tlinkpar=spectraB.kpar;
    xlinkpar=spectraB.P1Dpar/rms(spectraB.P1Dpar,'all');
    xlinkper_2=interp1(tlinkper, xlinkper, tlinkpar, 'linear','extrap');% This is the linear interpolation on kper based
    xlinkpar_2=interp1(tlinkpar, xlinkpar, tlinkper, 'linear','extrap');% This is the linear interpolation
         
    f7=figure(7);
    loglog(tlinkpar,xlinkpar, 'b',tlinkpar, xlinkper_2, 'r')  
    set(gca,'FontSize',20)
    x33=1;    x44=10;   x55=14.2857;
    x2=log10([1 10]); 
    y2=x2.*(-3)+0.4; %perpendicular
    x3=log10([0.099 1.85]); 
    y3=x3.*(-2)-1; %parallel
    x2p=log10([1.85 3]);
    y2p=x2p.*(-5)+0.29;
    xlim([0.05 53.36]); %ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    loglog(10.^x2,10.^y2,'--r')
    loglog(10.^x3,10.^y3,'--c')
    %loglog(10.^x2p,10.^y2p,'--m')
    %yline(0.0015,'--c')
    %xline(0.05,'--b')
    xline(1.85,'--c')
    xline(3.14,'--r')
    loglog([x33 x33],y1,'--k')
    loglog([x44 x44],y1,'--k')
    loglog([x55 x55],y1,'--k')
    txt3 = '$$ k_{\perp}d_e=1$$'; 
    txt4 = '$$ k_{\perp}d_i=1$$';  
    txt8 = '$$ k_{\perp}\lambda_D=1$$';
    txt9 = '$$ k_{\parallel}d_e=1$$';
    txt5 = '$$ k_{\perp}^{-3}$$';
    txt7 = '$$ k_{||}^{-2}$$';   
    text(3.5 , 2,txt3,'Interpreter','latex','FontSize',20); 
    text(1.1 , 30,txt4,'Interpreter','latex','FontSize',20)
    text(16 , 0.3,txt8,'Interpreter','latex','FontSize',20)
    text(7.5 , 2,txt9,'Interpreter','latex','FontSize',20);
    text(1.1 , 0.2,txt5,'Interpreter','latex','FontSize',20);
    text(0.3 , 2,txt7,'Interpreter','latex','FontSize',20); 
    xlabel('$k_{\perp}d_i$','Interpreter','latex','FontSize',20)
    ylabel('$P_{\tilde{B}}1D$','Interpreter','latex','FontSize',20)
    legend({'$P_{\tilde{B}}1D_\|$','$P_{\tilde{B}}1D_{\perp}$'},'Interpreter','latex','FontSize',20)
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % This part it to plot the ratio between Bper and Bpar and vice-versa
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bperoverpar=xlinkper_2.'./xlinkpar;
    bparoverper=xlinkpar_2./xlinkper;
    
    f8=figure(8);
    loglog(tlinkpar,bperoverpar)
 
    f9=figure(9);
    loglog(tlinkper,bparoverper)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    
    %----------------------------------------------------------------------
    %cd /disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics/2D_1D_PSD
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images';
    %cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/Images';
    % Save the plots
    cd '/Volumes/PSC_DiRAC_DATA';
    %-------------------------------------------------------------------------
    saveas(f3,'log10_120_Bperf.png');
    saveas(f666,'timesnol_12.png');
    saveas(f667,'timelin_12.png');
    
    saveas(f1,strcat(S(i).name,'2Dspectrum_Bpd_Isotropic_ppc50.png'));
    saveas(f5,strcat(S(i).name,'1DperBvi.png'));
    saveas(f55,strcat(S(i).name,'1DperBvini_120.png'));
    saveas(f555,strcat(S(i).name,'1DperBvini_ratios.png'));
    saveas(f6,strcat(S(i).name,'1DparBvi.png'));
    saveas(f66,strcat(S(i).name,'1DparBvini_120.png'));
    saveas(f7,strcat(S(i).name,'Bparper.png'));
    
    saveas(f5555,strcat(S(i).name,'1DperAll_120.png'));
    saveas(f6666,strcat(S(i).name,'1DparAll_120.png'));
    
    %saveas(f1,strcat(S(i).name,'2Dspectrum_Bpd_CB04_not.','epsc'));
    %saveas(f12,strcat(S(i).name,'2Dspectrum_B_trans_CB04.png'));
    %saveas(f2,strcat(S(i).name,'2Dspectrum_E_CB04.png'));
    %saveas(f3,strcat(S(i).name,'2Dspectrum_vi_CB04.png'));
    %saveas(f4,strcat(S(i).name,'2Dspectrum_ve_CB04.png'));
    %saveas(f41,strcat(S(i).name,'2Dspectrum_ni_CB04.png'));
    %saveas(f5,strcat(S(i).name,'1_perDspectrum_CB04_Bvipd_2.png'));
    %saveas(f6,strcat(S(i).name,'1_parDspectrum_CB04_Bvipd_2.png'));
    %saveas(f7,strcat(S(i).name,'1_parDspectrum_CB04_Bparper.png'));    
    %cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'; 
    cd '/Volumes/PSC_DiRAC_DATA/DATACB104';
   
