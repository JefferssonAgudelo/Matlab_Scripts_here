%#-----------------------------------------------------------------------
%C= 10^{-3}; beta =1
%#----------------------------------------------------------------------
%diary myDiaryFile


cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

%load ('/Volumes/PSC_DiRAC_DATA/nonlineartime_images/Vi_spectrum_fourier.mat')
load ('/Volumes/PSC_DiRAC_DATA/nonlineartime_images/Vi_cell_fourier_spectrum.mat')
S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
i=1;%19; %3, 22, 19
for i=2:40
    disp(strcat('Computing step ...',S(i).name)) %6 is 2000
    
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    V_1st = info.Groups(25).Name;
    %n_1st = info.Groups(23).Name;
    %diary off
    %s2='/hx/p0/3d';
    %sx = strcat(B_1st,s2);
    B_1st = info.Groups(18).Name;
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));

    %vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    %viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    %viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));

    
    res=0.06;
    threshold=52; %100; % Check if by changing this the extension in the kper will be reduced
    
    %cd '/disk/plasma2/jaa/CB8WAVES/Processing_data_Scripts/Matlab_Scripts/';
    cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
    spectraBi=CalcSpectra_CB(Bx,By,Bz,res,threshold);
    %spectravi=CalcSpectra_CB(vix,viy,viz,res,threshold);
    
    %This is not the best practise
    %predir=string(i);
    %assignin('base',strcat('spectravi_',predir), spectravi)

    %spectravis{i} =spectravi;
    spectraBis{i} =spectraBi;
    cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
    %path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
    

    
    %-------------------------------------------------------
    % To calculate the nonlinear time at sub-proton scales
    %$1/k^{2}_{\perp}B_{k}$
    %spectravi=spectravi_20;
    
    spectravi=spectravis{20};
    
    % one dimensional plots
    %----------------------------------------------------------------------
    %{
    tl1=1./(spectravi.kpar * 0.1);
    tnl_2 = 1./((spectravi.kper * 0.1) .* sqrt((spectravi.P1Dper / rms(spectravi.P1Dper)) .* (spectravi.kper) ));
    tnl_3 = 1./(spectravi.kper * 0.1 );
    tnl_4 = ((10^(4/3))/0.1 ).*((spectravi.kper).^(-2/3));
    
    f666=figure(666);
    loglog(spectravi.kper,tnl_2, '*-b')
    set(gca,'FontSize',18)
    y1=get(gca,'ylim');
    hold on
%    loglog(spectravi.kper, tnl_1, '*-k')
    loglog(spectravi.kper, tnl_4, '*-r')
    xline(1,'--k');
    yline(2.29,'--k');
    xlim([0.1 53.36]); ylim([10^(-1) 1000]);
%    legend({'$1/k^{2}_{\perp}v_{i,k}$','$1/k_{\perp}v_{i,k}$','$1/k_{\perp}$'},'Interpreter','latex')
    legend({'$1/k_{\perp}v_{i,k}$','$1/k_{\perp}$'},'Interpreter','latex')
    xlabel('$k_\perp d_i$','Interpreter','latex')
    ylabel('$\tau_{nl}\omega_{pi}$','Interpreter','latex')
    hold off
    
    f667=figure(667);
    loglog(spectravi.kpar,tl1, '*-b')
    set(gca,'FontSize',18)
    y1=get(gca,'ylim');
    hold on
    loglog(spectravi.kper, tnl_4, '*-r')
    loglog(spectravi.kper, tnl_2, '*-k')
    xline(1,'--k');
    xline(0.0001,'--k');
    yline(2.29,'--k');
    yline(208,'--k');
    %xlim([0.1 53.36]); ylim([10^(-1) 1000]);
%    legend({'$1/k^{2}_{\perp}v_{i,k}$','$1/k_{\perp}v_{i,k}$','$1/k_{\perp}$'},'Interpreter','latex')
    legend({'$1/k_{\parallel}v_{i,A}$','$10^{4/3}(k_{\perp}d_{i})^{-2/3}$','$1/k_{\perp}v_{i,k}$'},'Interpreter','latex')
    xlabel('$k d_i$','Interpreter','latex')
    ylabel('$\tau\omega_{pi}$','Interpreter','latex')
    hold off

    
    f668=figure(668);
    loglog(spectravi.kper,spectravi.P1Dper / rms(spectravi.P1Dper), '*-b')
    hold on
    x2=log10([1 10]); y2=x2.*(-3.4); x3=log10([0.099 1]); y3=x3.*(-1.7);
    x33=1;    x44=10;    x55=14.2857;
    y1=get(gca,'ylim');
    loglog(10.^x2,10.^y2,'--k')
    plot(spectravi.kper, -3*spectravi.kper)
    hold off
    %{
    f667=figure(667);
    loglog(spectraB.kpar, tl1, '*-r')
    set(gca,'FontSize',18)
    y1=get(gca,'ylim');
    hold on
    %loglog(squeeze( Bz(1,1,:)-Bz2(1,1,:)), 'r')
%    loglog(spectravi.kper, spectravi.P1Dper)
    xlabel('$k_{\parallel} d_i$','Interpreter','latex')
    ylabel('$\tau_{l}\omega_{pi}$','Interpreter','latex')
    xline(1,'--k');
    yline(10.13,'--k');
    legend({'$1/k_{\parallel}V_{A}$'},'Interpreter','latex')
    %xlim([0.3 52.36]); %ylim([10^(-5) 100])
    hold off
    %}
    
    cd '/Volumes/PSC_DiRAC_DATA';
    %-------------------------------------------------------------------------
    saveas(f666,strcat(S(i).name,'time_nolin.png'));
    %saveas(f667,'timelin_12.png');
    cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
    
    %-------------------------------------------------------    
    %}

    %----------------------------------------------------------------------
    % Two dimensional plots
    %----------------------------------------------------------------------
    tl_12 = (1./(spectravi.kpar * 0.1))'.* ones(size(spectravi.P2D));
    tnl_k3 = ones(size(spectravi.P2D)) .* (1./(sqrt(spectravi.P2D .* spectravi.kper) .* (spectravi.kper.^3)));
    
    tnl_22 = ones(size(spectravi.P2D)) .* (1./((spectravi.kper * 0.1) .* sqrt((spectravi.P1Dper / rms(spectravi.P1Dper)) .* (spectravi.kper) )));
    tnl_52 = ones(size(spectravi.P2D)) .* (1./(sqrt((spectravi.P2D / rms(spectravi.P2D)) .* (spectravi.kper) ) .* (spectravi.kper * 0.1)));

    tnlk3_52 = tnl_k3 - tnl_52;
    tnlk3_tl_12 = tnl_k3 ./ tl_12;
    
    tnl_42 = ones(size(spectravi.P2D)) .* (((10^(4/3))/0.1 ).*((spectravi.kper).^(-2/3)));
    tnl_62 = ones(size(spectravi.P2D)) .* (((10^(4/3))/0.1 ).*((spectravi.kper).^(-1/3)));
    
    tnl_22_over_tl_12= tnl_22 ./ tl_12;
    tnl_52_over_tl_12= tnl_52 ./ tl_12;
    
    tnl_42_over_tl_12= tnl_42 ./ tl_12;
    tnl_62_over_tl_12= tnl_62 ./ tl_12;
    
    k23=(10^(-4/3))*(10)^(2/3);
    k13=(10^(-4/3))*(10)^(1/3);

    cd /Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here
    f1=figure(1);
    ax(1)=subplot(2,2,1);
    %[C,h]=contourf(spectravi.kpar,spectravi.kper,log10(tnl_42_over_tl_12'),'LevelStep',0.2);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(tnlk3_tl_12'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([-3 3]);
    colormap(ax(1),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{nl}/\tau_{l} \ theoric_{2/3}$','Interpreter','latex')
    hold on
    x=log10([0.1 50]);
    y=(3/2)*x;
    loglog(10.^x,10.^y,'--k')
    ylim([0.28 52])
    txt = '\bf{3/2 \rightarrow}';
    text(0.1,0.5,txt)
    ylim([0.28 52])
    hold off

    ax(2)=subplot(2,2,2);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(tnl_62_over_tl_12'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([-3 3]);
    colormap(ax(2),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{nl}/\tau_{l} \ theoric_{1/3}$','Interpreter','latex')
    hold on
    x3=log10([0.1 70]);
    y3=(3)*x3;
    loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt2 = '\bf{\leftarrow 3}';
    text(5,28,txt2)
    ylim([0.28 52])
    hold off

%
    subplot(2,2,3)
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(tnl_22_over_tl_12'),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([0 13]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{nl}/\tau_{l} \ 1D$','Interpreter','latex')
    hold on
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 52])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    ylim([0.28 52])
    yline(10,'--k');
    %xline(1,'--k'); %xline(10,'--k'); xline(k23,'-.k'); xline(k13,'-k');
    hold off
%}
    ax(3)=subplot(2,2,3);
    yyy=tnl_52_over_tl_12';
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(yyy),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    yyy(isinf(yyy)|isnan(yyy)) = 1e-6;
    yyy=log10(yyy);
    maxy=max(max(yyy));
    caxis([-maxy maxy]);
    colormap(ax(3),'redblue')
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}\tau_{nl}/\tau_{l} \ 2D $','Interpreter','latex')
    hold on
    %x=log10([0.06 k23]);
    x=log10([0.06 30]);
    x3=log10([k23 1]);
    y=(3/2)*(x)+1;
    y3=(3)*(x3)+2;
    loglog(10.^x,10.^y,'--k'); loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt = '\bf{\leftarrow 3/2 }';
    txt2 = '\bf{3 \rightarrow}';
    text(0.25,0.5,txt)
    text(0.3,28,txt2)
    ylim([0.28 52])
    yline(10,'--k');
    %xline(1,'--k'); %xline(10,'--k'); xline(k23,'-.k'); xline(k13,'-k');
    hold off
    
    ax(4)=subplot(2,2,4);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(spectravi.P2D'),'LevelStep',0.5);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([0 13]);
    colormap(ax(4),jet)
    colorbar;
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$\log_{10}{P_{\mathbf{v}_{i}}2D}, \ \omega_{pi}t=12$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    hold on
    ylim([0.28 52])
    yline(10,'--k');
    xline(1,'--k'); %xline(10,'--k');
    %xline(k23,'-.k'); xline(k13,'-k');   
    hold off
    
    
    cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
    path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
    load ('/Volumes/PSC_DiRAC_DATA/nonlineartime_images/Vi_cell_fourier_spectrum.mat')
    
    LL=40;
    SP2Dn_nl=zeros(LL,1);
    SP2Dn_l=zeros(LL,1);
    SP2Dn_CB=zeros(LL,1);
    
    % Multiplying the two plots 
    for i=1:LL
        %i=20;
        spectravi=spectravis{i};
    
        %----------------------------------------------------------------------
        % Two dimensional plots
        %----------------------------------------------------------------------
        tl_12 = (1./(spectravi.kpar * 0.1))'.* ones(size(spectravi.P2D));
        tnl_22 = ones(size(spectravi.P2D)) .* (1./((spectravi.kper * 0.1) .* sqrt((spectravi.P1Dper / rms(spectravi.P1Dper)) .* (spectravi.kper) )));
        tnl_52 = ones(size(spectravi.P2D)) .* (1./(sqrt((spectravi.P2D / rms(spectravi.P2D)) .* (spectravi.kper) ) .* (spectravi.kper * 0.1)));
        tnl_42 = ones(size(spectravi.P2D)) .* (((10^(4/3))/0.1 ).*((spectravi.kper).^(-2/3)));
        tnl_62 = ones(size(spectravi.P2D)) .* (((10^(4/3))/0.1 ).*((spectravi.kper).^(-1/3)));    
        tnl_22_over_tl_12= tnl_22 ./ tl_12;
        tnl_52_over_tl_12= tnl_52 ./ tl_12;
        tnl_42_over_tl_12= tnl_42 ./ tl_12;
        tnl_62_over_tl_12= tnl_62 ./ tl_12;
        
        k23=(10^(-4/3))*(10)^(2/3);
        k13=(10^(-4/3))*(10)^(1/3);
        %{
    f2=figure(2);
    %[C,h]=contourf(spectravi.kpar,spectravi.kper,log10(tnl_52_over_tl_12'.*spectravi.P2D'),'LevelStep',0.2);
    [C,h]=contourf(spectravi.kpar,spectravi.kper,abs(log10(tnl_52_over_tl_12').*log10(spectravi.P2D')),'LevelStep',0.2);
    set(h,'LineColor','none')
    lim = caxis; 
    %caxis([-3 3]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    %title('$\log_{10} \left(\tau_{nl}/\tau_{l} \cdot P_{\mathbf{v}_{i}} \right)$','Interpreter','latex')
    title('$|\log_{10} \tau_{nl}/\tau_{l} \cdot \log_{10} P_{\mathbf{v}_{i}} |$','Interpreter','latex')
    hold on
    %x=log10([0.06 k23]);
    x=log10([0.06 30]);
    x3=log10([k23 1]);
    y=(3/2)*(x)+1;
    y3=(3)*(x3)+2;
    loglog(10.^x,10.^y,'--k'); loglog(10.^x3,10.^y3,'--k')
    ylim([0.28 52])
    txt = '\bf{\leftarrow 3/2 }';
    txt2 = '\bf{3 \rightarrow}';
    text(0.25,0.5,txt)
    text(0.3,28,txt2)
    ylim([0.28 52])
    yline(10,'--k');
    %xline(1,'--k'); %xline(10,'--k'); xline(k23,'-.k'); xline(k13,'-k');
    hold off
    %}
        
        %------------------------------------------------------------------
        % This is how to define the linear and non-linear times
        %------------------------------------------------------------------
        % The linear time is tl = 1/w where omega is the frequency of the
        % waves that are solutions to the disspersion relation. There is a
        % different definition for each type of wave. Alfven or fast waves.
        
        
       % For Alfven waves in the inertial range.
       
       tl_aw = 1 / k_par  
        
        
        %------------------------------------------------------------------
        
        %Integrating to see the amount of power that is above and below the
        %zero line
        %------------------------------------------------------------------
        tnl_52_over_tl_12=tnlk3_tl_12; % This doesn't work straight
        %-------------------------------------------------
        P2D=spectravi.P2D'/max(max(spectravi.P2D)); 
        tnlLtl = tnl_52_over_tl_12 > 1;
        P2D_nlLl=P2D.*tnlLtl';
        
        tlLtnl = tnl_52_over_tl_12 < 1;
        P2D_lLnl=P2D.*tlLtnl';
        
        tnlCBl = (tnl_52_over_tl_12 > 0.9) & (tnl_52_over_tl_12 < 1.1);
        P2D_nlCBl = P2D.*tnlCBl';
        %-------------------------------------------------
        %-------------------------------------------------  
        %Now integrating the power spectra
        SumTP2D=sum(P2D,'all');
        SumP2DnlLl=sum(P2D_nlLl,'all');
        SumP2DlLnl=sum(P2D_lLnl,'all');
        SumP2DnlCBl=sum(P2D_nlCBl,'all');
        % PCB as substracting T- (nl + l)
        %SumP2DnlCBl=SumTP2D-(SumP2DnlLl+SumP2DlLnl);
        
        SP2Dn_nl_i=SumP2DnlLl/SumTP2D;
        SP2Dn_l_i=SumP2DlLnl/SumTP2D;
        SP2Dn_CB_i=SumP2DnlCBl/SumTP2D;
        
        SP2Dn_nl(i)=SP2Dn_nl_i;
        SP2Dn_l(i)=SP2Dn_l_i;
        SP2Dn_CB(i)=SP2Dn_CB_i; 
        
        %-------------------------------------------------
        %{
        cd /Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here
        f3=figure(3);
        yyy=log10(tnl_52_over_tl_12').*log10(spectravi.P2D');
        [C,h]=contourf(spectravi.kpar,spectravi.kper,yyy,'LevelStep',2);
        yyy(isinf(yyy)|isnan(yyy)) = 0;
        set(h,'LineColor','none')
        maxy=max(max(abs(yyy)));
        lim = caxis; 
        caxis([-maxy maxy]);
        colormap('redblue')
        colorbar
        set(gca,'XScale','log','YScale','log','FontSize',18)
        xlabel('$k_\| d_i$','Interpreter','latex')
        ylabel('$k_\perp d_i$','Interpreter','latex')
        title('$\log_{10} \tau_{nl}/\tau_{l} \cdot \log_{10} P_{\mathbf{v}_{i}} \ \omega_{pi}t=$'+string((i-1)*(6)),'Interpreter','latex')
        hold on
        %x=log10([0.06 k23]);
        x=log10([0.06 30]);
        x3=log10([k23 1]);
        y=(3/2)*(x)+1;
        y3=(3)*(x3)+2;
        loglog(10.^x,10.^y,'--k'); loglog(10.^x3,10.^y3,'--k')
        ylim([0.28 52])
        txt = '\bf{\leftarrow 3/2 }';
        txt2 = '\bf{3 \rightarrow}';
        text(0.25,0.5,txt)
        text(0.3,28,txt2)
        ylim([0.28 52])
        yline(10,'--k');
        %xline(1,'--k'); %xline(10,'--k'); xline(k23,'-.k'); xline(k13,'-k');
        hold off
        
        
        f4=figure(4);
        yyy=log10(spectravi.P2D');
        [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(P2D),'LevelStep',0.5);
        yyy(isinf(yyy)|isnan(yyy)) = 0;
        set(h,'LineColor','none')
        maxy=max(max(abs(yyy)));
        lim = caxis; 
        %caxis([-maxy maxy]);
        colormap('redblue')
        colorbar
        set(gca,'XScale','log','YScale','log','FontSize',18)
        xlabel('$k_\| d_i$','Interpreter','latex')
        ylabel('$k_\perp d_i$','Interpreter','latex')
        title('$\log_{10} P_{\mathbf{v}_{i}} \ \omega_{pi}t=$'+string((i-1)*(6)),'Interpreter','latex')
        ylim([0.28 52])
        txt = '\bf{\leftarrow \tau_{nl} = \tau_{l}}';
        text(0.25,0.5,txt,'FontSize',18)
        hold on
        contour(spectravi.kpar,spectravi.kper,log10(tnl_52_over_tl_12'),[-0.01 0 -0.01],'-k');
        hold off
        %}
        %--------------------------------------------------------------------
        
      
        % These are the plots showing what is above and below CB
        %
        f55=figure(55);
        [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(P2D_nlLl),'LevelStep',0.5);
        set(h,'LineColor','none')
        lim = caxis; 
        %caxis([-maxy maxy]);
        colormap(jet)
        colorbar
        set(gca,'XScale','log','YScale','log','FontSize',18)
        xlabel('$k_\| d_i$','Interpreter','latex')
        ylabel('$k_\perp d_i$','Interpreter','latex')
        title('$\log_{10} P_{\mathbf{v}_{i}} \ \omega_{pi}t=$'+string((i-1)*(6)),'Interpreter','latex')
        ylim([0.28 52])
        txt = '\bf{\leftarrow \tau_{nl} = \tau_{l}}';
        text(0.25,0.5,txt,'FontSize',18)
        hold on
        contour(spectravi.kpar,spectravi.kper,log10(tnl_52_over_tl_12'),[-0.1 0 0.1],'-k');
        hold off
        
        
        f6=figure(6);
        [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(P2D_lLnl),'LevelStep',0.5);
        set(h,'LineColor','none')
        lim = caxis; 
        %caxis([-maxy maxy]);
        colormap(jet)
        colorbar
        set(gca,'XScale','log','YScale','log','FontSize',18)
        xlabel('$k_\| d_i$','Interpreter','latex')
        ylabel('$k_\perp d_i$','Interpreter','latex')
        title('$\log_{10} P_{\mathbf{v}_{i}} \ \omega_{pi}t=$'+string((i-1)*(6)),'Interpreter','latex')
        ylim([0.28 52])
        txt = '\bf{\leftarrow \tau_{nl} = \tau_{l}}';
        text(0.25,0.5,txt,'FontSize',18)
        hold on
        contour(spectravi.kpar,spectravi.kper,log10(tnl_52_over_tl_12'),[-0.1 0 0.1],'-k');
        hold off
        
        f7=figure(7);
        [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(P2D_nlCBl),'LevelStep',0.5);
        set(h,'LineColor','none')
        lim = caxis; 
        %caxis([-maxy maxy]);
        colormap(jet)
        colorbar
        set(gca,'XScale','log','YScale','log','FontSize',18)
        xlabel('$k_\| d_i$','Interpreter','latex')
        ylabel('$k_\perp d_i$','Interpreter','latex')
        title('$\log_{10} P_{\mathbf{v}_{i}} \ \omega_{pi}t=$'+string((i-1)*(6)),'Interpreter','latex')
        ylim([0.28 52])
        txt = '\bf{\leftarrow \tau_{nl} = \tau_{l}}';
        text(0.25,0.5,txt,'FontSize',18)
        hold on
        contour(spectravi.kpar,spectravi.kper,log10(tnl_52_over_tl_12'),[-0.1 0 0.1],'-k');
        hold off
        
        f8=figure(8);
        [C,h]=contourf(spectravi.kpar,spectravi.kper,log10(P2D_nlLl),'LevelStep',0.5);
        set(h,'LineColor','none')
        lim = caxis; 
        %caxis([-maxy maxy]);
        colormap(jet)
        colorbar
        set(gca,'XScale','log','YScale','log','FontSize',18)
        xlabel('$k_\| d_i$','Interpreter','latex')
        ylabel('$k_\perp d_i$','Interpreter','latex')
        title('$\log_{10} P_{\mathbf{v}_{i}} \ \omega_{pi}t=$'+string((i-1)*(6)),'Interpreter','latex')
        ylim([0.28 52])
        txt = '\bf{\leftarrow \tau_{nl} = \tau_{l}}';
        text(0.25,0.5,txt,'FontSize',18)
        hold on
        [C2,h2]=contourf(spectravi.kpar,spectravi.kper,log10(P2D_lLnl),'LevelStep',0.5);
        set(h2,'LineColor','none')
        [C3,h3]=contourf(spectravi.kpar,spectravi.kper,log10(P2D_nlCBl),'LevelStep',0.5);
        contour(spectravi.kpar,spectravi.kper,log10(tnl_52_over_tl_12'),[-0.1 0 0.1],'-k');
        hold off
        %}
        
        %cd '/Volumes/PSC_DiRAC_DATA/nonlineartime_images/2021_nlinear_time/';
        %-------------------------------------------------------------------------
        %saveas(f3,'2021_timelin_2D_' + string(i)+'.png');
    end

     time=[0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,96,102,108,114,120,126,132,156,168,180,...
            192,204,216,228,240,264,288,312,336,360,384,408,432,456,480];
        %-------------------------------------------------
        f11=figure(11);
        plot(time,SP2Dn_nl,'-r','LineWidth',1.5)
        set(gca,'XScale','lin','YScale','log','FontSize',20)
        xlabel('$t \omega_{pi}$','Interpreter','latex')
        ylabel('$\log_{10}P_{i}/P_{T}$','Interpreter','latex')
        %title('$\log_{10} \tau_{nl}/\tau_{l} \cdot \log_{10} P_{\mathbf{v}_{i}} \ \omega_{pi}t=$'+string((i-1)*(6)),'Interpreter','latex')
        hold on
        plot(time,SP2Dn_l,'-b','LineWidth',1.5)
        plot(time,SP2Dn_CB,'-k','LineWidth',1.5)
        plot(time,SP2Dn_CB+SP2Dn_CB+SP2Dn_l,'*m','LineWidth',1.5)
        legend('$P_{\tau_{nl}>\tau_{l}}$','$P_{\tau_{nl}<\tau_{l}}$','$P_{\tau_{nl}=\tau_{l}}$','Interpreter','latex','FontSize',20)
        hold off
       
    
end
