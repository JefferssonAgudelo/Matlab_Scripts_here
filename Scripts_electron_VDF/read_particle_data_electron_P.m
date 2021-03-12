%#-----------------------------------------------------------------------
%This is the script to read the particle infoprmation
% Currently the program does not create a xdmf but only a .h5
%#----------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/test_electron_project/run_030820';
% This is the info for 50 particles
%cd '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/test_electron_project/run_150820';
%path = '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/test_electron_project/run_150820';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_trillian/electron_project_19_uns_50';
cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_trillian/electron_project19_uns_trillian';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_trillian/electron_project19_uns_trillian';

cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/ISO_3D/ISO_3D_50_80';
path = '/Volumes/PSC_DiRAC_DATA/TEST2020_1/ISO_3D/ISO_3D_50_80';

cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Extended';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Extended';

cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/prt_outputs'
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/prt_outputs'

cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/short';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/short';

cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_test_190121';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_test_190121';

cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296_03';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296_03';

cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296/nostrahl';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296/nostrahl';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mime=200;
%mime=50;
mime=300;
S = dir(fullfile(path,'prt.*'));
%S = dir(fullfile(path,'pfd.*'));
N=numel(S); % number of files to use
H=zeros(N);
i=1; %prt 800
for i=1:28
    disp(strcat('Computing step ...',S(i).name)) 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID); %this is the command to show the information in the file 
    info = h5info(fileID);    
    %info = h5info(fileID,'TextEncoding','UTF-8');
    
    prt = info.Groups(1).Name;
    %read the information from the fields of the struct
    prtfields=h5read(fileID,strcat(prt,'/p0/1d'));
    %x = prtfields.x;
    %y = prtfields.y;
    %z = prtfields.z;
    px = prtfields.px;
    py = prtfields.py;
    pz = prtfields.pz;
    q = prtfields.q;
    %m = prtfields.m;
    %w = prtfields.w;
    
    %separate the popullations
    %{
    q_i_index = find(q==1);
    %Coordinates of ions
    %xi=x(q_i_index); yi=y(q_i_index); zi=z(q_i_index);
    pxi=px(q_i_index); pyi=py(q_i_index); pzi=pz(q_i_index);
    qi=q(q_i_index);
    clearvars px py pz
    %Calculate perpendicular
    pxyi=sqrt((pxi.*pxi) + (pyi.*pyi));
    clearvars pxi
    [Ni,Xedgesi,Redgesi] = histcounts2(pzi,pxyi,80);
    Ni=Ni(1:79,1:79)';    
    Xedges2i=Xedgesi(2:80)';
    Redges2i=Redgesi(2:80)';
    Xedges3i=Xedges2i-(Xedges2i(2)-Xedges2i(1))/2;    
    Redges3i=Redges2i-(Redges2i(2)-Redges2i(1))/2;
    [Ci,Di]=meshgrid(Xedges3i,Redges3i);
   %}
    
    q_e_index = find(q==-1);
    %Coordinates of electrons
    %xe=x(q_e_index); ye=y(q_e_index); ze=z(q_e_index);
    pxe=px(q_e_index); pye=py(q_e_index); pze=pz(q_e_index);
    %qe=q(q_e_index);
    
    clearvars px py pz q
    %calculate perpendicular
    pxye=sqrt((pxe.*pxe) + (pye.*pye));
    clearvars pxe pye
   %{
    [Ne,Xedgese,Redgese] = histcounts2(pze,pxye,80);
    Ne=Ne(1:79,1:79)';    
    Xedges2e=Xedgese(2:80)';
    Redges2e=Redgese(2:80)';
    Xedges3e=Xedges2e-(Xedges2e(2)-Xedges2e(1))/2;    
    Redges3e=Redges2e-(Redges2e(2)-Redges2e(1))/2;
    [Ce,De]=meshgrid(Xedges3e,Redges3e);
    %}
    
    Xedges2D = linspace(-1.2,1.2);
    Yedges2D = linspace(0,2.4);
    [Ne,Xedgese,Redgese] = histcounts2(pze,pxye,Xedges2D,Yedges2D);
    Ne=Ne(1:99,1:99)';    
    Xedges2e=Xedgese(2:100)';
    Redges2e=Redgese(2:100)';
    Xedges3e=Xedges2e-(Xedges2e(2)-Xedges2e(1))/2;    
    Redges3e=Redges2e-(Redges2e(2)-Redges2e(1))/2;
    [Ce,De]=meshgrid(Xedges3e,Redges3e);
    
    x = linspace(-pi,pi);
    y = linspace(0,2*pi);
    [X,Y] = meshgrid(x,y);
    Z = sqrt(X.*X + Y.*Y);
    v = [1,1];
    
    %%Parameters for background plasma
    %{
    kb = 1.; mu0 =1.; mi_over_me_ =50.; vA_over_c_ = 0.06; B0_ = vA_over_c_;
    mi = 1.; me = 1. / mi_over_me_;
    omega_pi=1;
    Omega_e=vA_over_c_*sqrt(mi_over_me_)*omega_pi;
    Omega_i=Omega_e / mi_over_me_;
    TOmega_e=1/Omega_e; 
    TOmega_i=1/Omega_i;
    Ttotal=300*TOmega_e;
    dt=0.0078;
    N_tsteps=Ttotal/dt;
    N_Oe=TOmega_e/dt; % this is the number of time steps required for 1 TOmega_e
    t_final=00000;
    time = (t_final/N_tsteps);
    time_Oe = t_final/N_Oe;
    %}
    
    
    
    %-----------------------------------------------------------------------
    %{
    f1=figure(1);
    %h = histogram2(pzi,pxyi,'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    %c = h.NumBins;
    %h.NumBins = [300 300];
    pcolor(Xedges3i,Redges3i, Ni./(Di))
    lim = caxis;  
    caxis([0.01 1e4]);
    colormap(bone(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-0.2 0.2]); ylim([0 0.015]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{i\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{i\perp}$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{i}$$'+string(i),'Interpreter','latex','FontSize',20) 
    
    clearvars h
    
    %f12=figure(12);
    %h = histogram2(pzi,pyi,'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    %c = h.NumBins;
    %h.NumBins = [300 300];
    %lim = caxis;  
    %caxis([0.01 1e4]);
    %colormap(bone(25))
    %set(gca,'ColorScale','log','FontSize',20)
    %colorbar
    %xlim([-0.2 0.2]); ylim([-0.2 0.2]);
    %set(gca,'XScale','lin','YScale','lin','FontSize',20)
    %xlabel('$$v_{i\|}$$','Interpreter','latex','FontSize',20)
    %ylabel('$$v_{iy}$$','Interpreter','latex','FontSize',20)
    %title('$$VDF_{i}$$'+string(i),'Interpreter','latex','FontSize',20)
    
    clearvars h
    %}
    %
    Nenan=Ne./(De);
    Nenan(Nenan==0)=NaN;
    Nenanes19{i}=Nenan;
end
    
    %----------------------------------------------------------------------
    %This is to store the initial time
    %Nenan0=Nenan;
    
    %Nenan_190121_1=Nenan;
    %Nenan=Nenan_16296_03_0;
    %Nenan_i21=Nenan-Nenan0;
    
    %Nenan_16296_0, Nenan_16296_1, Nenan_16296_03_0, Nenan_16296_03_1,
    %Nenan_190121_0, Nenan_190121_1
    
    %Nenan0 = Nenan_16296_03_0;
    %Nenan1 = Nenan_16296_03_1;
    
    d_step_out=25000; %5000
    %i=2;
    %d_step_out=150000;
    
    %i=28; %(0,1,9,27) 
    time_Oe_i=d_step_out*dt*vA_over_c_*sqrt(mime)*(i-1); %This is in terms of omega_e
    %---------------------------------------------------------------------- 
    
    f2=figure(2);
    %h = histogram2(pze,pxye,'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    %c = h.NumBins; %h.NumBins = [300 300];
    h=pcolor(Xedges3e,Redges3e, Nenan1);
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    caxis([1 1e7]); %caxis([0.01 1e4]); for the other 2 and 5
    %colormap(cool(25));
    colormap(jet(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 2.4]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{e} \ t = $$'+string(round(time_Oe_i,2)) + '$$\Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    [M,c]=contour(X,Y,Z,v);
    c.LineWidth = 3;
    grid off
    pbaspect([1 1 1])
    hold off
    saveas(f2,strcat('veparper_16296_03_'+string(round(time_Oe_i,2)),'.png'));
    
    
    
    %{
    [Nper,edgesper,nbinsper] = histcounts(pxye,80);
    Nper=Nper';
    N2per=Nper(1:79);
    edges2per=edgesper(2:80)';
    edges3per=edges2per-(edges2per(2)-edges2per(1))/2;
    f42=figure(42);
    plot(edges3per,N2per./(edges3per*2*pi))
    xlabel('$v_{s,\perp}$','Interpreter','latex','FontSize',20)
    ylabel('Counts','Interpreter','latex','FontSize',20)
    title('$$v_{\perp}$$'+string(i),'Interpreter','latex','FontSize',20)
    xline(0,'--k')
    legend('$v_{e,\perp,real}$','Interpreter','latex','FontSize',20)
    %}
    
    %edges = linspace(-2,2);
    edges = linspace(-1,1);
    [Nr,edges] = histcounts(pze,edges);
    %{
    f67=figure(67);
    semilogy(edges(2:100),Nr,'-+r')
    hold on
    semilogy(edges(2:100),Ncs,'k')
    xline(0,'k')
    xlabel('$v_{s,\parallel}$','Interpreter','latex','FontSize',20)
    ylabel('Counts','Interpreter','latex','FontSize',20)
    title('$$Distribution$$','Interpreter','latex','FontSize',20)
    legend('$v_{e,\parallel,real}$','$v_{e,\parallel,fake}$','Interpreter','latex','FontSize',20)
    hold off
    %}
    [Nrpdf,edges] = histcounts(pze,edges,'Normalization', 'probability');
    f670=figure(670);
    %semilogy(edges(2:100),Nrpdf,'-+r')
    semilogy(edges(2:100),Nr/sum(Nr),'-+r')
    hold on
    %semilogy(edges(2:100),Ncspdf,'k')
    semilogy(edges(2:100),Ncs/sum(Ncs),'k')
    xline(0,'k')
    ylim([1e-7, 0.5])
    xlabel('$v_{s,\parallel}$','Interpreter','latex','FontSize',20)
    ylabel('Probability','Interpreter','latex','FontSize',20)
    title('$$Distribution \ t =$$' + string(round(time_Oe_i,2)) + '$$\Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    legend('$v_{e,\parallel,real}$','$v_{e,\parallel,fake}$','Interpreter','latex','FontSize',20)
    hold off
    
    %cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/Images_this_run';
    %saveas(f67,strcat('vpar_RF_19uns_50_counts_'+string(i),'.png'));
    saveas(f670,strcat('pdf_RF_16296_03_'+string(round(time_Oe_i,2)),'.png'));
    
   
    %cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/differences/deltas';
    cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/differences/deltas_init'
    
    cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/short/images/deltas';
    Nenan00=Nenanes19_short{1};
    for i=2:21
    %----------------------------------------------------------------------
    Nenan0 = Nenanes19_short{i-1}; 
    Nenan1 = Nenanes19_short{i};
    d_step_out=500;
    time_Oe_im1=d_step_out*dt*vA_over_c_*sqrt(mime)*(i-2);
    time_Oe_i=d_step_out*dt*vA_over_c_*sqrt(mime)*(i-1);
    %----------------------------------------------------------------------
    %Nenan00=Nenanes19{1};
    %
    f6=figure(6);
    h=pcolor(Xedges3e,Redges3e, abs(Nenan1 - Nenan00)./Nenan00);
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    caxis([0.01 50]); %caxis([0.01 1e4]); for the other 2 and 5
    colormap(jet(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 2.4]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    %title('$$\Delta VDF_{e}($$'+string(round(time_Oe_1,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    %title('$$|\Delta VDF_{e}($$'+string(round(time_Oe_27,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}|$$','Interpreter','latex','FontSize',20)
    title('$$|\Delta VDF_{e}| / VDF_{e,0}, \ t_{i}=$$'+string(round(time_Oe_i,2)) + '$$ \ \Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    %title('$$VDF_{e}(t_{i}) - VDF_{e}(t_{0}), \ t_{i}=$$'+string(round(time_Oe_i,2)) + '$$ \ \Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    %title('$$VDF_{e,03}(t_{1}) - VDF_{e}(t_{1})$$','Interpreter','latex','FontSize',20)
    [M,c]=contour(X,Y,Z,v);
    c.LineWidth = 3;
    grid off
    pbaspect([1 1 1])
    hold off
    
    saveas(f6,strcat('diff_init_19_' + string(i) +'.png'));
    end
    
    f66=figure(66);
    h=pcolor(Xedges3e,Redges3e, (Nenan0 - Nenan1));
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    %caxis([1 1e7]); %caxis([0.01 1e4]); for the other 2 and 5
    colormap(jet(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 2.4]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    %title('$$\Delta VDF_{e}($$'+string(round(time_Oe_1,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    %title('$$|\Delta VDF_{e}($$'+string(round(time_Oe_27,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}|$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{e}(t_{0}) - VDF_{e}(t_{1})$$','Interpreter','latex','FontSize',20)
    %title('$$VDF_{e}(t_{1}) - VDF_{e,03}(t_{1})$$','Interpreter','latex','FontSize',20)
    [M,c]=contour(X,Y,Z,v);
    c.LineWidth = 3;
    grid off
    pbaspect([1 1 1])
    hold off
    
    f77=figure(77);
    h=pcolor(Xedges3e,Redges3e, abs(Nenan1 - Nenan0));
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    %caxis([1 1e7]); %caxis([0.01 1e4]); for the other 2 and 5
    colormap(jet(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 2.4]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    %title('$$\Delta VDF_{e}($$'+string(round(time_Oe_1,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    %title('$$|\Delta VDF_{e}($$'+string(round(time_Oe_27,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}|$$','Interpreter','latex','FontSize',20)
    title('$$|VDF_{e}(t_{1}) - VDF_{e}(t_{0})|$$','Interpreter','latex','FontSize',20)
    %title('$$|VDF_{e,03}(t_{1}) - VDF_{e}(t_{1})|$$','Interpreter','latex','FontSize',20)
    [M,c]=contour(X,Y,Z,v);
    c.LineWidth = 3;
    grid off
    pbaspect([1 1 1])
    hold off
    
    f7=figure(7);
    h=pcolor(Xedges3e,Redges3e, (Nenan1 ./ Nenan0));
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    %caxis([1 1e7]); %caxis([0.01 1e4]); for the other 2 and 5
    colormap(jet(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 2.4]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    %title('$$\Delta VDF_{e}($$'+string(round(time_Oe_1,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    %title('$$|\Delta VDF_{e}($$'+string(round(time_Oe_27,2)) +'/'+string(round(time_Oe_0,2)) + '$$ )\Omega_{e}^{-1}|$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{e}(t_{1})/VDF_{e}(t_{0})$$','Interpreter','latex','FontSize',20)
    %title('$$VDF_{e,03}(t_{1})/VDF_{e}(t_{1})$$','Interpreter','latex','FontSize',20)
    [M,c]=contour(X,Y,Z,v);
    c.LineWidth = 3;
    grid off
    pbaspect([1 1 1])
    hold off
    
    saveas(f6,strcat('veparperp_m1_0.png'));
    saveas(f66,strcat('veparperp_m0_1.png'));
    saveas(f77,strcat('veparperp_absm0_1.png'));
    saveas(f7,strcat('veparperp_div_1_0.png'));
    %}
    
    end
    
    %saveas(f42,strcat('vpar_RF_19_uns_'+string(i),'.png'));
    %{
    me=m(q_e_index); we=w(q_e_index);
    mi=m(q_i_index); wi=w(q_i_index);    
    f23=figure(23);
    histogram(me,90)
    hold on
    histogram(mi,90)
    title('mass')
    hold off
    f24=figure(24);
    histogram(we,90)
    hold on
    histogram(wi,90)
    title('weigth')
    hold off
  
    f25=figure(25);
    histogram(qe,90)
    hold on
    histogram(qi,90)
    title('charge')
    hold off
    %}

    %{
    f779=figure(779);
    histogram(pzecs,90,'Normalization','probability');
    hold on
    histogram(pze,90,'Normalization','probability');
    xlabel('$v_{s,\parallel}$','Interpreter','latex','FontSize',20)
    ylabel('Counts','Interpreter','latex','FontSize',20)
    title('$$Comparing$$'+string(i),'Interpreter','latex','FontSize',20)
    xline(0,'--k')
    legend('$v_{e_f0,\parallel}$','$v_{e,\parallel,real}$','Interpreter','latex','FontSize',20)
    hold off
    %}
    %clearvars h
    
    %{
    %f22=figure(22);
    %h = histogram2(pze,pye,'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    %c = h.NumBins;
    %h.NumBins = [300 300];
    %lim = caxis;  
    %caxis([0.01 1e4]);
    %colormap(bone(25))
    %set(gca,'ColorScale','log','FontSize',20)
    %colorbar
    %xlim([-2 2]); ylim([-2 2]);
    %set(gca,'XScale','lin','YScale','lin','FontSize',20)
    %xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    %ylabel('$$v_{ey}$$','Interpreter','latex','FontSize',20)
    %title('$$VDF_{e}$$'+string(i),'Interpreter','latex','FontSize',20)
    %}
    clearvars h
    
 %   
end
    

 map2 = [0.2 0.1 0.5
    0.1 0.5 0.8
    0.2 0.7 0.6
    0.8 0.7 0.3
    0.9 1 0];


    
    
        %TF = [x y z px py pz q m w];
    %bin the particles to see how many fill a condition
    % for example lets begin with the velocity 
    %pm=sqrt(px.*px + py.*py + pz.*pz);
    %pmi=sqrt(pxi.*pxi + pyi.*pyi + pzi.*pzi);
    %pme=sqrt(pxe.*pxe + pye.*pye + pze.*pze);    
    %{
    signos_i =sign(pyi./pxi);
    signos2_i = signos_i(~isnan(signos_i));
    if length(signos_i) == length(signos2_i)
        pmi_per = pmixy.* signos_i;
    end
    
    signos_e =sign(pye./pxe);
    signos2_e = signos_e(~isnan(signos_e));
    if length(signos_e) == length(signos2_e)
        pme_per = pmexy.* signos_e;
    end
        
    %}
    % Calculate the histcounts
    %[N,pziedges,pmixyedges] = histcounts2(pzi,pmixy,50);
    %[Vzi,Vxyi] = meshgrid(pziedges(1:18),pmixyedges(1:14));
    %f1=figure(1);
    %pcolor(N);
    
    %f2=figure(2);
    %%set(gca,'XScale','log','YScale','log','FontSize',18)
    %histogram(pxyi)
    %ylabel('$v_{i\perp}$','Interpreter','latex','FontSize',18)
    %xlabel('$v_{i\parallel}$','Interpreter','latex','FontSize',18)
    
    
    
    f123=figure(123);
    Xedges = [-Inf -4:0.2:4 Inf];
    Yedges = [-Inf -4:0.2:4 Inf];
    h = histogram2(pze,pme_per,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off');
    %h.NumBins = [50 50];
    lim = caxis; 
    caxis([0 50]);
    colormap(hot)
    colorbar
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlim([-4 4]); ylim([-4 4]);
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    title('$$Electron \ distribution$$','Interpreter','latex','FontSize',20)
    

    f124=figure(124);
    Xedges = [-Inf -4:0.2:4 Inf];
    Yedges = [-Inf -4:0.2:4 Inf];
    h = histogram2(pze,pye,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off');
    %h.NumBins = [50 50];
    lim = caxis; 
    caxis([0 50]);
    colormap(hot)
    colorbar
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlim([-4 4]); ylim([-4 4]);
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{ey}$$','Interpreter','latex','FontSize',20)
    title('$$Electron \ distribution$$','Interpreter','latex','FontSize',20)

    
    f126=figure(126);
    Xedges = [-Inf -4:0.2:4 Inf];
    Yedges = [-Inf -4:0.2:4 Inf];
    h = histogram2(pze,pxe,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off');
    %h.NumBins = [50 50];
    lim = caxis; 
    caxis([0 50]);
    colormap(hot)
    colorbar
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlim([-4 4]); ylim([-4 4]);
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{ex}$$','Interpreter','latex','FontSize',20)
    title('$$Electron \ distribution$$','Interpreter','latex','FontSize',20)

    
    f125=figure(125);
    Xedges = [-Inf -0.4:0.01:0.4 Inf];
    Yedges = [-Inf 0:0.01:0.4 Inf];
    h = histogram2(pzi,pmixy,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off');
    %h.NumBins = [50 50];
    lim = caxis; 
    caxis([0 50]);
    colormap(hot)
    colorbar
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlim([-0.3 0.3]); ylim([0 0.3]);
    xlabel('$$v_{i\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{i\perp}$$','Interpreter','latex','FontSize',20)
    title('$$Ion \ distribution$$','Interpreter','latex','FontSize',20)
    %}
    

    
    
    
    % Making plots
    %----------------------------------------------------------------------
    
    f1=figure(14);
    set(gca,'FontSize',18)
    scatter(y, py)
    ylabel('$p_{y}$','Interpreter','latex','FontSize',18)
    xlabel('$y$','Interpreter','latex','FontSize',18)
    f2=figure(24);
    set(gca,'FontSize',18)
    scatter(yi, pyi)
    ylabel('$p_{yi}$','Interpreter','latex','FontSize',18)
    xlabel('$y_{i}$','Interpreter','latex','FontSize',18)
    f3=figure(34);
    set(gca,'FontSize',18)
    scatter(ye, pye)
    ylabel('$p_{ye}$','Interpreter','latex','FontSize',18)
    xlabel('$y_{e}$','Interpreter','latex','FontSize',18)
    f7=figure(144);
    set(gca,'FontSize',18)
    scatter(z, pz)
    ylabel('$p_{z}$','Interpreter','latex','FontSize',18)
    xlabel('$z$','Interpreter','latex','FontSize',18)
    f8=figure(244);
    set(gca,'FontSize',18)
    scatter(zi, pzi)
    ylabel('$p_{zi}$','Interpreter','latex','FontSize',18)
    xlabel('$z_{i}$','Interpreter','latex','FontSize',18)
    f9=figure(344);
    set(gca,'FontSize',18)
    scatter(ze, pze)
    ylabel('$p_{ze}$','Interpreter','latex','FontSize',18)
    xlabel('$z_{e}$','Interpreter','latex','FontSize',18)
    f4=figure(44); % for ions
    set(gca,'FontSize',18)
    hpmi=histogram(pmi,100);
    ylabel('Number of particles i','Interpreter','latex','FontSize',18)
    xlabel('$$|v_{i}|$$','Interpreter','latex','FontSize',18)
    f5=figure(54); % for electrons
    set(gca,'FontSize',18)
    hpme=histogram(pme,100);
    ylabel('Number of particles e','Interpreter','latex','FontSize',18)
    xlabel('$$|v_{e}|$$','Interpreter','latex','FontSize',18)
    f6=figure(64); % for electrons
    set(gca,'FontSize',18)
    hpm=histogram(pm,100);
    ylabel('Number of particles','Interpreter','latex','FontSize',18)
    xlabel('$$|v|$$','Interpreter','latex','FontSize',18)
    

    
    
    f91=figure(94);
    histfit(pye,100);
    
    %f7=figure(74);
    %histfit(pxe,100);
    %pd = fitdist(pxe,'Normal') %The intervals next to the parameter estimates are the 95% confidence intervals for the distribution parameters.
    
    %cd '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/whistler_1/Images_particles';
    saveas(f1,strcat(S(i).name,'scatteri.png'));
    saveas(f2,strcat(S(i).name,'scattere.png'));
    saveas(f3,strcat(S(i).name,'scatter.png'));
    saveas(f4,strcat(S(i).name,'histoi.png'));
    saveas(f5,strcat(S(i).name,'histoe.png'));
    saveas(f6,strcat(S(i).name,'histo.png'));
%end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract a subset of the data that also follow the same distribution



    %These are diffenrent ways to save the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ae = [xe ye ze pxe pye pze];  
    fileID = fopen('pxyze.txt','w');
    fprintf(fileID,'%6s %10s %10s %10s %10s %10s\n','xe','ye','ze','pxe','pye','pze'); % These print the header
    nbytese = fprintf(fileID,'%8.8f %8.8f %8.8f %8.8f %8.8f %8.8f\n',Ae'); %without the ' it changes the data
    %This is the last used.
    fclose(fileID);
    
    Ai = [xi yi zi pxi pyi pzi];  
    fileID = fopen('pxyzi.txt','w');
    fprintf(fileID,'%6s %10s %10s %10s %10s %10s\n','xi','yi','zi','pxi','pyi','pzi'); % These print the header
    nbytesi = fprintf(fileID,'%8.8f %8.8f %8.8f %8.8f %8.8f %8.8f\n',Ai'); %without the ' it changes the data
    %This is the last used.
    fclose(fileID);
    
    
    
    

    [Vye,Vze] = meshgrid(pye,pze);
    figure44=figure(44);
    pcolor(Vye,Vze,Vze)
    
    histogram2(Vye,Vze)
    
    
    
    Fe = Ye.*exp(-Ye.^2-Ze.^2);
    surf(Ye,Ze,Fe)
    
    
    
    
    
    pmei=[pme, pxe];
    fileIDe = fopen('pme.txt','w');
    fprintf(fileID,'%6s %10s\n','pme'); % These print the header
    nbytese = fprintf(fileIDe,'%f %f\n',pmei'); %This is the last used.
    fclose(fileIDe);
    
    fileID=fopen('pxyzi.txt','r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    %[filename,~,~,encoding] = fopen(fid)
    A = fread(fid,[2000 6],'float');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %Get the file name and character encoding for the open file. Use ~ in place of output arguments you want to omit.
    %[filename,~,~,encoding] = fopen(fileID)
    fileID=fopen('pxyz.txt','r');
    formatSpec = '%f';
    A = fscanf(fileID,formatSpec);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x=double(x);
    y=double(y);
    z=double(z);
    px=double(px);
    
    % Interpolation to see the distribution as a function of the position
    [qx,qy] = meshgrid(linspace(min(x),max(x),6),linspace(min(y),max(y),6));
    Fpx = TriScatteredInterp(x,y,px);
    qpx = Fpx(qx,qy);
    
    % Interpolation to see the distribution as a function of the position
    [qx,qy,qz] = meshgrid(linspace(min(x),max(x),5),linspace(min(y),max(y),10), linspace(min(x),max(x),100));
    %F = TriScatteredInterp(x,y,px);
    Fpx = ScatteredInterpolant(x,y,z,px);
    qpxz = Fpx(qx,qy,qz);
    
    
    % Making plots
    %----------------------------------------------------------------------
 
    f4=figure(4);
    pcolor(qx,qy,log(qpx))
    colorbar
   
   f5=figure(5);
   scatter3(x,y,z,x)
  
    f5=figure(5);
   scatter3(x(:), y(:), z(:),px.*px);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Z] = sphere(16);
x = [0.5*X(:); 0.75*X(:); X(:)];
y = [0.5*Y(:); 0.75*Y(:); Y(:)];
z = [0.5*Z(:); 0.75*Z(:); Z(:)];
 S = repmat([50,25,10],numel(X),1);
C = repmat([1,2,3],numel(X),1);
s = S(:);
c = C(:);

figure
scatter3(x,y,z,s,c)
view(40,35)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
edges = [0:0.03:3];
 f1=figure(1);
 subplot(3,2,1)
 h1 = histogram(x1,edges);
 subplot(3,2,2)
 h2 = histogram(x2,edges);
 subplot(3,2,3)
 h3 = histogram(x3,edges);
 subplot(3,2,4)
 h4 = histogram(x4,edges);
 subplot(3,2,5)
 h5 = histogram(x5,edges);
 %legend('x1','x2','x3','x4','x5')
 %hold off
 
 
 f2=figure(2);
 loglog(x1)
 hold on
 loglog(x2);
 loglog(x3);
 loglog(x4);
 loglog(x5);
 legend('x1','x2','x3','x4','x5')
 hold off
 
 
 f3=figure(3);
 subplot(3,2,1)
 plot(x1)
 subplot(3,2,2)
 plot(x2)
 subplot(3,2,3)
 plot(x3)
 subplot(3,2,4)
 plot(x4)
 subplot(3,2,5)
 plot(x5)
 
