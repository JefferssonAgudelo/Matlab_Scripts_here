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

cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/short';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/short';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296/nostrahl';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296/nostrahl';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mime=300;
S = dir(fullfile(path,'prt.*'));
N=numel(S); % number of files to use
H=zeros(N);
i=1;

% This is the loop to calculate the VDF of the particles
%
for i=1:N
    disp(strcat('Computing step ...',S(i).name)) 
    fileID =  S(i).name; %change the 1 per i
    h5disp(fileID); %this is the command to show the information in the file 
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
    
    Nenan=Ne./(De);
    Nenan(Nenan==0)=NaN;
    Nenanes19{i}=Nenan;
end
   

    %}
    
    %save('2DVDF.mat')
%}

% when loading the file:
load('2DVDF_16296_nostrahl_25000')

    x = linspace(-pi,pi);
    y = linspace(0,2*pi);
    [X,Y] = meshgrid(x,y);
    Z = sqrt(X.*X + Y.*Y);
    v = [1,1];
    v2 = [0.5,0.5];
    
    %%Parameters for background plasma
    %
    kb = 1.; mu0 =1.; mi_over_me_ =300.; vA_over_c_ = 0.008; B0_ = vA_over_c_;
    mi = 1.; me = 1. / mi_over_me_;
    omega_pi=1;
    Omega_e=vA_over_c_*sqrt(mi_over_me_)*omega_pi;
    Omega_i=Omega_e / mi_over_me_;
    TOmega_e=1/Omega_e; 
    TOmega_i=1/Omega_i;
    Ttotal=300*TOmega_e;
    dt=0.0033;
    N_tsteps=Ttotal/dt;
    N_Oe=TOmega_e/dt; % this is the number of time steps required for 1 TOmega_e
    t_final=00000;
    time = (t_final/N_tsteps);
    time_Oe = t_final/N_Oe;
    %--------------------------------------------------------------------
    %thermal speed
    %Vai = vae *(ne*me)/(ni*mi)
    % Change this accordingly to run parameters
    %--------------------------------------------------------------------
    % 16296
    %{
    ni_ = 1.;
    ne_core_ = 0.978075171;
    ne_strahl_ = 0.021924829;
    beta_i_par_ = 1.200545495;
    beta_e_core_par_ = 3.506865593;
    beta_e_strahl_par_ = 0.053475199;
    %These are the temperature anisotropies
    Ti_per_over_Ti_par_ = 1.;
    Te_c_per_over_Te_c_par_ = 1.118098593;
    Te_s_per_over_Te_s_par_ = 2.400797005;
    vthe_c_par_=sqrt(beta_e_core_par_)*sqrt(ni_/ne_core_)*sqrt(mi_over_me_)*vA_over_c_;
    %} 
    %--------------------------------------------------------------------
    %
    ni_ = 1.;
    ne_core_ = 1.;
    ne_strahl_ = 0.;
    beta_i_par_ = 1.;%1.200545495;
    beta_e_core_par_ = 1.;%3.506865593;
    beta_e_strahl_par_ = 1.;%0.053475199;
    %These are the temperature anisotropies
    Ti_per_over_Ti_par_ = 1.;
    Te_c_per_over_Te_c_par_ = 1.;%1.118098593;
    Te_s_per_over_Te_s_par_ = 1.;%2.400797005;
    vthe_c_par_=sqrt(beta_e_core_par_)*sqrt(ni_/ne_core_)*sqrt(mi_over_me_)*vA_over_c_;
    %}
    %--------------------------------------------------------------------
    % 190121
    %{
    ni_ = 1.;
    ne_core_ = 0.967666899;
    ne_strahl_ = 0.032333101;
    beta_i_par_ = 0.417714;
    beta_e_core_par_ = 0.853699832;
    beta_e_strahl_par_ = 0.05450854;
    %These are the temperature anisotropies
    Ti_per_over_Ti_par_ = 1.;
    Te_c_per_over_Te_c_par_ = 0.747716263;
    Te_s_per_over_Te_s_par_ = 0.370285197;
    vthe_c_par_=sqrt(beta_e_core_par_)*sqrt(ni_/ne_core_)*sqrt(mi_over_me_)*vA_over_c_;
    %}
    %--------------------------------------------------------------------

    %----------------------------------------------------------------------
    %d_step_out=25000;
    d_step_out=25000;
    Nenan00=Nenanes19{1};
    max00=max(max(Nenan00));
    B=Nenan00; B(isnan(B))=0; Bmean=mean(B,'all');
    Brange = 0.05*Bmean/max00;
    %----------------------------------------------------------------------

    
    %i=3;
    %i=7
    for i=2:7
    %----------------------------------------------------------------------
    Nenan0 = Nenanes19{i-1}; 
    Nenan1 = Nenanes19{i};
    time_Oe_im1=d_step_out*dt*vA_over_c_*sqrt(mime)*(i-2);
    time_Oe_i=d_step_out*dt*vA_over_c_*sqrt(mime)*(i-1);
    %----------------------------------------------------------------------
    %Nenan00=Nenanes19{1};
    %
    %cd /Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here
    
    length = 20;
    red = [1, 0, 0]; white = [1, 1, 1]; blue = [0, 0, 1];
    pink = [255, 192, 203]/255;
    colors_p1 = [linspace(white(1),red(1),length)', linspace(white(2),red(2),length)', linspace(white(3),red(3),length)'];
    colors_p2 = [linspace(white(1),blue(1),length)', linspace(white(2),blue(2),length)', linspace(white(3),blue(3),length)'];
    
    
    %{
    f6=figure(6);
    h(1)=pcolor(Xedges3e,Redges3e, (Nenan1 - Nenan00)./Nenan00);
    hold on
    %colormap(colors_p1)
    h(2)=pcolor(Xedges3e,Redges3e, (- Nenan1 + Nenan00)./Nenan00);
    colormap([colors_p1;colors_p2])
    %colormap([hot(25);winter(25)])
    %colormap('red')
    set(h(1), 'EdgeColor', 'none');
    lim = caxis;  
    caxis([0.01 500]); %caxis([0.01 1e4]); for the other 2 and 5
    %colormap('redblue')
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 2.4]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    title('$$|\Delta VDF_{e}| / VDF_{e,0}, \ t_{i}=$$'+string(round(time_Oe_i,2)) + '$$ \ \Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    [M,c]=contour(X,Y,Z,v);
    c.LineWidth = 3;
    grid off
    pbaspect([1 1 1])
    hold off
    %}
    
    f7=figure(7);
    ax(1)=subplot(2,1,1);
    h=pcolor(Xedges3e,Redges3e, (Nenan1 - Nenan00)./max00);
    hold on
    set(h, 'EdgeColor', 'none');
    lim = caxis;  
    caxis([Brange 1]); 
    colormap(ax(1),colors_p1)
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 1.2]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    title('$$\Delta VDF_{e} / VDF_{e,max0}$$','Interpreter','latex','FontSize',20)
    [M,c1]=contour(X,Y,Z,v); c1.LineWidth = 3;
    xline(vthe_c_par_,'--k','LineWidth',1);
    xline(-vthe_c_par_,'--k','LineWidth',1);
    %[M,c2]=contour(X,Y,Z,v2); c2.LineWidth = 3;
    grid off
    pbaspect([1 0.5 0.5])
    hold off
    
    ax(2)=subplot(2,1,2);
    h=pcolor(Xedges3e,Redges3e, -(Nenan1 - Nenan00)./max00);
    hold on
    colormap(ax(2),colors_p2)
    set(h, 'EdgeColor', 'none');
    lim = caxis;  
    caxis([Brange 1]); 
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-1.2 1.2]); ylim([0 1.2]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    title('$$-\Delta VDF_{e} / VDF_{e,max0}, \ t_{i}=$$'+string(round(time_Oe_i,2)) + '$$ \ \Omega_{e}^{-1}$$','Interpreter','latex','FontSize',20)
    [M,c]=contour(X,Y,Z,v);
    c.LineWidth = 3;
    xline(vthe_c_par_,'--k','LineWidth',1);
    xline(-vthe_c_par_,'--k','LineWidth',1);
    grid off
    pbaspect([1 0.5 0.5])
    hold off
    
    
    %{
    f8=figure(8);
    ax3 = getCustomAxesPos(3,2,1);
    pcolor(ax3(1,2),Xedges3e,Redges3e, -(Nenan1 - Nenan00)./Nenan00);
    %}
    saveas(f7,strcat('diff_init_16296_nostrahl_25000_' + string(i) +'.png'));
    end
    

    
