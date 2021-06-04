%#-----------------------------------------------------------------------
%This is the script to read the particle information to compare the output
%moments from the fields data with the computed moments using the partcile
%data
%#----------------------------------------------------------------------


% For the fields data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/field_outputs'
pathf = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/field_outputs';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
i_f = 6;
Sf = dir(fullfile(pathf,'*.h5'));
Nf=numel(Sf); % number of files to use
Hf=zeros(Nf);
disp(strcat('Computing step ...',Sf(i_f).name)) 
fileIDf =  Sf(i_f).name; %change the 1 per i
%h5disp(fileIDf); %this is the command to show the information in the file 
infof = h5info(fileIDf);  

coord_1st_x=infof.Groups(2).Name;
coord_1st_y=infof.Groups(3).Name;
coord_1st_z=infof.Groups(4).Name;

coord_x=h5read(fileIDf,strcat(coord_1st_x,'/crd[0]/p0/1d'));
coord_y=h5read(fileIDf,strcat(coord_1st_y,'/crd[1]/p0/1d'));
coord_z=h5read(fileIDf,strcat(coord_1st_z,'/crd[2]/p0/1d'));

%T_1st = infof.Groups(1).Name;
%E_1st = infof.Groups(17).Name;
%B_1st = infof.Groups(18).Name;
%J_1st = infof.Groups(19).Name;
V_1st = infof.Groups(25).Name;
n_1st = infof.Groups(23).Name;

%Bx=h5read(fileIDf,strcat(B_1st,'/hx/p0/3d'));
%By=h5read(fileIDf,strcat(B_1st,'/hy/p0/3d'));
%Bz=h5read(fileIDf,strcat(B_1st,'/hz/p0/3d'));
%Ex=h5read(fileIDf,strcat(E_1st,'/ex/p0/3d'));
%Ey=h5read(fileIDf,strcat(E_1st,'/ey/p0/3d'));
%Ez=h5read(fileIDf,strcat(E_1st,'/ez/p0/3d'));
vix=h5read(fileIDf,strcat(V_1st,'/vx_i/p0/3d'));
viy=h5read(fileIDf,strcat(V_1st,'/vy_i/p0/3d'));
viz=h5read(fileIDf,strcat(V_1st,'/vz_i/p0/3d'));
%vex=h5read(fileIDf,strcat(V_1st,'/vx_e/p0/3d'));
%vey=h5read(fileIDf,strcat(V_1st,'/vy_e/p0/3d'));
%vez=h5read(fileIDf,strcat(V_1st,'/vz_e/p0/3d'));
ni = h5read(fileIDf,strcat(n_1st,'/n_i/p0/3d'));
%ne = h5read(fileIDf,strcat(n_1st,'/n_e/p0/3d'));
%--------------------------------------------------------------------------

vi_index_y = find( ((1<=coord_y) & (coord_y<=1.002)) & ((1<=coord_z) & (coord_z<=1.002)) );
vix2=vix(vi_index_y);


%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% For the particle data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/prt_outputs'
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/prt_outputs';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mime=300;
S = dir(fullfile(path,'prt.*'));
N=numel(S); % number of files to use
H=zeros(N);
i=6; %prt 800
    

%for i=1:28
    disp(strcat('Computing step ...',S(i).name)) 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID); %this is the command to show the information in the file 
    info = h5info(fileID);    
    %info = h5info(fileID,'TextEncoding','UTF-8');
    
    prt = info.Groups(1).Name;
    %read the information from the fields of the struct
    prtfields=h5read(fileID,strcat(prt,'/p0/1d'));
    x = prtfields.x;
    y = prtfields.y;
    z = prtfields.z;
    px = prtfields.px;
    py = prtfields.py;
    pz = prtfields.pz;
    q = prtfields.q;
    %m = prtfields.m;
    %w = prtfields.w;
    
    %separate the popullations
    %
    q_i_index = find(q==1);
    %Coordinates of ions
    xi=x(q_i_index); yi=y(q_i_index); zi=z(q_i_index);
    pxi=px(q_i_index); pyi=py(q_i_index); pzi=pz(q_i_index);
    qi=q(q_i_index);
    
    %create the grid to see the values in 2D
    nbins=[600 600];
    [Nyz,yedges,zedges] = histcounts2(yi,zi,nbins);
    Nyzn=Nyz./80; %This is equivalent to the density 
    
    % The index works edges(n) <= index(n) <= edges(n+1) !
    j=1; jj=j+1;
    k=1; kk=k+1;
    index2 = find(((yedges(1,j) <= yi) & (yi <= yedges(1,jj))) & ((zedges(1,k) <= zi) & (zi <= zedges(1,kk))));
    Nyzn(j,k)
            
    % Get the indices for a specific pixel
    %index2 = find(((1 <= yi) & (yi <= 1.002)) & ((1 <= zi) & (zi <= 1.002)));
    
    % Get the velocities 
    pxi2=px(index2); pyi2=py(index2); pzi2=pz(index2);
    
    %----------------------------------------------------------------------
    vvx = mean(pxi2);
    vvvx=mean(pxi2 - vvx)
    %----------------------------------------------------------------------
    
    % to get the delta value and the limits to compute the PMF
    pxi2_s=sort(pxi2);
    pxi2_s1=pxi2_s(1:end-1);    pxi2_s2=pxi2_s(2:end);
    deltapxi2s = pxi2_s2-pxi2_s1; deltapxi2_m=mean(deltapxi2s);
    pixa=min(pxi2); pixb=max(pxi2); 
    edpix=pixa:deltapxi2_m:pixb;
    
    pyi2_s=sort(pyi2);
    pyi2_s1=pyi2_s(1:end-1);    pyi2_s2=pyi2_s(2:end);
    deltapyi2s = pyi2_s2-pyi2_s1; deltapyi2_m=mean(deltapyi2s);
    piya=min(pyi2); piyb=max(pyi2);
    edpiy=piya:deltapyi2_m:piyb;
    
    pzi2_s=sort(pzi2);
    pzi2_s1=pzi2_s(1:end-1);    pzi2_s2=pzi2_s(2:end);
    deltapzi2s = pzi2_s2-pzi2_s1; deltapzi2_m=mean(deltapzi2s);
    piza=min(pzi2); pizb=max(pzi2);
    edpiz=piza:deltapzi2_m:pizb;    
    
    [Npxi,edges] = histcounts(pxi2_s,edpix);
    
    %Computing the moments    
    
    %pdf for the specific pixel 
    sm=sum(Npxi,'all');
    Npxid=Npxi*(Nyzn(j,k)/sm);
    sum(Npxid) % this is equivalent to the zeroth moment of f
    
    pxi2_s_mv=(pxi2_s(1:end-1) + pxi2_s(2:end))/2;
    Upix_1 = sum(pxi2_s_mv.*Npxid')    
    dUpix_1 = sum( (pxi2_s_mv-Upix_1).*Npxid' )
        
    
    vix(1,j,k)
    Nyzn(j,k)
    ni(1,j,k)
    
    
    figure(3212)
    plot((pxi2(1:end-1).*Npxi')/(sm*80))
   
    Uix = sum(pxi2)/
    Nyz(2,2)    
   

    f567=figure(567)
    plot(pxi2_s(1:end-1),Npxid,'*')
    
    mean(pxi2_s)
    
    hold on
    scatter(pxi2_s(2:end),Npxid)
    hold off
    
    f897=figure(897);
    plot(edges(1:end-1),Npxi)

    f789=figure(789)
    semilogy(deltapxi2s)
    hold on
    semilogy(deltapyi2s)
    semilogy(deltapzi2s)
    hold off
    

    
    %get the values 
    %clearvars px py pz
    
    %Calculate perpendicular
    pxyi=sqrt((pxi.*pxi) + (pyi.*pyi));
    %clearvars pxi pyi
    
    %----------------------------------------------------------------------
    % Testing plots
    
    f675 = figure(675);
    histogram(pxi2,100)
    hold on
    histogram(pyi2,100)
    histogram(pzi2,100)
    hold off
    
    f435=figure(435);
    subplot(1,3,1)
%    hc=pcolor(yedges(1:end-1),zedges(1:end-1),Nyz./80);
    hc=pcolor(Nyzn(1:10,1:10));
    set(hc,'edgecolor','none')
    colorbar
    caxis([0 2]);
    subplot(1,3,2)
    nis = squeeze(ni);
    hc2=pcolor(nis(1:10,1:10));
    set(hc2,'edgecolor','none')
    colorbar
    caxis([0 2]);
    subplot(1,3,3)
    hc2=pcolor((nis(1:10,1:10) - Nyzn(1:10,1:10))./nis(1:10,1:10));
    set(hc2,'edgecolor','none')
    colorbar
    
    
    %----------------------------------------------------------------------
    
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