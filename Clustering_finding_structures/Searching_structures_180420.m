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
for i=5:6
    disp(strcat('Computing step ...',S(i).name))
    
    % Read files 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    E_1st = info.Groups(17).Name;
    B_1st = info.Groups(18).Name;
    T_1st = info.Groups(1).Name;
    J_1st = info.Groups(19).Name;
    V_1st = info.Groups(25).Name;
    %n_1st = info.Groups(23).Name;
    %diary off
    %s2='/hx/p0/3d';
    %sx = strcat(B_1st,s2);
    
    %ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    %ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));    
    Ex=h5read(fileID,strcat(E_1st,'/ex/p0/3d'));
    Ey=h5read(fileID,strcat(E_1st,'/ey/p0/3d'));
    Ez=h5read(fileID,strcat(E_1st,'/ez/p0/3d'));
    B = sqrt(Bx.^2 + By.^2 + Bz.^2) ;
    E = sqrt(Ex.^2 + Ey.^2 + Ez.^2) ; 
    E_par = (Ex.*Bx + Ey.*By + Ez.*Bz)./(E.*B);
    E_par_22 = E_par.*E_par;
    %E_par2= ((Ex.*Bx).^2 + (Ey.*By).^2 + (Ez.*Bz).^2)./(E.*B).^2; %This is the one that must be different from zero
    clearvars Bx By Bz Ex Ey Ez

    
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    %Vi2 = vix.^2 + viy.^2 + viz.^2 ;
    Vi_m = sqrt(vix.^2 + viy.^2 + viz.^2) ;
    clearvars vix viy viz
    
    vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'));
    %Ve2 = vex.^2 + vey.^2 + vez.^2 ;  
    Vem = sqrt(vex.^2 + vey.^2 + vez.^2) ;
    clearvars vex vey vez
        
    Jx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    Jy=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    Jz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));
    J2 = Jx.^2 + Jy.^2 + Jz.^2 ; 
    J_m = sqrt(Jx.^2 + Jy.^2 + Jz.^2);
    clearvars Jx Jy Jz
    
    Txx_i=h5read(fileID,strcat(T_1st,'/Txx_i/p0/3d'));
    Tyy_i=h5read(fileID,strcat(T_1st,'/Tyy_i/p0/3d'));
    Tzz_i=h5read(fileID,strcat(T_1st,'/Tzz_i/p0/3d'));
    Ti = (Txx_i + Tyy_i + Tzz_i).*(1/3) ;
    clearvars Txx_i Tyy_i Tzz_i
    
    Txx_e=h5read(fileID,strcat(T_1st,'/Txx_e/p0/3d'));
    Tyy_e=h5read(fileID,strcat(T_1st,'/Tyy_e/p0/3d'));
    Tzz_e=h5read(fileID,strcat(T_1st,'/Tzz_e/p0/3d'));    
    Te = (Txx_e + Tyy_e + Tzz_e).*(1/3) ;
    clearvars Txx_e Tyy_e Tzz_e
   
    J2_mean=mean(J2,'all');
    J4_mean=mean(J2.*J2, 'all');
    J2_std=sqrt(J4_mean - J2_mean^2);
    J2th_3=J2_mean+3*J2_std;
    J2th_4=J2_mean+4*J2_std;    

    Vi2_mean=mean(Vi2,'all');
    Vi4_mean=mean(Vi2.*Vi2, 'all');
    Vi2_std=sqrt(Vi4_mean - Vi2_mean^2);
    Vi2th_3=Vi2_mean+3*Vi2_std;
    Vi2th_4=Vi2_mean+4*Vi2_std;
    
    Ve2_mean=mean(Ve2,'all');
    Ve4_mean=mean(Ve2.*Ve2, 'all');
    Ve2_std=sqrt(Ve4_mean - Ve2_mean^2);
    Ve2th_3=Ve2_mean+3*Ve2_std;
    
    Ti_mean=mean(Ti,'all');
    Ti2_mean=mean(Ti.*Ti, 'all');
    Ti_std=sqrt(Ti2_mean - Ti_mean^2);
    Tith_3=Ti_mean+3*Ti_std;
    
    Te_mean=mean(Te,'all');
    Te2_mean=mean(Te.*Te, 'all');
    Te_std=sqrt(Te2_mean - Te_mean^2);
    Teth_3=Te_mean+3*Te_std;
    
    B_mean=mean(B,'all');
    B2_mean=mean(B.*B, 'all');
    B_std=sqrt(B2_mean - B_mean^2);
    Bth_m3=B_mean-3*B_std;    
    E_par_22_mean=mean(E_par_22, 'all');
    E_par_44_mean=mean(E_par_22.^2, 'all');
    E_par_22_std=sqrt(E_par_44_mean - E_par_22_mean^2);
    %E_par_22_th_m3=E_par_22_mean-3*E_par_22_std;    
    E_par_22_th_32=E_par_22_mean+3.2*E_par_22_std;        
    %E_par2_mean=mean(E_par2,'all');
    %E_par4_mean=mean(E_par2.*E_par2, 'all');
    %E_par2_std=sqrt(E_par4_mean - E_par2_mean^2);
    %E_par2th_3=E_par2_mean+3*E_par2_std;        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First to calculate and save the individual statistics of the variables along
    % the box    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Current, velocities and Temperatures
    J2_above = J2 > J2th_3;    %J2_above_values=J2.*J2_above;
    Vi2_above = Vi2 > Vi2th_3;
    Ve2_above = Ve2 > Ve2th_3;
    Ti_above = Ti > Tith_3;
    Te_above = Te > Teth_3;
    
    J2_stats = regionprops3(J2_above,'all');
    Vi2_stats = regionprops3(Vi2_above,'all');
    Ve2_stats = regionprops3(Ve2_above,'all');
    Ti_stats = regionprops3(Ti_above,'all');
    Te_stats = regionprops3(Te_above,'all');
    
    clear J2_above Vi2_above Ve2_above Ti_above Te_above
    
    cd '/Volumes/PSC_DiRAC_DATA/Geometric_stats';
    filename_1 = 'J2Vie2Tie2_geometry';
    save(filename_1,'J2_stats');
    save(filename_1,'Ve2_stats','-append','-nocompression')
    save(filename_1,'Te_stats','-append','-nocompression')  
    
    %%% Magnetic field and parallel componentts of the electric field
    B_zero = B < Bth_m3;
    %E_par_zero = E_par < E_par_22_th_m3; %This seens to be not well define. The code takes too much
    E_par_above = E_par_22 > E_par_22_th_32;
    %A=find(E_par_22 > (E_par_22_mean + 3.2*E_par_22_std));
    
    Bze_stats = regionprops3(B_zero,'all');
    %Epar_zero_stats = regionprops3(E_par_zero,'all');
    Epar_above_stats = regionprops3(E_par_above,'all');
    %Epar2_above_stats = regionprops3(E_par2_above,'all');    
    
    filename_3 = 'BzeEpa_geometry';
    save(filename_3,'Bze_stats');
    save(filename_3,'Epar_above_stats','-append','-nocompression')
    
    clear B_zero 
    clear Epar_above_stats
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now calculate and save the statistics of the intersections
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Current, velocities and Temperatures
    J2Vie2Tie = J2 > J2th_3 & Vi2 > Vi2th_3 & Ti > Tith_3 & Ve2 > Ve2th_3 & Te > Teth_3; 
    %[rJVieTie,cJVieTie,vJVieTie] = ind2sub(size(J2),find(J2 > J2th_3 & Vi2 > vi2th_3 & Ti2 > Ti2th_3 & Ve2 > ve2th_3 & Te2 > Te2th_3)); AJVieTie = [rJVieTie,cJVieTie,vJVieTie];
    J2Vie2Tie_stats = regionprops3(J2Vie2Tie,'all');
    J2Vie2 = J2 > J2th_3 & Vi2 > Vi2th_3 & Ve2 > Ve2th_3; 
    J2Vie2_stats = regionprops3(J2Vie2,'all');
    J2Tie = J2 > J2th_3 & Ti > Tith_3 & Te > Teth_3; 
    J2Tie_stats = regionprops3(J2Tie,'all');
    
    filename_2 = 'Intersection_J2Vie2Tie2_geometry';
    save(filename_2,'J2Vie2Tie_stats');
    save(filename_2,'J2Vie2_stats','-append','-nocompression')
    save(filename_2,'J2Tie_stats','-append','-nocompression')
    
    %%% Magnetic field and parallel componentts of the electric field    
    BzeEpa = B < Bth_m3 & E_par_22 > E_par_22_th_32;
    BzeEpa_stats = regionprops3(BzeEpa,'all');
     
    %%% Magnetic field, Epar, J %%%%%%% 
    BzeJ2 = J2 > J2th_3 & B < Bth_m3;
    BzeJ2_stats = regionprops3(BzeJ2,'all');  
    EparJ2 = J2 > J2th_3 & E_par_22 > E_par_22_th_32;
    EparJ2_stats = regionprops3(EparJ2,'all');
    BzeEparJ2 = J2 > J2th_3 & E_par_22 > E_par_22_th_32 & B < Bth_m3;
    BzeEparJ2_stats = regionprops3(BzeEparJ2,'all');
    
    filename_4 = 'Intersection_BzeEpa_geometry';
    save(filename_4,'BzeEpa_stats');
    save(filename_4,'BzeJ2_stats','-append','-nocompression')
    save(filename_4,'EparJ2_stats','-append','-nocompression')
    save(filename_4,'BzeEparJ2_stats','-append','-nocompression')
    
    clear BzeEpa_stats 
    
    % To make the statistics about the geometrical properties
    %Lstats = regionprops3(B_zero,'all'); L=bwlabel(B_zero(:,:,1)); CC = bwconncomp(B_zero); Cstats = regionprops3(CC);    
    %%% Load the file 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Contents of workspace before loading file:')
    whos
    
    disp('Contents of J2Vie2Tie2_geometry.mat:')
    whos('-file','J2Vie2Tie2_geometry.mat')
    
    disp('Contents of Intersection_J2Vie2Tie2_geometry.mat:')
    whos('-file','Intersection_J2Vie2Tie2_geometry.mat')
    
    disp('Contents of BzeEpa_geometry.mat:')
    whos('-file','BzeEpa_geometry.mat')
    
    disp('Contents of Intersection_BzeEpa_geometry.mat:')
    whos('-file','Intersection_BzeEpa_geometry.mat')
    
    load('J2Vie2Tie2_geometry.mat')
    disp('Contents of workspace after loading file:')
    whos    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Multiplying the indexes of J2, B2, I can find the comun regions
    % Making plots
    %----------------------------------------------------------------------
    
    f1=figure(1);
    contourf(J2_above(:,:,10))
    colorbar
    f2=figure(2);
    contourf(B2_above(:,:,10))
    colormap(hot)
    colorbar
    
    f3=figure(3);
    contourf(intersect_JB(:,:,10))
    colorbar
    
    f4=figure(4);
    pcolor(rJB2a,cJB2a,vJB2a)
    
   
    f5=figure(5);
    pcolor(E_par_zero(:,:,1000))
    colormap('gray');
    
    f6=figure(6);
    pcolor(E_par_above(:,:,1000))
    colormap('gray');
    
    f7=figure(7);
    pcolor(E_par_above(:,:,1000))
    colormap('gray');
    
    f8=figure(8);
    pcolor(B_zero_Epar_above(:,:,1000))
     
    f9=figure(9);
    pcolor(B_zero_Epar2_above(:,:,1000))
    
    pcolor(B_zero(:,:,1000))
    
   
    f11=figure(11);
    contourf(B_zero(:,:,183))
    hold on
    scatter(62.4613423373760,392.601640022051)
    
    f12=figure(12);
    contourf(squeeze(B_zero(:,:,74)))
    hold on
    scatter(191.599143698154,394.816519489787)
    
  
    f13=figure(13);
    contourf(squeeze(B_zero(190,:,:)))
    hold on
    scatter(73.8454494098059, 394.816519489787)
    %{191.599143698154,394.816519489787,73.8454494098059}
    
    
    
    
    f1=figure(1);
    [C,h]=contourf(spectraB.kpar,spectraB.kper,log10(spectraB.P2D'),'LevelStep',0.5);
    %[C]=contourf(spectraB.kper,spectraB.kpar,log10(spectraE.P2D),'LevelStep',0.5);
    set(h,'LineColor','none')
    lim = caxis; 
    caxis([0 13]);
    colormap(jet)
    colorbar
    set(gca,'XScale','log','YScale','log','FontSize',18)
    xlabel('$k_\| d_i$','Interpreter','latex')
    ylabel('$k_\perp d_i$','Interpreter','latex')
    title('$P_{\mathbf{B}}2D$','Interpreter','latex')
    x=log10([0.1 1]);
    x3=log10([1 70]);
    y=(3/2)*x;
    y3=(3)*x3;
    hold on
    loglog(10.^x,10.^y,'--k')
    loglog(10.^x3,10.^y3,'--b')
    ylim([0.28 52])
    txt = '\bf{3/2 \rightarrow}';
    txt2 = '\bf{\leftarrow 3}';
    text(0.2,0.5,txt)
    text(5,28,txt2)
    hold off
    
    %{
    
    %}
    
    % this is to plot parallel and perpendicualr in the same plot even they
    % have diffenret dimensions
    tlin=spectraB.kpar;
    xlin=spectraB.P1Dpar/rms(spectraB.P1Dpar);
    ttlin=spectraB.kper;
    ylin=spectraB.P1Dper/rms(spectraB.P1Dper);
    xxlin=interp1(tlin, xlin, ttlin, 'linear','extrap');	
    
    
    %----------------------------------------------------------------------
    %cd /disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics/2D_1D_PSD
    cd '/Volumes/PSC_DiRAC_DATA/tmp_images';
    % Save the plots
    %-------------------------------------------------------------------------
    saveas(f1,strcat(S(i).name,'2Dspectrum_Bpd_CB04.png'));


     
    %cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'; 
    cd '/Volumes/PSC_DiRAC_DATA/DATACB104';
    
end
