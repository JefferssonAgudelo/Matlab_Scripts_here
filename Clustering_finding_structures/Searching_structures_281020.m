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
i=1;

%for i=5:6
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
    
    %
    %ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    %ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));    
    Ex=h5read(fileID,strcat(E_1st,'/ex/p0/3d'));
    Ey=h5read(fileID,strcat(E_1st,'/ey/p0/3d'));
    Ez=h5read(fileID,strcat(E_1st,'/ez/p0/3d'));
    B = sqrt(Bx.^2 + By.^2 + Bz.^2) ;
    %E = sqrt(Ex.^2 + Ey.^2 + Ez.^2) ; 
    E_par = (Ex.*Bx + Ey.*By + Ez.*Bz)./(B);
    %E_par_mag=sqrt(E_par.*E_par);
    
    %E_par2= ((Ex.*Bx).^2 + (Ey.*By).^2 + (Ez.*Bz).^2)./(E.*B).^2; %This is the one that must be different from zero
    clearvars Bx By Bz Ex Ey Ez
    %}
    Jx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    Jy=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    Jz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));    
    Ex=h5read(fileID,strcat(E_1st,'/ex/p0/3d'));
    Ey=h5read(fileID,strcat(E_1st,'/ey/p0/3d'));
    Ez=h5read(fileID,strcat(E_1st,'/ez/p0/3d'));
    J = sqrt(Jx.^2 + Jy.^2 + Jz.^2) ;
    E = sqrt(Ex.^2 + Ey.^2 + Ez.^2) ; 
    EJ = (Ex.*Jx + Ey.*Jy + Ez.*Jz);
    clearvars Jx Jy Jz Ex Ey Ez

    
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    Vi_m = sqrt(vix.^2 + viy.^2 + viz.^2) ;
    clearvars vix viy viz
    
    vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d')); 
    Ve_m = sqrt(vex.^2 + vey.^2 + vez.^2) ;
    clearvars vex vey vez
        
    Jx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    Jy=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    Jz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));
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
   
    Jm_mean=mean(J_m,'all');
    J2_mean=mean(J_m.*J_m, 'all');
    Jm_std=sqrt(J2_mean - Jm_mean^2);
    Jmth_3=Jm_mean+3*Jm_std;
    Jmth_4=Jm_mean+4*Jm_std;    

    Vim_mean=mean(Vi_m,'all');
    Vi2_mean=mean(Vi_m.*Vi_m, 'all');
    Vi_std=sqrt(Vi2_mean - Vim_mean^2);
    Vimth_3=Vim_mean+3*Vi_std;
    Vimth_4=Vim_mean+4*Vi_std;
    
    Vem_mean=mean(Ve_m,'all');
    Ve2_mean=mean(Ve_m.*Ve_m, 'all');
    Ve_std=sqrt(Ve2_mean - Vem_mean^2);
    Vemth_3=Vem_mean+3*Ve_std;
    Vemth_4=Vem_mean+4*Ve_std;
    
    Ti_mean=mean(Ti,'all');
    Ti2_mean=mean(Ti.*Ti, 'all');
    Ti_std=sqrt(Ti2_mean - Ti_mean^2);
    Tith_3=Ti_mean+3*Ti_std;
    Tith_4=Ti_mean+4*Ti_std;
    
    Te_mean=mean(Te,'all');
    Te2_mean=mean(Te.*Te, 'all');
    Te_std=sqrt(Te2_mean - Te_mean^2);
    Teth_3=Te_mean+3*Te_std;
    Teth_4=Te_mean+4*Te_std;
    
    %
    B_mean=mean(B,'all');
    B2_mean=mean(B.*B, 'all');
    B_std=sqrt(B2_mean - B_mean^2);
    Bth_m3=B_mean-3*B_std; 
   
    E_par_mean=mean(E_par, 'all');
    E_par_2_mean=mean(E_par.^2, 'all');
    E_par_std=sqrt(E_par_2_mean - E_par_mean^2);
    
    E_par_mag_mean=mean(E_par_mag, 'all');
    E_par_mag_2_mean=mean(E_par_mag.^2, 'all');
    E_par_mag_std=sqrt(E_par_mag_2_mean - E_par_mag_mean^2);
    
    E_parth_4min=E_par_mean-4*E_par_std;    
    E_parth_4plu=E_par_mean+4*E_par_std;
    E_parth_3min=E_par_mean-3*E_par_std;    
    E_parth_3plu=E_par_mean+3*E_par_std;
    E_parth_2min=E_par_mean-2*E_par_std;    
    E_parth_2plu=E_par_mean+2*E_par_std;
    
    Epar_below4 = E_par < E_parth_4min;
    Epar_min4_stats = regionprops3(Epar_below4,'all');
    Epar_above4 = E_par > E_parth_4plu;
    Epar_plu4_stats = regionprops3(Epar_above4,'all');
    
    Epar_below3 = E_par < E_parth_3min;
    Epar_min3_stats = regionprops3(Epar_below3,'all');
    Epar_above3 = E_par > E_parth_3plu;
    Epar_plu3_stats = regionprops3(Epar_above3,'all');  
    
    Epar_below2 = E_par < E_parth_2min;
    Epar_min2_stats = regionprops3(Epar_below2,'all');
    Epar_above2 = E_par > E_parth_2plu;
    Epar_plu2_stats = regionprops3(Epar_above2,'all');  

    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    filename_1 = 'Epar_plu_min4_geometry';
    save(filename_1,'Epar_plu4_stats');
    save(filename_1,'Epar_min4_stats','-append','-nocompression') 

    
    EJ_mean=mean(EJ, 'all');
    EJ_2_mean=mean(EJ.^2, 'all');
    EJ_std=sqrt(EJ_2_mean - EJ_mean^2);
    
    EJth_4min=EJ_mean-4*EJ_std;    
    EJth_4plu=EJ_mean+4*EJ_std;
    EJth_3min=EJ_mean-3*EJ_std;    
    EJth_3plu=EJ_mean+3*EJ_std;
    EJth_2min=EJ_mean-2*EJ_std;    
    EJth_2plu=EJ_mean+2*EJ_std;
    
    EJ_below4 = EJ < EJth_4min;
    EJ_min4_stats = regionprops3(EJ_below4,'all');
    EJ_above4 = EJ > EJth_4plu;
    EJ_plu4_stats = regionprops3(EJ_above4,'all');
    
    EJ_below3 = EJ < EJth_3min;
    EJ_min3_stats = regionprops3(EJ_below3,'all');
    EJ_above3 = EJ > EJth_3plu;
    EJ_plu3_stats = regionprops3(EJ_above3,'all');   
    
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    filename_1 = 'JE_plu_min4_geometry';
    save(filename_1,'EJ_plu4_stats');
    save(filename_1,'EJ_min4_stats','-append','-nocompression')
    
    
    %EJ_below2 = EJ < EJth_2min;
    %EJ_min2_stats = regionprops3(EJ_below2,'all');
    %EJ_above2 = EJ > EJth_2plu;
    %EJ_plu2_stats = regionprops3(EJ_above2,'all');   
    
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First to calculate and save the individual statistics of the variables along
    % the box    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Current, velocities and Temperatures
    Jm_above = J_m > Jmth_4;    %J2_above_values=J2.*J2_above;
    Vim_above = Vi_m > Vimth_4;
    Vem_above = Ve_m > Vemth_4;
    Ti_above = Ti > Tith_4;
    Te_above = Te > Teth_4;
    
    Jm4_stats = regionprops3(Jm_above,'all');
    Vim4_stats = regionprops3(Vim_above,'all');
    Vem4_stats = regionprops3(Vem_above,'all');
    Ti4_stats = regionprops3(Ti_above,'all');
    Te4_stats = regionprops3(Te_above,'all');
    
    clear Jm_above Vim_above Vem_above Ti_above Te_above
    
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    filename_1 = 'JmViemTiem4_geometry';
    save(filename_1,'Jm4_stats');
    save(filename_1,'Vim4_stats','-append','-nocompression')
    save(filename_1,'Vem4_stats','-append','-nocompression')
    save(filename_1,'Ti4_stats','-append','-nocompression') 
    save(filename_1,'Te4_stats','-append','-nocompression') 
       
    %%% Current, velocities and Temperatures
    JmViemTie4 = J_m > Jmth_4 & Vi_m > Vimth_4 & Ti > Tith_4 & Ve_m > Vemth_4 & Te > Teth_4; 
    JmViemTie4_stats = regionprops3(JmViemTie4,'all');
    filename_2 = 'Int_JmViemTie4_geometry';
    save(filename_2,'JmViemTie4_stats');
    JmViemTie3 = J_m > Jmth_3 & Vi_m > Vimth_3 & Ti > Tith_3 & Ve_m > Vemth_3 & Te > Teth_3; 
    JmViemTie3_stats = regionprops3(JmViemTie3,'all');
    filename_2 = 'Int_JmViemTie3_geometry';
    save(filename_2,'JmViemTie3_stats');
    
    % Current and velocities and Current and temperatures  
    JmViem4 = J_m > Jmth_4 & Vi_m > Vimth_4 & Ve_m > Vemth_4; 
    JmViem4_stats = regionprops3(JmViem4,'all');
    JmTie4 = J_m > Jmth_4 & Ti > Tith_4 & Te > Teth_4; 
    JmTie4_stats = regionprops3(JmTie4,'all');
    %save(filename_2,'JmViem4_stats','-append','-nocompression')
    %save(filename_2,'JmTie4_stats','-append','-nocompression')
    JmViem3 = J_m > Jmth_3 & Vi_m > Vimth_3 & Ve_m > Vemth_3; 
    JmViem3_stats = regionprops3(JmViem3,'all');
    JmTie3 = J_m > Jmth_3 & Ti > Tith_3 & Te > Teth_3; 
    JmTie3_stats = regionprops3(JmTie3,'all');
    %save(filename_2,'JmViem3_stats','-append','-nocompression')
    %save(filename_2,'JmTie3_stats','-append','-nocompression')
    %}
    
    
    %Current, velocities, temperatures and Epar 
    JmViemTieEpar_plu4 = J_m > Jmth_4 & Vi_m > Vimth_4 & Ti > Tith_4 & ...
        Ve_m > Vemth_4 & Te > Teth_4 & E_par > E_parth_4plu; 
    JmViemTieEpar_plu4_stats = regionprops3(JmViemTieEpar_plu4,'all');
  
    JmViemTieEpar_min4 = J_m > Jmth_4 & Vi_m > Vimth_4 & Ti > Tith_4 & ...
        Ve_m > Vemth_4 & Te > Teth_4 & E_par < E_parth_4min; 
    JmViemTieEpar_min4_stats = regionprops3(JmViemTieEpar_min4,'all');
    
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    filename_1 = 'JmViemTieEpar_plu_min4_geometry';
    save(filename_1,'JmViemTieEpar_plu4_stats');
    save(filename_1,'JmViemTieEpar_min4_stats','-append','-nocompression')
    
    
    %Current, velocities, temperatures and EJ 
    JmViemTieEJ_plu4 = J_m > Jmth_4 & Vi_m > Vimth_4 & Ti > Tith_4 & ...
        Ve_m > Vemth_4 & Te > Teth_4 & EJ > EJth_4plu; 
    JmViemTieEJ_plu4_stats = regionprops3(JmViemTieEJ_plu4,'all');
  
    JmViemTieEJ_min4 = J_m > Jmth_4 & Vi_m > Vimth_4 & Ti > Tith_4 & ...
        Ve_m > Vemth_4 & Te > Teth_4 & EJ < EJth_4min; 
    JmViemTieEJ_min4_stats = regionprops3(JmViemTieEJ_min4,'all');
   
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    filename_1 = 'JmViemTieEJ_plu_min4_geometry';
    save(filename_1,'JmViemTieEJ_plu4_stats');
    save(filename_1,'JmViemTieEJ_min4_stats','-append','-nocompression')
  
    
    
    %Current, velocities, temperatures and EJ 3
    JmViemTieEJ_plu3 = J_m > Jmth_3 & Vi_m > Vimth_3 & Ti > Tith_3 & ...
        Ve_m > Vemth_3 & Te > Teth_3 & EJ > EJth_3plu; 
    JmViemTieEJ_plu3_stats = regionprops3(JmViemTieEJ_plu3,'all');
  
    JmViemTieEJ_min3 = J_m > Jmth_3 & Vi_m > Vimth_3 & Ti > Tith_3 & ...
        Ve_m > Vemth_3 & Te > Teth_3 & EJ < EJth_3min; 
    JmViemTieEJ_min3_stats = regionprops3(JmViemTieEJ_min3,'all');
   
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
    filename_1 = 'JmViemTieEJ_plu_min3_geometry';
    save(filename_1,'JmViemTieEJ_plu3_stats');
    save(filename_1,'JmViemTieEJ_min3_stats','-append','-nocompression')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Geometric_stats';
 
    %dummy_stats=JmViemTie4_stats; %select the variables to analyse
    dummy_stats=JmViemTieEJ_plu3_stats;
    dummy_stats(dummy_stats.Volume < 16*16*16, :) = [];   %filter
    
    dummy_pal=dummy_stats.PrincipalAxisLength;
    dummy_angle=dummy_stats.Orientation;
    dummy_positions=0.06*dummy_stats.Centroid;
    dummy_equidiameter=0.06*dummy_stats.EquivDiameter;
    
    dummy_image=dummy_stats.Image;
    dummy_image1=cell2mat(dummy_stats.Image(1));
    
    f11=figure(11);
    volshow(cell2mat(dummy_stats.Image(1)))
    f12=figure(12);
    volshow(cell2mat(dummy_stats.Image(2)))
    f13=figure(13);
    volshow(cell2mat(dummy_stats.Image(3)))
    f14=figure(14);
    volshow(cell2mat(dummy_stats.Image(4)))
    f15=figure(15);
    volshow(cell2mat(dummy_stats.Image(5)))
    f16=figure(16);
    volshow(cell2mat(dummy_stats.Image(6)))
    
    
    
    dum2 = permute(E_par,[2 1 3]);
    dum2 = dum2(:)';
    f12=figure(12);
    h=histogram(dum2,100);
    h.Normalization='pdf';
    h.FaceColor='b';
    hold on
    xline(0,'--k')
    hold off
 
    
 
    
%end
