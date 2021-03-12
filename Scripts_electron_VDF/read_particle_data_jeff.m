%#-----------------------------------------------------------------------
%This is the script to read the particle infoprmation
% Currently the program does not create a xdmf but only a .h5
%#----------------------------------------------------------------------

% This is the info for 10 particles
cd '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/whistler_1/particles';
path = '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/whistler_1/particles';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size of the arrays     256000 
% size i                   5325
% size e                   5325 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the info for 50 particles
cd '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/whistler_1/particles_50';
path = '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/whistler_1/particles_50';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size of the arrays    1280000
% size i                  23774
% size e                  23732 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% this to check the distribution of particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/my_case_particles';
path = '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/my_case_particles';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/test_electron_project/run_22_07_20';
path = '/Users/jeffersson_agudelo/Documents/PSC_here/PSC_This_mac/test_runs_here/test_electron_project/run_22_07_20';


cd '/Volumes/PSC_DiRAC_DATA/RUNS2020/checking_subset_particle/reduced_box';
path = '/Volumes/PSC_DiRAC_DATA/RUNS2020/checking_subset_particle/reduced_box';

cd '/Volumes/PSC_DiRAC_DATA/RUNS2020/checking_subset_particle';
path = '/Volumes/PSC_DiRAC_DATA/RUNS2020/checking_subset_particle';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296/nostrahl';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Kathleen/electron_16296/nostrahl';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = dir(fullfile(path,'prt.*'));
N=numel(S); % number of files to use
H=zeros(N);
i=1;
%for i=1:6
    disp(strcat('Computing step ...',S(i).name)) 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID); %this is the command to show the information in the file 
    info = h5info(fileID);
    prt = info.Groups(1).Name;
    %read the information from the fields of the struct
    prtfields=h5read(fileID,strcat(prt,'/p0/1d'));
    x8 = prtfields.x;
    y8 = prtfields.y;
    z8 = prtfields.z;
    px = prtfields.px;
    py = prtfields.py;
    pz = prtfields.pz;
    q = prtfields.q;
    m = prtfields.m;
    w = prtfields.w;
 
    %se parate the popullations
    q_i_index = find(q==1);
    q_e_index = find(q==-1);
    
    %Coordinates of ions
    xi=x(q_i_index); yi=y(q_i_index); zi=z(q_i_index);
    pxi=px(q_i_index); pyi=py(q_i_index); pzi=pz(q_i_index);
    qi=q(q_i_index); %mi=m(q_i_index); wi=w(q_i_index);
    
    %coordinates of electrons
    xe=x(q_e_index); ye=y(q_e_index); ze=z(q_e_index);
    pxe=px(q_e_index); pye=py(q_e_index); pze=pz(q_e_index);
    qe=q(q_e_index); %me=m(q_e_index); we=w(q_e_index);
    
    %TF = [x y z px py pz q m w];
    
    %bin the particles to see how many fill a condition
    % for example lets begin with the velocity 
    
    pm=sqrt(px.*px + py.*py + pz.*pz);
    pmi=sqrt(pxi.*pxi + pyi.*pyi + pzi.*pzi);
    pme=sqrt(pxe.*pxe + pye.*pye + pze.*pze);
    
    
    % Making plots
    %----------------------------------------------------------------------
    
    f1=figure(14);
    set(gca,'FontSize',18)
    scatter(x, px)
    ylabel('$p_{x}$','Interpreter','latex','FontSize',18)
    xlabel('$x$','Interpreter','latex','FontSize',18)
    f2=figure(24);
    set(gca,'FontSize',18)
    scatter(xi, pxi)
    ylabel('$p_{xi}$','Interpreter','latex','FontSize',18)
    xlabel('$x_{i}$','Interpreter','latex','FontSize',18)
    f3=figure(34);
    set(gca,'FontSize',18)
    scatter(xe, pxe)
    ylabel('$p_{xe}$','Interpreter','latex','FontSize',18)
    xlabel('$x_{e}$','Interpreter','latex','FontSize',18)
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
    
    f9=figure(94);
    histfit(pxe,100);
    
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
    

    
    %----------------------------------------------------------------------
    %cd /disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Statistics/2D_1D_PSD
    cd '/Volumes/PSC_DiRAC_DATA/tmp_images';
    % Save the plots
    %-------------------------------------------------------------------------
    saveas(f1,strcat(S(i).name,'2Dspectrum_Bpd_CB04.png'));
    %saveas(f12,strcat(S(i).name,'2Dspectrum_B_trans_CB04.png'));
    %saveas(f2,strcat(S(i).name,'2Dspectrum_E_CB04.png'));
    %saveas(f3,strcat(S(i).name,'2Dspectrum_vi_CB04.png'));
    %saveas(f4,strcat(S(i).name,'2Dspectrum_ve_CB04.png'));
    %saveas(f41,strcat(S(i).name,'2Dspectrum_ni_CB04.png'));
    saveas(f5,strcat(S(i).name,'1_perDspectrum_CB04_Bvipd_2.png'));
    saveas(f6,strcat(S(i).name,'1_parDspectrum_CB04_Bvipd_2.png'));
    %saveas(f7,strcat(S(i).name,'1_parDspectrum_CB04_Bparper.png'));    
    %cd '/disk/plasmaz/jaa/Critical_Balance/CB8waves_104_1/Raw_outputs'; 
    cd '/Volumes/PSC_DiRAC_DATA/DATACB104';
    
