%#-----------------------------------------------------------------------
%This script is to compute the shapes based on thresholds. Another way 
% to measure how is the distribution of scales by changing the intensity 
% in the magnetic field first and then in the current density 
%#----------------------------------------------------------------------

cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data'; %This is important because the xdmf files are in that directory
path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
i=19;

%for i=5:6
    disp(strcat('Computing step ...',S(i).name))
    
    % Read files 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    
    %Load the variables and make the arrays
    B_1st = info.Groups(18).Name;
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    B = sqrt(Bx.^2 + By.^2 + Bz.^2) ;
    clearvars Bx By Bz
    B_mean=mean(B,'all');
    B2_mean=mean(B.*B, 'all');
    B_std=sqrt(B2_mean - B_mean^2);
    Bmax=max(max(max(B)));
    Bmin=min(min(min(B)));
    
    % Choose the number of bins to make the histogram in intensity
    nbins=20; % number of bins
    deltai=(Bmax-Bmin)/nbins;
    
    % Make the bins
    binsB=zeros(nbins,1); %bin the intensity
    for i=1:nbins
        binsB(i)=i*deltai;
    end
    
    % Apply the filter by bin and perform the shape statistics 
    j=0;
    B_Inten_stats=cell(nbins,1); %Create the cell table to store the tables
    for i=1:nbins
        %B_Inten = j*deltai < B < i*deltai; This is binning in values 
        B_Inten = i*deltai < B ; %This is above
        j=i;
        B_Inten_stats_i = regionprops3(B_Inten,'Volume','Centroid','EquivDiameter',...
            'PrincipalAxisLength','Orientation','ConvexVolume','Solidity','SurfaceArea');
        dummy_stats=B_Inten_stats_i;
        dummy_stats(dummy_stats.Volume < 33, :) = []; %(de^3 4pi/3)
        B_Inten_stats{i}=dummy_stats; %B_Inten_stats_i;
        clear B_Inten_stats_i
    end
    
    % Look at the min and max of rpar, D and Asp to what edges should be
    % define for the histogram in distances and Asp
    for i=2:nbins
        Bxy_pal=B_Inten_stats{i}.PrincipalAxisLength;
        rpar=Bxy_pal(:,1);
        rparmax=0.06*max(rpar);
        rparmin=0.06*min(rpar);
        Diameter=sqrt(Bxy_pal(:,2).^2 + Bxy_pal(:,3).^2);
        Dmax=0.06*max(Diameter);
        Dmin=0.06*min(Diameter);
        Asp_ratio=rpar./Diameter;
        Aspmax=max(Asp_ratio);
        Aspmin=min(Asp_ratio);
        result = [rparmin rparmax Dmin Dmax Aspmin Aspmax];
        disp(result)       
    end
    
    % Define the edges
    edgespar = linspace(0,140,81);
    edgesD = linspace(0,40,41);
    edgesAsp = linspace(0,50,41);
    
    % Create the histograms in distances and Asp
    %A =(1:nbins)';
    %T=cell(nbins,1);
    AA_rpar=zeros(nbins,length(edgespar)-1);
    AA_D=zeros(nbins,length(edgesD)-1);
    AA_Asp=zeros(nbins,length(edgesAsp)-1);
    
    for i = 1:nbins-1 
        Bxy_pal=B_Inten_stats{i}.PrincipalAxisLength;
        rpar=Bxy_pal(:,1);
        Diameter=sqrt(Bxy_pal(:,2).^2 + Bxy_pal(:,3).^2);
        Asp_ratio=rpar./Diameter;
     
        Asp_mean=mean(Asp_ratio); Asp_std=std(Asp_ratio);
        rpar_mean=mean(rpar*0.06); rpar_std=std(rpar*0.06);
        Diameter_mean=mean(Diameter*0.06); Diameter_std=std(Diameter*0.06);
        
        [N_rpar,ed_rpar] = histcounts(rpar,edgespar);
        [N_D,ed_D] = histcounts(Diameter,edgesD);
        [N_Asp,ed_Asp] = histcounts(Asp_ratio,edgesAsp);
        
        AA_rpar(i,:)=N_rpar;
        AA_D(i,:)=N_D;
        AA_Asp(i,:)=N_Asp;
        
        %{
        T_i = table(1);
        T_i.N_rpar = {N_rpar}; 
        T_i.ed_rpar = {ed_rpar};
        T_i.N_D = {N_D}; 
        T_i.ed_D = {ed_D};
        T_i.N_Asp = {N_Asp}; 
        T_i.ed_Asp = {ed_Asp};
        T{i}=T_i;
        clearvars T_i
        %}
    end
    
    B0=0.1;
    AA_rpar(AA_rpar==0)=NaN;    
    f1=figure(1);
    h=pcolor(ed_rpar(1:80)',((binsB/(B0)).^2)',AA_rpar);
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    colormap(jet(20))
    set(gca,'ColorScale','lin','FontSize',18)
    colorbar
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$L_{\parallel}/d_{i}$$','Interpreter','latex','FontSize',20)
    ylabel('$$|B/B_{0}|^2$$','Interpreter','latex','FontSize',20)
    hold off

    AA_D(AA_D==0)=NaN;    
    f2=figure(2);
    h=pcolor(ed_D(1:40)',((binsB/(B0)).^2)',AA_D);
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    colormap(jet(20))
    set(gca,'ColorScale','lin','FontSize',18)
    colorbar
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$L_{D}/d_{i}$$','Interpreter','latex','FontSize',20)
    ylabel('$$|B/B_{0}|^2$$','Interpreter','latex','FontSize',20)
    hold off
    
    AA_Asp(AA_Asp==0)=NaN;
    f3=figure(3);
    h=pcolor(ed_Asp(1:40)',((binsB/(B0)).^2)',AA_Asp);
    set(h, 'EdgeColor', 'none');
    hold on
    lim = caxis;  
    colormap(jet(20))
    set(gca,'ColorScale','lin','FontSize',18)
    colorbar
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$L_{\parallel}/L_{D}$$','Interpreter','latex','FontSize',20)
    ylabel('$$|B/B_{0}|^2$$','Interpreter','latex','FontSize',20)
    hold off
     
    %zzz=cell2mat(T{i}.N_rpar);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First to calculate and save the individual statistics of the variables along
    % the box    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    cd '/Volumes/PSC_DiRAC_DATA';
    saveas(f1,'B2_Lpar_Above.png');
    saveas(f2,'B2_LD_Above.png');
    saveas(f3,'B2_Asp_Above.png');
      
%end
