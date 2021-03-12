cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Electron_project/electron_project5';
path = '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/Electron_project/electron_project5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mime=200;
%S = dir(fullfile(path,'prt.*'));
S = dir(fullfile(path,'pfd.*'));
N=numel(S); % number of files to use
H=zeros(N);
%i=2; %prt 800
for i=2:12
    disp(strcat('Computing step ...',S(i).name)) 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
    %T_1st = info.Groups(1).Name;
    %E_1st = info.Groups(17).Name;
    %B_1st = info.Groups(18).Name;
    %J_1st = info.Groups(19).Name;
    V_1st = info.Groups(25).Name;
    n_1st = info.Groups(23).Name;
    %diary off
    %s2='/hx/p0/3d';
    %sx = strcat(B_1st,s2);
    %Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    %By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    %Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    %Ex=h5read(fileID,strcat(E_1st,'/ex/p0/3d'));
    %Ey=h5read(fileID,strcat(E_1st,'/ey/p0/3d'));
    %Ez=h5read(fileID,strcat(E_1st,'/ez/p0/3d'));
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'));
    ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    
    
    
    %calculate perpendicular
    vxyi=sqrt((vix(1,:,:).*vix(1,:,:)).*(viy(1,:,:).*viy(1,:,:)));
%    clearvars pxi
    vxye=sqrt((vex(1,:,:).*vex(1,:,:)).*(vey(1,:,:).*vey(1,:,:)));
 %   clearvars pxe
    
    %
    f1=figure(1);
    h = histogram2(viz,vxyi,'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    c = h.NumBins;
    h.NumBins = [150 150];
    lim = caxis;  
    %caxis([0.01 1e4]);
    colormap(bone(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-0.006 0.006]); ylim([0 3.5e-5]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{i\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{i\perp}$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{i}$$'+string(i),'Interpreter','latex','FontSize',20) 
    
    clearvars h
    
    f12=figure(12);
    h = histogram2(viz(1,:,:),viy(1,:,:),'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    c = h.NumBins;
    h.NumBins = [150 150];
    lim = caxis;  
    %caxis([0.01 1e4]);
    colormap(bone(25))
    set(gca,'ColorScale','log','FontSize',20)
    colorbar
    xlim([-7e-3 7e-3]); ylim([-7e-3 7e-3]);
    set(gca,'XScale','lin','YScale','lin','FontSize',20)
    xlabel('$$v_{i\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{iy}$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{i}$$'+string(i),'Interpreter','latex','FontSize',20)
    
    clearvars h
    %}
    %
    
    f2=figure(2);
    h = histogram2(vez,vxye,'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    c = h.NumBins;
    h.NumBins = [150 150];
    lim = caxis;  
    %caxis([0.01 1e4]); %caxis([0.01 1e4]); for the other 2 and 5
    colormap(bone(25))
    set(gca,'ColorScale','log','FontSize',18)
    colorbar
    xlim([-0.08 0.08]); ylim([0 3e-3]);
    set(gca,'XScale','lin','YScale','lin','FontSize',18)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{e}$$'+string(i),'Interpreter','latex','FontSize',20)
   
    clearvars h
    
    f22=figure(22);
    h = histogram2(vez(1,:,:),vey(1,:,:),'DisplayStyle','tile','ShowEmptyBins','on');%,'Normalization','pdf');
    c = h.NumBins;
    h.NumBins = [150 150];
    lim = caxis;  
   % caxis([0.01 1e4]);
    colormap(bone(25))
    set(gca,'ColorScale','log','FontSize',20)
    colorbar
    xlim([-0.05 0.05]); ylim([-0.07 0.07]);
    set(gca,'XScale','lin','YScale','lin','FontSize',20)
    xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
    ylabel('$$v_{ey}$$','Interpreter','latex','FontSize',20)
    title('$$VDF_{e}$$'+string(i),'Interpreter','latex','FontSize',20)
 
    clearvars h
    %}
    saveas(f1,strcat('25_viparperp_'+string(i),'.png'));
    saveas(f2,strcat('25_veparperp_'+string(i),'.png'));
    saveas(f12,strcat('25_vipary_'+string(i),'.png'));
    saveas(f22,strcat('25_vepary_'+string(i),'.png'));
end
