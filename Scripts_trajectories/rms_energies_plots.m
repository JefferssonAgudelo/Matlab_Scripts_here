cd '/Volumes/PSC_DiRAC_DATA/RUNS2020/REAL_8_small_08_cor';
path = '/Volumes/PSC_DiRAC_DATA/RUNS2020/REAL_8_small_08_cor';

S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
Brms=zeros(N,1);
Jrms=zeros(N,1);
Vierms=zeros(N,1);
Bmean=zeros(N,1);
Jmean=zeros(N,1);
Viemean=zeros(N,1);


for i=1:N
    disp(strcat('Computing step ...',S(i).name))
    % Read files 
    fileID =  S(i).name; %change the 1 per i
    %h5disp(fileID);
    info = h5info(fileID);
  
    B_1st = info.Groups(18).Name;
    V_1st = info.Groups(25).Name;
    n_1st = info.Groups(23).Name;
    J_1st = info.Groups(19).Name;
    
    Bx=h5read(fileID,strcat(B_1st,'/hx/p0/3d'));
    By=h5read(fileID,strcat(B_1st,'/hy/p0/3d'));
    Bz=h5read(fileID,strcat(B_1st,'/hz/p0/3d'));
    B = sqrt(Bx.^2 + By.^2 + Bz.^2) ; % this is the one for the velocity    
    clearvars Bx By Bz 
    B_mean=mean(B,'all');
    B2_mean=mean(B.*B, 'all');
    Brms_i=sqrt(B2_mean - B_mean^2);
    clearvars B
    
    jx=h5read(fileID,strcat(J_1st,'/jx/p0/3d'));
    jy=h5read(fileID,strcat(J_1st,'/jy/p0/3d'));
    jz=h5read(fileID,strcat(J_1st,'/jz/p0/3d'));
    J = sqrt(jx.^2 + jy.^2 + jz.^2) ; % this is the one for the velocity    
    clearvars jx jy jz 
    J_mean=mean(J,'all');
    J2_mean=mean(J.*J, 'all');
    Jrms_i=sqrt(J2_mean - J_mean^2);
    clearvars J
    
    %
    vix=h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'));
    viy=h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'));
    viz=h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'));
    vex=h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'));
    vey=h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'));
    vez=h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'));
    ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'));
    ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'));
    
    viex=(vix.*ni + vex.*ne)./(ni+ne); 
    viey=(viy.*ni + vey.*ne)./(ni+ne);
    viez=(viz.*ni + vez.*ne)./(ni+ne);
    clearvars vix viy viz vex vey vez ni ne
    Vie = sqrt(viex.^2 + viey.^2 + viez.^2) ; % this is the one for the velocity    
    clearvars viex viey viez 
    Vie_mean=mean(Vie,'all');
    Vie2_mean=mean(Vie.*Vie, 'all');
    Vierms_i=sqrt(Vie2_mean - Vie_mean^2);
    clearvars Vie      
    
    Brms(i)=Brms_i;
    Jrms(i)=Jrms_i;
    Vierms(i)=Vierms_i;
    
    Bmean(i)=B_mean;
    Jmean(i)=J_mean;
    Viemean(i)=Vie_mean;
    
end

dt=0.0322;
out_t=500;
tt=out_t*dt;
t=[0,1,2,3,4,5,6]*tt;    

%-------------------------------------------------------------------------
% Load the file BviJE3 to make the plots. That file has the data already 
%-------------------------------------------------------------------------
t = table2array(stats(1:33,1)); 
t = t(~isnan(t));
Brms = table2array(stats(:,9));
Brms = Brms(~isnan(Brms));
Jrms = table2array(stats(:,10));
Jrms = Jrms(~isnan(Jrms));
virms = table2array(stats(:,11));
virms = virms(~isnan(virms));

%-------------------------------------------------------------------------
% Load the file BviJE3 to make the plots. That file has the data already 
%-------------------------------------------------------------------------
t=BviJE3(1:33,7);
Brms = BviJE3(1:33,1); % Magnitudes
Jrms = BviJE3(1:33,5);
virms = BviJE3(1:33,3);

Brms = BviJE3(1:33,2); % real rms
Jrms = BviJE3(1:33,6);
virms = BviJE3(1:33,4);




f1=figure(1);
    plot(t,Brms,'k', t, Jrms, 'b',t, virms, 'r','LineWidth',1.5) 
    %plot(t,Brms/rms(Brms),'k', t, Jrms/rms(Jrms), 'b',t, virms/rms(virms), 'r','LineWidth',1.5)
    set(gca,'FontSize',18)
    %xlim([0.03 52.36]); ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    x3=120;
    plot([x3 x3],y1,'--', 'LineWidth',1.5)
    % Text in the plot
    xlabel('$\omega_{pi}t$','Interpreter','latex','FontSize',28)
    ylabel('$\psi^{rms}$','Interpreter','latex','FontSize',28)  
    %ylabel('$\langle \psi \rangle$','Interpreter','latex')  
    txt8 = '$\omega_{pi}t=120$';
    %text(100 , 0.13,txt8,'Interpreter','latex')
    %title('$P1D_\|$','Interpreter','latex')
    legend({'$B^{rms}$','$J^{rms}$','$v_{i}^{rms}$'},'Interpreter','latex','FontSize',25)
    %legend({'$\langle B \rangle$','$\langle J \rangle$','$\langle v_{i} \rangle$'},'Interpreter','latex')
    xlim([0 480]);
    hold off
    
    f1=figure(1);
    plot(t,Brms,'-*k', t, Jrms, '-*b',t, Vierms, '-*r','LineWidth',1.5) 
    %plot(t,Brms/rms(Brms),'k', t, Jrms/rms(Jrms), 'b',t, virms/rms(virms), 'r','LineWidth',1.5)
    set(gca,'FontSize',18)
    %xlim([0.03 52.36]); ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    x3=120;
    plot([x3 x3],y1,'--', 'LineWidth',1.5)
    % Text in the plot
    xlabel('$\omega_{pi}t$','Interpreter','latex','FontSize',28)
    ylabel('$\psi^{rms}$','Interpreter','latex','FontSize',28)  
    %ylabel('$\langle \psi \rangle$','Interpreter','latex')  
    txt8 = '$\omega_{pi}t=120$';
    %text(100 , 0.13,txt8,'Interpreter','latex')
    %title('$P1D_\|$','Interpreter','latex')
    legend({'$B^{rms}$','$J^{rms}$','$v^{rms}$'},'Interpreter','latex','FontSize',25)
    %legend({'$\langle B \rangle$','$\langle J \rangle$','$\langle v_{i} \rangle$'},'Interpreter','latex')
    xlim([0 100]);
    hold off
    
    f2=figure(2);
    plot(t,Bmean,'-*k', t, Jmean, '-*b',t, Viemean, '-*r','LineWidth',1.5) 
    %plot(t,Brms/rms(Brms),'k', t, Jrms/rms(Jrms), 'b',t, virms/rms(virms), 'r','LineWidth',1.5)
    set(gca,'FontSize',18)
    %xlim([0.03 52.36]); ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    x3=120;
    plot([x3 x3],y1,'--', 'LineWidth',1.5)
    % Text in the plot
    xlabel('$\omega_{pi}t$','Interpreter','latex','FontSize',28)
    ylabel('$\langle \psi \rangle$','Interpreter','latex','FontSize',28)  
    %ylabel('$\langle \psi \rangle$','Interpreter','latex')  
    txt8 = '$\omega_{pi}t=120$';
    %text(100 , 0.13,txt8,'Interpreter','latex')
    %title('$P1D_\|$','Interpreter','latex')
    legend({'$\langle B \rangle$','$\langle J \rangle$','$\langle v \rangle$'},'Interpreter','latex','FontSize',25)
    %legend({'$\langle B \rangle$','$\langle J \rangle$','$\langle v_{i} \rangle$'},'Interpreter','latex')
    xlim([0 100]);
    hold off
    
    % Save the plots
    cd '/Volumes/PSC_DiRAC_DATA';
    %-------------------------------------------------------------------------
    saveas(f1,'rms_CB103_BJV.png');
    saveas(f2,'mean_CB103_BJV.png');
    
    f3=figure(3);
    semilogy(t(1:32),abs(Brms(2:33)-Brms(1:32)),'k', t(1:32),abs(Jrms(2:33)-Jrms(1:32)), 'b',t(1:32),abs(virms(2:33)-virms(1:32)), 'r','LineWidth',1.5) 
    %plot(t(1:32),abs(Brms(2:33)-Brms(1:32)),'k', t(1:32),abs(Jrms(2:33)-Jrms(1:32)), 'b',t(1:32),abs(virms(2:33)-virms(1:32)), 'r','LineWidth',1.5) 
    set(gca,'FontSize',18)
    xlim([0 480]); %ylim([10^(-5) 10^(-1)]);
    y1=get(gca,'ylim');
    hold on
    x3=120;
    plot([x3 x3],y1,'--', 'LineWidth',1.5)
    % Text in the plot
    xlabel('$\omega_{pi}t$','Interpreter','latex','FontSize',28)
    %ylabel('$|\Delta\psi^{rms}|$','Interpreter','latex','FontSize',22)  
    ylabel('$|\Delta \langle \psi \rangle|$','Interpreter','latex','FontSize',28)  
    txt8 = '$\omega_{pi}t=120$';
    %text(100 , 0.13,txt8,'Interpreter','latex')
    %title('$P1D_\|$','Interpreter','latex')
    %legend({'$|\Delta B^{rms}|$','$|\Delta J^{rms}|$','$|\Delta v_{i}^{rms}|$'},'Interpreter','latex','FontSize',22)
    legend({'$|\Delta \langle B \rangle|$','$|\Delta \langle J \rangle|$','$|\Delta \langle v_{i} \rangle|$'},'Interpreter','latex','FontSize',25)
    hold off
    
%-------------------------------------------------------------------------
%  
%-------------------------------------------------------------------------
    
    
    
 
 Bmag = table2array(stats(1:33,2));
 Bmag = Bmag(~isnan(Bmag));
 Emag = table2array(stats(1:33,3));
 Emag = Emag(~isnan(Emag));
 vimag = table2array(stats(1:33,5));
 vimag = vimag(~isnan(vimag));
 vemag = table2array(stats(1:33,6));
 vemag = vemag(~isnan(vemag));
 vithmag = table2array(stats(1:33,7));
 vithmag = vithmag(~isnan(vithmag));
 vethmag = table2array(stats(1:33,8));
 vethmag = vethmag(~isnan(vethmag));
 
 mi2me=100;
 
 UE=0.5*Emag.^2;
 UM=0.5*Bmag.^2;
 
 Ki = 0.5*vimag.^2;
 Ke = (1\mi2me)*0.5*vemag.^2;
 
 Kith = 0.5*vithmag.^2;
 Keth = (1\mi2me)*0.5*vethmag.^2;
 
 Utotal = UE + UM +Ki + Ke + Kith+ Keth;
 
 
 Bmag2 = Bx.^2 +By.^2 +Bz.^2;
    Bmag2_2 = Bmag2(:);
    Brms_1 = sqrt(mean(Bmag2_2));
 
 f2=figure(2);
    %semilogy(t,Utotal,'k', t, UE, 'b',t, UM, 'r','LineWidth',1.5) 
    semilogy(t,Utotal/rms(Utotal),'k', t, UE/rms(UE), 'b',t, UM/rms(UM), 'r','LineWidth',1.5) 
    set(gca,'FontSize',18)
    %xlim([0.03 52.36]); ylim([10^(-5) 100]);
    y1=get(gca,'ylim');
    hold on
    x3=96;
    plot([x3 x3],y1,'--', 'LineWidth',1.5)
    % Text in the plot
    xlabel('$\omega_{p,i}t$','Interpreter','latex')
    ylabel('Energy density','Interpreter','latex')  
    txt8 = '$\omega_{p,i}t=96$';
    %text(100 , 0.13,txt8,'Interpreter','latex')
    %title('$P1D_\|$','Interpreter','latex')
    legend({'$\langle U_{T} \rangle$','$$\langle U_{E} \rangle$$','$$\langle U_{B} \rangle$$'},'Interpreter','latex')
    hold off
    
    
    
 cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images';
    % Save the plots
    %-------------------------------------------------------------------------
    saveas(f1,'rms_this_m_bF.png');
    saveas(f3,'deltarms_this_m_bF.png');
    saveas(f2,'energy_density.png');