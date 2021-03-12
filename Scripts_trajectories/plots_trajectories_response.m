%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First load the files per trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently importing by hand


cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Trajectories'

load('8to11trajectories.mat');

d=exyz10(:,4);ex=exyz10(:,1);ey=exyz10(:,2);ez=exyz10(:,3);
hx=hxyz10(:,1);hy=hxyz10(:,2);hz=hxyz10(:,3);
jx=jxyz10(:,1);jy=jxyz10(:,2);jz=jxyz10(:,3);
vix=vixyz10(:,1);viy=vixyz10(:,2);viz=vixyz10(:,3);
vex=vexyz10(:,1);vey=vexyz10(:,2);vez=vexyz10(:,3);
Tix=Tixyz10(:,1);Tiy=Tixyz10(:,2);Tiz=Tixyz10(:,3);
Tex=Texyz10(:,1);Tey=Texyz10(:,2);Tez=Texyz10(:,3);
ne=nine10(:,1);ni=nine10(:,2);


d=exyz8(:,4);ex=exyz8(:,1);ey=exyz8(:,2);ez=exyz8(:,3);
hx=hxyz8(:,1);hy=hxyz8(:,2);hz=hxyz8(:,3);
jx=jxyz8(:,1);jy=jxyz8(:,2);jz=jxyz8(:,3);
vix=vixyz8(:,1);viy=vixyz8(:,2);viz=vixyz8(:,3);
vex=vexyz8(:,1);vey=vexyz8(:,2);vez=vexyz8(:,3);
Tix=Tixyz8(:,1);Tiy=Tixyz8(:,2);Tiz=Tixyz8(:,3);
Tex=Texyz8(:,1);Tey=Texyz8(:,2);Tez=Texyz8(:,3);
ne=nine8(:,1);ni=nine8(:,2);

d2=exyz11(:,4);ex2=exyz11(:,1);ey2=exyz11(:,2);ez2=exyz11(:,3);
hx2=hxyz11(:,1);hy2=hxyz11(:,2);hz2=hxyz11(:,3);
jx2=jxyz11(:,1);jy2=jxyz11(:,2);jz2=jxyz11(:,3);
vix2=vixyz11(:,1);viy2=vixyz11(:,2);viz2=vixyz11(:,3);
vex2=vexyz11(:,1);vey2=vexyz11(:,2);vez2=vexyz11(:,3);
Tix2=Tixyz11(:,1);Tiy2=Tixyz11(:,2);Tiz2=Tixyz11(:,3);
Tex2=Texyz11(:,1);Tey2=Texyz11(:,2);Tez2=Texyz11(:,3);
ne2=nine11(:,1);ni2=nine11(:,2);

d = d(~isnan(d)); 
ex = ex(~isnan(ex)); ey = ey(~isnan(ey)); ez = ez(~isnan(ez)); 
hx = hx(~isnan(hx)); hy = hy(~isnan(hy)); hz = hz(~isnan(hz));
jx = jx(~isnan(jx)); jy = jy(~isnan(jy)); jz = jz(~isnan(jz));
vix = vix(~isnan(vix)); viy = viy(~isnan(viy)); viz = viz(~isnan(viz));
vex = vex(~isnan(vex)); vey = vey(~isnan(vey)); vez = vez(~isnan(vez));
Tix = Tix(~isnan(Tix)); Tiy = Tiy(~isnan(Tiy)); Tiz = Tiz(~isnan(Tiz));
Tex = Tex(~isnan(Tex)); Tey = Tey(~isnan(Tey)); Tez = Tez(~isnan(Tez));
ne = ne(~isnan(ne)); ni = ni(~isnan(ni)); 
em=sqrt(ex.^2 + ey.^2 + ez.^2); 
hm=sqrt(hx.^2 + hy.^2 + hz.^2); 
jm=sqrt(jx.^2 + jy.^2 + jz.^2); 
vim=sqrt(vix.^2 + viy.^2 + viz.^2); 
vem=sqrt(vex.^2 + vey.^2 + vez.^2); 
Tim=(Tix + Tiy + Tiz)/3;
Tem=(Tex + Tey + Tez)/3;
Vailocal=hm./sqrt(ni);
csmi=sqrt(Tim + Vailocal.^2);
csme=sqrt(100*Tem + Vailocal.^2);

d2 = d2(~isnan(d2)); 
ex2 = ex2(~isnan(ex2)); ey2 = ey2(~isnan(ey2)); ez2 = ez2(~isnan(ez2)); 
hx2 = hx2(~isnan(hx2)); hy2 = hy2(~isnan(hy2)); hz2 = hz2(~isnan(hz2));
jx2 = jx2(~isnan(jx2)); jy2 = jy2(~isnan(jy2)); jz2 = jz2(~isnan(jz2));
vix2 = vix2(~isnan(vix2)); viy2 = viy2(~isnan(viy2)); viz2 = viz2(~isnan(viz2));
vex2 = vex2(~isnan(vex2)); vey2 = vey2(~isnan(vey2)); vez2 = vez2(~isnan(vez2));
Tix2 = Tix2(~isnan(Tix2)); Tiy2 = Tiy2(~isnan(Tiy2)); Tiz2 = Tiz2(~isnan(Tiz2));
Tex2 = Tex2(~isnan(Tex2)); Tey2 = Tey2(~isnan(Tey2)); Tez2 = Tez2(~isnan(Tez2));
ne2 = ne2(~isnan(ne2)); ni2 = ni2(~isnan(ni2)); 
em2=sqrt(ex2.^2 + ey2.^2 + ez2.^2); 
hm2=sqrt(hx2.^2 + hy2.^2 + hz2.^2); 
jm2=sqrt(jx2.^2 + jy2.^2 + jz2.^2); 
vim2=sqrt(vix2.^2 + viy2.^2 + viz2.^2); 
vem2=sqrt(vex2.^2 + vey2.^2 + vez2.^2); 
Tim2=(Tix2 + Tiy2 + Tiz2)/3;
Tem2=(Tex2 + Tey2 + Tez2)/3;
Vailocal2=hm2/sqrt(ni2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
%B are normalized to the background magnetic field 
B0=0.1;
%n is normalise to the proton density n0=1
n0=1;
%Vi normalised to the proton alfven speed V_a=0.1c The alfven speed of
%electrons is Vae=10Vai
Vai=0.1; Vae=1;
%E = vB, as B0=0.1 and v in units of c, E0= B0*Vai/c
E0=0.01;
%T=(beta/2)B^2, Thus, as beta=1 
T0=0.005;
%The current density J0=1?
J0=Vai;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 % This makes the plot of the magnetic field components
 start_p=18.4; delta=8.5; %square2 trajectory 10
 
 %start_p=12; delta=3; %square1_vertical trajectory 8
 %start_p=20; delta=3;%square1
 %start_p=10; delta=7; %square3 trajectory 11
 
 vim=vim2; ni=ni2; vem=vem2;  
 hm=hm2; Tim=Tim2; Tem=Tem2;
 
 % This is the original one
 %-------------------------------------------------------------------------
  % this one is the plot of Ti,Te,B,ni,ve,vi
 f7=figure(7);
% f7.PaperPositionMode = 'auto'; fig_pos = f7.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
 xpos=31.5;%25.5;
 xlpos=30;
 %-------------------------------------------------------------------------
 subplot(5,1,1),
 plot(d,vim/Vai,'r','LineWidth',2),set(gca,'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlabel('$r/d_{i}$','Position',[median(d)+0.2 -0.1],'Interpreter','latex','FontSize',18)
 xlim([0 xlpos]); 
 maxi=abs(max(vim/Vai));
 mini=abs(min(vim/Vai));
 ylimi=max(maxi,mini);
 ylim([0 ylimi]);
 yticks([0.2 0.4 0.6]);
 %yticks([1 2]);
 ylabel('$$\frac{|v_{i}|}{V_{A0}}$$','Rotation',0,'Position',[xpos 0],'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p 0 delta 1], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 subplot(5,1,2),
 plot(d,vem/Vai,'b','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 maxe=abs(max(vem/Vai));
 mine=abs(min(vem/Vai));
 ylime=max(maxe,mine);
 ylim([0 ylime]);
 %set(gca,'FontSize',15)
 yticks([2 4 6 8]);
 xlim([0 xlpos]); %yticks([1 2 5]);%6 8 10]);
 ylabel('$$\frac{|v_{e}|}{V_{A0}}$$','Rotation',0,'Position',[xpos 0],'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p 0 delta 7], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 subplot(5,1,3);
 plot(d,ni/n0,'c','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 xlpos]); 
 ylimn=max(ni/n0);
 ylim([0 ylimn]);  
 yticks([0.5 1 1.5]);
 ylabel('$$\frac{n_{i}}{n_{0}}$$','Rotation',0,'Position',[xpos 0],'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 subplot(5,1,4);
 plot(d,hm/B0,'k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 xlpos]); 
 ylimb=max(hm/B0);
 ylim([0 ylimb]); 
 %yticks([0.4 0.8 1.2]);
 yticks([1 2]);
 ylabel('$$\frac{|B|}{B_{0}}$$','Rotation',0,'Position',[xpos 0],'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p 0 delta 3], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 subplot(5,1,5);
 plot(d,Tim/T0,'r', d,Tem/T0,'b','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 maxTi=abs(max(Tim/T0));
 maxTe=abs(max(Tem/T0));
 ylimT=max(maxTi,maxTe);
 ylim([0 ylimT]);
 xlim([0 xlpos]); yticks([1 2 3]);
 ylabel('$$\frac{T_{s}}{T_{0}}$$','Rotation',0,'Position',[xpos 0],'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.1 .81 .8 .18]); set(ha(2),'position',[.1 .63 .8 .18]); 
 set(ha(3),'position',[.1 .45 .8 .18]); set(ha(4),'position',[.1 .27 .8 .18]); 
 set(ha(5),'position',[.1 .09 .8 .18]);
 lgd=legend({'$i$','$e$'},'Location','northwest','Interpreter','latex');
 lgd.FontSize = 14;
 %-------------------------------------------------------------------------
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % This is the one in order x,y,z,m trajectory 10
 %-------------------------------------------------------------------------
 f8=figure(8);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 yyaxis left
 line1 = plot(d,vim/Vai,'r','LineWidth',2); set(gca,'FontSize',16)
 xlim([0 xlpos]);
 maxi=abs(max(vim/Vai));
 mini=abs(min(vim/Vai));
 ylim1=max(maxi,mini);
 ylim([0 ylim1]);
 yticks([0.2 0.4 0.6]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{|v_{i}|}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 yyaxis right
 line2 = plot(d, hm/B0, 'k','LineWidth',2);set(gca,'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlim([0 xlpos]);
 ylim([0 max(hm/B0)]) 
 %%yticks([1 2]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{|B|}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 xlabel('$r/d_{i}$','Position',[median(d)+0.2 -0.1],'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 yyaxis left
 plot(d,viz/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--'); 
 xlim([0 xlpos]);
 ylim1=max(abs(max(viz/Vai)), abs(min(viz/Vai)));
 ylim([-ylim1 ylim1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iz}}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 yyaxis right
 plot(d,hz/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 xlpos]);  
 ylimh=max(abs(max(hz/B0)), abs(min(hz/B0)));
 ylim([-ylimh ylimh]);
 %%yticks([0 1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{z}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 yyaxis left
 plot(d,viy/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--');
 xlim([0 xlpos]); 
 ylimi=max(abs(max(viy/Vai)), abs(min(viy/Vai)));
 ylim([-ylimi ylimi]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iy}}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 yyaxis right
 plot(d,hy/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 xlpos]); 
 ylimh=max(abs(max(hy/B0)), abs(min(hy/B0)));
 ylim([-ylimh ylimh]);
 %%yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{y}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 yyaxis left
 plot(d,vix/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 yline(0,'--'); 
 %xline(12,'--'); 
 xlim([0 xlpos]); 
 ylimi=max(abs(max(vix/Vai)), abs(min(vix/Vai)));
 ylim([-ylimi ylimi]);
 %%yticks([-0.1 0.1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ix}}{V_{A0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 yyaxis right
 plot(d,hx/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',16)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 xlim([0 xlpos]); 
 %%yticks([-0.5 0.5]); 
 ylimh=max(abs(max(hx/B0)), abs(min(hx/B0)));
 ylim([-ylimh ylimh]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{x}}{B_{0}}$$','Rotation',0,'Interpreter','latex','FontSize',18)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .74 .8 .21]); set(ha(2),'position',[.12 .53 .8 .21]);
 set(ha(3),'position',[.12 .32 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 %legend({'$v_{i}$', '$B$' },'Location','northwest','Interpreter','latex')
 hL = legend([line1,line2],{'$v_{i}$', '$B$'},'Interpreter','latex');
 set(hL);
 
 %-------------------------------------------------------------------------
 
 
 % This is the one with x,y,z, T2 T3
 %-------------------------------------------------------------------------
 f83=figure(83);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 str = '#EDB120';
 color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
 yyaxis left
 %plot(d,vim/Vai,'-r',d,vim2/Vai,'-m' ,'LineWidth',2),set(gca,'FontSize',20)
 line1 = plot(d,vim/Vai,'-r','LineWidth',2);set(gca,'FontSize',20);
 hold on 
 line2 = plot(d,vim2/Vai,'-','Color', color ,'LineWidth',2);%,set(gca,'FontSize',20)
 hold off
 xlim([0 24]); 
 ymax=max(max(vim/Vai), max(vim2/Vai));
 ylim([0 ymax]);
 %yticks([0.3 2 3]);
 yticks([0.2 0.4 0.6]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{v_{i}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 yyaxis right
 line3 = plot(d, hm/B0, '-k','LineWidth',2); set(gca,'FontSize',20)
 hold on
 line4 = plot(d, hm2/B0, '-b','LineWidth',2); set(gca,'FontSize',20)
 hold off
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlim([0 24]); 
 ylim([0 max(hm/B0)]) 
 yticks([1 2]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{B}{B_{0}}\right|$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 yyaxis left
 plot(d,viz/Vai,'-r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 hold on
 plot(d,viz2/Vai,'-','Color', color,'LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 hold off
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); 
 xlim([0 24]);  %ylim([-0.5 2]); yticks([0.5 1.5]);
 ymax1=abs(max(max(viz/Vai), max(viz2/Vai)));
 ymin1=abs(min(min(viz/Vai), min(viz2/Vai)));
 ymax=max(ymax1, ymin1);
 ylim([-ymax ymax]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iz}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hz/B0, '-k',d, hz2/B0, '-b','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 yline(0,'--');
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); 
 xlim([0 24]);  
 ylimh=max(abs(max(hz/B0)), abs(min(hz/B0)));
 ylim([-ylimh ylimh]);
 yticks([-1 0 1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{z}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 yyaxis left
 plot(d,viy/Vai,'-r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 hold on
 plot(d,viy2/Vai,'-','Color', color,'LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 hold off
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--'); 
 xlim([0 24]); %ylim([-0.7 0.7]); yticks([-0.5 0.5]);%6 8 10]);
 yymax1=abs(max(max(viy/Vai), max(viy2/Vai)));
 ymin1=abs(min(min(viy/Vai), min(viy2/Vai)));
 ymax=max(ymax1, ymin1);
 ylim([-ymax ymax]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iy}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hy/B0, '-k',d, hy2/B0, '-b','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 24]); 
 ylimh=max(abs(max(hy/B0)), abs(min(hy/B0)));
 ylim([-ylimh ylimh]);
 yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{y}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 yyaxis left
 plot(d,vix/Vai,'-r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 hold on
 plot(d,vix2/Vai,'-','Color', color,'LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 hold off
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 yline(0,'--'); 
 %xline(12,'--'); 
 xlim([0 24]); 
 %ylim([-0.7 0.7]); 
 ymax1=abs(max(max(vix/Vai), max(vix2/Vai)));
 ymin1=abs(min(min(vix/Vai), min(vix2/Vai)));
 ymax=max(ymax1, ymin1);
 ylim([-ymax ymax]);
 yticks([-0.2 0.2]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ix}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hx/B0, '-k',d, hx2/B0, '-b','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 xlim([0 24]); 
 yticks([-0.5 0.5]); 
 ylimh=max(abs(max(hx/B0)), abs(min(hx/B0)));
 ylim([-ylimh ylimh]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{x}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21]);
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 %legend({'$v_{iT2}$','$v_{iT3}$','$B_{T2}$','$B_{T3}$'},'Location','northwest','Interpreter','latex')
 hL = legend([line1,line2, line3, line4],{'$v_{iT2}$','$v_{iT3}$','$B_{T2}$','$B_{T3}$'},'Interpreter','latex');
 set(hL);
 %-------------------------------------------------------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here/'
 % Correlation plots
 %
 %-------------------------------------------------------------------------
 f331=figure(331); 
 xlpos=30;
 subplot(4,1,1),
 Xx=vix/Vai; Yx=hx/B0; t2=d; Dd=10; Wd=20; dd=20; %this is to avoid the first data where there is an error
 colorbar1 = [0 0 0] ; coloredge1 = [0 0 0]; 
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);% plot(t,dXdY,'k','LineWidth',2);
 set(gca,'xtick',[],'FontSize',16)
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{x},B_{x})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBx}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,2);
  Xx=viy/Vai; Yx=hy/B0; t2=d;
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);% plot(t,dXdY,'k','LineWidth',2); 
 set(gca,'xtick',[],'FontSize',16)
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{y},B_{y})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBy}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 Xx=viz/Vai; Yx=hz/B0; t2=d;
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);% plot(t,dXdY,'k','LineWidth',2); 
 set(gca,'xtick',[],'FontSize',16)
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{z},B_{z})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBz}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 Xx=ni%vim/Vai; 
 Yx=hm/B0; t2=d;
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);% plot(t,dXdY,'k','LineWidth',2); 
 set(gca,'FontSize',16)
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(|v|,|B|)$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{|v||B|}$$','Interpreter','latex','FontSize',18) 
 xlabel('$r/d_{i}$','Position',[median(t) -(1.1*ylim1)],'Interpreter','latex')
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(4),'position',[.12 .73 .8 .21]); set(ha(3),'position',[.12 .52 .8 .21]);
 set(ha(2),'position',[.12 .31 .8 .21]); set(ha(1),'position',[.12 .11 .8 .21]);
 %legend({'$T_{1}$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 
 %-------------------------------------------------------------------------
 f332=figure(332); 
 subplot(4,1,1),
 Xx=vix/Vai; Yx=hx/B0; t2=d; Dd=20; Wd=10;
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); dd=1; %this is to avoid the first data where there is an error 
 %plot(t,dXdY,'-sk','LineWidth',2); 
 bar(t,dXdY,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
 set(gca,'xtick',[],'FontSize',16);
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{x},B_{x})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBx}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,2);
  Xx=viy/Vai; Yx=hy/B0; t2=d;
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 %plot(t,dXdY,'-sk','LineWidth',2); 
 bar(t,dXdY,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
 set(gca,'xtick',[],'FontSize',16)
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{y},B_{y})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBy}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 Xx=viz/Vai; Yx=hz/B0; t2=d;
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 %plot(t,dXdY,'-sk','LineWidth',2); 
 bar(t,dXdY,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
 set(gca,'xtick',[],'FontSize',16)
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{z},B_{z})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBz}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 Xx=vim/Vai; Yx=hm/B0; t2=d;
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 bar(t,dXdY,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
 %plot(t,dXdY,'-sk','LineWidth',2); 
 set(gca,'FontSize',16)
 yline(0,'--')
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 ylim([-ylim1 ylim1])
 xlim([0 xlpos])
 %yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(|v|,|B|)$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{|v||B|}$$','Interpreter','latex','FontSize',18) 
 xlabel('$r/d_{i}$','Position',[median(t) -(1.1*ylim1)],'Interpreter','latex')
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(4),'position',[.12 .73 .8 .21]); set(ha(3),'position',[.12 .52 .8 .21]);
 set(ha(2),'position',[.12 .31 .8 .21]); set(ha(1),'position',[.12 .11 .8 .21]);
 %legend({'$T_{1}$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 %}
 
 
 % T2T3
 %
 %------------------------------------------------------------------------- 
 f341=figure(341); 
 subplot(4,1,1),
 Xx=vix/Vai; Yx=hx/B0; t2=d; Dd=20; Wd=20; dd=20; %this is to avoid the first data where there is an error
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,20);
 colorbar1 = [0 0 0] ; coloredge1 = [0 0 0];
 colorbar2 = [0 .5 .5 ] ; coloredge2 = [0 .9 .9];
 dlen=length(t);  
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5; 
 set(gca,'xtick',[],'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=vix2/Vai; Yx=hx2/B0; t2=d;  
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line2 = bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line2.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{x},B_{x})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBx}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %hold off
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 Xx=viy/Vai; Yx=hy/B0; t2=d;  
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); 
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=viy2/Vai; Yx=hy2/B0; t2=d;  
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line2 =  bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line2.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{y},B_{y})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBy}$$','Interpreter','latex','FontSize',18) 
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 Xx=viz/Vai; Yx=hz/B0; t2=d; 
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=viz2/Vai; Yx=hz2/B0; t2=d;  
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); 
 line2 =  bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line2.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{z},B_{z})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBz}$$','Interpreter','latex','FontSize',18) 
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 Xx=vim/Vai; Yx=hm/B0; t2=d; 
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); 
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5;
 set(gca,'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=vim2/Vai; Yx=hm2/B0; t2=d;  
 [t,dXdY]=Local_corr(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line2 =  bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line2.FaceAlpha = 0.5;
 set(gca,'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(|v|,|B|)$$','Interpreter','latex','FontSize',18)
 ylabel('$$\rho_{|v||B|}$$','Interpreter','latex','FontSize',18) 
 xlabel('$r/d_{i}$','Position',[median(t) -(ylim3)],'Interpreter','latex','FontSize',18)
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(4),'position',[.12 .73 .8 .21]); set(ha(3),'position',[.12 .52 .8 .21]);
 set(ha(2),'position',[.12 .31 .8 .21]); set(ha(1),'position',[.12 .11 .8 .21]);
 hL = legend([line1,line2],{'$T_{2}$', '$T_{3}$'},'Interpreter','latex');
 set(hL);
 %legend({'$v_{iT2}$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 
 %Calculating covariance by intervals 
 %------------------------------------------------------------------------- 
 f342=figure(342); 
 subplot(4,1,1),
 colorbar1 = [0 0 0] ; coloredge1 = [0 0 0];
 colorbar2 = [0 .5 .5 ] ; coloredge2 = [0 .9 .9];
 Xx=vix/Vai; Yx=hx/B0; t2=d; Dd=20; Wd=20; dd=20; %this is to avoid the first data where there is an error 
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); 
 line1 =  bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=vix2/Vai; Yx=hx2/B0; t2=d;  
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); 
 line2 = bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'r','LineWidth',2); 
 set(gca,'xtick',[],'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{x},B_{x})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBx}$$','Interpreter','latex','FontSize',18) 
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %hold off
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 Xx=viy/Vai; Yx=hy/B0; t2=d; 
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line1=bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=viy2/Vai; Yx=hy2/B0; t2=d;  
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'r','LineWidth',2); 
 set(gca,'xtick',[],'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{y},B_{y})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBy}$$','Interpreter','latex','FontSize',18) 
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 Xx=viz/Vai; Yx=hz/B0; t2=d; 
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line1=bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);%plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5;
 set(gca,'xtick',[],'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=viz2/Vai; Yx=hz2/B0; t2=d;  
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); 
 bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'r','LineWidth',2); 
 set(gca,'xtick',[],'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(v_{z},B_{z})$$','Interpreter','latex','FontSize',20)
 ylabel('$$\rho_{vBz}$$','Interpreter','latex','FontSize',18) 
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 Xx=vim/Vai; Yx=hm/B0; t2=d;  
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t);  
 line1=bar(t,dXdY,'FaceColor',colorbar1,'EdgeColor',coloredge1,'LineWidth',2);% plot(t,dXdY,'k','LineWidth',2); 
 line1.FaceAlpha = 0.5;
 set(gca,'FontSize',16)
 ylim1 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 hold on
 Xx=vim2/Vai; Yx=hm2/B0; t2=d;  
 [t,dXdY]=Local_corr2(t2,Xx,Yx,Dd,Wd);
 dlen=length(t); 
 bar(t,dXdY,'FaceColor',colorbar2,'EdgeColor',coloredge2,'LineWidth',2);%plot(t,dXdY,'r','LineWidth',2); 
 set(gca,'FontSize',16)
 ylim2 =max(abs(min(dXdY(dd:dlen))),max(dXdY(dd:dlen)));
 yline(0,'--')
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 %yticks([-0.5 0.5])
 xlim([0 xlpos])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.4])
 %ylabel('$$\rho(|v|,|B|)$$','Interpreter','latex','FontSize',18)
 ylabel('$$\rho_{|v||B|}$$','Interpreter','latex','FontSize',18) 
 xlabel('$r/d_{i}$','Position',[median(t) -(ylim3+0.05)],'Interpreter','latex','FontSize',18)
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(4),'position',[.12 .73 .8 .21]); set(ha(3),'position',[.12 .52 .8 .21]);
 set(ha(2),'position',[.12 .31 .8 .21]); set(ha(1),'position',[.12 .11 .8 .21]);
 hL = legend([line1,line2],{'$T_{2}$', '$T_{3}$'},'Interpreter','latex');
 set(hL);
 %legend({'$v_{iT2}$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %}
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Save the plots
  %-------------------------------------------------------------------------
  cd '/Volumes/PSC_DiRAC_DATA/';

  saveas(f7,'trajectory_BnvieTie_T1_square.png');
  saveas(f8,'trajectory_viB_T1_xyz_square.png');
  saveas(f331,'trajectory_dev_rhovb_T1_xyz_square.png');
  saveas(f332,'trajectory_cov_rhovb_T1_xyz_square.png');
  
  saveas(f341,'trajectory_dev_rhovb_T2T3_xyz_square.png');
  saveas(f342,'trajectory_cov_rhovb_T2T3_xyz_square.png');
  
  saveas(f83,'trajectory_viB_T2T3_xyz.png');
  saveas(f777,'trajectory_10vive_csm.png');
  saveas(f33,'trajectory_rhovb_10_xyz_square.png');
  saveas(f34,'trajectory_rhovb_T1T2_xyz_square.png');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %-------------------------------------------------------------------------
 