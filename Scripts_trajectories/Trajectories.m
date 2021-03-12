%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First load the files per trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently importing by hand



cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/Trajectories'

load('8to11trajectories.mat');

d=exyz8(:,4);ex=exyz8(:,1);ey=exyz8(:,2);ez=exyz8(:,3);
hx=hxyz8(:,1);hy=hxyz8(:,2);hz=hxyz8(:,3);
jx=jxyz8(:,1);jy=jxyz8(:,2);jz=jxyz8(:,3);
vix=vixyz8(:,1);viy=vixyz8(:,2);viz=vixyz8(:,3);
vex=vexyz8(:,1);vey=vexyz8(:,2);vez=vexyz8(:,3);
Tix=Tixyz8(:,1);Tiy=Tixyz8(:,2);Tiz=Tixyz8(:,3);
Tex=Texyz8(:,1);Tey=Texyz8(:,2);Tez=Texyz8(:,3);
ne=nine8(:,1);ni=nine8(:,2);

d=exyz9(:,4);ex=exyz9(:,1);ey=exyz9(:,2);ez=exyz9(:,3);
hx=hxyz9(:,1);hy=hxyz9(:,2);hz=hxyz9(:,3);
jx=jxyz9(:,1);jy=jxyz9(:,2);jz=jxyz9(:,3);
vix=vixyz9(:,1);viy=vixyz9(:,2);viz=vixyz9(:,3);
vex=vexyz9(:,1);vey=vexyz9(:,2);vez=vexyz9(:,3);
Tix=Tixyz9(:,1);Tiy=Tixyz9(:,2);Tiz=Tixyz9(:,3);
Tex=Texyz9(:,1);Tey=Texyz9(:,2);Tez=Texyz9(:,3);
ne=nine9(:,1);ni=nine9(:,2);

d=exyz10(:,4);ex=exyz10(:,1);ey=exyz10(:,2);ez=exyz10(:,3);
hx=hxyz10(:,1);hy=hxyz10(:,2);hz=hxyz10(:,3);
jx=jxyz10(:,1);jy=jxyz10(:,2);jz=jxyz10(:,3);
vix=vixyz10(:,1);viy=vixyz10(:,2);viz=vixyz10(:,3);
vex=vexyz10(:,1);vey=vexyz10(:,2);vez=vexyz10(:,3);
Tix=Tixyz10(:,1);Tiy=Tixyz10(:,2);Tiz=Tixyz10(:,3);
Tex=Texyz10(:,1);Tey=Texyz10(:,2);Tez=Texyz10(:,3);
ne=nine10(:,1);ni=nine10(:,2);

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
 
 emjm=ex.*jx + ey.*jy + ez.*jz; 
 
 f3=figure(3);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 plot(d,hx/B0,'b'),set(gca,'FontSize',15)
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 24]); ylim([-1.3 1.3]); yticks([-0.5 0.5]); 
 %xlim([0 30]); ylim([-3.5 3.5]); yticks([-1.5 0 1.5]);
 yline(0,'--'); 
 %xline(jmax_pos,'--'); xline(hmin_pos,'--');
 ylabel('$$\tilde{B}_{x}$$','Interpreter','latex')
 rectangle('Position', [start_p -3.5 delta 7], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2),
 plot(d,hy/B0,'r'),set(gca,'xtick',[],'FontSize',15)
 set(gca,'FontSize',15)
 %y1=get(gca,'ylim');
 xlim([0 24]); ylim([-1 1]); yticks([-0.5 0.5]);
 %xlim([0 30]); ylim([-3.5 3.5]); yticks([-1.5 0 1.5]);
 yline(0,'--'); 
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 ylabel('$$\tilde{B}_{y}$$','Interpreter','latex')
 rectangle('Position', [start_p -3.5 delta 7], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 plot(d,hz/B0,'c'),set(gca,'xtick',[],'FontSize',15)
 xlim([0 24]); ylim([-1.5 1.5]); yticks([-1.0 1.0]);
 %xlim([0 30]); ylim([-3.5 3.5]); yticks([-1.5 0 1.5]);
 yline(0,'--'); 
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 ylabel('$$\tilde{B}_{z}$$','Interpreter','latex')
 rectangle('Position', [start_p -3.5 delta 7], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 plot(d,emjm/(E0*0.1),'k'),set(gca,'xtick',[],'FontSize',15)
 maxi=abs(max(emjm/(E0*0.1)));
 mini=abs(min(emjm/(E0*0.1)));
 ylim1=max(maxi,mini);
 ylim([-ylim1 ylim1]);
 xlim([0 24]); %ylim([0 0.3]); yticks([1.0]);
 xline(9,'--'); 
 %xlim([0 30]); ylim([0 3]); yticks([1 2]); 
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 ylabel('$$|\tilde{B}|$$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 3], 'FaceColor', [0.7 0.7 0.7 0.2])
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21])
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21])
 
 
 
 % this one is the plot of Ti,Te,B,ni,ne,ve,vi
 f777=figure(777);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 plot(d,vim./csmi,'k', d,vim./Vailocal,'r',d,Vailocal./csmi,'g'),set(gca,'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlabel('$r/d_{i}$','Interpreter','latex')
 xlim([0 30]); 
 ylim([0 1.5]); yticks([0.5]);
 ylabel('$$\left|\frac{v_{i}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 1.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2),
 plot(d,vem./csme,'k',d,vem./Vailocal,'b',d,Vailocal./csme,'g'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 set(gca,'FontSize',15)
 xlim([0 30]); ylim([0 4]); 
 yticks([1 2 3]);%6 8 10]);
 ylabel('$$\left|\frac{v_{e}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 6], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 plot(d,ni/n0,'c',d,ne/n0,'m'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 30]); ylim([0 1.8]); yticks([0.5 1]);
 ylabel('$$\frac{n_{i}}{n_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 plot(d,hm/B0,'k'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 30]);  ylim([0 2]); yticks([1 2]);
 ylabel('$$\left|\frac{B}{B_{0}}\right|$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 3], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21])
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21])
 %legend({'$i$','$e$'},'Location','northwest','Interpreter','latex')

 
 
 % This is the original one
 %-------------------------------------------------------------------------
  % this one is the plot of Ti,Te,B,ni,ve,vi
 f7=figure(7);
 %-------------------------------------------------------------------------
 subplot(5,1,1),
 plot(d,vim/Vai,'r','LineWidth',2),set(gca,'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlabel('$r/d_{i}$','Interpreter','latex')
 xlim([0 24]); 
 maxi=abs(max(vim/Vai));
 mini=abs(min(vim/Vai));
 ylimi=max(maxi,mini);
 ylim([0 ylimi]);
 yticks([0.2 0.4 0.6]);
 %yticks([1 2]);
 ylabel('$$\left|\frac{v_{i}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 1], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,2),
 plot(d,vem/Vai,'b','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 maxe=abs(max(vem/Vai));
 mine=abs(min(vem/Vai));
 ylime=max(maxe,mine);
 ylim([0 ylime]);
 %set(gca,'FontSize',15)
 yticks([2 4 6 8]);
 xlim([0 24]); %yticks([1 2 5]);%6 8 10]);
 ylabel('$$\left|\frac{v_{e}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 7], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,3);
 plot(d,ni/n0,'c','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 24]); 
  ylim([0 max(ni/n0)]);  
 yticks([0.5 1 1.5]);
 ylabel('$$\frac{n_{i}}{n_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,4);
 plot(d,hm/B0,'k','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 24]);  
 ylim([0 max(hm/B0)]); 
 %yticks([0.4 0.8 1.2]);
 yticks([1 2]);
 ylabel('$$\left|\frac{B}{B_{0}}\right|$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 3], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,5);
 plot(d,Tim/T0,'r', d,Tem/T0,'b','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 maxTi=abs(max(Tim/T0));
 minTi=abs(max(Tem/T0));
 ylimT=max(maxi,mini);
 ylim([0 4]);
 xlim([0 24]); yticks([1 2 3]);
 ylabel('$$\frac{T_{i}}{T_{0}}, \frac{T_{e}}{T_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.1 .81 .8 .18]); set(ha(2),'position',[.1 .63 .8 .18]); 
 set(ha(3),'position',[.1 .45 .8 .18]); set(ha(4),'position',[.1 .27 .8 .18]); 
 set(ha(5),'position',[.1 .09 .8 .18]);
 legend({'$i$','$e$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 
 %vim=vim2; vix=vix2; viy=viy2; viz=viz2; 
 %hm=hm2; hx=hx2; hy=hy2; hz=hz2;
 
 
 % subset of B and vi components
 %-------------------------------------------------------------------------
 f8=figure(8);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 yyaxis left
 line1 = plot(d,vim/Vai,'r','LineWidth',2); set(gca,'FontSize',20);
 %xlim([0 24]);
 xlim([0 30]);
 maxi=abs(max(vim/Vai)); mini=abs(min(vim/Vai)); %ylim1=max(maxi,mini);
 ylim([0 0.7]);
 yticks([0.2 0.4 0.6]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{v_{i}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 yyaxis right
 line2 = plot(d, hm/B0, 'k','LineWidth',2); set(gca,'FontSize',20);
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %xlim([0 24]);
 xlim([0 30]);
 ylim([0 1.6]);% ylim([0 max(hm/B0)]) 
 yticks([0.4 0.8 1.2]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{B}{B_{0}}\right|$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 yyaxis left
 plot(d,vix/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 yline(0,'--'); 
 %xline(12,'--'); 
 %xlim([0 24]); 
 xlim([0 30]);
 ylimi=max(abs(max(vix/Vai)), abs(min(vix/Vai)));
 ylim([-ylimi ylimi]);
 %%yticks([-0.1 0.1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ix}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hx/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 %xlim([0 24]); 
 xlim([0 30]);
 %%yticks([-0.5 0.5]); 
 ylimh=max(abs(max(hx/B0)), abs(min(hx/B0)));
 ylim([-ylimh ylimh]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{x}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 yyaxis left
 plot(d,viy/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--');
 xlim([0 30]);
 %xlim([0 24]); 
 ylimi=max(abs(max(viy/Vai)), abs(min(viy/Vai)));
 ylim([-ylimi ylimi]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iy}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hy/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 %xlim([0 24]); 
 xlim([0 30]);
 ylimh=max(abs(max(hy/B0)), abs(min(hy/B0)));
 ylim([-ylimh ylimh]);
 %%yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{y}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 yyaxis left
 plot(d,viz/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--'); 
 %xlim([0 24]);
 xlim([0 30]);
 ylim1=max(abs(max(viz/Vai)), abs(min(viz/Vai)));
 ylim([-ylim1 ylim1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iz}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hz/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 30]);
 %xlim([0 24]);  
 ylimh=max(abs(max(hz/B0)), abs(min(hz/B0)));
 ylim([-ylimh ylimh]);
 %%yticks([0 1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{z}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21]);
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 %hL=legend({'$v_{i}$', '$B$' },'Location','northwest','Interpreter','latex');
 hL = legend([line1,line2],{'$v_{i}$', '$B$'},'Interpreter','latex');
 set(hL);

 
 %-------------------------------------------------------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Calculating the covariance point by point
 t2=d; Xx=vim/Vai; Yx=hm/B0;
 [t,dXdY]=Local_corr(t2,Xx,Yx,10,10);
 f49=figure(49);
 plot(t,dXxdYx,'k')
 
 f444=figure(444);
 plot(t2,Xx,'--r',t2,Yx,'--g')
 hold on
 plot(t2,Xx,'-k',t2,Yx,'-k')
 hold off

 
 Yx=hm/B0;
 f44=figure(44)
 scatter3(d,Xx,Yx)
 xlabel('d')
 ylabel('vi')
 zlabel('B')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % This is the one in order x,y,z,m
 %-------------------------------------------------------------------------
 f8=figure(8);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 yyaxis left
 line1 = plot(d,vim/Vai,'r','LineWidth',2); set(gca,'FontSize',20)
 xlim([0 30]);
 %xlim([0 30]);
 maxi=abs(max(vim/Vai));
 mini=abs(min(vim/Vai));
 ylim1=max(maxi,mini);
 ylim([0 ylim1]);
 %%yticks([0.3 2 3]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{v_{i}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 yyaxis right
 line2 = plot(d, hm/B0, 'k','LineWidth',2);set(gca,'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlim([0 30]);
 %xlim([0 30]);
 ylim([0 max(hm/B0)]) 
 %%yticks([1 2]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{B}{B_{0}}\right|$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 yyaxis left
 plot(d,viz/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--'); 
 xlim([0 30]);
 %lim([0 30]);
 ylim1=max(abs(max(viz/Vai)), abs(min(viz/Vai)));
 ylim([-ylim1 ylim1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iz}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hz/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 30]);
 %xlim([0 24]);  
 ylimh=max(abs(max(hz/B0)), abs(min(hz/B0)));
 ylim([-ylimh ylimh]);
 %%yticks([0 1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{z}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 yyaxis left
 plot(d,viy/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--');
 xlim([0 30]);
 %xlim([0 24]); 
 ylimi=max(abs(max(viy/Vai)), abs(min(viy/Vai)));
 ylim([-ylimi ylimi]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iy}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hy/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 %xlim([0 24]); 
 xlim([0 30]);
 ylimh=max(abs(max(hy/B0)), abs(min(hy/B0)));
 ylim([-ylimh ylimh]);
 %%yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{y}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 yyaxis left
 plot(d,vix/Vai,'r','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 yline(0,'--'); 
 %xline(12,'--'); 
 %xlim([0 24]); 
 xlim([0 30]);
 ylimi=max(abs(max(vix/Vai)), abs(min(vix/Vai)));
 ylim([-ylimi ylimi]);
 %%yticks([-0.1 0.1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ix}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hx/B0, 'k','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 %xlim([0 24]); 
 xlim([0 30]);
 %%yticks([-0.5 0.5]); 
 ylimh=max(abs(max(hx/B0)), abs(min(hx/B0)));
 ylim([-ylimh ylimh]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{x}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21]);
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 %legend({'$v_{i}$', '$B$' },'Location','northwest','Interpreter','latex')
 hL = legend([line1,line2],{'$v_{i}$', '$B$'},'Interpreter','latex');
 set(hL);
 %-------------------------------------------------------------------------
 
 
 %-------------------------------------------------------------------------
 f33=figure(33); 
 subplot(4,1,1),
 Xx=vix/Vai;
 Yx=hx/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 yline(0,'--')
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim([-ylim1 ylim1])
 xlim([0 30])
 yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(v_{x},B_{x})$$','Interpreter','latex','FontSize',20)
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 Xx=viy/Vai;
 Yx=hy/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 yline(0,'--')
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim([-ylim1 ylim1])
 xlim([0 30])
 yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(v_{y},B_{y})$$','Interpreter','latex','FontSize',20)
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 Xx=viz/Vai;
 Yx=hz/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 yline(0,'--')
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim([-ylim1 ylim1])
 xlim([0 30])
 yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(v_{z},B_{z})$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 Xx=vim/Vai;
 Yx=hm/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'FontSize',20);
 yline(0,'--')
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim([-ylim1 ylim1])
 xlim([0 30])
 yticks([-0.5 0.5])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(|v|,|B|)$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(4),'position',[.12 .73 .8 .21]); set(ha(3),'position',[.12 .52 .8 .21]);
 set(ha(2),'position',[.12 .31 .8 .21]); set(ha(1),'position',[.12 .11 .8 .21]);
 legend({'$T_{1}$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 
 
 
 f34=figure(34); 
 subplot(4,1,1),
 Xx=vix/Vai; Yx=hx/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 line1 = plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 hold on
 Xx=vix2/Vai; Yx=hx2/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 line2 = plot(d,covXxYx_i,'r','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 hold off
 yline(0,'--')
 ylim2 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 yticks([-0.5 0.5])
 xlim([0 24])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(v_{x},B_{x})$$','Interpreter','latex','FontSize',20)
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %hold off
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 Xx=viy/Vai; Yx=hy/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 hold on
 Xx=viy2/Vai; Yx=hy2/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'r','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 hold off
 yline(0,'--')
 ylim2 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 yticks([-0.5 0.5])
 xlim([0 24])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(v_{y},B_{y})$$','Interpreter','latex','FontSize',20)
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 Xx=viz/Vai; Yx=hz/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 hold on
 Xx=viz2/Vai; Yx=hz2/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'r','LineWidth',2); set(gca,'xtick',[],'FontSize',20)
 hold off
 yline(0,'--')
 ylim2 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 yticks([-0.5 0.5])
 xlim([0 24])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(v_{z},B_{z})$$','Interpreter','latex','FontSize',20)
 %xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 Xx=vim/Vai; Yx=hm/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'k','LineWidth',2); set(gca,'FontSize',20);
 ylim1 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 hold on
 Xx=vim2/Vai; Yx=hm2/B0;
 covXxYx_i=((Xx-mean(Xx)).*(Yx-mean(Yx)))./((length(Xx)-1)*std(Xx)*std(Yx));
 covXxYx_i=covXxYx_i/max(covXxYx_i);
 plot(d,covXxYx_i,'r','LineWidth',2); set(gca,'FontSize',20);
 hold off
 yline(0,'--')
 ylim2 =max(abs(min(covXxYx_i)),max(covXxYx_i));
 ylim3=max(ylim1,ylim2);
 ylim([-ylim3 ylim3])
 yticks([-0.5 0.5])
 xlim([0 24])
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$\rho(|v|,|B|)$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(4),'position',[.12 .73 .8 .21]); set(ha(3),'position',[.12 .52 .8 .21]);
 set(ha(2),'position',[.12 .31 .8 .21]); set(ha(1),'position',[.12 .11 .8 .21]);
 hL = legend([line1,line2],{'$T_{2}$', '$T_{3}$'},'Interpreter','latex');
 set(hL);
 %legend({'$v_{iT2}$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 
 
 
 
 %-------------------------------------------------------------------------
 f83=figure(83);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 yyaxis left
 plot(d,vim/Vai,'-r',d,vim2/Vai,'-m' ,'LineWidth',2),set(gca,'FontSize',20)
 xlim([0 24]); 
 ymax=max(max(vim/Vai), max(vim2/Vai));
 ylim([0 ymax]);
 yticks([0.3 2 3]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{v_{i}}{V_{A0}}\right|$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d, hm/B0, '-k',d, hm2/B0, '-c','LineWidth',2),set(gca,'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlim([0 24]); 
 ylim([0 max(hm/B0)]) 
 %yticks([1 2]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{B}{B_{0}}\right|$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 yyaxis left
 plot(d,vix/Vai,'-r',d,vix2/Vai,'-m','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 yline(0,'--'); 
 xline(12,'--'); 
 xlim([0 24]); 
 %ylim([-0.7 0.7]); 
 ymax1=abs(max(max(vix/Vai), max(vix2/Vai)));
 ymin1=abs(min(min(vix/Vai), min(vix2/Vai)));
 ymax=max(ymax1, ymin1);
 ylim([-ymax ymax]);
 yticks([-0.1 0.1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ix}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hx/B0, '-k',d, hx2/B0, '-c','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
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
 subplot(4,1,3);
 yyaxis left
 plot(d,viy/Vai,'-r',d,viy2/Vai,'-m','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
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
 plot(d,hy/B0, '-k',d, hy2/B0, '-c','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
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
 plot(d,viz/Vai,'-r',d,viz2/Vai,'-m','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); 
 xlim([0 24]);  %ylim([-0.5 2]); yticks([0.5 1.5]);
 ymax1=abs(max(max(viz/Vai), max(viz2/Vai)));
 ymin1=abs(min(min(viz/Vai), min(viz2/Vai)));
 ymax=max(ymax1, ymin1);
 ylim([-ymax ymax]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{iz}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hz/B0, '-k',d, hz2/B0, '-c','LineWidth',2),set(gca,'xtick',[],'FontSize',20)
 yline(0,'--');
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); 
 xlim([0 24]);  
 ylimh=max(abs(max(hz/B0)), abs(min(hz/B0)));
 ylim([-ylimh ylimh]);
 yticks([0 1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{z}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21]);
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 legend({'$v_{iT2}$','$v_{iT3}$','$B_{T2}$','$B_{T3}$'},'Location','northwest','Interpreter','latex')
 %-------------------------------------------------------------------------
 
 
 
 % This is the one with x,y,z
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
 
 
 
 x=ni;
 y=ne;
 
 x1 = x(1:1000);
 x2 = x(2:1001);
 dx=x2-x1; dx2=dx.*dx;
 normdx=sqrt(sum(dx2));
 dxn=dx/normdx;
 
 y1 = y(1:1000);
 y2 = y(2:1001);
 dy=y2-y1;
 dy2=dy.*dy;
 normdy=sqrt(sum(dy2));
 dyn=dy/normdy;

 %error=(dyn-dxn)./dxn;

 
figure333=figure(333)
plot(d(1:1000),dxn,d(1:1000),dyn)
legend({'$v_{ix}$', '$Bx$' },'Location','northwest','Interpreter','latex')

figure334=figure(334)
 [c,lags] = xcorr(dxn,dyn);
 stem(lags,c)
 
 figure335=figure(335)
 [c,lags] = xcorr(x,y);
 stem(lags,c)
 
% t=
 figure336=figure(336)
 scatter(x,y,'b')
 hold on
 plot(x,x,'r')
 hold off
 
  % subset of B and vi components
 f9=figure(9);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 yyaxis left
 plot(d,vem/Vai,'b'),set(gca,'FontSize',20)
 xlim([0 30]); 
 %ylim([0 0.5]); 
 yticks([1 3]);
 set(gca,'ycolor','k') 
 ylabel('$$\left| \frac{v_{e}}{V_{A0}} \right|$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d, hm/B0, 'k'),set(gca,'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlim([0 30]); 
 %ylim([0 2.5]); 
 yticks([1 2]);
 set(gca,'ycolor','k') 
 ylabel('$$\left|\frac{B}{B_{0}}\right|$$','Interpreter','latex','FontSize',20)
 xlabel('$r/d_{i}$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 yyaxis left
 plot(d,vex/Vai,'b'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 xlim([0 30]);
 yline(0,'--'); 
 xline(12,'--'); 
 %ylim([-0.3 0.3]); yticks([-0.1 0.1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ex}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hx/B0, 'k'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 xlim([0 30]); 
 yticks([-0.5 0.5]); %ylim([-0.3 0.3]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{x}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 yyaxis left
 plot(d,vey/Vai,'b'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 30]); %ylim([-0.7 0.7]); yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ey}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hy/B0, 'k'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 30]); 
 %ylim([-1 1]); 
 yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{y}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -2 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 yyaxis left
 plot(d,vez/Vai,'b'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 30]);  %ylim([-0.5 2]); 
 yticks([-1.5 1.5]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{v_{ez}}{V_{A0}}$$','Interpreter','latex','FontSize',20)
 yyaxis right
 plot(d,hz/B0, 'k'),set(gca,'xtick',[],'FontSize',20)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 xlim([0 30]);  
 %ylim([0 2]); 
 yticks([0 1]);
 set(gca,'ycolor','k') 
 ylabel('$$\frac{B_{z}}{B_{0}}$$','Interpreter','latex','FontSize',20)
 rectangle('Position', [start_p -0.5 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21]);
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 legend({'$v_{e}$', '$B$' },'Location','northwest','Interpreter','latex')

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Save the plots
  %-------------------------------------------------------------------------
  cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images';
  cd '/Volumes/PSC_DiRAC_DATA/';

  %saveas(f1,'trajectory_JEBxyz_1.png');
  %saveas(f2,'trajectory_BnvieTie_6.png');
  %saveas(f3,'trajectory_Bmxyz_11_square.png');
  %saveas(f4,'trajectory_Emxyz_6.png');
  %saveas(f5,'trajectory_BnJE_6.png');
  %saveas(f6,'trajectory_BnJETie_6.png');
  saveas(f7,'trajectory_BnvieTie_8_square.png');
  saveas(f8,'trajectory_viB_10_xyz_square.png');
  saveas(f83,'trajectory_viB_T2T3_xyz.png');
  saveas(f9,'trajectory_veB_10_square.png');
  saveas(f777,'trajectory_10vive_csm.png');
  saveas(f33,'trajectory_rhovb_10_xyz_square.png');
  saveas(f34,'trajectory_rhovb_T1T2_xyz_square.png');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  
  
  %{
d=exyz1(:,4);ex=exyz1(:,1);ey=exyz1(:,2);ez=exyz1(:,3);
hx=hxyz1(:,1);hy=hxyz1(:,2);hz=hxyz1(:,3);
jx=jxyz1(:,1);jy=jxyz1(:,2);jz=jxyz1(:,3);
vix=vixyz1(:,1);viy=vixyz1(:,2);viz=vixyz1(:,3);
vex=vexyz1(:,1);vey=vexyz1(:,2);vez=vexyz1(:,3);
Tix=Tixyz1(:,1);Tiy=Tixyz1(:,2);Tiz=Tixyz1(:,3);
Tex=Texyz1(:,1);Tey=Texyz1(:,2);Tez=Texyz1(:,3);
ne=nine1(:,1);ni=nine1(:,2);

d=exyz4(:,4);ex=exyz4(:,1);ey=exyz4(:,2);ez=exyz4(:,3);
hx=hxyz4(:,1);hy=hxyz4(:,2);hz=hxyz4(:,3);
jx=jxyz4(:,1);jy=jxyz4(:,2);jz=jxyz4(:,3);
vix=vixyz4(:,1);viy=vixyz4(:,2);viz=vixyz4(:,3);
vex=vexyz4(:,1);vey=vexyz4(:,2);vez=vexyz4(:,3);
Tix=Tixyz4(:,1);Tiy=Tixyz4(:,2);Tiz=Tixyz4(:,3);
Tex=Texyz4(:,1);Tey=Texyz4(:,2);Tez=Texyz4(:,3);
ne=nine4(:,1);ni=nine4(:,2);

d=exyz5(:,4);ex=exyz5(:,1);ey=exyz5(:,2);ez=exyz5(:,3);
hx=hxyz5(:,1);hy=hxyz5(:,2);hz=hxyz5(:,3);
jx=jxyz5(:,1);jy=jxyz5(:,2);jz=jxyz5(:,3);
vix=vixyz5(:,1);viy=vixyz5(:,2);viz=vixyz5(:,3);
vex=vexyz5(:,1);vey=vexyz5(:,2);vez=vexyz5(:,3);
Tix=Tixyz5(:,1);Tiy=Tixyz5(:,2);Tiz=Tixyz5(:,3);
Tex=Texyz5(:,1);Tey=Texyz5(:,2);Tez=Texyz5(:,3);
ne=nine5(:,1);ni=nine5(:,2);

d=exyz6(:,4);ex=exyz6(:,1);ey=exyz6(:,2);ez=exyz6(:,3);
hx=hxyz6(:,1);hy=hxyz6(:,2);hz=hxyz6(:,3);
jx=jxyz6(:,1);jy=jxyz6(:,2);jz=jxyz6(:,3);
vix=vixyz6(:,1);viy=vixyz6(:,2);viz=vixyz6(:,3);
vex=vexyz6(:,1);vey=vexyz6(:,2);vez=vexyz6(:,3);
Tix=Tixyz6(:,1);Tiy=Tixyz6(:,2);Tiz=Tixyz6(:,3);
Tex=Texyz6(:,1);Tey=Texyz6(:,2);Tez=Texyz6(:,3);
ne=nine6(:,1);ni=nine6(:,2);

d=exyz7(:,4);ex=exyz7(:,1);ey=exyz7(:,2);ez=exyz7(:,3);
hx=hxyz7(:,1);hy=hxyz7(:,2);hz=hxyz7(:,3);
jx=jxyz7(:,1);jy=jxyz7(:,2);jz=jxyz7(:,3);
vix=vixyz7(:,1);viy=vixyz7(:,2);viz=vixyz7(:,3);
vex=vexyz7(:,1);vey=vexyz7(:,2);vez=vexyz7(:,3);
Tix=Tixyz7(:,1);Tiy=Tixyz7(:,2);Tiz=Tixyz7(:,3);
Tex=Texyz7(:,1);Tey=Texyz7(:,2);Tez=Texyz7(:,3);
ne=nine7(:,1);ni=nine7(:,2);
%}
  
 %-------------------------------------------------------------------------
 figure666=figure(666);
 plot(d,Vailocal/Vai,'b'),set(gca,'FontSize',15)
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 24]); ylim([0 2]); yticks([0.5 1.5]); 
 %xlim([0 30]); ylim([-3.5 3.5]); yticks([-1.5 0 1.5]);
 yline(0,'--'); 
 %xline(jmax_pos,'--'); xline(hmin_pos,'--');
 ylabel('$$\tilde{B}_{x}$$','Interpreter','latex')
 rectangle('Position', [start_p -3.5 delta 7], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------

 
 
 
 f9=figure(9);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 plot( d,vem/Vai,'b', d, hm/B0, 'k'),set(gca,'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlim([10 20]); ylim([0 8]); yticks([2 4 6]);
 xlabel('$r (d_{i})$','Interpreter','latex')
 ylabel('$$|\tilde{\psi}|$$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 8], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 plot(d,vex/Vai,'b', d,hx/B0, 'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 yline(0,'--');
 xlim([10 20]); ylim([-3.5 2.5]); yticks([-1.5 1.5]);
 ylabel('$$\tilde{\psi}_{x}$$','Interpreter','latex')
 rectangle('Position', [start_p -3.5 delta 6], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 plot(d,vey/Vai,'b', d,hy/B0, 'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--');
 xlim([10 20]); ylim([-2.5 2.5]); yticks([-1.5 1.5]);
 ylabel('$$\tilde{\psi}_{y}$$','Interpreter','latex')
 rectangle('Position', [start_p -2.5 delta 5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 plot(d,vez/Vai,'b', d,hz/B0, 'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 yline(0,'--');
 xlim([10 20]); ylim([-2.5 5.5]); yticks([-1.5 1.5]);
 ylabel('$$\tilde{\psi}_{z}$$','Interpreter','latex')
 rectangle('Position', [start_p -2.5 delta 8], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21]);
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 legend({'$v_{e}$', '$B$' },'Location','northwest','Interpreter','latex')

 
 f71=figure(71);
 %-------------------------------------------------------------------------
 subplot(5,1,1),
 plot(d,vim./Vailocal,'r'),set(gca,'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 24]); ylim([0 0.8]); yticks([0.5]);
 ylabel('$$|\tilde{v}_{i}|$$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 1], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,2),
 plot(d,vem./Vailocal,'b'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 set(gca,'FontSize',15)
 xlim([0 24]); ylim([0 3]); yticks([1 2 5]);%6 8 10]);
 ylabel('$$|\tilde{v}_{e}|$$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 6], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,3);
 plot(d,ni/n0,'c'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 24]); ylim([0 1.5]); yticks([0.5 1]);
 ylabel('$$\tilde{n}_{i}$$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,4);
 plot(d,hm/B0,'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 24]);  ylim([0 2]); yticks([1 2]);
 ylabel('$$|\tilde{B}|$$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 3], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(5,1,5);
 plot(d,Tim/T0,'r', d,Tem/T0,'b'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 xlim([0 24]); ylim([0 2]); yticks([1 2 3]);
 ylabel('$$\tilde{T}_{i}, \tilde{T}_{e}$$','Interpreter','latex')
 rectangle('Position', [start_p 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.1 .81 .8 .18]); set(ha(2),'position',[.1 .63 .8 .18]); 
 set(ha(3),'position',[.1 .45 .8 .18]); set(ha(4),'position',[.1 .27 .8 .18]); 
 set(ha(5),'position',[.1 .09 .8 .18]);
 legend({'$i$','$e$'},'Location','northwest','Interpreter','latex')

 
 
  % subset of B and vi components
 f81=figure(81);
 %-------------------------------------------------------------------------
 subplot(4,1,1),
 yyaxis left
 plot(d,vim./Vailocal,'r'),set(gca,'FontSize',15)
 %xlim([10 17]); 
 ylim([0 0.5]); yticks([0.3 2 3]);
 set(gca,'ycolor','k') 
 ylabel('$$|\tilde{v}_{i}|$$','Interpreter','latex')
 yyaxis right
 plot(d, hm/B0, 'k'),set(gca,'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %xlim([10 17]); 
 ylim([0 2.5]); yticks([1 2]);
 set(gca,'ycolor','k') 
 ylabel('$$|\tilde{B}|$$','Interpreter','latex')
 xlabel('$r (d_{i})$','Interpreter','latex')
 %rectangle('Position', [start_p 0 delta 2.5], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,2);
 yyaxis left
 plot(d,vix./Vailocal,'r'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 %xlim([10 17]); 
 ylim([-0.3 0.3]); yticks([-0.1 0.1]);
 set(gca,'ycolor','k') 
 ylabel('$$\tilde{v}_{ix}$$','Interpreter','latex')
 yyaxis right
 plot(d,hx/B0, 'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--'); %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %yline(0,'--'); 
 %xlim([10 17]); 
 yticks([0.5]); %ylim([-0.3 0.3]);
 set(gca,'ycolor','k') 
 ylabel('$$\tilde{B}_{x}$$','Interpreter','latex')
 %rectangle('Position', [start_p -1.8 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,3);
 yyaxis left
 plot(d,viy./Vailocal,'r'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 %xlim([10 17]); %ylim([-0.7 0.7]); yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\tilde{v}_{iy}$$','Interpreter','latex')
 yyaxis right
 plot(d,hy/B0, 'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 %xlim([10 17]); 
 ylim([-0.7 0.7]); yticks([-0.5 0.5]);%6 8 10]);
 set(gca,'ycolor','k') 
 ylabel('$$\tilde{B}_{y}$$','Interpreter','latex')
 %rectangle('Position', [start_p -1.8 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 subplot(4,1,4);
 yyaxis left
 plot(d,viz./Vailocal,'r'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 %xlim([10 17]);  %ylim([-0.5 2]); yticks([0.5 1.5]);
 set(gca,'ycolor','k') 
 ylabel('$$\tilde{v}_{iz}$$','Interpreter','latex')
 yyaxis right
 plot(d,hz/B0, 'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %yline(0,'--'); 
 %xlim([10 17]);  
 ylim([0 1.5]); yticks([1]);
 set(gca,'ycolor','k') 
 ylabel('$$\tilde{B}_{z}$$','Interpreter','latex')
 %rectangle('Position', [start_p -0.5 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 %-------------------------------------------------------------------------
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21]);
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21]);
 %legend({'$v_{i}$', '$B$' },'Location','northwest','Interpreter','latex')

  
  
  
  
  

 subplot(2,2,1),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 subplot(2,2,2),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 subplot(2,2,3),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 subplot(2,2,4),plot(rand(3)),set(gca,'xtick',[],'ytick',[])
 ha=get(gcf,'children');
 set(ha(1),'position',[.5 .1 .4 .4])
 set(ha(2),'position',[.1 .1 .4 .4])
 set(ha(3),'position',[.5 .5 .4 .4])
 set(ha(4),'position',[.1 .5 .4 .4])
 

 
 x = linspace(0, 10, 10);
y = randi(9, 1, 10);
figure(1)
hp = plot(x, y, 'bp');
hold on
ybars = [2 6];
patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
plot(x, y, 'bp')
hold off
axis([0  10    0  10])


 f1=figure(1);
 subplot(3,1,1),
 plot(d,hx/rms(hx),'b', d,hy/rms(hy),'r', d,hz/rms(hz),'c'),set(gca,'FontSize',15)
 y1=get(gca,'ylim'); %xline(40);
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 30]); %yticks([-1 1]);
 ylabel('$$\tilde{B}_{i}$$','Interpreter','latex')
 subplot(3,1,2),
 plot(d,ex/rms(ex),'b', d,ey/rms(ey),'r', d,ez/rms(ez),'c'),set(gca,'xtick',[], 'FontSize',15)
 xlim([0 30]); yticks([-2 2]); %xline(40);
 ylabel('$$\tilde{E}_{i}$$','Interpreter','latex')
 subplot(3,1,3),
 plot(d,jx/rms(jx),'b', d,jy/rms(jy),'r', d,jz/rms(jz),'c'),set(gca,'xtick',[], 'FontSize',15) 
 xlim([0 30]); %xline(40);
 ylabel('$$\tilde{J}_{i}$$','Interpreter','latex')
 ha=get(gcf,'children');
 set(ha(1),'position',[.1 .68 .8 .29]); set(ha(2),'position',[.1 .4 .8 .29])
 set(ha(3),'position',[.1 .11 .8 .29])
 legend({'x','y','z'},'Interpreter','latex')
 %yticks([-1 1]);
 %hold off

 f5=figure(5);
 subplot(4,1,1),
 %plot(d,em/rms(em),'b'),set(gca,'FontSize',15)
 plot(d,Tim/rms(Tim),'b', d,Tem/rms(Tem),'r'),set(gca,'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 30]); ylim([0 4]); yticks([1 2 3]);
 %ylabel('$$|\tilde{E}|$$','Interpreter','latex')
 ylabel('$$\tilde{T}_{i,e}$$','Interpreter','latex')
 subplot(4,1,2),
 plot(d,jm/rms(jm),'r'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 set(gca,'FontSize',15)
 xlim([0 30]); ylim([0 4]); yticks([1 2 3]);
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$|\tilde{J}|$$','Interpreter','latex')
 subplot(4,1,3);
 plot(d,ni/rms(ni),'c'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 30]); ylim([0 2]); yticks([0.5 1.5]);
 ylabel('$$\tilde{n}_{i}$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(4,1,4);
 plot(d,hm/rms(hm),'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 30]);  ylim([0 2]); yticks([0.5 1.5]);
 ylabel('$$|\tilde{B}|$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ha=get(gcf,'children');
 set(ha(1),'position',[.1 .73 .8 .21]); set(ha(2),'position',[.1 .52 .8 .21])
 set(ha(3),'position',[.1 .31 .8 .21]); set(ha(4),'position',[.1 .11 .8 .21])

 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epar=(ex.*hx+ey.*hy+ez.*hz)./hm;
y_hxmax = max(abs(hx)); y_hymax = max(abs(hy)); y_hzmax = max(abs(hz));
y_hxyz = max(max(y_hxmax, y_hymax), y_hzmax);

%position of the maximum jm
aa=find(jm==max(jm));
jmax_pos=aa(1,1)*jxyz4(3,4); %The position of the first value

aa=find(hm==min(hm));
hmin_pos=aa(1,1)*hxyz4(3,4);
%delta = jmax_pos - hmin_pos;
delta=9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make the plots

 % For the files1 the xlime is up to 130

 
 f2=figure(2);
 subplot(4,1,1),
 plot(d,hm/B0,'k'),set(gca,'FontSize',15);%plot(d,hm/rms(hm),'k'),set(gca,'FontSize',15)
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 45]); ylim([0 3]); yticks([1 2]);
 ylabel('$$|\tilde{B}|$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(4,1,2),
 plot(d,ni/n0,'c'),set(gca,'xtick',[],'FontSize',15);%plot(d,ni/rms(ni),'c'),set(gca,'xtick',[],'FontSize',15)
 set(gca,'FontSize',15)
 xlim([0 45]); ylim([0 2.5]);  yticks([0.5 1 1.5]);
 ylabel('$$\tilde{n}_{i}$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(4,1,3),
 plot(d,Tim/T0,'b', d,Tem/T0,'r'),set(gca,'xtick',[],'FontSize',15)
 set(gca,'FontSize',15)
 xlim([0 45]); ylim([0 4]);  yticks([1 2 3]);
 ylabel('$$\tilde{T}_{i}, \tilde{T}_{e}$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(4,1,4);
 semilogy(d,vim/Vai,'b',d,vem/Vai,'r'),set(gca,'xtick',[],'FontSize',15)
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 xlim([0 45]); 
 ylabel('$$|\tilde{v}_{i}|, |\tilde{v}_{e}|$$','Interpreter','latex')
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21])
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21])
 legend({'$i$','$e$'},'Location','northwest','Interpreter','latex')
 %legend({'$|\tilde{v}_{i}|$','$|\tilde{v}_{e}|$'},'Location','northwest','Interpreter','latex')
 
  
 % This makes the plot of the electric field components
 f4=figure(4);
 subplot(4,1,1),
 plot(d,ex/E0,'b'),set(gca,'FontSize',15)
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 24]); ylim([-1 1]);yticks([-0.3 0 0.3]);yline(0,':');%ylim([-5 5]);yticks([-2 0 2]);yline(0,':');
 ylabel('$$\tilde{E}_{x}$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos -5 delta 10], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(4,1,2),
 plot(d,ey/E0,'r'),set(gca,'xtick',[],'FontSize',15)
 set(gca,'FontSize',15)
 xlim([0 45]); ylim([-1 1]);yticks([-0.3 0 0.3]);yline(0,':');
 ylabel('$$\tilde{E}_{y}$$','Interpreter','latex')
  %rectangle('Position', [hmin_pos -5 delta 10], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(4,1,3);
 plot(d,ez/E0,'c'),set(gca,'xtick',[],'FontSize',15)
 %plot(d,epar/rms(epar),'c'),set(gca,'xtick',[],'FontSize',15)
 xlim([0 45]); ylim([-1 1]);yticks([-0.3 0 0.3]);yline(0,':');
 ylabel('$$\tilde{E}_{z}$$','Interpreter','latex')
  %rectangle('Position', [hmin_pos -5 delta 10], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(4,1,4);
 plot(d,em/E0,'k'),set(gca,'xtick',[],'FontSize',15)
 xlim([0 45]); ylim([0 0.5]); yticks([0.3]);
 ylabel('$$|\tilde{E}|$$','Interpreter','latex')
  %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ha=get(gcf,'children');
 set(ha(1),'position',[.12 .73 .8 .21]); set(ha(2),'position',[.12 .52 .8 .21])
 set(ha(3),'position',[.12 .31 .8 .21]); set(ha(4),'position',[.12 .11 .8 .21])
 
 f6=figure(6);
 subplot(5,1,1),
 plot(d,em/E0,'b'),set(gca,'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 xlabel('$r (d_{i})$','Interpreter','latex')
 xlim([0 45]); ylim([0 0.5]); yticks([0.3]);
 ylabel('$$|\tilde{E}|$$','Interpreter','latex')
 subplot(5,1,2),
 plot(d,jm/J0,'r'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 set(gca,'FontSize',15)
 xlim([0 45]); ylim([0 6]); yticks([2 4 ]);%6 8 10]);
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ylabel('$$|\tilde{J}|$$','Interpreter','latex')
 subplot(5,1,3);
 plot(d,ni/n0,'c'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 45]); ylim([0 2.5]); yticks([0.5 1.5]);
 ylabel('$$\tilde{n}_{i}$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 subplot(5,1,4);
 plot(d,hm/B0,'k'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 xlim([0 45]);  ylim([0 3]); yticks([1 2]);
 ylabel('$$|\tilde{B}|$$','Interpreter','latex')
 subplot(5,1,5);
 plot(d,Tim/T0,'b', d,Tem/T0,'r'),set(gca,'xtick',[],'FontSize',15)
 %xline(jmax_pos,'--');  xline(hmin_pos,'--');
 %ybars = [0 4]; patch([hmin_pos jmax_pos jmax_pos hmin_pos], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 xlim([0 45]); ylim([0 4]); yticks([1 2 3]);
 %ylabel('$$|\tilde{E}|$$','Interpreter','latex')
 ylabel('$$\tilde{T}_{i}, \tilde{T}_{e}$$','Interpreter','latex')
 %rectangle('Position', [hmin_pos 0 delta 4], 'FaceColor', [0.7 0.7 0.7 0.2])
 ha=get(gcf,'children');
 set(ha(1),'position',[.1 .81 .8 .18]); set(ha(2),'position',[.1 .63 .8 .18]); 
 set(ha(3),'position',[.1 .45 .8 .18]); set(ha(4),'position',[.1 .27 .8 .18]); 
 set(ha(5),'position',[.1 .09 .8 .18]);
 legend({'$i$','$e$'},'Location','northwest','Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Defining functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,dXdY]=Local_corr(t2,Xx,Yx,dd,window)
 % time/distance, var1, var2, increment, window  
 %Define the arrays
 %Ex: t2=d; Xx=vim/Vai; Yx=hm/B0;
 %Define the increment/window
 %dd=10; window=10; % 10 is the optimun value that do not shift the plot
 %Filter data
 b = (1/window)*ones(1,window);
 a = 1;
 Xx = filter(b,a,Xx);
 Yx = filter(b,a,Yx);  
 %Calculate the increments (derivatives) 
 Xx2=Xx(1+dd:length(Xx)); Xx1=Xx(1:length(Xx)-dd);
 dXx=(Xx2 - Xx1)/dd;
 Yx2=Yx(1+dd:length(Yx)); Yx1=Yx(1:length(Yx)-dd);
 dYx=(Yx2 - Yx1)/dd;
 %outputs
 dXdY=dXx.*dYx;
 t=t2(1:length(t2)-dd);
 end
 