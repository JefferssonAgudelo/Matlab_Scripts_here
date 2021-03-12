%%%%%%%%%%% Checking the initial distribution parallel distribution %%%%%%%
%%----------------------------------------------------------------------
%// ----------------------------------------------------------------------
%// Parameters for background plasma
 kb = 1.;
 mu0 =1.;
 mi_over_me_ =300.;
 vA_over_c_ = .008;
 B0_ = vA_over_c_;
 mi = 1.;
 me = 1. / mi_over_me_;
 mes = 1. / mi_over_me_;

 omega_pi=1;
 Omega_e=vA_over_c_*sqrt(mi_over_me_)*omega_pi;
 Omega_i=Omega_e / mi_over_me_;
 TOmega_e=1/Omega_e; %2*pi/Omega_e;
 TOmega_i=1/Omega_i;
 Ttotal=300*TOmega_e;
 dt=0.0033;
 %dt=3e-5;
 %dt=0.001;
 %dt=0.04;
 N_tsteps=Ttotal/dt;
 N_Oe=TOmega_e/dt; % this is the number of time steps required for 1 TOmega_e
 
 t_final=80000;
 time = (t_final/N_tsteps);
 time_Oe = t_final/N_Oe;
 
 %----------------------------------------------------------------------
 
 %%These values are in units of c from the fitting   

    %// These values are in units of c from the fitting    
    ni_ = 1.;
    ne_core_ = 0.978075171;
    ne_strahl_ = 0.021924829;
    beta_i_par_ = 1.200545495;
    beta_e_core_par_ = 3.506865593;
    beta_e_strahl_par_ = 0.053475199;
    %//These are the temperature anisotropies
    Ti_per_over_Ti_par_ = 1.;
    Te_c_per_over_Te_c_par_ = 1.118098593;
    Te_s_per_over_Te_s_par_ = 2.400797005;
    %//This are the parallel values of the temperature used in this simulation 
    Ti_par_ = beta_i_par_ * (B0_*B0_) / ( 2. * kb * mu0 *ni_ );
    Te_core_par_ = beta_e_core_par_ * (B0_*B0_) / ( 2. * kb * mu0 *ne_core_ );
    Te_strahl_par_ = beta_e_strahl_par_ * (B0_*B0_) / (2. * kb * mu0 *ne_strahl_);
    Ti_per_ = Ti_per_over_Ti_par_ * Ti_par_;
    Te_core_per_ = Te_c_per_over_Te_c_par_ * Te_core_par_;
    Te_strahl_per_ = Te_s_per_over_Te_s_par_ * Te_strahl_par_;
    vth_s_par = 2730000;
    VB_s_par = -5910000;
    %vth_s_par_over_VB_s = 1.0*vth_s_par / VB_s_par; 
    vth_s_par_over_VB_s = 0.3*vth_s_par / VB_s_par; 
    
    vA_i_over_c_ = vA_over_c_;
    vA_e_core_over_c_ = vA_over_c_ * ((mi * ni_) / (me * ne_core_));
    vA_e_strahl_over_c_ = vA_over_c_ * ((mi * ni_) / (me * ne_strahl_));
    
    %//double vth_s_per_over_VB_s = 0.38
    %// Bulk velocities. These values are in units of c from the fitting. Equivalent to (U/VA)*(VA/c)
    Vi_per_ = 0. * vA_over_c_;
    Vi_par_ = 0. * vA_over_c_ ;
    Ve_strahl_per_ = 0.* vA_over_c_;
    Ve_strahl_par_ = vA_over_c_ * (beta_e_strahl_par_ / vth_s_par_over_VB_s) * sqrt(mi_over_me_ / ne_strahl_);
    Ve_core_per_ = 0. * vA_over_c_;
    Ve_core_par_ = -Ve_strahl_par_ * ne_strahl_ / ne_core_; 
 
 
 %----------------------------------------------------------------------  
 npt.ni  = ni_; npt.Tix = Ti_per_; npt.Tiy = Ti_per_; npt.Tiz = Ti_par_;
 npt.pix = Vi_per_; npt.piy = Vi_per_; npt.piz = Vi_par_;

 npt.nec = ne_core_; npt.Tecx = Te_core_per_; npt.Tecy = Te_core_per_;
 npt.Tecz = Te_core_par_; npt.pecx = Ve_core_per_; npt.pecy = Ve_core_per_;
 npt.pecz = Ve_core_par_;

 npt.nes = ne_strahl_; npt.Tesx = Te_strahl_per_; npt.Tesy = Te_strahl_per_;
 npt.Tesz = Te_strahl_par_; npt.pesx = Ve_strahl_per_;
 npt.pesy = Ve_strahl_per_; npt.pesz = Ve_strahl_par_;
 %----------------------------------------------------------------------  

 rng('default')
 %NN=1000000;
 %NN=28800000;
 NN=7200000;
 %NN=108000000;
 %pzi_f=npt.piz * randn(NN,1);
 %Tzi_f=npt.Tiz * ones(NN,1);
 %mi_f=mi * ones(NN,1);
 %betai_f=beta_i_par_ * ones(NN,1);
 
 pzec_f=npt.pecz ;
 Tzec_f=npt.Tecz ;
 
 me_f=me ;
 betaec_f=beta_e_core_par_ ;
 
 pzes_f=npt.pesz ;
 Tzes_f=npt.Tesz ;
 %me_f=me * ones(NN,1);
 betaes_f=beta_e_strahl_par_ ;
 
 
 a=0; b=NN;
 %r = a + (b-a).*rand(NN,1);
 k=5;
 rng(k);
 NNi=1*NN;
 NNc=fix(ne_core_*NN);
 NNs=fix(ne_strahl_*NN);
 %ran1=rand(NN,1)/(NN+1); %randn doesnt allow to use pd.
 %ran1=randi([0 NN],NN,1)/(NN+1);
 %ran2=randi([0 NN],NN,1)/(NN+1);
 ran3=randi([0 NNc],NNc,1)/(NNc+1); %randn doesnt allow to use pd.
 ran4=randi([0 NNc],NNc,1)/(NNc+1);
 ran5=randi([0 NNs],NNs,1)/(NNs+1); %randn doesnt allow to use pd.
 ran6=randi([0 NNs],NNs,1)/(NNs+1);
 
 %ran1 = random() / ((float) RAND_MAX + 1);
% double pxi = npt->p[0] +
%     sqrtf(-2.f*npt->T[0]/npt->m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
 
 %betaec_f=1; betaes_f=1;
 %pzec_f_1 = pzec_f+sqrt(-2*(Tzec_f/me_f)*(betaec_f*betaec_f).*log(1-ran3)).* cos(2*pi*ran4);
 %pzes_f_1 = pzes_f+sqrt(-2*(Tzes_f/me_f)*(betaes_f*betaes_f).*log(1-ran5)).* cos(2*pi*ran6);
 %pzecs = pzec_f_1 + pzes_f_1; 
 %pxec_f_1 = pxec_f+sqrt(-2*(Txec_f/me_f)*(betaec_f*betaec_f).*log(1-ran3x)).* cos(2*pi*ran4x);
 %pxes_f_1 = pxes_f+sqrt(-2*(Txes_f/me_f)*(betaes_f*betaes_f).*log(1-ran5x)).* cos(2*pi*ran6x);
 %pxecs = pxec_f_1 + pxes_f_1;
 
 pzec_f_1 = pzec_f+sqrt(-(2*Tzec_f/me_f).*log(1 - ran3)).* cos(2*pi*ran4);
 pzes_f_1 = pzes_f+sqrt(-(2*Tzes_f/me_f).*log(1 - ran5)).* cos(2*pi*ran6);
 %pzecs = pzec_f_1 + pzes_f_1;
 


 %edges = linspace(-0.06,0.06);
 edges = linspace(-1,1);
 %linspace(-3,3);
 %[-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
 [Nc,edges] = histcounts(pzec_f_1,edges);
 [Ns,edges] = histcounts(pzes_f_1,edges);
 Ncs=Nc+Ns;
 f666=figure(666);
 semilogy(edges(2:100),Nc,'r') %plot(edges(2:100),Nc,'r')
 hold on
 semilogy(edges(2:100),Ns,'b')
 semilogy(edges(2:100),Ncs,'k')
 xline(0,'k')
 xlabel('$v_{s,\parallel}$','Interpreter','latex','FontSize',20)
 ylabel('Counts','Interpreter','latex','FontSize',20)
 title('$$Fake \ distribution$$','Interpreter','latex','FontSize',20)
 legend('$v_{ec,\parallel}$','$v_{es,\parallel}$','$v_{ecs,\parallel}$','Interpreter','latex','FontSize',20)
 hold off
 %saveas(f666,strcat('vepar_fake_counts_19_uns_50_0','.png'));
 
 [Ncpdf,edges] = histcounts(pzec_f_1,edges,'Normalization', 'probability');
 [Nspdf,edges] = histcounts(pzes_f_1,edges,'Normalization', 'probability');
 Ncspdf=Ncpdf+Nspdf;
 
 
 f667=figure(667);
 semilogy(edges(2:100),Nc/sum(Ncs),'-*r') % Everything is normalised to the total number of counts
 hold on
 %semilogy(edges(2:100),Ncpdf,'r')
 %semilogy(edges(2:100),Nspdf,'b')
 %semilogy(edges(2:100),Ncspdf,'k')
 semilogy(edges(2:100),Ns/sum(Ncs),'-*b')
 semilogy(edges(2:100),Ncs/sum(Ncs),'-*k')
 xline(0,'k')
  xlabel('$v_{s,\parallel}$','Interpreter','latex','FontSize',20)
 ylabel('Probability','Interpreter','latex','FontSize',20)
 title('$$Fake \ distribution$$','Interpreter','latex','FontSize',20)
 legend('$v_{ec,\parallel}$','$v_{es,\parallel}$','$v_{ecs,\parallel}$','Interpreter','latex','FontSize',20)
 hold off
 
 
 %saveas(f667,strcat('vepar_fake_prob_19_uns_50ext_0','.png'));
 %cd '/Volumes/PSC_DiRAC_DATA';
 %cd '/Volumes/PSC_DiRAC_DATA/ELECTRON_PROJECT/electron_project_mac/electron_test/Images_this_run';
 %saveas(f667,strcat('vepar_fake_prob_16296_03','.png'));
 %saveas(f667,strcat('vepar_fake_prob_19_01','.png'));
 
 
 