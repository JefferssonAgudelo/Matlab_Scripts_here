%------------------------------------------------------------------
% This is how to define the linear and non-linear times
%------------------------------------------------------------------
% The linear time is tl = 1/w where omega is the frequency of the
% waves that are solutions to the disspersion relation. There is a
% different definition for each type of wave. Alfven or fast waves.

% Plasma Constants in cgs to make it constant with most of the papers
% notation
%-------------------------------------------------------------------
mu0=1;
eps0=1;
kb=1;
c=1;
mi=1;
mime=100;
me= mi/mime;
qi=1;
qe=-qi;
di=1;
 

% Plasma parameters
%-------------------------------------------------------------------
Va= 0.1;
B0 = Va;

%Ions 
ni=1;
Ti=1;
Ti_per=1;
Ti_par_per = 1; 
Ti_par = Ti_per / Ti_par_per; 
Vth_i_per = sqrt(2 * Ti_per / mi);
Vth_i_par = sqrt(2 * Ti_par / mi);
beta_i= (8*pi*ni*kb*Ti) / (B0*B0);
beta_i_per= (8*pi*ni*kb*Ti_per) / (B0*B0);
beta_i_par= (8*pi*ni*kb*Ti_par) / (B0*B0);
Omega_i = (qi * B0) / (mi * c);
omega_pi = sqrt( (4*pi*ni*qi*qi) /( mi ));
lambda_Di = sqrt( (kb * Ti_per) / (4*pi*ni*qi*qi) ) ;
rho_i = Vth_i_per / Omega_i ;

%Va = B0/sqrt(4*pi*ni*mi);

% Electrons
ne=1;
Te=1;
Te_per=1;
Te_par_per = 1; 
Te_par = Te_per / Te_par_per; 
Vth_e_per = sqrt(2 * Te_per / me);
Vth_e_par = sqrt(2 * Te_par / me);
neni = ne/ni;
memi = 1/mime;
beta_e= (8*pi*ne*kb*Te) / (B0*B0);
beta_e_per= (8*pi*ne*kb*Te_per) / (B0*B0);
beta_e_par= (8*pi*ne*kb*Te_par) / (B0*B0);
Omega_e = (qe * B0) / (me * c) ;
omega_pe = sqrt( (4*pi*ne*qe*qe) /( me ));
lambda_De = sqrt( (kb * Te) / (4*pi*ne*qe*qe) ) ;
rho_e = Vth_e_per / Omega_e ;
d_i=rho_i ./ sqrt(beta_i);

%Ion acustic scaling
Vs = sqrt(Te / mi); 
Rs = Vs / Omega_i;
VsVa = Vs./Va; 
%-------------------------------------------------------------------

%Here is where I am going to include the input for k_per and k_par
%-------------------------------------------------------------------
k_per = 1;

k_par = 1;
%-------------------------------------------------------------------
% Using the notation of Daniel's Multiscale solar wind paper 2019
%-------------------------------------------------------------------
% For Alfven waves in the inertial range. (eq. 163)

w_aw = k_par * Va;
v_aw_phase = w_aw / k_par ;
tl_aw = 1 / (k_par * v_aw_phase);   
%-------------------------------------------------------------------

%-------------------------------------------------------------------
% For Kinetic Alfven waves (eq. 165)
TeTi = Te./Ti;
w_kaw = ( k_par * Va * k_per * rho_i ) / sqrt( beta_i + ( 2 / ( 1 + ( TeTi ) ) ));
v_kaw_phase = w_kaw / k_par;
k = sqrt( k_per.*k_per + k_par.*k_par);
tl_kaw = 1 / (k * v_kaw_phase);

% For Kinetic Alfven waves (eq. 41) From the Boldyrev et al., 2013
TiTe = Ti./Te;
w_kaw_b = ( Va * Rs * k_par * k_per  ) * sqrt( (1 + TiTe) / ( 1  + (VsVa.^2)*( 1 + TeTi) ));
v_kaw_phase_b = w_kaw_b / k_par;
tl_kaw_b = 1 / (k * v_kaw_phase_b);

%-------------------------------------------------------------------

%-------------------------------------------------------------------
% For Alfve/Ion-Cyclotron waves (eq. 168)

w_IC = ( (Omega_i * k.^2 * d_i.^2)/2 )* sqrt( 1 + 4/(k.^2 * d_i.^2) - 1 ); 
v_IC_phase = w_IC / k_par;
tl_IC = 1 / (k * v_kaw_phase);
%-------------------------------------------------------------------


%-------------------------------------------------------------------
% For Ion-acustic waves (eq. 165)

w_IA = ( k_par * Va * k_per * rho_i ) / sqrt( beta_i + ( 2 / ( 1 + ( Te/Ti ) ) ));
v_IA_phase = w_IA / k_par;
tl_IA = 1 / (k * v_kaw_phase);
%-------------------------------------------------------------------