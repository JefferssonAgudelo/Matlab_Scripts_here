%--------------------------------------------------------------------------
%This script is to compute the energy budget
% Most of the plotss are in the Computing_energy_reconnection2.m file

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data'; %This is important because the xdmf files are in that directory
path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

% useful parameters used in the simulation
%--------------------------------------------------------------------------
kb = 1.; mu0 =1.; mi_over_me_ =100.; vA_over_c_ = 0.1; B0_ = vA_over_c_;
mi = 1.; me = 1. / mi_over_me_; eps0=1.;
omega_pi=1;
Omega_e=vA_over_c_*sqrt(mi_over_me_)*omega_pi;
Omega_i=Omega_e / mi_over_me_;
TOmega_e=1/Omega_e; 
TOmega_i=1/Omega_i;
dt=0.06;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use <---------------- This is the N 
H=zeros(N);
%in=30;
%i=2;
%--------------------------------------------------------------------------

% To save some variables in a file for all the times. Uncomment this when
% running again the program to generate the data
%--------------------------------------------------------------------------
%{
Vectorsc=cell(6,N);
Tensorsc=cell(4,N);
Densitiesc=cell(2,N);
Poyn_theoc=cell(3,N);
kin_energy_terms_ic=cell(4,N);
kin_energy_terms_ec=cell(4,N);
Int_energy_terms_ic=cell(7,N);
Int_energy_terms_ec=cell(7,N);
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
nx1=166;    nx2=333;
ny1=232;    ny2=400;
nz1=1200;   nz2=1367;

px2 = 168; py2 = 168; pz2 = 168;
start = [nx1 ny1 nz1];
count = [px2 py2 pz2];
%--------------------------------------------------------------------------

resx=0.06; %0.06081
resy=0.06;
resz=0.06; %0.06152

Lx=px2*resx; k0x=1/Lx;
Ly=px2*resy; k0y=1/Ly;
Lz=px2*resz; k0z=1/Lz;
k00=sqrt(k0x*k0x + k0y*k0y + k0z*k0z);
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the whole core of the analysis and the data generation !! 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for in=1:N
%--------------------------------------------------------------------------
disp(strcat('Computing step ...',S(in).name))
%--------------------------------------------------------------------------
%N_steps = 2000;
%time_t = dt*N_steps; %(In terms of 1/omega_pi)
str = string(S(in).name);
newStr = extractBetween(str,"pfd.","_p000000.h5");
N_steps = str2double(newStr);
time_t = dt*N_steps;
%--------------------------------------------------------------------------

fileID =  S(in).name; %change the 1 per i
%h5disp(fileID);
info = h5info(fileID);
T_1st = info.Groups(1).Name;
E_1st = info.Groups(17).Name;
B_1st = info.Groups(18).Name;
J_1st = info.Groups(19).Name;
V_1st = info.Groups(25).Name;
n_1st = info.Groups(23).Name;
%--------------------------------------------------------------------------


%------------------------------------------------------------------
B{1,1} =h5read(fileID,strcat(B_1st,'/hx/p0/3d'),start,count); %Bx=Bx(nx1:nx2,ny1:ny2,nz1:nz2);
B{1,2} =h5read(fileID,strcat(B_1st,'/hy/p0/3d'),start,count); %By=By(nx1:nx2,ny1:ny2,nz1:nz2);
B{1,3} =h5read(fileID,strcat(B_1st,'/hz/p0/3d'),start,count); %Bz=Bz(nx1:nx2,ny1:ny2,nz1:nz2);
E{1,1} =h5read(fileID,strcat(E_1st,'/ex/p0/3d'),start,count); %Ex=Ex(nx1:nx2,ny1:ny2,nz1:nz2);
E{1,2} =h5read(fileID,strcat(E_1st,'/ey/p0/3d'),start,count); %Ey=Ey(nx1:nx2,ny1:ny2,nz1:nz2);
E{1,3} =h5read(fileID,strcat(E_1st,'/ez/p0/3d'),start,count); %Ez=Ez(nx1:nx2,ny1:ny2,nz1:nz2);
J{1,1} =h5read(fileID,strcat(J_1st,'/jx/p0/3d'),start,count); %Jx=Jx(nx1:nx2,ny1:ny2,nz1:nz2);
J{1,2} =h5read(fileID,strcat(J_1st,'/jy/p0/3d'),start,count); %Jy=Jy(nx1:nx2,ny1:ny2,nz1:nz2);
J{1,3} =h5read(fileID,strcat(J_1st,'/jz/p0/3d'),start,count); %Jz=Jz(nx1:nx2,ny1:ny2,nz1:nz2);
vi{1,1} =h5read(fileID,strcat(V_1st,'/vx_i/p0/3d'),start,count); %vix=vix(nx1:nx2,ny1:ny2,nz1:nz2);
vi{1,2} =h5read(fileID,strcat(V_1st,'/vy_i/p0/3d'),start,count); %viy=viy(nx1:nx2,ny1:ny2,nz1:nz2);
vi{1,3} =h5read(fileID,strcat(V_1st,'/vz_i/p0/3d'),start,count); %viz=viz(nx1:nx2,ny1:ny2,nz1:nz2);
ve{1,1} =h5read(fileID,strcat(V_1st,'/vx_e/p0/3d'),start,count); %vex=vex(nx1:nx2,ny1:ny2,nz1:nz2);
ve{1,2} =h5read(fileID,strcat(V_1st,'/vy_e/p0/3d'),start,count); %vey=vey(nx1:nx2,ny1:ny2,nz1:nz2);
ve{1,3} =h5read(fileID,strcat(V_1st,'/vz_e/p0/3d'),start,count); %vez=vez(nx1:nx2,ny1:ny2,nz1:nz2);
nic1{1,1} = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'),start,count); %ni=ni(nx1:nx2,ny1:ny2,nz1:nz2);
nec1{1,1} = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'),start,count); %ne=ne(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_i{1,1} =h5read(fileID,strcat(T_1st,'/Txx_i/p0/3d'),start,count); %Pxxi=Pxxi(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_i{1,2} =h5read(fileID,strcat(T_1st,'/Txy_i/p0/3d'),start,count); %Pxyi=Pxyi(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_i{1,3} =h5read(fileID,strcat(T_1st,'/Txz_i/p0/3d'),start,count); %Pxzi=Pxzi(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_i{2,2} =h5read(fileID,strcat(T_1st,'/Tyy_i/p0/3d'),start,count); %Pyyi=Pyyi(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_i{2,3} =h5read(fileID,strcat(T_1st,'/Tyz_i/p0/3d'),start,count); %Pyzi=Pyzi(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_i{3,3} =h5read(fileID,strcat(T_1st,'/Tzz_i/p0/3d'),start,count); %Pzzi=Pzzi(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_e{1,1} =h5read(fileID,strcat(T_1st,'/Txx_e/p0/3d'),start,count); %Pxxe=Pxxe(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_e{1,2} =h5read(fileID,strcat(T_1st,'/Txy_e/p0/3d'),start,count); %Pxye=Pxye(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_e{1,3} =h5read(fileID,strcat(T_1st,'/Txz_e/p0/3d'),start,count); %Pxze=Pxze(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_e{2,2} =h5read(fileID,strcat(T_1st,'/Tyy_e/p0/3d'),start,count); %Pyye=Pyye(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_e{2,3} =h5read(fileID,strcat(T_1st,'/Tyz_e/p0/3d'),start,count); %Pyze=Pyze(nx1:nx2,ny1:ny2,nz1:nz2);
Pij_e{3,3} =h5read(fileID,strcat(T_1st,'/Tzz_e/p0/3d'),start,count); %Pzze=Pzze(nx1:nx2,ny1:ny2,nz1:nz2);
%------------------------------------------------------------------

%------------------------------------------------------------------
%vectors as 1 arrow 3 columns
%mi=1;
mime=100;
%me=mi/mime;
qi=1;
qe=-1;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Now lets do the reference frame transformation
%--------------------------------------------------------------------------
% This is to get the quantities in the slide that is related to the correct plane
% or it seems that way.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
xl=linspace(nx1,nx1+px2,px2)*resx;
yl=linspace(ny1,ny1+py2,py2)*resy;
zl=linspace(nz1,nz1+pz2,pz2)*resz;
%Pick the slide number along z-axes to get the plots
slide_n =83;
%--------------------------------------------------------------------------
[xm,ym] = meshgrid(xl,yl);
[xm3,ym3, zm3] = meshgrid(xl,yl,zl);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Now it's time to define the values for the slides coordinates
x_pi=60; y_pi=80; z_pi = 46; radi=100;
%y_pi=60;

% Extract the slide. This is to define the plane using a point a a normal
% vector

% To define the plane lets do the following: 
% 1) Pick the major axis along the current structure (i1p)
% 2) Define the plane perpendicular to the major axis. To do this lets pick
% a generic perpendicular vector 

% This is the point
pt = [x_pi y_pi z_pi]; %pt = [0 0 0];

% This is the normal vector to the plane perpendicular to the current (Vec adjacent)
vec_a = [1 0 -3.3333];
%vec_a = [-3.3333 0 1]; %This doesn't work as it is. There is a problem.
%See comment below
mag_va = sqrt(vec_a(1)*vec_a(1) + vec_a(2)*vec_a(2) + vec_a(3)*vec_a(3));
% This is the normal vector to the plane along the current (vec perpendicular) 
vec_p = [vec_a(1) vec_a(2) -(vec_a(1) + vec_a(2))/vec_a(3)]; 
mag_vp = sqrt(vec_p(1)*vec_p(1) + vec_p(2)*vec_p(2) + vec_p(3)*vec_p(3));
% This is the vector from the right hand rule
vec_rh = cross(vec_a, vec_p);
mag_vrh = sqrt(vec_rh(1)*vec_rh(1) + vec_rh(2)*vec_rh(2) + vec_rh(3)*vec_rh(3));
uip = vec_rh/mag_vrh; ujp = vec_a/mag_va; ukp = vec_p/mag_vp;

%uip=[1 0 0]; ujp=[0 1 0]; ukp=[0 0 1];
%--------------------------------------------------------------------------
%To get the slides using the normal vector along the current filament
vec=vec_a; % <----------
%To get the slides using the normal vector perpendicular direction (z'')
%vec=vec_p;

% ------> NOTE: 
% When I use uip=[1 0 0]; ujp=[0 1 0]; ukp=[0 0 1]; and vec=ujp or vec=ukp
% it doesnt work =[] lots of nans
%vec=ukp;
%--------------------------------------------------------------------------
  
% Now lets compute the quantities in the new reference frame  (rh,a,p) 
%--------------------------------------------------------------------------
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
% For vectors apply the function directly
%--------------------------------------------------------------------------
% Magnetic field in the prime reference frame
Qb = cell(1,3);
Qb{1,1} = B{1,1}; Qb{1,2} = B{1,2}; Qb{1,3} = B{1,3};% -0.1; % it is very important to remember this substraction
B = Ref_Frame_2(uip,ujp,ukp,Qb); clearvars Qb;
% Electric field
Qb=E;
E = Ref_Frame_2(uip,ujp,ukp,Qb); clearvars Qb
% Current
Qb=J;
J = Ref_Frame_2(uip,ujp,ukp,Qb); clearvars Qb
% Ion velocity
Qb=vi;
vi = Ref_Frame_2(uip,ujp,ukp,Qb); clearvars Qb
% Electron velocity
Qb=ve;
ve = Ref_Frame_2(uip,ujp,ukp,Qb); clearvars Qb
%--------------------------------------------------------------------------

% The densities which are scalars are invariant uner coordinates transformations. 

%--------------------------------------------------------------------------
%Now the proper transformation of the tensor.
%--------------------------------------------------------------------------
%For ions
GPij=cell(3,3);
%--------------------------------------------------------------------------
% Make the tensor [3, 3] for the ions
GPij{1,1}=Pij_i{1,1}; GPij{1,2}=Pij_i{1,2}; GPij{1,3}=Pij_i{1,3};
GPij{2,1}=Pij_i{1,2}; GPij{2,2}=Pij_i{2,2}; GPij{2,3}=Pij_i{2,3};
GPij{3,1}=Pij_i{1,3}; GPij{3,2}=Pij_i{2,3}; GPij{3,3}=Pij_i{3,3};
Qb=GPij;
Pij_i = Ref_Frame_2(uip,ujp,ukp,Qb); clearvars Qb
% Make the tensor [3, 3] for the electrons
GPij{1,1}=Pij_e{1,1}; GPij{1,2}=Pij_e{1,2}; GPij{1,3}=Pij_e{1,3};
GPij{2,1}=Pij_e{1,2}; GPij{2,2}=Pij_e{2,2}; GPij{2,3}=Pij_e{2,3};
GPij{3,1}=Pij_e{1,3}; GPij{3,2}=Pij_e{2,3}; GPij{3,3}=Pij_e{3,3};
Qb=GPij;
Pij_e = Ref_Frame_2(uip,ujp,ukp,Qb); clearvars Qb
%--------------------------------------------------------------------------



%In the new reference frame the (rh,a,p) calculate the energy terms
%--------------------------------------------------------------------------
% This is to unlock. These are the lines to calculate the energy terms
%{\  
% Get again the densities right
ni=nic1{1,1};
ne=nec1{1,1};
% Kinetic energy densities 
%------------------------------------------------------------------
Ki = 0.5.*ni.* mi.*(vi{1,1}.*vi{1,1} + vi{1,2}.*vi{1,2} + vi{1,3}.*vi{1,3});
Ke = 0.5.*ne.* me.*(ve{1,1}.*ve{1,1} + ve{1,2}.*ve{1,2} + ve{1,3}.*vi{1,3});
%------------------------------------------------------------------

% Internal energy densities from Fadanelli et al, 2020 
%------------------------------------------------------------------
Ui = 0.5.*(Pij_i{1,1} + Pij_i{2,2} + Pij_i{3,3});
Ue = 0.5.*(Pij_e{1,1} + Pij_e{2,2} + Pij_e{3,3});
%------------------------------------------------------------------

%For the partial time derivative
%------------------------------------------------------------------
%{
Ki_1 = Ki;
Ke_1 = Ke;
Ui_1 = Ui;
Ue_1 = Ue;
%clearvars -except Ki_0 Ke_0 Ui_0 Ue_0 
partialKi = (Ki_1 - Ki_0)/6;
partialKe = (Ke_1 - Ke_0)/6;
partialUi = (Ui_1 - Ui_0)/6;
partialUe = (Ue_1 - Ue_0)/6;
%}
%------------------------------------------------------------------

%Derivates%
%------------------------------------------------------------------
%Grad_B
[Bx_x,Bx_y,Bx_z] = gradient(B{1,1});
[By_x,By_y,By_z] = gradient(B{1,2});
[Bz_x,Bz_y,Bz_z] = gradient(B{1,3});
grad_B{1,1} = Bx_x/resx; grad_B{1,2} = By_x/resx; grad_B{1,3} = Bz_x/resx;
grad_B{2,1} = Bx_y/resy; grad_B{2,2} = By_y/resy; grad_B{2,3} = Bz_y/resy;
grad_B{3,1} = Bx_z/resz; grad_B{3,2} = By_z/resz; grad_B{3,3} = Bz_z/resz;

clearvars Bx_x Bx_y Bx_z By_x By_y By_z Bz_x Bz_y Bz_z 
Div_B=grad_B{1,1} + grad_B{2,2} + grad_B{3,3};

%Grad_vi
[vix_x,vix_y,vix_z] = gradient(vi{1,1});
[viy_x,viy_y,viy_z] = gradient(vi{1,2});
[viz_x,viz_y,viz_z] = gradient(vi{1,3});
grad_vi{1,1} = vix_x/resx; grad_vi{1,2} = viy_x/resx; grad_vi{1,3} = viz_x/resx;
grad_vi{2,1} = vix_y/resy; grad_vi{2,2} = viy_y/resy; grad_vi{2,3} = viz_y/resy;
grad_vi{3,1} = vix_z/resz; grad_vi{3,2} = viy_z/resz; grad_vi{3,3} = viz_z/resz;

clearvars vix_x vix_y vix_z viy_x viy_y viy_z viz_x viz_y viz_z 
Div_vi=grad_vi{1,1} + grad_vi{2,2} + grad_vi{3,3};

%Grad_ve
[vex_x,vex_y,vex_z] = gradient(ve{1,1});
[vey_x,vey_y,vey_z] = gradient(ve{1,2});
[vez_x,vez_y,vez_z] = gradient(ve{1,3});
grad_ve{1,1} = vex_x/resx; grad_ve{1,2} = vey_x/resx; grad_ve{1,3} = vez_x/resx;
grad_ve{2,1} = vex_y/resy; grad_ve{2,2} = vey_y/resy; grad_ve{2,3} = vez_y/resy;
grad_ve{3,1} = vex_z/resz; grad_ve{3,2} = vey_z/resz; grad_ve{3,3} = vez_z/resz;

clearvars vex_x vex_y vex_z vey_x vey_y vey_z vez_x vez_y vez_z 
Div_ve=grad_ve{1,1} + grad_ve{2,2} + grad_ve{3,3};

%Grad_Pi
[Pxxi_x, Pxxi_y, Pxxi_z] = gradient(Pij_i{1,1}); 
[Pyyi_x, Pyyi_y, Pyyi_z] = gradient(Pij_i{2,2}); 
[Pzzi_x, Pzzi_y, Pzzi_z] = gradient(Pij_i{3,3}); 
[Pxyi_x, Pxyi_y, Pxyi_z] = gradient(Pij_i{1,2}); 
[Pxzi_x, Pxzi_y, Pxzi_z] = gradient(Pij_i{1,3}); 
[Pyzi_x, Pyzi_y, Pyzi_z] = gradient(Pij_i{2,3});
Div_Pi_t{1,1} = Pxxi_x/resx + Pxyi_y/resy + Pxzi_z/resz;  
Div_Pi_t{1,2} = Pxyi_x/resx + Pyyi_y/resy + Pyzi_z/resz; 
Div_Pi_t{1,3} = Pxzi_x/resx + Pyzi_y/resy + Pzzi_z/resz;
clearvars Pxxi_y Pxxi_z Pyyi_x Pyyi_z Pzzi_x Pzzi_z Pxyi_z Pxzi_y Pyzi_x
clearvars Pxxi_x Pxyi_x Pxyi_y Pyyi_y Pyzi_y Pyzi_z Pxzi_x Pxzi_z Pzzi_z

%Grad_Pe
[Pxxe_x, Pxxe_y, Pxxe_z] = gradient(Pij_e{1,1}); 
[Pyye_x, Pyye_y, Pyye_z] = gradient(Pij_e{2,2}); 
[Pzze_x, Pzze_y, Pzze_z] = gradient(Pij_e{3,3}); 
[Pxye_x, Pxye_y, Pxye_z] = gradient(Pij_e{1,2}); 
[Pxze_x, Pxze_y, Pxze_z] = gradient(Pij_e{1,3}); 
[Pyze_x, Pyze_y, Pyze_z] = gradient(Pij_e{2,3});
Div_Pe_t{1,1} = Pxxe_x/resx + Pxye_y/resy + Pxze_z/resz;  
Div_Pe_t{1,2} = Pxye_x/resx + Pyye_y/resy + Pyze_z/resz; 
Div_Pe_t{1,3} = Pxze_x/resx + Pyze_y/resy + Pzze_z/resz;
clearvars Pxxe_y Pxxe_z Pyye_x Pyye_z Pzze_x Pzze_z Pxye_z Pxze_y Pyze_x
clearvars Pxxe_x Pxye_x Pxye_y Pyye_y Pyze_y Pyze_z Pxze_x Pxze_z Pzze_z

%Grad_Ki
[Ki_x,Ki_y,Ki_z] = gradient(Ki);
grad_Ki{1,1} = Ki_x/resx; 
grad_Ki{1,2} = Ki_y/resy; 
grad_Ki{1,3} = Ki_z/resz;
clearvars Ki_x Ki_y Ki_z

%Grad_Ki
[Ke_x,Ke_y,Ke_z] = gradient(Ke);
grad_Ke{1,1} = Ke_x/resx; 
grad_Ke{1,2} = Ke_y/resy; 
grad_Ke{1,3} = Ke_z/resz;
clearvars Ke_x Ke_y Ke_z

%Grad_Ui
[Ui_x,Ui_y,Ui_z] = gradient(Ui);
grad_Ui{1,1} = Ui_x/resx; 
grad_Ui{1,2} = Ui_y/resy; 
grad_Ui{1,3} = Ui_z/resz;
clearvars Ui_x Ui_y Ui_z

%Grad_Ue
[Ue_x,Ue_y,Ue_z] = gradient(Ue);
grad_Ue{1,1} = Ue_x/resx; 
grad_Ue{1,2} = Ue_y/resy; 
grad_Ue{1,3} = Ue_z/resz;
clearvars Ue_x Ue_y Ue_z
%------------------------------------------------------------------


% Dyadic products
%------------------------------------------------------------------
%vivi and %veve
vivi=cell(3,3);
veve=cell(3,3);
for i=1:3
    for j=1:3
        vivi{i,j} = vi{1,i}.*vi{1,j};
        veve{i,j} = ve{1,i}.*ve{1,j};
    end
end

%Pivi
Pit_vi{1,1} = Pij_i{1,1}.*vi{1,1} + Pij_i{1,2}.*vi{1,2} + Pij_i{1,3}.*vi{1,3};
Pit_vi{1,2} = Pij_i{1,2}.*vi{1,1} + Pij_i{2,2}.*vi{1,2} + Pij_i{2,3}.*vi{1,3};
Pit_vi{1,3} = Pij_i{1,3}.*vi{1,1} + Pij_i{2,3}.*vi{1,2} + Pij_i{3,3}.*vi{1,3};

%Peve
Pet_ve{1,1} = Pij_e{1,1}.*ve{1,1} + Pij_e{1,2}.*ve{1,2} + Pij_e{1,3}.*ve{1,3};
Pet_ve{1,2} = Pij_e{1,2}.*ve{1,1} + Pij_e{2,2}.*ve{1,2} + Pij_e{2,3}.*ve{1,3};
Pet_ve{1,3} = Pij_e{1,3}.*ve{1,1} + Pij_e{2,3}.*ve{1,2} + Pij_e{3,3}.*ve{1,3};
%------------------------------------------------------------------


% Scalar products
%------------------------------------------------------------------
% Pit_grad_vi. this the the same double dot procut
Pit_grad_vi = Pij_i{1,1}.*grad_vi{1,1} + Pij_i{1,2}.*grad_vi{1,2} + Pij_i{1,3}.*grad_vi{1,3} +...
              Pij_i{1,2}.*grad_vi{2,1} + Pij_i{2,2}.*grad_vi{2,2} + Pij_i{2,3}.*grad_vi{2,3} +...
              Pij_i{1,3}.*grad_vi{3,1} + Pij_i{2,3}.*grad_vi{3,2} + Pij_i{3,3}.*grad_vi{3,3};

% Pet_grad_ve
Pet_grad_ve = Pij_e{1,1}.*grad_ve{1,1} + Pij_e{1,2}.*grad_ve{1,2} + Pij_e{1,3}.*grad_ve{1,3} +...
              Pij_e{1,2}.*grad_ve{2,1} + Pij_e{2,2}.*grad_ve{2,2} + Pij_e{2,3}.*grad_ve{2,3} +...
              Pij_e{1,3}.*grad_ve{3,1} + Pij_e{2,3}.*grad_ve{3,2} + Pij_e{3,3}.*grad_ve{3,3};

%--------------------------------------------------------------------------          
% Diagonal          
% Pit_grad_vi_ii
Pit_grad_vi_ii = Pij_i{1,1}.*grad_vi{1,1} + Pij_i{2,2}.*grad_vi{2,2} + Pij_i{3,3}.*grad_vi{3,3};

% Pet_grad_ve_ii
Pet_grad_ve_ii = Pij_e{1,1}.*grad_ve{1,1} + Pij_e{2,2}.*grad_ve{2,2} + Pij_e{3,3}.*grad_ve{3,3};

% Off diagonal
% Pit_grad_vi_ij
Pit_grad_vi_ij = Pij_i{1,2}.*grad_vi{1,2} + Pij_i{1,3}.*grad_vi{1,3} + Pij_i{1,2}.*grad_vi{2,1} +...
                 Pij_i{2,3}.*grad_vi{2,3} + Pij_i{1,3}.*grad_vi{3,1} + Pij_i{2,3}.*grad_vi{3,2};

% Pet_grad_ve_ij
Pet_grad_ve_ij =  Pij_e{1,2}.*grad_ve{1,2} + Pij_e{1,3}.*grad_ve{1,3} + Pij_e{1,2}.*grad_ve{2,1} +...
                  Pij_e{2,3}.*grad_ve{2,3} + Pij_e{1,3}.*grad_ve{3,1} + Pij_e{2,3}.*grad_ve{3,2};
%--------------------------------------------------------------------------          
              
% vi_Div_Pit
vi_Div_Pit = vi{1,1}.*Div_Pi_t{1,1} + vi{1,2}.*Div_Pi_t{1,2} + vi{1,3}.*Div_Pi_t{1,3}; 

% ve_Div_Pet
ve_Div_Pet = ve{1,1}.*Div_Pe_t{1,1} + ve{1,2}.*Div_Pe_t{1,2} + ve{1,3}.*Div_Pe_t{1,3}; 

%--------------------------------------------------------------------------
% Div_Pi_dot_vi
Div_Pis_dot_vi = Pit_grad_vi + vi_Div_Pit;
% Div_Pe_dot_e
Div_Pes_dot_ve = Pet_grad_ve + ve_Div_Pet;
%--------------------------------------------------------------------------

% Ki_Div_vi
Ki_Div_vi = Ki.*Div_vi; 

% Ke_Div_ve
Ke_Div_ve = Ke.*Div_ve; 

% Ui_Div_vi
Ui_Div_vi = Ui.*Div_vi; 

% Ke_Div_ve
Ue_Div_ve = Ue.*Div_ve; 

%qi_ni_vi_E
qi_ni_vi_E = qi.*ni.*(vi{1,1}.*E{1,1} + vi{1,2}.*E{1,2} + vi{1,3}.*E{1,3});

%qe_ne_ve_E
qe_ne_ve_E = qe.*ne.*(ve{1,1}.*E{1,1} + ve{1,2}.*E{1,2} + ve{1,1}.*E{1,3});

%vi_grad_Ki
vi_grad_Ki = vi{1,1}.*grad_Ki{1,1} + vi{1,2}.*grad_Ki{1,2} + vi{1,3}.*grad_Ki{1,3}; 

%ve_grad_Ke
ve_grad_Ke = ve{1,1}.*grad_Ke{1,1} + ve{1,2}.*grad_Ke{1,2} + ve{1,3}.*grad_Ke{1,3}; 

%vi_grad_Ui
vi_grad_Ui = vi{1,1}.*grad_Ui{1,1} + vi{1,2}.*grad_Ui{1,2} + vi{1,3}.*grad_Ui{1,3}; 

%ve_grad_Ue
ve_grad_Ue = ve{1,1}.*grad_Ue{1,1} + ve{1,2}.*grad_Ue{1,2} + ve{1,3}.*grad_Ue{1,3}; 
%------------------------------------------------------------------


% From the definition in the code I think Tij is actually Pij as
% Pij = m * int f(vi - ui)(vj - uj) dv^3 and not int( fvivj )dv^3 which is
% PPij = Pij + nmuiuj 
% (HERE PPij is for the total energy and Pij is the thermal, the other way aroaund from the notes)

% The internal energy can be compute using (1/2)Pij_i/ni?
%------------------------------------------------------------------
PPij_i=cell(3,3);
PPij_e=cell(3,3);
for i=1:3
    for j=1:3
        if(isempty(Pij_i{i,j}) == 1)
            continue
        end
        PPij_i{i,j} = Pij_i{i,j} + ni.*mi.*vi{1,i}.*vi{1,j};
        PPij_e{i,j} = Pij_e{i,j} + ne.*me.*ve{1,i}.*ve{1,j};
    end
end
%------------------------------------------------------------------


% Calculating the temperatures
%------------------------------------------------------------------
Tij_i=cell(3,3);
Tij_e=cell(3,3);
for i=1:3
    for j=1:3
        if(isempty(Pij_i{i,j}) == 1)
            continue
        end
        Tij_i{i,j} =Pij_i{i,j}./(3.*ni);
        Tij_e{i,j} =Pij_e{i,j}./(3.*ne);
    end
end
%------------------------------------------------------------------


%Lets compute the averages to see what is what
%------------------------------------------------------------------
ni_bar=mean(ni,'all');
ne_bar=mean(ne,'all');

PPij_i_bar=zeros(3);
PPij_e_bar=zeros(3);
Tij_i_bar=zeros(3);
Tij_e_bar=zeros(3);
for i=1:3
    for j=1:3
        PPij_i_bar(i,j)=mean(PPij_i{i,j}, 'all');
        PPij_e_bar(i,j)=mean(PPij_e{i,j}, 'all');
        Tij_i_bar(i,j)=mean(Tij_i{i,j}, 'all');
        Tij_e_bar(i,j)=mean(Tij_e{i,j}, 'all');
    end 
end

%{
Pij_i{1,1} = Pxxi; Pij_i{1,2} = Pxyi; Pij_i{1,3} = Pxzi; 
Pij_i{2,1} = Pxyi; Pij_i{2,2} = Pyyi; Pij_i{2,3} = Pyzi; 
Pij_i{3,1} = Pxzi; Pij_i{3,2} = Pyzi; Pij_i{3,3} = Pzzi; 
clearvars Pxxi Pxyi Pxzi Pyyi Pyzi Pzzi
Pij_e{1,1} = Pxxe; Pij_e{1,2} = Pxye; Pij_e{1,3} = Pxze; 
Pij_e{2,1} = Pxye; Pij_e{2,2} = Pyye; Pij_e{2,3} = Pyze; 
Pij_e{3,1} = Pxze; Pij_e{3,2} = Pyze; Pij_e{3,3} = Pzze; 
clearvars Pxxe Pxye Pxze Pyye Pyzi Pzze
%}
Pij_i_bar=zeros(3);
Pij_e_bar=zeros(3);
for i=1:3
    for j=1:3
        Pij_i_bar(i,j)=mean(Pij_i{i,j}, 'all');
        Pij_e_bar(i,j)=mean(Pij_e{i,j}, 'all');
    end 
end
%------------------------------------------------------------------


% To calculate the grad of the heat flux vector wang et al, 2015
%------------------------------------------------------------------
p_i = (Pij_i{1,1} + Pij_i{2,2} + Pij_i{3,3})./3;  
p_e = (Pij_e{1,1} + Pij_e{2,2} + Pij_e{3,3})./3;  

vth_i = sqrt((Tij_i{1,1} + Tij_i{2,2} + Tij_i{3,3})./(3*mi));
vth_e = sqrt((Tij_e{1,1} + Tij_e{2,2} + Tij_e{3,3})./(3*me));


%(HERE PPij is for the total energy and Qij is the (v-u)..., the other way aroaund from the notes)
%dkQijk_i{1,1} = vth_i.*k00.*( Pij_i{1,1} - Pij_i_bar(1,1) - ((ni - ni_bar).*Tij_i_bar(1,1)));
dk_Qijk_i=cell(3,3);
dk_Qijk_e=cell(3,3);
for i=1:3
    for j=1:3
        if(isempty(PPij_i{i,j}) == 1)
            continue
        end
        dk_Qijk_i{i,j} = vth_i.*k00.*( Pij_i{i,j} - Pij_i_bar(i,j) - ((i == j).*(ni - ni_bar).*Tij_i_bar(i,j)));
        dk_Qijk_e{i,j} = vth_e.*k00.*( Pij_e{i,j} - Pij_e_bar(i,j) - ((i == j).*(ne - ne_bar).*Tij_e_bar(i,j)));        
    end
end

%------------------------------------------------------------------

%Now lets compute the poynting terms
%------------------------------------------------------------------
mu0=1;
eps0=1;
% The poynting vector
Poyn_v=cell(1,3);
Poyn_v{1,1}=(1/mu0).*(E{1,2}.*B{1,3} - E{1,3}.*B{1,2}); 
Poyn_v{1,2}=(1/mu0).*(E{1,3}.*B{1,1} - E{1,1}.*B{1,3});
Poyn_v{1,3}=(1/mu0).*(E{1,1}.*B{1,2} - E{1,2}.*B{1,1});

% Now the divergence of the poyting vector
[Poyn_vx_x,Poyn_vx_y,Poyn_vx_z] = gradient(Poyn_v{1,1});
[Poyn_vy_x,Poyn_vy_y,Poyn_vy_z] = gradient(Poyn_v{1,2});
[Poyn_vz_x,Poyn_vz_y,Poyn_vz_z] = gradient(Poyn_v{1,3});
grad_Poyn_v{1,1} = Poyn_vx_x/resx; grad_Poyn_v{1,2} = Poyn_vy_x/resx; grad_Poyn_v{1,3} = Poyn_vz_x/resx;
grad_Poyn_v{2,1} = Poyn_vx_y/resy; grad_Poyn_v{2,2} = Poyn_vy_y/resy; grad_Poyn_v{2,3} = Poyn_vz_y/resy;
grad_Poyn_v{3,1} = Poyn_vx_z/resz; grad_Poyn_v{3,2} = Poyn_vy_z/resz; grad_Poyn_v{3,3} = Poyn_vz_z/resz;

clearvars Poyn_vx_x Poyn_vx_y Poyn_vx_z Poyn_vy_x Poyn_vy_y Poyn_vy_z Poyn_vz_x Poyn_vz_y Poyn_vz_z 
Div_Poyn_v=grad_Poyn_v{1,1} + grad_Poyn_v{2,2} + grad_Poyn_v{3,3};


% J dot E
JdotE = J{1,1}.*E{1,1} + J{1,2}.*E{1,2} + J{1,3}.*E{1,3};  

% Electromagnetic energy 
E2 = E{1,1}.^2 + E{1,2}.^2 + E{1,3}.^2; 
B2 = B{1,1}.^2 + B{1,2}.^2 + B{1,3}.^2; 
E_em = 0.5.*( eps0.*E2 + B2./mu0 );
%------------------------------------------------------------------

% Now lets compute the time changes of the energy quantities
%------------------------------------------------------------------
% Time change of the electromagnetic energy 
dt_E_em = - JdotE - Div_Poyn_v;  

% Time change of the electromagnetic energy 
%------------------------------------------------------------------
%}
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

% From here onwards is where using the other vec vector doesn't work. The following is the error 
%---------------------------------------------
%Index in position 2 exceeds array bounds.
%Error in Vars_on_Slice (line 21)
%X(isnan(X(:,1)),:) = []; X=X';
%---------------------------------------------

%--------------------------------------------------------------------------
% calculate the quantities on the slides for scalars
%--------------------------------------------------------------------------
nic{1,1}=ni; 
nec{1,1}=ne;
Gp=nic;
[ni_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Gp=nec;
[ne_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
%--------------------------------------------------------------------------
% Scalar terms from the energy analysis
%--------------------------------------------------------------------------
% Poything theorem terms
%--------------------------------------------------------------------------
Gp=Poyn_v;
[Poyn_v_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

clearvars JdotEc Div_Poyn_vc dt_E_emc
JdotEc{1,1}=JdotE; 
Gp=JdotEc;
[JdotE_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Div_Poyn_vc{1,1}=Div_Poyn_v;
Gp=Div_Poyn_vc;
[Div_Poyn_v_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
dt_E_emc{1,1}=dt_E_em;
Gp=dt_E_emc;
[dt_E_em_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

clearvars Gdum
Gdum{1,1}=Ki_Div_vi./Ki_Div_vi;%dt_E_em./dt_E_em; %dt_E_em+JdotE+Div_Poyn_v; %Ki_Div_vi + vi_Div_Pit
clearvars Gdum_Gp;
Gp=Gdum;
[Gdum_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP


% Kinetic energy equation Fadanelli
%--------------------------------------------------------------------------
Ki_Div_vic{1,1}=Ki_Div_vi;
Gp=Ki_Div_vic;
[Ki_Div_vi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Ke_Div_vec{1,1}=Ke_Div_ve;
Gp=Ke_Div_vec;
[Ke_Div_ve_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

vi_Div_Pitc{1,1}=vi_Div_Pit; 
Gp=vi_Div_Pitc;
[vi_Div_Pit_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
ve_Div_Petc{1,1}=ve_Div_Pet;
Gp=ve_Div_Petc;
[ve_Div_Pet_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

qi_ni_vi_Ec{1,1}=qi_ni_vi_E;
Gp=qi_ni_vi_Ec;
[qi_ni_vi_E_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
qe_ne_ve_Ec{1,1}=qe_ne_ve_E;
Gp=qe_ne_ve_Ec;
[qe_ne_ve_E_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

%{
partialKic{1,1}=partialKi;
Gp=partialKic;
[partialKi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
partialKec{1,1}=partialKe;
Gp=partialKec;
[partialKe_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
%}

vi_grad_Kic{1,1}=vi_grad_Ki;
Gp=vi_grad_Kic;
[vi_grad_Ki_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
ve_grad_Kec{1,1}=ve_grad_Ke;
Gp=ve_grad_Kec;
[ve_grad_Ke_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

% Potential energy equation Fadanelli
%--------------------------------------------------------------------------
Ui_Div_vic{1,1}=Ui_Div_vi;
Gp=Ui_Div_vic;
[Ui_Div_vi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Ue_Div_vec{1,1}=Ue_Div_ve;
Gp=Ue_Div_vec;
[Ue_Div_ve_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

Pit_grad_vic{1,1}=Pit_grad_vi;
Gp=Pit_grad_vic;
[Pit_grad_vi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Pet_grad_vec{1,1}=Pet_grad_ve;
Gp=Pet_grad_vec;
[Pet_grad_ve_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

% Diagonal
Pit_grad_vi_iic{1,1}=Pit_grad_vi_ii;
Gp=Pit_grad_vi_iic;
[Pit_grad_vi_ii_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Pet_grad_ve_iic{1,1}=Pet_grad_ve_ii;
Gp=Pet_grad_ve_iic;
[Pet_grad_ve_ii_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

%Off diagonal
Pit_grad_vi_ijc{1,1}=Pit_grad_vi_ij;
Gp=Pit_grad_vi_ijc;
[Pit_grad_vi_ij_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Pet_grad_ve_ijc{1,1}=Pet_grad_ve_ij;
Gp=Pet_grad_ve_ijc;
[Pet_grad_ve_ij_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

% Div_Pi_dot_vi_Gp
Div_Pis_dot_vic{1,1}=Div_Pis_dot_vi;
Gp=Div_Pis_dot_vic;
[Div_Pis_dot_vi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
% Div_Pe_dot_ve_Gp
Div_Pes_dot_vec{1,1}=Div_Pes_dot_ve;
Gp=Div_Pes_dot_vec;
[Div_Pes_dot_ve_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

%{
partialUic{1,1}=partialUi;
Gp=partialUic;
[partialUi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
partialUec{1,1}=partialUe;
Gp=partialUec;
[partialUe_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
%}

vi_grad_Uic{1,1}=vi_grad_Ui;
Gp=vi_grad_Uic;
[vi_grad_Ui_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
ve_grad_Uec{1,1}=ve_grad_Ue;
Gp=ve_grad_Uec;
[ve_grad_Ue_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
%--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
% calculate the quantities on the slides for vectors
%--------------------------------------------------------------------------
Gp=B;
[B_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Gp=E;
[E_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Gp=J;
[J_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Gp=vi;
[vi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Gp=ve;
[ve_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
%--------------------------------------------------------------------------
%
% calculate the quantities on the slides for tensors
%--------------------------------------------------------------------------
Gp=Pij_i;
[Pij_i_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
Gp=Pij_e;
[Pij_e_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

Gp=dk_Qijk_i;
[dk_Qijk_i_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

Gp=dk_Qijk_e;
[dk_Qijk_e_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

Gp=grad_vi;
[grad_vi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

Gp=grad_ve;
[grad_ve_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
%--------------------------------------------------------------------------------------

% Total derivatives Fadanelli
%--------------------------------------------------------------------------------------
dkQiik_05_i_Gp{1,1} = (dk_Qijk_i_Gp{1,1} + dk_Qijk_i_Gp{2,2} + dk_Qijk_i_Gp{3,3})/2; 
dkQiik_05_e_Gp{1,1} = (dk_Qijk_e_Gp{1,1} + dk_Qijk_e_Gp{2,2} + dk_Qijk_e_Gp{3,3})/2;

dtUi_Gp{1,1} = - Ui_Div_vi_Gp{1,1} - Pit_grad_vi_Gp{1,1} - dkQiik_05_i_Gp{1,1};
dtUe_Gp{1,1} = - Ue_Div_ve_Gp{1,1} - Pet_grad_ve_Gp{1,1} - dkQiik_05_e_Gp{1,1};

dtKi_Gp{1,1} = - Ki_Div_vi_Gp{1,1} - vi_Div_Pit_Gp{1,1} + qi_ni_vi_E_Gp{1,1};
dtKe_Gp{1,1} = - Ke_Div_ve_Gp{1,1} - ve_Div_Pet_Gp{1,1} + qe_ne_ve_E_Gp{1,1};
%--------------------------------------------------------------------------
%}


% To Make the plots while saving the *_Gp data. <--------------------------
%This might be useful to keep checking. That is why I don't remove from here yet 
%--------------------------------------------------------------------------
%{

%--------------------------------------------------------------------------
% The field lines need to be calculated on the new reference frame 
%--------------------------------------------------------------------------
% Plot the field over the slide and the vectors in the NEW frame. 
%--------------------------------------------------------------------------
X = B_Gp{1,3};
%--------------------------------------------------------------------------
% built the grid to do the right plots
% this migth ruin the plots as it was used for the untranformed quantities
size_X=size(X);
xll=linspace(1,size_X(2),size_X(2))*0.06;
yll=linspace(1,size_X(1),size_X(1))*0.06;
[XLL, YLL] = meshgrid(xll,yll);
%--------------------------------------------------------------------------

% The point where the reconnection is happening is around 
xp_mr = xll(99); % this is in the new RF
yp_mr = yll(127);
%--------------------------------------------------------------------------

% This are the projections for the magnetic field
%ubx = slice_Gp{1,1};   vby = slice_Gp{1,2}; wbz = slice_Gp{1,3};
ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};

%ubx = B_Gp{1,2};   vby = B_Gp{1,1}; wbz = B_Gp{1,3}; %Is it this way? it doesn't look that way

bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

% This is for the electron velocities vectors
%Vev_x = slice_Gp{2,1}; Vev_y = slice_Gp{2,2}; Vev_z = slice_Gp{2,3}; 
Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,2}; 

%Vev_x = ve_Gp{1,2}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,3}; 

Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

Vivm=sqrt(Viv_x.*Viv_x + Viv_y.*Viv_y + Viv_z.*Viv_z);
Viv_x=Viv_x.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));
Viv_y=Viv_y.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));

%NN=50; xstart = (max(xll)/1.2)*rand(NN,1); ystart = (max(yll)/1.2)*rand(NN,1); %Do this just once

NN=40; 
xstart1 = (resx)*ones(1,NN); ystart1 = linspace(0,yll(end),NN); 
xstart2 = linspace(0,xll(end),NN); ystart2 = (resy)*ones(1,NN); 
xstart3 = (xll(end)-resx)*ones(1,NN); ystart3 = linspace(0,yll(end),NN); 
xstart4 = linspace(0,xll(end),NN); ystart4 = (yll(end)-resy)*ones(1,NN); 
xstart = horzcat(xstart1,xstart2,xstart3,xstart4);
ystart = horzcat(ystart1,ystart2,ystart3,ystart4);

sv=7;
%sv=5;
%--------------------------------------------------------------------------

%PLOTS!
%--------------------------------------------------------------------------
%Plots on the slide. The velocities ions, electrons, current and magnetic field components
%--------------------------------------------------------------------------
%colorFVBl=[0.5 0.5 0.5 0.5];
colorFVBl=[0 0 0 0.5];

%{
f1=figure(1);
for i=1:12
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(4,3,i,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=vi_Gp{1,i};
    titl = 'vi';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=ve_Gp{1,i2};
    titl = 've';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=J_Gp{1,i3};
    titl = 'J';
    s=i3;
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=B_Gp{1,i4};
    if (i4==2)
        dum_paux=B_Gp{1,i4};
        dum_p = dum_paux ;%+ 0.1;
    end
    titl = 'B';
    s=i4;    
end
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
hcb=colorbar;
if(i==11)
    colormap(h1,hot);
    %caxis([-lim_yp 0])
else
    colormap(h1,BWR);
    caxis([-lim_yp lim_yp])
end
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
%hcb.Title.String='';
%set(get(hcb,'Title'),'String',titl + string(s) ,'FontSize',14)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = '$p / d_{i}$';
    ax.YTick = [2 4 6 8 10];
end
if (9<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];    
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
%set(h2,'AutoScale','on', 'AutoScaleFactor', 2, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
%}

%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f1,strcat('vi_ve_J_B_vev_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

%end


%-------------------------------------------------------------------------
aaa=stream2(XLL,YLL,ubx,vby,xstart,ystart);
bbb=cell2mat(aaa');
ccc=round(100.*bbb)./100;

cuni1=unique(ccc(:,1),'first');
cuni2=unique(ccc(:,2),'first');
[CU1,CU2]=meshgrid(cuni1,cuni2);
CUZ=cos(CU1);

%f79=figure(79)
%contourf(CU1,CU2,CUZ)
%hold on
%plot(ccc(:,1) , ccc(:,2), '*r')
%k=5;
%eps=1e-1;
%result =find(xll(10)-ccc(k,1) < eps);

sv=3;
f77=figure(77);
hlines=streamline(aaa);
hold on
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%set(h2,'AutoScale','on', 'AutoScaleFactor', 5, 'Color', 'b');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
hold off

f78=figure(78);
hlines=streamline(aaa);
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', 5, 'Color', 'b');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$v_{e} \ vectors \ and \ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
hold off

cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f77,strcat('Bl_slideRF2_per'+ string(N_steps) +'.png'));
saveas(f78,strcat('ve_vec_Bl_slideRF2_per'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

%cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
%scatter
%view(3);


%f78=figure(78)
%scatter(ccc(:,1) , ccc(:,2))
%--------------------------------------------------------------------------


%Remove this to see the rest. This is important!!
%
sv=3;
f13=figure(13);
for i=1:2
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(1,2,i,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i==1)
    i1=i;
    dum_p=B_Gp{1,2};
    titl = 'B_{a}';
    s=i1;
elseif (i==2)
    i2=i-3;
    %dum_p=J_Gp{1,2};
    %titl = 'J_{a}';
    dum_p=ve_Gp{1,3};
    titl = 've_{r}';
    s=i2;
%elseif (i==3)
%    i3=i-6;
%    dum_p=ve_Gp{1,3};
%    titl = 've_{r}';
%    s=i3;
end
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
if(i==1) 
    colormap(h1,cool);
    %caxis([0 lim_yp])
    %caxis([-lim_yp 0])
    caxis([-0.14 0])
elseif(i==2)
    %caxis([-lim_yp lim_yp])
    caxis([-0.28 0.28])
    colormap(h1,BWR);
end
set(hc,'edgecolor','none')
hcb=colorbar;
%colormap(jet);
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.XTick = [];
ax.YTick = [];
hcb.Ruler.Exponent = -1;
if (mod(i,3)==1)
    ax.YLabel.String = '$p / d_{i}$';
    ax.YTick = [2 4 6 8 10];
end
%if (1<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];    
%end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
if (i==2)
%%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%%set(h2,'AutoScale','on', 'AutoScaleFactor', 6, 'Color', 'k');
%colorbar('off')
end
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%

%
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f13,strcat('Ba_ver_slideRF2_per_2'+ string(N_steps) +'.png'));
%saveas(f1,strcat('vi_ve_J_B_vev_no_Bl_slideRF2_per'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%

clf(f77)
clf(f78)
clf(f13)

%}
%--------------------------------------------------------------------------

% Get the 1D cuts to see how the terms releates each other (keep just the energy terms)
%{
%--------------------------------------------------------------------------
ypc_i=127;
cut_1d=ve_Gp{1,1}';
ve_1_1d=cut_1d(:,ypc_i); 
cut_1d=ve_Gp{1,2}';
ve_2_1d=cut_1d(:,ypc_i); 
cut_1d=ve_Gp{1,3}';
ve_3_1d=cut_1d(:,ypc_i); 

% For electrons
% Kinetic
%cut_1d=partialKe_Gp{1,1}'; %partialKe_1d = cut_1d(:,ypc_i);
cut_1d=ve_grad_Ke_Gp{1,1}'; ve_grad_Ke_1d = cut_1d(:,ypc_i);
cut_1d=ve_Div_Pet_Gp{1,1}'; ve_Div_Pet_1d = cut_1d(:,ypc_i);
cut_1d=Ke_Div_ve_Gp{1,1}'; Ke_Div_ve_1d = cut_1d(:,ypc_i);
cut_1d=qe_ne_ve_E_Gp{1,1}'; qe_ne_ve_E_Gp_1d = cut_1d(:,ypc_i);

% Internal
%cut_1d=partialUe_Gp{1,1}'; partialUe_1d = cut_1d(:,ypc_i);
cut_1d=ve_grad_Ue_Gp{1,1}'; ve_grad_Ue_1d = cut_1d(:,ypc_i);
cut_1d=dkQiik_05_e_Gp{1,1}'; dkQiik_05_e_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Pet_grad_ve_Gp{1,1}'; Pet_grad_ve_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Pet_grad_ve_ii_Gp{1,1}'; Pet_grad_ve_ii_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Pet_grad_ve_ij_Gp{1,1}'; Pet_grad_ve_ij_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Ue_Div_ve_Gp{1,1}'; Ue_Div_ve_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Div_Pes_dot_ve_Gp{1,1}'; Div_Pes_dot_ve_Gp_1d = cut_1d(:,ypc_i);

% These are the plots of the 1D profiles
%--------------------------------------------------------------------------
%
% <--- interal {} (not internal energy) one
%{
f213 = figure(213);
%subplot(1,2,1)
%plot(xll,partialKe_1d + ve_grad_Ke_1d,'k')
%hold on
plot(xll,ve_Div_Pet_1d,'-*k')
hold on
plot(xll,Ke_Div_ve_1d,'b')
plot(xll,-qe_ne_ve_E_Gp_1d,'m')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('u_{e} \cdot \nabla \cdot P_{e}','K_{e}\nabla \cdot u_{e}','-qE\cdot v_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('Energy rate')
hold off
%subplot(1,2,2)
%dum_p=ve_Gp{1,2};
%hc = pcolor(XLL,YLL,dum_p);
%set(hc,'edgecolor','none')
%}

lw=1.5;
f214 = figure(214);
%subplot(1,2,1)
plot(xll,ve_1_1d,'-b','LineWidth',lw)
hold on
plot(xll,ve_2_1d,'-r','LineWidth',lw)
plot(xll,ve_3_1d,'-k','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
legend('u_{p}','u_{a}','u_{r}')
xlabel('r / d_{i}')
ylabel('components u_{e}')
xlim([4,8])
hold off

f215 = figure(215);
% Kinetic
subplot(1,2,1)
plot(xll,ve_grad_Ke_1d,'-','LineWidth',lw)
hold on
plot(xll,ve_Div_Pet_1d,'-','LineWidth',lw)
plot(xll,Ke_Div_ve_1d,'-','LineWidth',lw)
plot(xll,-qe_ne_ve_E_Gp_1d,'LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
legend('u_{e} \cdot \nabla K','u_{e} \cdot \nabla \cdot P_{e}','K_{e}\nabla \cdot u_{e}','-qE\cdot v_{e}')
xlabel('r / d_{i}')
%xlim([4,8])
ylabel('kinetic Energy rate terms')
hold off
% Internal
%------------------------------------
subplot(1,2,2)
plot(xll,ve_grad_Ue_1d,'LineWidth',lw)
hold on
plot(xll,dkQiik_05_e_Gp_1d,'LineWidth',lw)
plot(xll,Pet_grad_ve_Gp_1d,'LineWidth',lw)
plot(xll,Pet_grad_ve_ii_Gp_1d,'LineWidth',lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'LineWidth',lw)
plot(xll,Ue_Div_ve_Gp_1d,'LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
legend('u_{e} \cdot \nabla U','\nabla \cdot q','\nabla u_{e}:P_{e}','\nabla u_{e}:P_{iie}','\nabla u_{e}:P_{ije}','U \nabla \cdot u_{e} ')
xlabel('r / d_{i}')
%xlim([4,8])
ylabel('Internal Energy rate terms')
hold off

f216 = figure(216);
plot(xll,Div_Pes_dot_ve_Gp_1d,'-k','LineWidth',lw)
hold on
plot(xll,Pet_grad_ve_Gp_1d,'b','LineWidth',lw)
plot(xll,ve_Div_Pet_1d,'m','LineWidth',lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('\nabla \cdot( P_{e} \cdot u_{e}) ','P_{e}:\nabla u_{e}','u_{e} \cdot(\nabla \cdot P_{e})')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('Kinetic Energy rate terms')
xlim([4,8])
hold off

f217 = figure(217);
plot(xll,Pet_grad_ve_Gp_1d,'-k','LineWidth',lw)
hold on
plot(xll,Pet_grad_ve_ii_Gp_1d,'b','LineWidth',lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'m','LineWidth',lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('P_{e}:\nabla u_{e}','(P_{e}:\nabla u_{e})_{ii}','(P_{e}:\nabla u_{e})_{ij}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('Energy rate terms')
%xlim([4,8])
hold off

f220 = figure(220);
dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum1,'b','LineWidth',lw)
hold on
dum2=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum2,'m','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('10+(P_{e}:\nabla u_{e}-(P_{e}:\nabla u_{e})_{ii})/P_{e}:\nabla u_{e}',...
    '10+(P_{e}:\nabla u_{e}-(P_{e}:\nabla u_{e})_{ij})/P_{e}:\nabla u_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('Energy rate terms')%xlim([4,8])
hold off

%{
f218 = figure(218);
dum1=sign(Pet_grad_ve_ii_Gp_1d).*log(Pet_grad_ve_ii_Gp_1d./Pet_grad_ve_Gp_1d);
plot(xll,dum1,'b','LineWidth',lw)
hold on
dum2=sign(Pet_grad_ve_ij_Gp_1d).*log(Pet_grad_ve_ij_Gp_1d./Pet_grad_ve_Gp_1d);
plot(xll,dum2,'m','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('(P_{e}:\nabla u_{e})_{ii}/P_{e}:\nabla u_{e}','(P_{e}:\nabla u_{e})_{ij}/P_{e}:\nabla u_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('log(P_{e}:\nabla u_{e}_{ij})/P_{e}:\nabla u_{e})')
%xlim([4,8])
hold off


f219 = figure(219);
dum1=(Pet_grad_ve_ii_Gp_1d./Pet_grad_ve_Gp_1d);
semilogy(xll,dum1,'b','LineWidth',lw)
hold on
dum2=1+(Pet_grad_ve_ij_Gp_1d./Pet_grad_ve_Gp_1d);
semilogy(xll,dum2,'m','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('(P_{e}:\nabla u_{e})_{ii}/P_{e}:\nabla u_{e}','1+(P_{e}:\nabla u_{e})_{ij}/P_{e}:\nabla u_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('(P_{e}:\nabla u_{e}_{ij})/P_{e}:\nabla u_{e})')
%xlim([4,8])
hold off
%}

cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f220,strcat('PGu_PGuij_over_PGu_10_'+ string(N_steps) +'.png'));
%saveas(f217,strcat('PGu_PGuii_PGuij'+ string(N_steps) +'.png'));
%saveas(f216,strcat('divPu_PGu_udivP'+ string(N_steps) +'.png'));
%saveas(f215,strcat('kinetic_internal_subset'+ string(N_steps) +'.png'));
saveas(f214,strcat('upuaur_sub'+ string(N_steps) +'.png'));
%--------------------------------------------------------------------------

cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

%}
%--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
% This is to save the variables *_Gp to store them for all the times and
% just load that file allowing to pick different times
%--------------------------------------------------------------------------


% This also shows which are the variables that are currently available in
% the new reference frame over the slide
%--------------------------------------------------------------------------
%{
% Vectors
B_Gp_j =B_Gp;
Poyn_v_Gp_j = Poyn_v_Gp; 

%Tensors
Pij_e_Gp_j = Pij_e_Gp;
Pij_i_Gp_j = Pij_i_Gp;

% Scalars
dt_E_em_Gp_j = dt_E_em_Gp; %This is the one calculated adding the terms 
Div_Poyn_v_Gp_j = Div_Poyn_v_Gp;

ne_Gp_j = ne_Gp;
ni_Gp_j = ni_Gp;

% Kinetic energy terms
ve_grad_Ke_Gp_j = ve_grad_Ke_Gp;
vi_grad_Ki_Gp_j = vi_grad_Ki_Gp;
ve_Div_Pet_Gp_j = ve_Div_Pet_Gp;
vi_Div_Pit_Gp_j = vi_Div_Pit_Gp;

Ke_Div_ve_Gp_j = Ke_Div_ve_Gp;
Ki_Div_vi_Gp_j = Ki_Div_vi_Gp;
qe_ne_ve_E_Gp_j = qe_ne_ve_E_Gp;
qi_ni_vi_E_Gp_j = qi_ni_vi_E_Gp;

% Computed by adding the other terms 
%-------------------------------
dtKe_Gp_j = dtKe_Gp;
dtKi_Gp_j = dtKi_Gp;
%-------------------------------

% Internal energy terms
ve_grad_Ue_Gp_j = ve_grad_Ue_Gp;
vi_grad_Ui_Gp_j = vi_grad_Ui_Gp;

dk_Qijk_e_Gp_j = dk_Qijk_e_Gp; 
dk_Qijk_i_Gp_j = dk_Qijk_i_Gp;
dkQiik_05_e_Gp_j = dkQiik_05_e_Gp;
dkQiik_05_i_Gp_j = dkQiik_05_i_Gp;

Pet_grad_ve_Gp_j = Pet_grad_ve_Gp;
Pit_grad_vi_Gp_j = Pit_grad_vi_Gp;

Pet_grad_ve_ii_Gp_j = Pet_grad_ve_ii_Gp;
Pit_grad_vi_ii_Gp_j = Pit_grad_vi_ii_Gp;

Pet_grad_ve_ij_Gp_j = Pet_grad_ve_ij_Gp;
Pit_grad_vi_ij_Gp_j = Pit_grad_vi_ij_Gp;

Ue_Div_ve_Gp_j = Ue_Div_ve_Gp; 
Ui_Div_vi_Gp_j = Ui_Div_vi_Gp; 

% Computed by adding the other terms 
%-------------------------------
dtUe_Gp_j = dtUe_Gp;
dtUi_Gp_j = dtUi_Gp;
%-------------------------------
%}
%--------------------------------------------------------------------------

%Save the data on the plane.  <----------------This is to save the data
%--------------------------------------------------------------------------

% Vectors
Vectorsc{1,in} =B_Gp;
Vectorsc{2,in} =E_Gp;
Vectorsc{3,in} =vi_Gp;
Vectorsc{4,in} =ve_Gp;
Vectorsc{5,in} =J_Gp;
Vectorsc{6,in} = Poyn_v_Gp; 

%Tensors
Tensorsc{1,in} = Pij_i_Gp;
Tensorsc{2,in} = Pij_e_Gp;
Tensorsc{3,in} = grad_vi_Gp;
Tensorsc{4,in} = grad_ve_Gp;

% Densities
Densitiesc{1,in} = ni_Gp;
Densitiesc{2,in} = ne_Gp;

% Poyn theorem terms
Poyn_theoc{1,in} = dt_E_em_Gp; %This is the one calculated adding the terms 
Poyn_theoc{2,in} = Div_Poyn_v_Gp;
Poyn_theoc{3,in} = JdotE_Gp;

%Kinetic Energy
kin_energy_terms_ic{1,in} = vi_grad_Ki_Gp;
kin_energy_terms_ic{2,in} = vi_Div_Pit_Gp;
kin_energy_terms_ic{3,in} = Ki_Div_vi_Gp;
kin_energy_terms_ic{4,in} = qi_ni_vi_E_Gp;

kin_energy_terms_ec{1,in} = ve_grad_Ke_Gp;
kin_energy_terms_ec{2,in} = ve_Div_Pet_Gp;
kin_energy_terms_ec{3,in} = Ke_Div_ve_Gp;
kin_energy_terms_ec{4,in} = qe_ne_ve_E_Gp;

% Internal energy terms
Int_energy_terms_ic{1,in} = vi_grad_Ui_Gp;
Int_energy_terms_ic{2,in} = dk_Qijk_i_Gp;
Int_energy_terms_ic{3,in} = dkQiik_05_i_Gp;
Int_energy_terms_ic{4,in} = Pit_grad_vi_Gp;
Int_energy_terms_ic{5,in} = Pit_grad_vi_ii_Gp;
Int_energy_terms_ic{6,in} = Pit_grad_vi_ij_Gp;
Int_energy_terms_ic{7,in} = Ui_Div_vi_Gp; 

Int_energy_terms_ec{1,in} = ve_grad_Ue_Gp;
Int_energy_terms_ec{2,in} = dk_Qijk_e_Gp; 
Int_energy_terms_ec{3,in} = dkQiik_05_e_Gp;
Int_energy_terms_ec{4,in} = Pet_grad_ve_Gp;
Int_energy_terms_ec{5,in} = Pet_grad_ve_ii_Gp;
Int_energy_terms_ec{6,in} = Pet_grad_ve_ij_Gp;
Int_energy_terms_ec{7,in} = Ue_Div_ve_Gp; 

%--------------------------------------------------------------------------


cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Save in a single file 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%{
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'

save('variabletimefile.mat','Vectorsc','Tensorsc','Densitiesc','Poyn_theoc',...
    'kin_energy_terms_ic','kin_energy_terms_ec',...
    'Int_energy_terms_ic','Int_energy_terms_ec')
%}
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Once the data are stored in a file  <------------------------------------
%--------------------------------------------------------------------------
load ('variabletimefile.mat') %<----- This is the one that actually has the energy terms in time
%dt_energydensities_ic=cell(5,N);

% This was a one time workspace thing.
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
load ('energy_things.mat') %   <------------------------


in=1;
in=16; %<--- 96
in=19;
for in=1:N-1
%--------------------------------------------------------------------------
disp(strcat('Computing step ...',S(in).name))
%--------------------------------------------------------------------------
%N_steps = 2000;
%time_t = dt*N_steps; %(In terms of 1/omega_pi)
dt = 0.06;
str = string(S(in).name);
newStr = extractBetween(str,"pfd.","_p000000.h5");
N_steps = str2double(newStr); time_t = dt*N_steps;

%--------------------------------------------------------------------------

% Vectors
B_Gp = Vectorsc{1,in} ;
%{
E_Gp = Vectorsc{2,in} ;
vi_Gp = Vectorsc{3,in} ;
ve_Gp = Vectorsc{4,in} ;
J_Gp = Vectorsc{5,in} ;
Poyn_v_Gp = Vectorsc{6,in} ; 

%Tensors
Pij_i_Gp = Tensorsc{1,in} ;
Pij_e_Gp = Tensorsc{2,in} ;
grad_vi_Gp = Tensorsc{3,in} ;
grad_ve_Gp = Tensorsc{4,in} ;

% Densities
ni_Gp = Densitiesc{1,in} ;
ne_Gp = Densitiesc{2,in} ;

% Poyn theorem terms
dt_E_em_Gp = Poyn_theoc{1,in} ; %This is the one calculated adding the terms 
Div_Poyn_v_Gp = Poyn_theoc{2,in} ;
JdotE_Gp = Poyn_theoc{3,in} ;

%Kinetic Energy
vi_grad_Ki_Gp = kin_energy_terms_ic{1,in} ;
vi_Div_Pit_Gp = kin_energy_terms_ic{2,in} ;
Ki_Div_vi_Gp = kin_energy_terms_ic{3,in} ;
qi_ni_vi_E_Gp = kin_energy_terms_ic{4,in} ;

ve_grad_Ke_Gp = kin_energy_terms_ec{1,in} ;
ve_Div_Pet_Gp = kin_energy_terms_ec{2,in} ;
Ke_Div_ve_Gp = kin_energy_terms_ec{3,in} ;
qe_ne_ve_E_Gp = kin_energy_terms_ec{4,in} ;

% Internal energy terms
vi_grad_Ui_Gp = Int_energy_terms_ic{1,in} ;
dk_Qijk_i_Gp = Int_energy_terms_ic{2,in} ;
dkQiik_05_i_Gp = Int_energy_terms_ic{3,in} ;
Pit_grad_vi_Gp = Int_energy_terms_ic{4,in} ;
Pit_grad_vi_ii_Gp = Int_energy_terms_ic{5,in} ;
Pit_grad_vi_ij_Gp = Int_energy_terms_ic{6,in} ;
Ui_Div_vi_Gp = Int_energy_terms_ic{7,in} ; 

ve_grad_Ue_Gp = Int_energy_terms_ec{1,in} ;
dk_Qijk_e_Gp = Int_energy_terms_ec{2,in} ; 
dkQiik_05_e_Gp = Int_energy_terms_ec{3,in} ;
Pet_grad_ve_Gp =Int_energy_terms_ec{4,in} ;
Pet_grad_ve_ii_Gp = Int_energy_terms_ec{5,in} ;
Pet_grad_ve_ij_Gp = Int_energy_terms_ec{6,in} ;
Ue_Div_ve_Gp = Int_energy_terms_ec{7,in} ; 

%--------------------------------------------------------------------------
% Div_Pi_dot_vi
Div_Pis_dot_vi_Gp{1,1} = Pit_grad_vi_Gp{1,1} + vi_Div_Pit_Gp{1,1};
% Div_Pe_dot_e
Div_Pes_dot_ve_Gp{1,1} = Pet_grad_ve_Gp{1,1} + ve_Div_Pet_Gp{1,1};
%--------------------------------------------------------------------------
%}
%--------------------------------------------------------------------------
%calculating the time derivatives
%{
%--------------------------------------------------------------------------
str_1 = string(S(in-1).name); newStr_1 = extractBetween(str_1,"pfd.","_p000000.h5");
N_steps_1 = str2double(newStr_1); time_1=dt*N_steps_1;
str_2 = string(S(in+0).name); newStr_2 = extractBetween(str_2,"pfd.","_p000000.h5");
N_steps_2 = str2double(newStr_2); time_2=dt*N_steps_2;
str_3 = string(S(in+1).name); newStr_3 = extractBetween(str_3,"pfd.","_p000000.h5");
N_steps_3 = str2double(newStr_3); time_3=dt*N_steps_3;

%--------------------------------------------------------------------------
densities_time_d=cell(5,3);
for k=in-1:in+1
    %disp(k)
B_Gp = Vectorsc{1,k} ; E_Gp = Vectorsc{2,k} ;
vi_Gp = Vectorsc{3,k} ; ve_Gp = Vectorsc{4,k} ;
ni_Gp = Densitiesc{1,k} ; ne_Gp = Densitiesc{2,k} ;

Pij_i_Gp = Tensorsc{1,k}; 
Pij_e_Gp = Tensorsc{1,k};

% Kinetic energy densities 
%--------------------------------------------------------------------------
Ki_Gp_1 = 0.5.*ni_Gp{1,1}.* mi.*(vi_Gp{1,1}.*vi_Gp{1,1} + vi_Gp{1,2}.*vi_Gp{1,2} + vi_Gp{1,3}.*vi_Gp{1,3});
Ke_Gp_1 = 0.5.*ne_Gp{1,1}.* me.*(ve_Gp{1,1}.*ve_Gp{1,1} + ve_Gp{1,2}.*ve_Gp{1,2} + ve_Gp{1,3}.*vi_Gp{1,3});
%--------------------------------------------------------------------------
% Internal energy densities 
%--------------------------------------------------------------------------
Ui_Gp_1 = 0.5.*(Pij_i_Gp{1,1} + Pij_i_Gp{2,2} + Pij_i_Gp{3,3});
Ue_Gp_1 = 0.5.*(Pij_e_Gp{1,1} + Pij_e_Gp{2,2} + Pij_e_Gp{3,3});
%--------------------------------------------------------------------------
% Electromagnetic energy density
%--------------------------------------------------------------------------
E2_Gp = E_Gp{1,1}.^2 + E_Gp{1,2}.^2 + E_Gp{1,3}.^2; 
B2_Gp = B_Gp{1,1}.^2 + B_Gp{1,2}.^2 + B_Gp{1,3}.^2; 
E_em_1 = 0.5.*( eps0.*E2_Gp + B2_Gp./mu0 );
%--------------------------------------------------------------------------
% The Three steps.
%--------------------------------------------------------------------------
densities_time_d{1,k-in+2} = Ki_Gp_1;
densities_time_d{2,k-in+2} = Ke_Gp_1;
densities_time_d{3,k-in+2} = Ui_Gp_1;
densities_time_d{4,k-in+2} = Ue_Gp_1;
densities_time_d{5,k-in+2} = E_em_1;
end

% These are the actual time derivatives (t32 + t21)/2
dtKi_21 = (densities_time_d{1,2} - densities_time_d{1,1})./(time_2 -time_1);  
dtKi_32 = (densities_time_d{1,3} - densities_time_d{1,2})./(time_3 -time_2);
dtKi_in = 0.5.*(dtKi_21 + dtKi_32);

dtKe_21 = (densities_time_d{2,2} - densities_time_d{2,1})./(time_2 -time_1);  
dtKe_32 = (densities_time_d{2,3} - densities_time_d{2,2})./(time_3 -time_2);
dtKe_in = 0.5.*(dtKe_21 + dtKe_32);

dtUi_21 = (densities_time_d{3,2} - densities_time_d{3,1})./(time_2 -time_1);  
dtUi_32 = (densities_time_d{3,3} - densities_time_d{3,2})./(time_3 -time_2);
dtUi_in = 0.5.*(dtUi_21 + dtUi_32);

dtUe_21 = (densities_time_d{4,2} - densities_time_d{4,1})./(time_2 -time_1);  
dtUe_32 = (densities_time_d{4,3} - densities_time_d{4,2})./(time_3 -time_2);
dtUe_in = 0.5.*(dtUe_21 + dtUe_32);

dtE_em_21 = (densities_time_d{5,2} - densities_time_d{5,1})./(time_2 -time_1);  
dtE_em_32 = (densities_time_d{5,3} - densities_time_d{5,2})./(time_3 -time_2);
dtE_em_in = 0.5.*(dtE_em_21 + dtE_em_32);
%}

%--------------------------------------------------------------------------
%Save the time derivatives at each time to plot them along with the others 
%{
dt_energydensities_ic{1,in} = dtKi_in;
dt_energydensities_ic{2,in} = dtKe_in;
dt_energydensities_ic{3,in} = dtUi_in;
dt_energydensities_ic{4,in} = dtUe_in;
dt_energydensities_ic{5,in} = dtE_em_in;
%}
%--------------------------------------------------------------------------

% This only works once dt_energydensities_ic is loaded
%{
%--------------------------------------------------------------------------
dtKi_in_Gp_2{1,1} = dt_energydensities_ic{1,in} ;
dtKe_in_Gp_2{1,1} = dt_energydensities_ic{2,in} ;
dtUi_in_Gp_2{1,1} = dt_energydensities_ic{3,in} ;
dtUe_in_Gp_2{1,1} = dt_energydensities_ic{4,in} ;
dtE_em_in_Gp_2{1,1} = dt_energydensities_ic{5,in} ;
%}

%--------------------------------------------------------------------------


% Get the 1D cuts to see how the terms releates each other (keep just the energy terms)
%--------------------------------------------------------------------------
%{
%----------------------------------------------------------------------
ypc_i=127; % <---- Here is where I choose the 1D line that crosses the reconnection event 
cut_1d=ve_Gp{1,1}';
ve_1_1d=cut_1d(:,ypc_i); 
cut_1d=ve_Gp{1,2}';
ve_2_1d=cut_1d(:,ypc_i); 
cut_1d=ve_Gp{1,3}';
ve_3_1d=cut_1d(:,ypc_i); 

% For electrons

% Time derivatives
cut_1d=dtKe_in_Gp_2{1,1}'; dtKe_in_Gp_1d = cut_1d(:,ypc_i);
cut_1d=dtUe_in_Gp_2{1,1}'; dtUe_in_Gp_1d = cut_1d(:,ypc_i);

% Kinetic
%cut_1d=partialKe_Gp{1,1}'; %partialKe_1d = cut_1d(:,ypc_i);
cut_1d=ve_grad_Ke_Gp{1,1}'; ve_grad_Ke_1d = cut_1d(:,ypc_i);
cut_1d=ve_Div_Pet_Gp{1,1}'; ve_Div_Pet_1d = cut_1d(:,ypc_i);
cut_1d=Ke_Div_ve_Gp{1,1}'; Ke_Div_ve_1d = cut_1d(:,ypc_i);
cut_1d=qe_ne_ve_E_Gp{1,1}'; qe_ne_ve_E_Gp_1d = cut_1d(:,ypc_i);

% Internal
%cut_1d=partialUe_Gp{1,1}'; partialUe_1d = cut_1d(:,ypc_i);
cut_1d=ve_grad_Ue_Gp{1,1}'; ve_grad_Ue_1d = cut_1d(:,ypc_i);
cut_1d=dkQiik_05_e_Gp{1,1}'; dkQiik_05_e_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Pet_grad_ve_Gp{1,1}'; Pet_grad_ve_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Pet_grad_ve_ii_Gp{1,1}'; Pet_grad_ve_ii_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Pet_grad_ve_ij_Gp{1,1}'; Pet_grad_ve_ij_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Ue_Div_ve_Gp{1,1}'; Ue_Div_ve_Gp_1d = cut_1d(:,ypc_i);
cut_1d=Div_Pes_dot_ve_Gp{1,1}'; Div_Pes_dot_ve_Gp_1d = cut_1d(:,ypc_i);
%}
%--------------------------------------------------------------------------

%  <--- interal {} (not internal energy) one
%{
f213 = figure(213);
%subplot(1,2,1)
%plot(xll,partialKe_1d + ve_grad_Ke_1d,'k')
%hold on
plot(xll,ve_Div_Pet_1d,'-*k')
hold on
plot(xll,Ke_Div_ve_1d,'b')
plot(xll,-qe_ne_ve_E_Gp_1d,'m')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('u_{e} \cdot \nabla \cdot P_{e}','K_{e}\nabla \cdot u_{e}','-qE\cdot v_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('Energy rate')
hold off
%subplot(1,2,2)
%dum_p=ve_Gp{1,2};
%hc = pcolor(XLL,YLL,dum_p);
%set(hc,'edgecolor','none')
%}

fs=18; lw=1.5;

% 1D plots for electrons no smoothing 
%--------------------------------------------------------------------------
%{
f214 = figure(214);
plot(xll,ve_1_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_2_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_3_1d,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}$','$u_{a}$','$u_{r}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$u_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f2145 = figure(2145);
plot(xll,ve_1_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_2_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_3_1d,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}$','$u_{a}$','$u_{r}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$u_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f215 = figure(215);
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ke_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_Div_Pet_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Ke_Div_ve_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,-qe_ne_ve_E_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{k}_{e}$','$u_{e} \cdot \nabla \cdot \overline{P}_{e}$',...
    '$\varepsilon^{k}_{e}\nabla \cdot u_{e}$','$-qE\cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Internal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ue_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,dkQiik_05_e_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,Ue_Div_ve_Gp_1d,'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{th}_{e}$','$\nabla \cdot h_{e}$',...
    '$\nabla u_{e}:\overline{P}_{e}$','$p\theta_{e}$','$PiD_{e}$',...
    '$\varepsilon^{th}_{e} \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylabel('Internal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f2155 = figure(2155);
% Kinetic
%subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ke_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_Div_Pet_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Ke_Div_ve_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,-qe_ne_ve_E_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{k}_{e}$','$u_{e} \cdot \nabla \cdot \overline{P}_{e}$',...
    '$\varepsilon^{k}_{e}\nabla \cdot u_{e}$','$-qE\cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Internal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ue_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,dkQiik_05_e_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,Ue_Div_ve_Gp_1d,'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{th}_{e}$','$\nabla \cdot h_{e}$',...
    '$\nabla u_{e}:\overline{P}_{e}$','$p\theta_{e}$','$PiD_{e}$',...
    '$\varepsilon^{th}_{e} \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylabel('Internal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f216 = figure(216);
plot(xll,Div_Pes_dot_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_Div_Pet_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\nabla \cdot( \overline{P}_{e} \cdot u_{e})$ ','$\overline{P}_{e}:\nabla u_{e}$',...
    '$u_{e} \cdot(\nabla \cdot \overline{P}_{e})$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
hold off

f2165 = figure(2165);
plot(xll,Div_Pes_dot_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_Div_Pet_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\nabla \cdot( \overline{P}_{e} \cdot u_{e})$ ','$\overline{P}_{e}:\nabla u_{e}$',...
    '$u_{e} \cdot(\nabla \cdot \overline{P}_{e})$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f217 = figure(217);
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\overline{P}_{e}:\nabla u_{e}$','$p\theta_{e}$','$PiD_{e}$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
hold off

f2175 = figure(2175);
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\overline{P}_{e}:\nabla u_{e}$','$p\theta_{e}$','$PiD_{e}$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f218 = figure(218);
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=(Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=(Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
hold off

f2185 = figure(2185);
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=(Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=(Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off
%}
%--------------------------------------------------------------------------


% 1D Plots for electros smoothed   <------------------------------------
%--------------------------------------------------------------------------
%{
f214 = figure(214);
plot(xll,smoothdata(ve_1_1d),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_2_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_3_1d),'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}$','$u_{a}$','$u_{r}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$u_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-0.5, 0.5])
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
hold off

f2145 = figure(2145);
plot(xll,smoothdata(ve_1_1d),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_2_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_3_1d),'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}$','$u_{a}$','$u_{r}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$u_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-0.5, 0.5])
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
hold off

f215 = figure(215);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(dtKe_in_Gp_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_grad_Ke_1d),'Color',colorblind(1,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_Div_Pet_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Ke_Div_ve_1d),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(-qe_ne_ve_E_Gp_1d),'Color',colorblind(6,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\partial \varepsilon^{k}_{e} /\partial t $',...
    '$u_{e} \cdot \nabla \varepsilon^{k}_{e}$',...
    '$u_{e} \cdot \nabla \cdot \overline{P}_{e}$',...
    '$\varepsilon^{k}_{e}\nabla \cdot u_{e}$',...
    '$-q_{e}n_{e}E\cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([-0.005, 0.005])
%ylim([-0.03, 0.03]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Internal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(dtUe_in_Gp_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_grad_Ue_1d),'Color',colorblind(1,:), 'LineWidth', lw)
plot(xll,smoothdata(dkQiik_05_e_Gp_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ii_Gp_1d),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ij_Gp_1d),'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,smoothdata(Ue_Div_ve_Gp_1d),'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\partial \varepsilon^{th}_{e} / \partial t $','$u_{e} \cdot \nabla \varepsilon^{th}_{e}$','$\nabla \cdot h_{e}$',...
    '$\nabla u_{e}:\overline{P}_{e}$','$p\theta_{e}$','$PiD_{e}$',...
    '$\varepsilon^{th}_{e} \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
%ylim([-0.03, 0.03]) %<--------------
ylim([-0.015, 0.015])
ylabel('Internal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f2155 = figure(2155);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(dtKe_in_Gp_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_grad_Ke_1d),'Color',colorblind(1,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_Div_Pet_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Ke_Div_ve_1d),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(-qe_ne_ve_E_Gp_1d),'Color',colorblind(6,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\partial \varepsilon^{k}_{e} /\partial t $',...
    '$u_{e} \cdot \nabla \varepsilon^{k}_{e}$',...
    '$u_{e} \cdot \nabla \cdot \overline{P}_{e}$',...
    '$\varepsilon^{k}_{e}\nabla \cdot u_{e}$',...
    '$-q_{e}n_{e}E\cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-0.005, 0.005])
%ylim([-0.03, 0.03]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Internal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(dtUe_in_Gp_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_grad_Ue_1d),'Color',colorblind(1,:), 'LineWidth', lw)
plot(xll,smoothdata(dkQiik_05_e_Gp_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ii_Gp_1d),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ij_Gp_1d),'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,smoothdata(Ue_Div_ve_Gp_1d),'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\partial \varepsilon^{th}_{e} / \partial t $','$u_{e} \cdot \nabla \varepsilon^{th}_{e}$','$\nabla \cdot h_{e}$',...
    '$\nabla u_{e}:\overline{P}_{e}$','$p\theta_{e}$','$PiD_{e}$',...
    '$\varepsilon^{th}_{e} \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
%ylim([-0.03, 0.03]) %<--------------
ylim([-0.015, 0.015])
ylabel('Internal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f216 = figure(216);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Div_Pes_dot_ve_Gp_1d),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pet_grad_ve_Gp_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_Div_Pet_1d),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\nabla \cdot( \overline{P}_{e} \cdot u_{e})$ ','$\overline{P}_{e}:\nabla u_{e}$',...
    '$u_{e} \cdot(\nabla \cdot \overline{P}_{e})$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
ylim([-0.04, 0.04])
%xlim([4,8])
hold off

f2165 = figure(2165);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Div_Pes_dot_ve_Gp_1d),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pet_grad_ve_Gp_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_Div_Pet_1d),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\nabla \cdot( \overline{P}_{e} \cdot u_{e})$ ','$\overline{P}_{e}:\nabla u_{e}$',...
    '$u_{e} \cdot(\nabla \cdot \overline{P}_{e})$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-0.04, 0.04])
hold off

f217 = figure(217);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pet_grad_ve_ii_Gp_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ij_Gp_1d),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\overline{P}_{e}:\nabla u_{e}$','$p\theta_{e}$','$PiD_{e}$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([-0.03, 0.03])
hold off

f2175 = figure(2175);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pet_grad_ve_ii_Gp_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ij_Gp_1d),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\overline{P}_{e}:\nabla u_{e}$','$p\theta_{e}$','$PiD_{e}$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-0.03, 0.03])
hold off

f218 = figure(218);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=smoothdata((Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=smoothdata((Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([0.1, 100])
hold off

f2185 = figure(2185);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=smoothdata((Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=smoothdata((Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([0.1, 100])
hold off

% This is aboutsome olf plotd of the diagonal and off diagonal terms
%--------------------------------------------------------------------------
%{
f218 = figure(218);
dum1=sign(Pet_grad_ve_ii_Gp_1d).*log(Pet_grad_ve_ii_Gp_1d./Pet_grad_ve_Gp_1d);
plot(xll,dum1,'b','LineWidth',lw)
hold on
dum2=sign(Pet_grad_ve_ij_Gp_1d).*log(Pet_grad_ve_ij_Gp_1d./Pet_grad_ve_Gp_1d);
plot(xll,dum2,'m','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('(P_{e}:\nabla u_{e})_{ii}/P_{e}:\nabla u_{e}','(P_{e}:\nabla u_{e})_{ij}/P_{e}:\nabla u_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('log(P_{e}:\nabla u_{e}_{ij})/P_{e}:\nabla u_{e})')
%xlim([4,8])
hold off

f219 = figure(219);
dum1=(Pet_grad_ve_ii_Gp_1d./Pet_grad_ve_Gp_1d);
semilogy(xll,dum1,'b','LineWidth',lw)
hold on
dum2=1+(Pet_grad_ve_ij_Gp_1d./Pet_grad_ve_Gp_1d);
semilogy(xll,dum2,'m','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('(P_{e}:\nabla u_{e})_{ii}/P_{e}:\nabla u_{e}','1+(P_{e}:\nabla u_{e})_{ij}/P_{e}:\nabla u_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('(P_{e}:\nabla u_{e}_{ij})/P_{e}:\nabla u_{e})')
%xlim([4,8])
hold off
%}
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f214,strcat('s_upuaur_'+ string(N_steps) +'.png'));
saveas(f2145,strcat('s_upuaur_subset_'+ string(N_steps) +'.png'));
saveas(f215,strcat('s_kinetic_internal_'+ string(N_steps) +'.png'));
saveas(f2155,strcat('s_kinetic_internal_subset_'+ string(N_steps) +'.png'));
saveas(f216,strcat('s_divPu_PGu_udivP_'+ string(N_steps) +'.png'));
saveas(f2165,strcat('s_divPu_PGu_udivP_subset_'+ string(N_steps) +'.png'));
saveas(f217,strcat('s_PGu_PGuii_PGuij_'+ string(N_steps) +'.png'));
saveas(f2175,strcat('s_PGu_PGuii_PGuij_subset_'+ string(N_steps) +'.png'));
saveas(f218,strcat('s_PGu_PGuij_over_PGu_10_'+ string(N_steps) +'.png'));
saveas(f2185,strcat('s_PGu_PGuij_over_PGu_10_subset'+ string(N_steps) +'.png'));

clf(f214); clf(f2145); 
clf(f215); clf(f2155);
clf(f216); clf(f2165);
clf(f217); clf(f2175);
clf(f218); clf(f2185);
%}
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here'

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Plots 2D on the slide last used

%--------------------------------------------------------------------------
%PLOTS!
%--------------------------------------------------------------------------

%repeated control
%--------------------------------------------------------------------------
%{
%--------------------------------------------------------------------------
X = B_Gp{1,3};
%--------------------------------------------------------------------------
% built the grid to do the right plots
% this migth ruin the plots as it was used for the untranformed quantities
size_X=size(X);
xll=linspace(1,size_X(2),size_X(2))*0.06;
yll=linspace(1,size_X(1),size_X(1))*0.06;
[XLL, YLL] = meshgrid(xll,yll);
%--------------------------------------------------------------------------
% The point where the reconnection is happening is around 
xp_mr = xll(99); % this is in the new RF
yp_mr = yll(127);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% This are the projections for the magnetic field
ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};
bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));
%{
% This is for the electron velocities vectors
Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,2}; % <----- this seems to be the one that works
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Viv_x = vi_Gp{1,3}; Viv_y = vi_Gp{1,1}; Viv_z = vi_Gp{1,2}; 
Vivm=sqrt(Viv_x.*Viv_x + Viv_y.*Viv_y + Viv_z.*Viv_z);
Viv_x=Viv_x.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));
Viv_y=Viv_y.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));
%NN=50; xstart = (max(xll)/1.2)*rand(NN,1); ystart = (max(yll)/1.2)*rand(NN,1); %Do this just once
%}
NN=40; 
xstart1 = (resx)*ones(1,NN); ystart1 = linspace(0,yll(end),NN); 
xstart2 = linspace(0,xll(end),NN); ystart2 = (resy)*ones(1,NN); 
xstart3 = (xll(end)-resx)*ones(1,NN); ystart3 = linspace(0,yll(end),NN); 
xstart4 = linspace(0,xll(end),NN); ystart4 = (yll(end)-resy)*ones(1,NN); 
xstart = horzcat(xstart1,xstart2,xstart3,xstart4);
ystart = horzcat(ystart1,ystart2,ystart3,ystart4);
%--------------------------------------------------------------------------
% Plots
colorFVBl=[0.5 0.5 0.5 0.5];
colorFVBl=[0 0 0 0.5];
%--------------------------------------------------------------------------
aaa=stream2(XLL,YLL,ubx,vby,xstart,ystart);
%--------------------------------------------------------------------------
%}
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% These are the plots including the streamlines
%--------------------------------------------------------------------------
X = B_Gp{1,3};
%--------------------------------------------------------------------------
% built the grid to do the right plots
% this migth ruin the plots as it was used for the untranformed quantities
size_X=size(X);
xll=linspace(1,size_X(2),size_X(2))*0.06; %Jeff flag it works with (1,size_x(2))
yll=linspace(1,size_X(1),size_X(1))*0.06;
[XLL, YLL] = meshgrid(xll,yll);

resx = abs(xll(2)-xll(1)); 
resy = abs(yll(2)-yll(1)); 
%--------------------------------------------------------------------------

% The point where the reconnection is happening is around 
xp_mr = xll(99); % this is in the new RF
yp_mr = yll(127);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% This are the projections for the magnetic field
ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};
%--------------------------------------------------------------------------

% This is to check a dummy field. The plot needs to be run after the ccc
% streamlines have been calculated.
%--------------------------------------------------------------------------
%{
N = 100 ;
xl=linspace(1,6,N);%-2:4/(N-1) :2 ;
yl=linspace(1,6,N);% -2:4/(N-1 ) : 2;
[ x , y]=meshgrid(xl, yl) ;
z = (x -3 ).* exp(-(x -3).^ 2 - (y-3) .^ 2 ) ;
[ px , py ] = gradient ( z , 4 / (N-1) ,4/(N-1) ) ;
ubx = px; vby= py; wbz =px./px;
XLL=x ; YLL=y; xll=xl; yll=yl; 
resx = abs(xll(2)-xll(1)); 
resy = abs(yll(2)-yll(1)); 
NN=20; 
xx = linspace(xl(1),xl(end),NN);
yy = linspace(yl(1),yl(end),NN);
[xx yy] = meshgrid(xx,yy);
aaa51=stream2(x,y,ubx,vby,xx,yy);

%-------------------------------------------------------------
f780=figure(780)
h1=subplot(2,2,1);
h2=quiver (px , py ) ;
set(h2,'AutoScale','on', 'AutoScaleFactor', 2, 'Color', 'k');
h3=subplot(2,2,2);
dum_p = abs(real(EIGV1).*saddle_cond12./real(EIGV1));
hc=pcolor(dum_p');
%hc=pcolor(ubx);
colorbar
set(hc,'edgecolor','none'); colormap(jet(2))
h3=subplot(2,2,3);
hlines1=streamline(aaa51);
set(hlines1,'LineWidth',1,'Color', 'k');
h3=subplot(2,2,4);
hlines1=streamline(ccc12);
set(hlines1,'LineWidth',1.5,'Color', 'm');
%-------------------------------------------------------------
%}
%--------------------------------------------------------------------------

bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
%ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
%vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

% This is for the electron velocities vectors
Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,2}; % <----- this seems to be the one that works
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
%Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
%Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

Viv_x = vi_Gp{1,3}; Viv_y = vi_Gp{1,1}; Viv_z = vi_Gp{1,2}; 
Vivm=sqrt(Viv_x.*Viv_x + Viv_y.*Viv_y + Viv_z.*Viv_z);
%Viv_x=Viv_x.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));
%Viv_y=Viv_y.*sqrt(1-((Viv_z./Vivm).*(Viv_z./Vivm)));


% In this par reference frame there is no projection
%ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};
%Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,2};
%Viv_x = vi_Gp{1,3}; Viv_y = vi_Gp{1,1}; Viv_z = vi_Gp{1,2}; 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


%NN=50; xstart = (max(xll)/1.2)*rand(NN,1); ystart = (max(yll)/1.2)*rand(NN,1); %Do this just once

NN=40; 
xstart1 = (resx)*ones(1,NN); ystart1 = linspace(0,yll(end),NN); 
xstart2 = linspace(0,xll(end),NN); ystart2 = (resy)*ones(1,NN); 
xstart3 = (xll(end)-resx)*ones(1,NN); ystart3 = linspace(0,yll(end),NN); 
xstart4 = linspace(0,xll(end),NN); ystart4 = (yll(end)-resy)*ones(1,NN); 
xstart = horzcat(xstart1,xstart2,xstart3,xstart4);
ystart = horzcat(ystart1,ystart2,ystart3,ystart4);
aaa=stream2(XLL,YLL,ubx,vby,xstart,ystart);

NN4=10; 
xx = linspace(xll(10),xll(end-10),NN4);
yy = linspace(yll(10),yll(end-10),NN4);
[xx, yy] = meshgrid(xx,yy);
aaa5=stream2(XLL,YLL,ubx,vby,xx,yy);

NN5=10; 
xx5 = linspace(xll(10),xll(end-10),NN5);
yy5 = linspace(yll(10),yll(end-10),NN5);
[xx5, yy5] = meshgrid(xx5,yy5);

aaa=stream2(XLL,YLL,ubx,vby,xstart,ystart);
aaam=stream2(XLL,YLL,-ubx,-vby,xstart,ystart);

%veve=stream2(XLL,YLL,Vev_x,Vev_y,xstart,ystart);
%vivi=stream2(XLL,YLL,Viv_x,Viv_y,xstart,ystart);
veve5=stream2(XLL,YLL,Vev_x,Vev_y,xx5,yy5);
vivi5=stream2(XLL,YLL,Viv_x,Viv_y,xx5,yy5);
bbbb5=stream2(XLL,YLL,ubx,vby,xx5,yy5);

veve5m=stream2(XLL,YLL,-Vev_x,-Vev_y,xx5,yy5);
vivi5m=stream2(XLL,YLL,-Viv_x,-Viv_y,xx5,yy5);
bbbb5m=stream2(XLL,YLL,-ubx,-vby,xx5,yy5);

%--------------------------------------------------------------------------
%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%--------------------------------------------------------------------------

% To get the critical points and separatrix curves in the plane
%%%%
%--------------------------------------------------------------------------
%{
% Get the critical points (ubx = 0;   vby = 0; wbz =0;)
% Define the zero
%zero_d=0.01;
zero_d=0.007;
%zero_d=0.008; %<---
ubx_cr_i=(abs(ubx) < zero_d);vby_cr_i=(abs(vby) < zero_d);wbz_cr_i=(abs(wbz) < zero_d);

ubx_cr_nan=isnan(ubx);
ubvb_cr_i=((abs(ubx) < zero_d) & (abs(vby) < zero_d));
%ubvbwb_cr_i=((abs(ubx) < zero_d) & (abs(vby) < zero_d) & (abs(wbz) <
%zero_d)); % This one is hardly fulfilled

% Now lets get the Jacobian
[ubx_x,ubx_y] = gradient(ubx);
[vby_x,vby_y] = gradient(vby);

%--------------------------------------------------------------------------
ubx_x2 = ubx_x(ubvb_cr_i); ubx_y2 = ubx_y(ubvb_cr_i);  % this in array representation
vby_x2 = vby_x(ubvb_cr_i); vby_y2 = vby_y(ubvb_cr_i);
XLc = XLL(ubvb_cr_i); YLc = YLL(ubvb_cr_i); % This is the first filter to select zeros

% Now get the eigen values
% [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V 
% whose columns are the corresponding right eigenvectors, so that A*V = V*D.

l2e=length(ubx_x2);
ubx1 = zeros(1,l2e); ubx2 = zeros(1,l2e); uby1 = zeros(1,l2e); uby2 = zeros(1,l2e); 
EIG=cell(2,l2e);

%UBC1 = zeros(size(ubx)); VBC1 = zeros(size(ubx));
%UBC2 = zeros(size(ubx)); VBC2 = zeros(size(ubx));
%EIGV1= zeros(size(ubx)); EIGV2= zeros(size(ubx));

UBC1 = NaN(size(ubx)); VBC1 = NaN(size(ubx));
UBC2 = NaN(size(ubx)); VBC2 = NaN(size(ubx));
EIGV1= NaN(size(ubx)); EIGV2= NaN(size(ubx));
ix = floor(XLc/resx); iy = floor(YLc/resy);

% Get the eigen values and vectors form the Jacobian
%--------------------------------------------------------------------------
for i=1:l2e
    %
    Jac = [ubx_x2(i), ubx_y2(i);vby_x2(i), vby_y2(i)];
    [Vast,Dast] = eig(Jac);
    EIG{1,i} = Vast;
    EIG{2,i} = Dast;
    vec1=EIG{1,i};
    eigenv=EIG{2,i};
    %
    ubx1_i=vec1(1,1); ubx2_i=vec1(1,2); 
    uby1_i=vec1(2,1); uby2_i=vec1(2,2);
    eigv1_i=eigenv(1,1);  eigv2_i=eigenv(2,2);
    
    %ubx1(1,i) = ubx1_i;     ubx2(1,i) = ubx2_i;
    %uby1(1,i) = uby1_i;     uby2(1,i) = uby2_i;
    
    % This is to get the position on the plane that satisfies that
    % condition
    
    %
    UBC1(ix(i),iy(i)) = ubx1_i;    VBC1(ix(i),iy(i)) = uby1_i;
    UBC2(ix(i),iy(i)) = ubx2_i;    VBC2(ix(i),iy(i)) = uby2_i;     
    EIGV1(ix(i),iy(i)) = eigv1_i;  EIGV2(ix(i),iy(i)) = eigv2_i;   
    %
end   
%--------------------------------------------------------------------------

% Now put the condition for saddle points with the eigen values:
% Real, opposite signs
%--------------------------------------------------------------------------
saddle_cond1 =  ((imag(EIGV1)==0) & (imag(EIGV2)==0)) ; % complex part equal zero
saddle_cond2 =  ( sign(real(EIGV1))~=sign(real(EIGV2)) ) ; % opposite sign

% This is spiral source 
%saddle_cond1 =  ((imag(EIGV1)~=0) & (imag(EIGV2)~=0)) ; % Complex
%saddle_cond2 =  ( (sign(real(EIGV1)) > 0) & (sign(real(EIGV2)) > 0 ) ) ; % real part positive sign

%saddle_cond12 = (saddle_cond1 & saddle_cond2); 
saddle_cond12 = EIGV1.*(saddle_cond1).* (saddle_cond2); % this allows to use the NaN conversion
A12=saddle_cond12;  A12(A12 == 0) = NaN;  %rrr=ubx(~isnan(A12));size(rrr);

% Get the coordinates of the points where the conditions are fullfilled
XLc2 = XLL(~isnan(A12)); YLc2 = YLL(~isnan(A12));  

% Get the streamlines of the B startting only at the saddle points
%ccc1=stream2(XLL,YLL,ubx,vby,XLc2,YLc2); %This is transported for some reason<------------
ccc12=stream2(XLL,YLL,ubx,vby,YLc2,XLc2); % This is the one that matchs the points <--------------
ccc12m=stream2(XLL,YLL,-ubx,-vby,YLc2,XLc2); 

% This doesn't work. There is no formation of streamlines 
%-------------------------------------------------
% Using the vectors of eigenvectors 
%ubx11=UBC2; vby11=VBC2;
%ccc12=stream2(XLL,YLL,ubx11,vby11,YLc2,XLc2);
%-------------------------------------------------

%--------------------------------------------------------------------------

% Plots
%colorFVBl=[0.5 0.5 0.5 0.5];
colorFVBl=[0 0 0 0.5];

% This is the best plot (This is not working and I don't know why)
f299=figure(299);
dum_p = (real(EIGV1).*saddle_cond12)./real(EIGV1); dum_p(dum_p == 0) = NaN; % This is to erase the unimportant points
%dum_p = (imag(EIGV1).*saddle_cond12)./imag(EIGV1); dum_p(dum_p == 0) = NaN; % This is to erase the unimportant points
dum_p2 = sqrt(ubx.^2 + vby.^2); %vcounl = [zero_d];
[M,c] =contour(xll,yll,dum_p2,zero_d); c.LineWidth = 1.5;
hold on 
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), ubx(1:sv:end,1:sv:end), vby(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'm');
%hc=pcolor(YLL,XLL,dum_p);
hc=pcolor(YLL,XLL,real(dum_p));
set(hc,'edgecolor','none'); colormap(jet(1))
hlines1=streamline(ccc12);set(hlines1,'LineWidth',1.0,'Color', 'g');% [0.2,0.1,0.5]);
hlines1=streamline(ccc12m);set(hlines1,'LineWidth',1.0,'Color', 'g');
hlines2=streamline(aaa);set(hlines2,'LineWidth',1.0, 'Color', colorFVBl);%colorFVBl);
hlines3=streamline(aaam);set(hlines3,'LineWidth',1.0, 'Color', colorFVBl);
hlines2=streamline(bbbb5);set(hlines2,'LineWidth',1.0, 'Color', colorFVBl);%colorFVBl);
hlines3=streamline(bbbb5m);set(hlines3,'LineWidth',1.0, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
%title('$$ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
title('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Separatrices$$','Interpreter','latex','FontSize',20)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
hold off

%clearvars ccc12 ccc12m aaa aaam bbbb5 bbbb5m 

cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f299,strcat('Bl_separatrices_2_'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
clf(f299)


end

%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% This plot is to see how the eigenvalues are distributed along the 2D plane
%--------------------------------------------------------------------------
%{
%--------------------------------------------------------------------------
f99=figure(99)
subplot(2,2,1)
hc=pcolor(real(EIGV1)'); set(hc,'edgecolor','none'); colorbar; title('real1')
subplot(2,2,2); dum_p = real(EIGV1).*imag(EIGV1);
hc=pcolor(dum_p'); set(hc,'edgecolor','none'); colorbar;  title('imag1')
%hc=pcolor(imag(EIGV1)'); set(hc,'edgecolor','none'); colorbar;  title('imag1')
subplot(2,2,3)
hc=pcolor(real(EIGV2)'); set(hc,'edgecolor','none'); colorbar; title('real2')
subplot(2,2,4); dum_p = real(EIGV2).*imag(EIGV2);
hc=pcolor(dum_p'); set(hc,'edgecolor','none'); colorbar;  title('imag1')
%hc=pcolor(imag(EIGV2)'); set(hc,'edgecolor','none'); colorbar;  title('imag2')

f199=figure(199)
subplot(1,3,1)
dum_p = (real(EIGV1).*saddle_cond1)./real(EIGV1); % This is to erase the unimportant points
%dum_p(dum_p == 0) = NaN;
hc=pcolor(dum_p'); set(hc,'edgecolor','none'); colorbar; title('cond1'); colormap(jet(2))
subplot(1,3,2)
dum_p = (real(EIGV1).*saddle_cond2)./real(EIGV1); % This is to erase the unimportant points
%dum_p(dum_p == 0) = NaN;
hc=pcolor(dum_p'); set(hc,'edgecolor','none'); colorbar; title('cond2')
subplot(1,3,3)
dum_p = (real(EIGV1).*saddle_cond12)./real(EIGV1); % This is to erase the unimportant points
dum_p(dum_p == 0) = NaN;
%dum_p = (real(EIGV1).*A)./real(EIGV1); % This is to erase the unimportant points
hc=pcolor(dum_p'); set(hc,'edgecolor','none'); colorbar; title('cond12')

% This is the best plot
f299=figure(299)
dum_p = (real(EIGV1).*saddle_cond12)./real(EIGV1); dum_p(dum_p == 0) = NaN; % This is to erase the unimportant points
dum_p2 = sqrt(ubx.^2 + vby.^2); %vcounl = [zero_d];
[M,c] =contour(xll,yll,dum_p2,zero_d); 
c.LineWidth = 1.5;
hold on 
hc=pcolor(YLL,XLL,dum_p);
set(hc,'edgecolor','none'); colormap(jet(1))
hlines1=streamline(ccc12);
set(hlines1,'LineWidth',1,'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
%title('$$ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
title('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Separatrices$$','Interpreter','latex','FontSize',20)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
hold off
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%ubx_x1 = ubx_x.*ubvb_cr_i; ubx_y1 = ubx_y.*ubvb_cr_i;  % this in matrix representation
%vby_x1 = vby_x.*ubvb_cr_i; vby_y1 = vby_y.*ubvb_cr_i;
%XLc = XLL.*ubvb_cr_i; YLc = YLL.*ubvb_cr_i;
%--------------------------------------------------------------------------

bbb1=stream2(XLL,YLL,real(UBC1),real(VBC1),XLc,YLc);
bbb2=stream2(XLL,YLL,imag(UBC1),imag(VBC1),XLc,YLc);
bbb3=stream2(XLL,YLL,real(UBC2),real(VBC2),XLc,YLc);
bbb4=stream2(XLL,YLL,imag(UBC2),imag(VBC2),XLc,YLc);

%--------------------------------------------------------------------------
figure213=figure(213);
%hold on
hlines1=streamline(bbb1);
set(hlines1,'LineWidth',1,'Color', 'b');
hold on
hlines2=streamline(bbb2);
set(hlines2,'LineWidth',1,'Color', 'r');
hlines3=streamline(bbb3);
set(hlines3,'LineWidth',1,'Color', 'c');
hlines4=streamline(bbb4);
set(hlines4,'LineWidth',1,'Color', 'm');
%hlines5=streamline(aaa5);
hlines6=streamline(aaa);
set(hlines6,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
%title('$$ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
title('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Blines$$','Interpreter','latex','FontSize',20)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
hold off
%--------------------------------------------------------------------------

%bbb5=stream2(XLc,YLc,ubx1,uby1,XLc,YLc);
% This works so far   

% Now the question is how to plot the tangent lines relates to the
% eigenvectors
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
figure212=figure(212);
h1=subplot(1,3,1);
dum_p=ubx;
hc=pcolor(XLL,YLL,dum_p); set(hc,'edgecolor','none')
lim_yp=max(max(max(abs(dum_p))));  caxis([-lim_yp lim_yp]);
colorbar
hold on
hlines=streamline(aaa);
hlines2=streamline(bbb5);
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
set(hlines2,'LineWidth',1.3, 'Color', 'm');
xlim([xll(1) xll(end)]);ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
hold off; colormap(h1,BWR);
h1=subplot(1,3,2);
%hc=pcolor(ubx_cr_nan); 
hc=pcolor(XLL,YLL,ones(size(ubx)).*(ubx_cr_i) -2*ones(size(ubx)).*vby_cr_i );
%hc=pcolor(XLL,YLL,vby.*vby_cr_i);
%hc=pcolor(XLL,YLL,ubx.*ubvbwb_cr_i);
%lim_yp=max(max(max(abs(dum_p))));  
%caxis(-2, -1, 0, 1);
set(hc,'edgecolor','none')
colorbar
hold on
hlines=streamline(aaa);
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
colormap(h1,jet(4));
hold off

h1=subplot(1,3,3);
%f123=figure(123)
%hc=pcolor(XLL,YLL,ones(size(ubx)).* ubx_cr_nan); 
%hc=pcolor(XLL,YLL,ones(size(ubx)).*(ubx_cr_i) -2*ones(size(ubx)).*vby_cr_i );
%hc=pcolor(XLL,YLL,ones(size(ubx)).*vby_cr_i);
hc=pcolor(XLL,YLL,ones(size(ubx)).*ubvb_cr_i);
%lim_yp=max(max(max(abs(dum_p))));  
%caxis(-2, -1, 0, 1);
set(hc,'edgecolor','none')
colorbar
hold on
hlines=streamline(aaa);
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
colormap(h1,jet(2));
hold off
%--------------------------------------------------------------------------
min(min(abs(ubx)))
mean(abs(ubx),'all','omitnan')
mean(std(abs(ubx),'omitnan'))
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%Plots on the slide. The velocities ions, electrons, current and magnetic
%field components
%--------------------------------------------------------------------------
%{
f1=figure(1);
for i=1:12
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(4,3,i,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=vi_Gp{1,i};
    titl = 'ui';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=ve_Gp{1,i2};
    titl = 'ue';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=J_Gp{1,i3};
    titl = 'J';
    s=i3;
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=B_Gp{1,i4};
    if (i4==2)
        dum_paux=B_Gp{1,i4};
        dum_p = dum_paux ;%+ 0.1;
    end
    titl = 'B';
    s=i4;    
end
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
hcb=colorbar;
if(i==11)
    colormap(h1,cool);
    %caxis([-lim_yp 0])
else
    colormap(h1,BWR);
    %caxis([-lim_yp lim_yp])
end
if (i==1 || i==2 || i==3)
    caxis([-0.07 0.07])
elseif (i==4 || i==6 || i==7 || i==9)
    caxis([-0.3 0.3])
elseif (i==5 || i==8)
    caxis([-1 1])  
elseif (i==10 || i==12)
    caxis([-0.15 0.15]) 
elseif (i==11)
    caxis([-0.17 0])
end
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex'; %hcb.Ruler.Exponent = -2;
%hcb.Title.String='';
%set(get(hcb,'Title'),'String',titl + string(s) ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    if (i==1)
        ax.YLabel.String = '$p / d_{i}$, $\mathbf{u}_{i}$';
    elseif (i==4)
        ax.YLabel.String = '$p / d_{i}$, $\mathbf{u}_{e}$';
    elseif (i==7)
        ax.YLabel.String = '$p / d_{i}$, $\mathbf{J}$';
    elseif (i==10)
        ax.YLabel.String = '$p / d_{i}$, $\mathbf{B}$';
    end
    ax.YTick = [2 4 6 8 10];
end
if (9<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];    
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
%set(h2,'AutoScale','on', 'AutoScaleFactor', 2, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ (p,a,r)$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
%
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f1,strcat('ui_ue_J_B_Bl_'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --->   Ba, uep uer  <-----
%{
f13=figure(13);
h1=subaxis(1,3,1,'SV',0.01,'SH',0.004,'MR',0.01,'ML',0.05,'PL',0.01,'PR',0.006);
dum_p=B_Gp{1,2};  
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h1,cool); caxis([-0.17 0])
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; 
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$B_{a}$','Interpreter','latex','FontSize',fs);
hold off
h2=subaxis(1,3,2,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=ve_Gp{1,1};    
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h2,BWR); caxis([-0.28 0.28]) 
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$u_{e,p}$','Interpreter','latex','FontSize',fs);
hold off
h3=subaxis(1,3,3,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=ve_Gp{1,3};    
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h3,BWR); caxis([-0.28 0.28]) 
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$u_{e,r}$','Interpreter','latex','FontSize',fs);
hold off
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',fs)

% --->   Ba, uip uir  <-----
%--------------------------------------------------------------------------
f1310=figure(1310);
h1=subaxis(1,3,1,'SV',0.01,'SH',0.004,'MR',0.01,'ML',0.05,'PL',0.01,'PR',0.006);
%dum_p=B_Gp{1,2};
dum_p=vi_Gp{1,2};
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h1,BWR); caxis([-0.1 0.1]) ; % colormap(h1,cool); caxis([-0.17 0])
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; 
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
%title('$B_{a}$','Interpreter','latex','FontSize',fs);
title('$u_{i,a}$','Interpreter','latex','FontSize',fs);
hold off
h2=subaxis(1,3,2,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=vi_Gp{1,1};    
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h2,BWR); caxis([-0.1 0.1]) 
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$u_{i,p}$','Interpreter','latex','FontSize',fs);
hold off
h3=subaxis(1,3,3,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=vi_Gp{1,3};    
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h3,BWR); caxis([-0.1 0.1]) 
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$u_{i,r}$','Interpreter','latex','FontSize',fs);
hold off
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',fs)
%--------------------------------------------------------------------------




%}
%--------------------------------------------------------------------------

% This is a plot that shows the distribution of the poynting theorem terms 
%--------------------------------------------------------------------------
%{
f144=figure(144);
subplot(1,4,1)
histogram(dtE_em_in_Gp_2{1,1},100)
xlim([-0.02, 0.02])
subplot(1,4,2)
histogram(Div_Poyn_v_Gp{1,1},100)
xlim([-0.02, 0.02])
subplot(1,4,3)
histogram(JdotE_Gp{1,1},100)
xlim([-0.02, 0.02])
subplot(1,4,4)
histogram(dtE_em_in_Gp_2{1,1}+Div_Poyn_v_Gp{1,1}+JdotE_Gp{1,1},100)
xlim([-0.02, 0.02])
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --->  dtE_em_in_Gp_2  Div_Poyn_v_Gp, JdotE    <----- plot of
%--------------------------------------------------------------------------
%{
f14=figure(14);
h1=subaxis(1,3,1,'SV',0.01,'SH',0.004,'MR',0.01,'ML',0.05,'PL',0.01,'PR',0.006);
dum_p=dtE_em_in_Gp_2{1,1};  
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h1,BWR); 
caxis([-6e-4 6e-4])
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -4;
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$d \varepsilon^{em}/dt$','Interpreter','latex','FontSize',fs);
hold off
h2=subaxis(1,3,2,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=Div_Poyn_v_Gp{1,1};    
Kd = medfilt2(dum_p,[5,5]); %filter
hc = pcolor(XLL,YLL,Kd);
%hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h2,BWR); 
caxis([-0.01 0.01]) 
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -2;
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$\nabla \cdot S$','Interpreter','latex','FontSize',fs);
hold off
h3=subaxis(1,3,3,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=JdotE_Gp{1,1};    
Kd = medfilt2(dum_p,[5,5]); %filter
hc = pcolor(XLL,YLL,Kd);
%hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=max(max(max(abs(dum_p)))); colormap(h3,BWR); 
caxis([-0.001 0.001]) 
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -2;
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$J \cdot E$','Interpreter','latex','FontSize',fs);
hold off
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',fs)
%}
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% --->  E components  <----- plot of
%--------------------------------------------------------------------------
%{
f1499=figure(1499);
h1=subaxis(1,3,1,'SV',0.01,'SH',0.004,'MR',0.01,'ML',0.05,'PL',0.01,'PR',0.006);
dum_p=E_Gp{1,1};  
hc = pcolor(XLL,YLL,dum_p);
%Kd = medfilt2(dum_p,[5,5]); %filter
%hc = pcolor(XLL,YLL,Kd);
hold on
hlinesbb5=streamline(bbbb5); set(hlinesbb5,'LineWidth',1, 'Color', colorFVBl);
hlinesbb5=streamline(bbbb5m); set(hlinesbb5,'LineWidth',1, 'Color', colorFVBl);
lim_yp=max(max(max(abs(dum_p)))); colormap(h1,BWR); 
%caxis([-6e-4 6e-4])
caxis([-lim_yp lim_yp]);
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -2;
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
%hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$E_{p}$','Interpreter','latex','FontSize',fs);
hold off

h2=subaxis(1,3,2,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=E_Gp{1,2};    
Kd = medfilt2(dum_p,[5,5]); %filter
%hc = pcolor(XLL,YLL,Kd);
hc = pcolor(XLL,YLL,dum_p);
hold on
hlinesbb5=streamline(bbbb5); set(hlinesbb5,'LineWidth',1, 'Color', colorFVBl);
hlinesbb5=streamline(bbbb5m); set(hlinesbb5,'LineWidth',1, 'Color', colorFVBl);
lim_yp=max(max(max(abs(dum_p)))); colormap(h2,BWR); 
%caxis([-0.01 0.01]) 
caxis([-lim_yp lim_yp]);
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -2;
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
%hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$E_{a}$','Interpreter','latex','FontSize',fs);
hold off

h3=subaxis(1,3,3,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
dum_p=E_Gp{1,3};    
Kd = medfilt2(dum_p,[5,5]); %filter
%hc = pcolor(XLL,YLL,Kd);
hc = pcolor(XLL,YLL,dum_p);
hold on
hlinesbb5=streamline(bbbb5); set(hlinesbb5,'LineWidth',1, 'Color', colorFVBl);
hlinesbb5=streamline(bbbb5m); set(hlinesbb5,'LineWidth',1, 'Color', colorFVBl);
lim_yp=max(max(max(abs(dum_p)))); colormap(h3,BWR); 
%caxis([-0.001 0.001]) 
caxis([-lim_yp lim_yp]);
set(hc,'edgecolor','none')
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -2;
hcb.Label.FontSize=fs; hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=fs; ax.XTick = []; ax.YTick = [];
%hcb.Ruler.Exponent = -1;
%ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
%hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart); set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
title('$E_{r}$','Interpreter','latex','FontSize',fs);
hold off
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',fs)
%}



% to save the plots of Ba, uep uer and dtE_em_in_Gp_2  Div_Poyn_v_Gp, JdotE
%--------------------------------------------------------------------------
%{
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f13,strcat('Ba_uep_uer_'+ string(N_steps) +'.png'));
saveas(f14,strcat('Poynting_'+ string(N_steps) +'.png'));
saveas(f1499,strcat('E_components_'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

clf(f1); clf(f13); clf(f14)
%--------------------------------------------------------------------------
%}
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
dCp=25;
%--------------------------------------------------------------------------
% This is the plot for the energy terms in the kinetic equation
%--------------------------------------------------------------------------

% 2D plots for the kinetic energy rate terms
%{
%--------------------------------------------------------------------------
GPijf{1,1}=dtKe_in_Gp_2{1,1}; GPijf{1,2}=ve_grad_Ke_Gp{1,1};
GPijf{1,3}=Ke_Div_ve_Gp{1,1}; GPijf{1,4}=ve_Div_Pet_Gp{1,1}; 
GPijf{1,5}=-qe_ne_ve_E_Gp{1,1};  
f4=figure(4);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
h1=subaxis(1,5,1,'SV',0,'SH',0.0,'MR',0,'ML',0,'PL',0.025,'PR',0);
dum_p=GPijf{1,1}; %dum_p2 = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-5 5e-5]); %caxis([-lim_yp lim_yp]); 
colormap(h1,BWR); set(hc,'edgecolor','none'); 
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -2;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs);
ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs);
xlh.Position(2) = xlh.Position(2) + 0.4; ylh.Position(1) = ylh.Position(1) + 0.4; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%xlh.Position(2) = xlh.Position(2) + abs(xlh.Position(2) * 0.1);
title('$\partial \varepsilon^{k}_{e} / \partial t $','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h2=subaxis(1,5,2,'SV',0,'SH',0.000,'MR',0,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,2}; 
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-4 5e-4]); %caxis([-lim_yp lim_yp]); 
colormap(h2,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex'; ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$u_{e} \cdot \nabla \varepsilon^{k}_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h3=subaxis(1,5,3,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,3}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-4 5e-4]); %caxis([-lim_yp lim_yp]); 
colormap(h3,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$\varepsilon^{k}_{e}\nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h4=subaxis(1,5,4,'SV',0,'SH',0.00,'MR',0,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,4}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h4,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$u_{e} \cdot \nabla \cdot \overline{P}_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h5=subaxis(1,5,5,'SV',0,'SH',0.00,'MR',0.001,'ML',0,'PL',0.02,'PR',0.005);
dum_p=GPijf{1,5}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h5,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$-q_{e}n_{e}E\cdot u_{e}$','Interpreter','latex','FontSize',fs);
hold off
clearvars GPijf
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 2D plots for the Internal energy rate terms
%--------------------------------------------------------------------------
GPijf{1,1}=dtUe_in_Gp_2{1,1}; GPijf{1,2}=ve_grad_Ue_Gp{1,1};
GPijf{1,3}=Ue_Div_ve_Gp{1,1}; GPijf{1,4}=Pet_grad_ve_Gp{1,1}; 
GPijf{1,5}=dkQiik_05_e_Gp{1,1};  

f5=figure(5);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
h1=subaxis(1,5,1,'SV',0,'SH',0.0,'MR',0,'ML',0,'PL',0.025,'PR',0);
dum_p=GPijf{1,1}; %dum_p2 = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-4 5e-4]); %caxis([-lim_yp lim_yp]); 
colormap(h1,BWR); set(hc,'edgecolor','none'); 
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -2;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs);
ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs);
xlh.Position(2) = xlh.Position(2) + 0.4; ylh.Position(1) = ylh.Position(1) + 0.4; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%xlh.Position(2) = xlh.Position(2) + abs(xlh.Position(2) * 0.1);
title('$\partial \varepsilon^{th}_{e} / \partial t $','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h2=subaxis(1,5,2,'SV',0,'SH',0.000,'MR',0,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,2}; 
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h2,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex'; ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$u_{e} \cdot \nabla \varepsilon^{th}_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h3=subaxis(1,5,3,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,3}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-2 5e-2]);  %caxis([-lim_yp lim_yp]); 
colormap(h3,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -2;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$\varepsilon^{th}_{e}\nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h4=subaxis(1,5,4,'SV',0,'SH',0.00,'MR',0,'ML',0,'PL',0.02,'PR',0);
dum_p=GPijf{1,4}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-2 5e-2]);  %caxis([-lim_yp lim_yp]); 
colormap(h4,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -2;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$\overline{P}_{e} : \nabla u_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h5=subaxis(1,5,5,'SV',0,'SH',0.00,'MR',0.001,'ML',0,'PL',0.02,'PR',0.005);
dum_p=GPijf{1,5}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-4 5e-4]); %caxis([-lim_yp lim_yp]); 
colormap(h5,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$\nabla \cdot h_{e}$','Interpreter','latex','FontSize',fs);
hold off
clearvars GPijf
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 2D plots for the diagonal and off diagonal
%--------------------------------------------------------------------------
%
GPijf{1,1}=Pet_grad_ve_Gp{1,1}; 
GPijf{1,2}=Pet_grad_ve_ii_Gp{1,1};
GPijf{1,3}=Pet_grad_ve_ij_Gp{1,1}; 
GPijf{1,4}=Pij_e_Gp{3,3} + Pij_e_Gp{1,1} + Pij_e_Gp{2,2}; 
GPijf{1,5}=Pij_e_Gp{1,3} + Pij_e_Gp{1,2} + Pij_e_Gp{2,3} ;  
GPijf{1,6}=grad_ve_Gp{3,3} + grad_ve_Gp{1,1} + grad_ve_Gp{2,2}; 
GPijf{1,7}=grad_ve_Gp{1,3} + grad_ve_Gp{1,2} + grad_ve_Gp{2,3} ;

f6=figure(6);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
h1=subaxis(1,7,1,'SV',0,'SH',0.0,'MR',0,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,1}; %dum_p2 = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h1,BWR); set(hc,'edgecolor','none'); 
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -2;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs);
ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs);
xlh.Position(2) = xlh.Position(2) + 0.4; ylh.Position(1) = ylh.Position(1) + 0.4; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%xlh.Position(2) = xlh.Position(2) + abs(xlh.Position(2) * 0.1);
title('$\overline{P}_{e} : \nabla u_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h2=subaxis(1,7,2,'SV',0,'SH',0.000,'MR',0,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,2}; 
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); 
%caxis([-lim_yp lim_yp]); 
colormap(h2,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex'; ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$p\theta_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h3=subaxis(1,7,3,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,3}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-4 5e-4]);  %caxis([-lim_yp lim_yp]); 
colormap(h3,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -4;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$Pi-D_{e}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h4=subaxis(1,7,4,'SV',0,'SH',0.00,'MR',0,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,4}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([0 2e-2]);  %caxis([-lim_yp lim_yp]); 
colormap(h4,jet); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -3;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${P}_{e,ii}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h5=subaxis(1,7,5,'SV',0,'SH',0.00,'MR',0.00,'ML',0,'PL',0.0,'PR',0.00);
dum_p=GPijf{1,5}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h5,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${P}_{e,ij}|_{i \neq j}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h4=subaxis(1,7,6,'SV',0,'SH',0.00,'MR',0,'ML',0,'PL',0.00,'PR',0);
dum_p=GPijf{1,6}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-2 2]);  %caxis([-lim_yp lim_yp]); 
colormap(h4,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -3;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${\nabla u}_{e,ii}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h5=subaxis(1,7,7,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.00,'PR',0.000);
dum_p=GPijf{1,7}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-2 2]); %caxis([-lim_yp lim_yp]); 
colormap(h5,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${\nabla u}_{e,ij}|_{i \neq j}$','Interpreter','latex','FontSize',fs);
hold off

clearvars GPijf

%--------------------------------------------------------------------------
%
% 2D plots for the diagonal and off diagonal for the ions
%--------------------------------------------------------------------------
%
GPijf{1,1}=Pit_grad_vi_Gp{1,1}; 
GPijf{1,2}=Pit_grad_vi_ii_Gp{1,1};
GPijf{1,3}=Pit_grad_vi_ij_Gp{1,1}; 
GPijf{1,4}=Pij_i_Gp{3,3} + Pij_i_Gp{1,1} + Pij_i_Gp{2,2}; 
GPijf{1,5}=Pij_i_Gp{1,3} + Pij_i_Gp{1,2} + Pij_i_Gp{2,3} ;  
GPijf{1,6}=grad_vi_Gp{3,3} + grad_vi_Gp{1,1} + grad_vi_Gp{2,2}; 
GPijf{1,7}=grad_vi_Gp{1,3} + grad_vi_Gp{1,2} + grad_vi_Gp{2,3} ;

f61=figure(61);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
h1=subaxis(1,7,1,'SV',0,'SH',0.0,'MR',0,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,1}; %dum_p2 = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h1,BWR); set(hc,'edgecolor','none'); 
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -3;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs);
ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs);
xlh.Position(2) = xlh.Position(2) + 0.4; ylh.Position(1) = ylh.Position(1) + 0.4; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%xlh.Position(2) = xlh.Position(2) + abs(xlh.Position(2) * 0.1);
title('$\overline{P}_{i} : \nabla u_{i}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h2=subaxis(1,7,2,'SV',0,'SH',0.000,'MR',0,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,2}; 
hc = pcolor(XLL,YLL,dum_p); 
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); 
%caxis([-lim_yp lim_yp]); 
colormap(h2,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -3;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex'; ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$p\theta_{i}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h3=subaxis(1,7,3,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,3}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-4 5e-4]);  %caxis([-lim_yp lim_yp]); 
colormap(h3,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; hcb.Ruler.Exponent = -4;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('$Pi-D_{i}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h4=subaxis(1,7,4,'SV',0,'SH',0.00,'MR',0,'ML',0,'PL',0.0,'PR',0);
dum_p=GPijf{1,4}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
%caxis([0 1e-2]);  %caxis([-lim_yp lim_yp]); 
colormap(h4,jet); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -3;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${P}_{i,ii}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h5=subaxis(1,7,5,'SV',0,'SH',0.00,'MR',0.00,'ML',0,'PL',0.0,'PR',0.00);
dum_p=GPijf{1,5}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-5e-3 5e-3]); %caxis([-lim_yp lim_yp]); 
colormap(h5,BWR); set(hc,'edgecolor','none'); hcb.Ruler.Exponent = -3;
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${P}_{i,ij}|_{i \neq j}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h4=subaxis(1,7,6,'SV',0,'SH',0.00,'MR',0,'ML',0,'PL',0.00,'PR',0);
dum_p=GPijf{1,6}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-1 1]);  %caxis([-lim_yp lim_yp]); 
colormap(h4,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside'; %hcb.Ruler.Exponent = -3;
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${\nabla u}_{i,ii}$','Interpreter','latex','FontSize',fs);
hold off
%--------------------------------------------------------------------------
h5=subaxis(1,7,7,'SV',0,'SH',0.00,'MR',0.000,'ML',0,'PL',0.00,'PR',0.000);
dum_p=GPijf{1,7}; 
hc = pcolor(XLL,YLL,dum_p);
hold on
lim_yp=0.7*max(max(max(abs(dum_p)))); clearvars dum_p;
caxis([-1 1]); %caxis([-lim_yp lim_yp]); 
colormap(h5,BWR); set(hc,'edgecolor','none');
hcb=colorbar; hcb.Location = 'northoutside';
hcb.Label.Interpreter = 'latex';hcb.TickLabelInterpreter = 'latex';
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';ax.YTick = [];
xlh=xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs); xlh.Position(2) = xlh.Position(2) + 0.4;
%ylh=ylabel('$p / d_{i}$','Interpreter','latex','FontSize',fs); ylh.Position(1) = ylh.Position(1) + 0.3; 
hlines=streamline(aaa);set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
title('${\nabla u}_{i,ij}|_{i \neq j}$','Interpreter','latex','FontSize',fs);
hold off

clearvars GPijf
%--------------------------------------------------------------------------

cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(4,strcat('2D_kinetic_energy_rates_'+ string(N_steps) +'.png'));
%saveas(5,strcat('2D_Internal_energy_rates_'+ string(N_steps) +'.png'));
saveas(6,strcat('2D_diag_offdiag_'+ string(N_steps) +'.png'));
clf(f4); clf(f5); 
clf(f6)
%}
%--------------------------------------------------------------------------

cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here'

%--------------------------------------------------------------------------
%}
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%To make the plots once I have loaded the data. This has the field lines
%and vectors plots
%--------------------------------------------------------------------------
%
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%--------------------------------------------------------------------------
% Plot the field over the slide and the vectors in the NEW frame. 
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here'
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%}

%f79=figure(79)
%contourf(CU1,CU2,CUZ)
%hold on
%plot(ccc(:,1) , ccc(:,2), '*r')
%k=5;
%eps=1e-1;
%result =find(xll(10)-ccc(k,1) < eps);

sv=3;
% This are plots to check that the transformation is (x,y,z) --> (p,a,r)
%--------------------------------------------------------------------------
%{

%--------------------------------------------------------------------------
Vev_x = ve_Gp{1,1}; Vev_y = ve_Gp{1,2}; Vev_z = ve_Gp{1,3};
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x_1=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y_1=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = ve_Gp{1,2}; % <----- this seems to be the one that works
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x_2=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y_2=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

Vev_x = ve_Gp{1,3}; Vev_y = ve_Gp{1,1}; Vev_z = -ve_Gp{1,2}; 
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x_3=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y_3=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

%-------------------------------
ubx = B_Gp{1,1};   vby = B_Gp{1,2}; wbz = B_Gp{1,3};
bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));
aaa_1=stream2(XLL,YLL,ubx,vby,xstart,ystart);

ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};
bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));
aaa_2=stream2(XLL,YLL,ubx,vby,xstart,ystart);

ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2}; %Is it this way? it doesn't look that way
bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));
aaa_3=stream2(XLL,YLL,ubx,vby,xstart,ystart);
%--------------------------------------------------------------------------

f110=figure(110);
h1=subplot(1,3,1);
dum_p=B_Gp{1,3}; hc=pcolor(dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none');
title('v1=vr')
colorbar; colormap(h1,BWR)
h2=subplot(1,3,2);
dum_p=B_Gp{1,1}; hc=pcolor(dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none');
title('v2=vp')
colorbar; colormap(h2,BWR);
h3=subplot(1,3,3);
dum_p=-B_Gp{1,2};
hc=pcolor(dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none');
title('v3=va')
colorbar; colormap(h3,BWR);


f111=figure(111);
h1=subplot(1,3,1);
dum_p=ve_Gp{1,3}; hc=pcolor(dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none');
title('v1=vr')
colorbar; colormap(h1,BWR)
h2=subplot(1,3,2);
dum_p=ve_Gp{1,1}; hc=pcolor(dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none');
title('v2=vp')
colorbar; colormap(h2,BWR);
h3=subplot(1,3,3);
dum_p=ve_Gp{1,2};
hc=pcolor(dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none');
title('v3=va')
colorbar; colormap(h3,BWR);

asf=3;
sv=3;

f112=figure(112);
h1=subplot(1,3,1);
hlines=streamline(aaa_1);
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x_1(1:sv:end,1:sv:end), Vev_y_1(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'b');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$v_{e,1} \ vectors \ and \ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
hold off
h2=subplot(1,3,2);
hlines=streamline(aaa_2);
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x_2(1:sv:end,1:sv:end), Vev_y_2(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'b');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$v_{e,2} \ vectors \ and \ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
hold off
h3=subplot(1,3,3);
hlines=streamline(aaa_3);
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x_3(1:sv:end,1:sv:end), Vev_y_3(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'b');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$v_{e,3} \ vectors \ and \ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
hold off
%}
%--------------------------------------------------------------------------


% These are the plots of the magnetic field lines and ve, vi and B vectors
% with time
%--------------------------------------------------------------------------
%
asf=3;
sv=3;
f1131=figure(1131);
%h1=subplot(1,3,1);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Blines$$','Interpreter','latex','FontSize',20)
h1=subaxis(1,3,1,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
hlines=streamline(aaa);
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'b');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$u_{e} \ vectors$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hold off
%h2=subplot(1,3,2);
h2=subaxis(1,3,2,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
hlines=streamline(aaa);
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Viv_x(1:sv:end,1:sv:end), Viv_y(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'r');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$u_{i} \ vectors$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.YLabel.String = '$''$'; ax.YTick = [];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hold off
%h3=subplot(1,3,3);
h3=subaxis(1,3,3,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
hlines=streamline(aaa);
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), ubx(1:sv:end,1:sv:end), vby(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'm');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$B \ vectors$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.YLabel.String = '$''$'; ax.YTick = [];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hold off
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f1131,strcat('Blines_ve_vi_B'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
clf(f1131)


%--------------------------------------------------------------------------
%}
%--------------------------------------------------------------------------

% These are the plots of the magnetic field lines, ve, vi, B lines 
% with time
%--------------------------------------------------------------------------
%
NN5=30; 
xx5 = linspace(xll(10),xll(end-10),NN5);
yy5 = linspace(yll(10),yll(end-10),NN5);
[xx5, yy5] = meshgrid(xx5,yy5);

aaa=stream2(XLL,YLL,ubx,vby,xstart,ystart);
aaam=stream2(XLL,YLL,-ubx,-vby,xstart,ystart);

%veve=stream2(XLL,YLL,Vev_x,Vev_y,xstart,ystart);
veve5=stream2(XLL,YLL,Vev_x,Vev_y,xx5,yy5);
%vivi=stream2(XLL,YLL,Viv_x,Viv_y,xstart,ystart);
vivi5=stream2(XLL,YLL,Viv_x,Viv_y,xx5,yy5);
bbbb5=stream2(XLL,YLL,ubx,vby,xx,yy);

veve5m=stream2(XLL,YLL,-Vev_x,-Vev_y,xx5,yy5);
vivi5m=stream2(XLL,YLL,-Viv_x,-Viv_y,xx5,yy5);
bbbb5m=stream2(XLL,YLL,-ubx,-vby,xx,yy);


asf=5;
sv=3;
f1141=figure(1141);
%h1=subplot(1,3,1);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Blines$$','Interpreter','latex','FontSize',20)
h1=subaxis(1,3,1,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
hlines=streamline(aaa); set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
hold on
%hlines=streamline(aaam); set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%hlinesve=streamline(veve); set(hlinesve,'LineWidth',1, 'Color', 'b');
hlinesve5=streamline(veve5); set(hlinesve5,'LineWidth',1, 'Color', 'r');
hlinesve5=streamline(veve5m); set(hlinesve5,'LineWidth',1, 'Color', 'b');
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'k');
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$u_{e} \ lines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.YLabel.String = '$p / d_{i}$'; ax.YTick = [2 4 6 8 10];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hold off
%h2=subplot(1,3,2);
h2=subaxis(1,3,2,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
hlines=streamline(aaa);
hold on
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%hlinesvi=streamline(vivi); set(hlinesvi,'LineWidth',1, 'Color', 'r');
hlinesvi5=streamline(vivi5); set(hlinesvi5,'LineWidth',1, 'Color', 'r');
hlinesvi5=streamline(vivi5m); set(hlinesvi5,'LineWidth',1, 'Color', 'b');
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Viv_x(1:sv:end,1:sv:end), Viv_y(1:sv:end,1:sv:end), 0);
%set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'k');
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$u_{i} \ streamlines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.YLabel.String = '$''$'; ax.YTick = [];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hold off
%h3=subplot(1,3,3);
h3=subaxis(1,3,3,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
hlines=streamline(aaa);
hold on
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), ubx(1:sv:end,1:sv:end), vby(1:sv:end,1:sv:end), 0);
%set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'm');
hlinesbb5=streamline(bbbb5); set(hlinesbb5,'LineWidth',1, 'Color', 'm');
hlinesbb5=streamline(bbbb5m); set(hlinesbb5,'LineWidth',1, 'Color', 'm');
xlim([xll(1) xll(end)]); ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$B \ lines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.YLabel.String = '$''$'; ax.YTick = [];
ax.XLabel.String = '$r / d_{i}$'; ax.XTick = [2 4 6 8];    
hold off

cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f1141,strcat('lines_B_ve_vi_2_'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
clf(f1141)


%--------------------------------------------------------------------------
%}
%--------------------------------------------------------------------------



% ve vectors and field lines 
%--------------------------------------------------------------------------
% overploting the positive and negative vector streamlines I got the
% backward propagation lines
ubx = B_Gp{1,3};   vby = B_Gp{1,1}; wbz = B_Gp{1,2};

bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

bbbb5=stream2(XLL,YLL,ubx,vby,xx,yy);
bbbb5m=stream2(XLL,YLL,-ubx,-vby,xx,yy);

ubx33 = B_Gp{1,3};
vby33 = B_Gp{1,1};
bbbb533=stream2(XLL,YLL,ubx33,vby33,xx,yy);
bbbb5m33=stream2(XLL,YLL,-ubx33,-vby33,xx,yy);

f79=figure(79);
%hlines=streamline(aaa);
subplot(1,2,1)
hc=contourf(XLL,YLL,vby33);
subplot(1,2,2)
hc=contourf(XLL,YLL,vby);


f78=figure(78);
%hlines=streamline(aaa);
%hc=contourf(XLL,YLL,wbz);
%hold on
hlines=streamline(bbbb5); set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
hold on
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'b');
hlines=streamline(bbbb5m); set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$u_{e} \ vectors \ and \ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
hold off
%}
%--------------------------------------------------------------------------




% Magnetic field lines
%--------------------------------------------------------------------------
%
f77=figure(77);
hlines2=streamline(aaa);
%contourf(sqrt(ubx.^2+vby.^2 + wbz.^2))
hold on
hlines=streamline(aaa);
%contour(wbz); %<------ trying to get the separatrix
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%set(h2,'AutoScale','on', 'AutoScaleFactor', asf, 'Color', 'b');
set(hlines,'LineWidth',1.3, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
%title('$$ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
title('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Blines$$','Interpreter','latex','FontSize',20)
ax = gca; ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
hold off

cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f77,strcat('Bl_slideRF2_per'+ string(N_steps) +'.png'));
%saveas(f78,strcat('ve_vec_Bl_slideRF2_per'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
clf(f77)
%}
%--------------------------------------------------------------------------

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%save('variabletimefile.mat','Vectorsc','Tensorsc','Densitiesc','Poyn_theoc',...
%    'kin_energy_terms_ic','kin_energy_terms_ec',...
%    'Int_energy_terms_ic','Int_energy_terms_ec','dt_energydensities_ic')
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% this part is exclusive to compute the averages in time.
% Once loaded the variabletimefile.mat
%--------------------------------------------------------------------------
dtKe_in_Gp_m=zeros(1,N-2); dtKe_in_Gp_std=zeros(1,N-2);
ve_grad_Ke_Gp_m=zeros(1,N-2); ve_grad_Ke_Gp_std=zeros(1,N-2);
ve_Div_Pet_Gp_m=zeros(1,N-2); ve_Div_Pet_Gp_std=zeros(1,N-2);
Ke_Div_ve_Gp_m=zeros(1,N-2); Ke_Div_ve_Gp_std=zeros(1,N-2);
qe_ne_ve_E_Gp_m=zeros(1,N-2); qe_ne_ve_E_Gp_std=zeros(1,N-2);

dtKi_in_Gp_m=zeros(1,N-2); dtKi_in_Gp_std=zeros(1,N-2);
vi_grad_Ki_Gp_m=zeros(1,N-2); vi_grad_Ki_Gp_std=zeros(1,N-2);
vi_Div_Pit_Gp_m=zeros(1,N-2); vi_Div_Pit_Gp_std=zeros(1,N-2);
Ki_Div_vi_Gp_m=zeros(1,N-2); Ki_Div_vi_Gp_std=zeros(1,N-2);
qi_ni_vi_E_Gp_m=zeros(1,N-2); qi_ni_vi_E_Gp_std=zeros(1,N-2);

dtUe_in_Gp_m = zeros(1,N-2); dtUe_in_Gp_std = zeros(1,N-2);
ve_grad_Ue_Gp_m =  zeros(1,N-2); ve_grad_Ue_Gp_std = zeros(1,N-2);
dkQiik_05_e_Gp_m = zeros(1,N-2); dkQiik_05_e_Gp_std = zeros(1,N-2);
Pet_grad_ve_Gp_m = zeros(1,N-2); Pet_grad_ve_Gp_std = zeros(1,N-2);
Pet_grad_ve_ii_Gp_m = zeros(1,N-2); Pet_grad_ve_ii_Gp_std = zeros(1,N-2);
Pet_grad_ve_ij_Gp_m = zeros(1,N-2); Pet_grad_ve_ij_Gp_std = zeros(1,N-2);
Ue_Div_ve_Gp_m = zeros(1,N-2); Ue_Div_ve_Gp_std = zeros(1,N-2);

dtUi_in_Gp_m = zeros(1,N-2); dtUi_in_Gp_std = zeros(1,N-2);
vi_grad_Ui_Gp_m =  zeros(1,N-2); vi_grad_Ui_Gp_std = zeros(1,N-2);
dkQiik_05_i_Gp_m = zeros(1,N-2); dkQiik_05_i_Gp_std= zeros(1,N-2);
Pit_grad_vi_Gp_m = zeros(1,N-2); Pit_grad_vi_Gp_std = zeros(1,N-2);
Pit_grad_vi_ii_Gp_m = zeros(1,N-2); Pit_grad_vi_ii_Gp_std = zeros(1,N-2);
Pit_grad_vi_ij_Gp_m = zeros(1,N-2); Pit_grad_vi_ij_Gp_std = zeros(1,N-2);
Ui_Div_vi_Gp_m = zeros(1,N-2); Ui_Div_vi_Gp_std = zeros(1,N-2);

for in=2:N-1
%--------------------------------------------------------------------------
disp(strcat('Computing step ...',S(in).name))
%--------------------------------------------------------------------------
%N_steps = 2000;
%time_t = dt*N_steps; %(In terms of 1/omega_pi)
dt = 0.06;
str = string(S(in).name);
newStr = extractBetween(str,"pfd.","_p000000.h5");
N_steps = str2double(newStr); time_t = dt*N_steps;

str_1 = string(S(in-1).name); newStr_1 = extractBetween(str_1,"pfd.","_p000000.h5");
N_steps_1 = str2double(newStr_1); time_1=dt*N_steps_1;
str_2 = string(S(in+0).name); newStr_2 = extractBetween(str_2,"pfd.","_p000000.h5");
N_steps_2 = str2double(newStr_2); time_2=dt*N_steps_2;
str_3 = string(S(in+1).name); newStr_3 = extractBetween(str_3,"pfd.","_p000000.h5");
N_steps_3 = str2double(newStr_3); time_3=dt*N_steps_3;
%--------------------------------------------------------------------------

% The are cell within cell that's why the {1,1} additional
% Vectors
B_Gp = Vectorsc{1,in}{1,1} ;
E_Gp = Vectorsc{2,in}{1,1} ;
vi_Gp = Vectorsc{3,in}{1,1} ;
ve_Gp = Vectorsc{4,in}{1,1} ;
J_Gp = Vectorsc{5,in}{1,1} ;
Poyn_v_Gp = Vectorsc{6,in}{1,1} ; 

%Tensors
Pij_i_Gp = Tensorsc{1,in}{1,1} ;
Pij_e_Gp = Tensorsc{2,in}{1,1} ;
grad_vi_Gp = Tensorsc{3,in}{1,1} ;
grad_ve_Gp = Tensorsc{4,in}{1,1} ;

% Densities
ni_Gp = Densitiesc{1,in}{1,1};
ne_Gp = Densitiesc{2,in}{1,1} ;

% Poyn theorem terms
dt_E_em_Gp = Poyn_theoc{1,in}{1,1} ; %This is the one calculated adding the terms 
Div_Poyn_v_Gp = Poyn_theoc{2,in}{1,1} ;
JdotE_Gp = Poyn_theoc{3,in}{1,1} ;

%Kinetic Energy
vi_grad_Ki_Gp = kin_energy_terms_ic{1,in}{1,1} ; 
vi_Div_Pit_Gp = kin_energy_terms_ic{2,in}{1,1} ;
Ki_Div_vi_Gp = kin_energy_terms_ic{3,in}{1,1} ;
qi_ni_vi_E_Gp = kin_energy_terms_ic{4,in}{1,1} ;

ve_grad_Ke_Gp = kin_energy_terms_ec{1,in}{1,1} ;
ve_Div_Pet_Gp = kin_energy_terms_ec{2,in}{1,1} ;
Ke_Div_ve_Gp = kin_energy_terms_ec{3,in}{1,1} ;
qe_ne_ve_E_Gp = kin_energy_terms_ec{4,in}{1,1} ;

% Internal energy terms
vi_grad_Ui_Gp = Int_energy_terms_ic{1,in}{1,1} ;
dk_Qijk_i_Gp = Int_energy_terms_ic{2,in}{1,1} ;
dkQiik_05_i_Gp = Int_energy_terms_ic{3,in}{1,1} ;
Pit_grad_vi_Gp = Int_energy_terms_ic{4,in}{1,1} ;
Pit_grad_vi_ii_Gp = Int_energy_terms_ic{5,in}{1,1} ;
Pit_grad_vi_ij_Gp = Int_energy_terms_ic{6,in}{1,1} ;
Ui_Div_vi_Gp = Int_energy_terms_ic{7,in}{1,1} ; 

ve_grad_Ue_Gp = Int_energy_terms_ec{1,in}{1,1} ;
dk_Qijk_e_Gp = Int_energy_terms_ec{2,in}{1,1} ; 
dkQiik_05_e_Gp = Int_energy_terms_ec{3,in}{1,1} ;
Pet_grad_ve_Gp =Int_energy_terms_ec{4,in}{1,1} ;
Pet_grad_ve_ii_Gp = Int_energy_terms_ec{5,in}{1,1} ;
Pet_grad_ve_ij_Gp = Int_energy_terms_ec{6,in}{1,1} ;
Ue_Div_ve_Gp = Int_energy_terms_ec{7,in}{1,1} ; 

dtKi_in_Gp = dt_energydensities_ic{1,in} ;
dtKe_in_Gp = dt_energydensities_ic{2,in} ;
dtUi_in_Gp = dt_energydensities_ic{3,in} ;
dtUe_in_Gp = dt_energydensities_ic{4,in} ;
dtE_em_in_Gp = dt_energydensities_ic{5,in} ;
%--------------------------------------------------------------------------

% Now lets plot the averages in time
%--------------------------------------------------------------------------

% First the terms of the kinetic energy
%--------------------------------------------------------------------------
dtKe_in_Gp_m(1,in) = mean(dtKe_in_Gp,'all'); dtKe_in_Gp_std(1,in) = std2(dtKe_in_Gp);
ve_grad_Ke_Gp_m(1,in) =  mean(ve_grad_Ke_Gp,'all'); ve_grad_Ke_Gp_std(1,in) = std2(ve_grad_Ke_Gp_m);
ve_Div_Pet_Gp_m(1,in) = mean(ve_Div_Pet_Gp,'all'); ve_Div_Pet_Gp_std(1,in) = std2(ve_Div_Pet_Gp); 
Ke_Div_ve_Gp_m(1,in) = mean(Ke_Div_ve_Gp,'all'); Ke_Div_ve_Gp_std(1,in) = std2(Ke_Div_ve_Gp);
qe_ne_ve_E_Gp_m(1,in) = mean(qe_ne_ve_E_Gp,'all'); qe_ne_ve_E_Gp_std(1,in) = std2(qe_ne_ve_E_Gp);

dtKi_in_Gp_m(1,in) = mean(dtKi_in_Gp,'all'); dtKi_in_Gp_std(1,in) = std2(dtKi_in_Gp);
vi_grad_Ki_Gp_m(1,in) =  mean(vi_grad_Ki_Gp,'all'); vi_grad_Ki_Gp_std(1,in) = std2(vi_grad_Ki_Gp_m);
vi_Div_Pit_Gp_m(1,in) = mean(vi_Div_Pit_Gp,'all'); vi_Div_Pit_Gp_std(1,in) = std2(vi_Div_Pit_Gp); 
Ki_Div_vi_Gp_m(1,in) = mean(Ki_Div_vi_Gp,'all'); Ki_Div_vi_Gp_std(1,in) = std2(Ki_Div_vi_Gp);
qi_ni_vi_E_Gp_m(1,in) = mean(qi_ni_vi_E_Gp,'all'); qi_ni_vi_E_Gp_std(1,in) = std2(qi_ni_vi_E_Gp);
%--------------------------------------------------------------------------


% Second the terms of the Internal energy
%--------------------------------------------------------------------------
dtUe_in_Gp_m(1,in) = mean(dtUe_in_Gp,'all'); dtUe_in_Gp_std(1,in) = std2(dtUe_in_Gp);
ve_grad_Ue_Gp_m(1,in) =  mean(ve_grad_Ue_Gp,'all'); ve_grad_Ue_Gp_std(1,in) = std2(ve_grad_Ue_Gp_m);
dkQiik_05_e_Gp_m(1,in) = mean(dkQiik_05_e_Gp,'all'); dkQiik_05_e_Gp_std(1,in) = std2(dkQiik_05_e_Gp);
Pet_grad_ve_Gp_m(1,in) = mean(Pet_grad_ve_Gp,'all'); Pet_grad_ve_Gp_std(1,in) = std2(Pet_grad_ve_Gp);
Pet_grad_ve_ii_Gp_m(1,in) = mean(Pet_grad_ve_ii_Gp,'all'); Pet_grad_ve_ii_Gp_std(1,in) = std2(Pet_grad_ve_ii_Gp);
Pet_grad_ve_ij_Gp_m(1,in) = mean(Pet_grad_ve_ij_Gp,'all'); Pet_grad_ve_ij_Gp_std(1,in) = std2(Pet_grad_ve_ij_Gp);
Ue_Div_ve_Gp_m(1,in) = mean(Ue_Div_ve_Gp,'all'); Ue_Div_ve_Gp_std(1,in) = std2(Ue_Div_ve_Gp);

dtUi_in_Gp_m(1,in) = mean(dtUi_in_Gp,'all'); dtUi_in_Gp_std(1,in) = std2(dtUi_in_Gp);
vi_grad_Ui_Gp_m(1,in) =  mean(vi_grad_Ui_Gp,'all'); vi_grad_Ui_Gp_std(1,in) = std2(vi_grad_Ui_Gp_m);
dkQiik_05_i_Gp_m(1,in) = mean(dkQiik_05_i_Gp,'all'); dkQiik_05_i_Gp_std(1,in) = std2(dkQiik_05_i_Gp);
Pit_grad_vi_Gp_m(1,in) = mean(Pit_grad_vi_Gp,'all'); Pit_grad_vi_Gp_std(1,in) = std2(Pit_grad_vi_Gp);
Pit_grad_vi_ii_Gp_m(1,in) = mean(Pit_grad_vi_ii_Gp,'all'); Pit_grad_vi_ii_Gp_std(1,in) = std2(Pit_grad_vi_ii_Gp);
Pit_grad_vi_ij_Gp_m(1,in) = mean(Pit_grad_vi_ij_Gp,'all'); Pit_grad_vi_ij_Gp_std(1,in) = std2(Pit_grad_vi_ij_Gp);
Ui_Div_vi_Gp_m(1,in) = mean(Ui_Div_vi_Gp,'all'); Ui_Div_vi_Gp_std(1,in) = std2(Ui_Div_vi_Gp); 
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now lets do the plots of these quantities in time
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%--------------------------------------------------------------------------

% Again, first the plot for the kinetic energy.

%Create the time array
time_line=zeros(1,N-2);
dt = 0.06;
for in=2:N-1
str = string(S(in).name);
newStr = extractBetween(str,"pfd.","_p000000.h5");
N_steps = str2double(newStr); time_t = dt*N_steps;
time_line(1,in) = time_t; 
end

% Plots in time along a 1D trajectory
%--------------------------------------------------------------------------
lw=1.5;
f901 = figure(901);
curve1 = dtUe_in_Gp_m + dtUe_in_Gp_std;
curve2 = dtUe_in_Gp_m - dtUe_in_Gp_std;
plot(time_line, dtUe_in_Gp_m, 'k', 'LineWidth', lw);
hold on;
[fillhandle,msg]=jbfill(time_line,curve1,curve2,rand(1,3),'none',0,rand(1,1));
%[fillhandle,msg]=jbfill(time_line,curve1,curve2,blue,edge,add,transparency);
hold off

load('colorblind_colormap.mat');

% These are the time plots without normalization
meueXi_e=dtKe_in_Gp_m + ve_grad_Ke_Gp_m + ve_Div_Pet_Gp_m + Ke_Div_ve_Gp_m -qe_ne_ve_E_Gp_m;
miuiXi_i=dtKi_in_Gp_m + vi_grad_Ki_Gp_m + vi_Div_Pit_Gp_m + Ki_Div_vi_Gp_m -qi_ni_vi_E_Gp_m;
%
fs=18;
f902 = figure(902);
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(time_line, dtKe_in_Gp_m, 'Color',colorblind(1,:), 'LineWidth', lw);
hold on
plot(time_line, ve_grad_Ke_Gp_m, 'Color',colorblind(2,:), 'LineWidth', lw);
plot(time_line, ve_Div_Pet_Gp_m, 'Color',colorblind(6,:), 'LineWidth', lw);
plot(time_line, Ke_Div_ve_Gp_m, 'Color',colorblind(7,:), 'LineWidth', lw);
plot(time_line, -qe_ne_ve_E_Gp_m, 'Color',colorblind(5,:), 'LineWidth', lw);
plot(time_line, meueXi_e, 'Color',colorblind(8,:), 'LineWidth', lw);
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
ylabel('Kinetic energy rate terms','Interpreter','latex','FontSize',fs)
xlabel('$t \omega_{pi}$','Interpreter','latex','FontSize',fs)
%legend('$dK_{e}/dt$','$u_{e}\cdot \nabla K_{e}$','$u_{e} \cdot \nabla \cdot P_{e}$',...
 %   '$K_{e} \nabla \cdot u_{e}$','$q_{e}n_{e}u_{e}\cdot E$','Interpreter','latex','FontSize',fs)
legend('$d\varepsilon^{k}_{e}/dt$','$u_{e}\cdot \nabla \varepsilon^{k}_{e}$','$u_{e} \cdot \nabla \cdot \overline{P}_{e}$',...
    '$\varepsilon^{k}_{e} \nabla \cdot u_{e}$','$-q_{e}n_{e}u_{e}\cdot E$','$m_{e}u_{e}\Xi_{e}^{1}$','Interpreter','latex','FontSize',fs)
xlim([time_line(1), time_line(end)])
hold off
h2=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(time_line, dtKi_in_Gp_m, 'Color',colorblind(1,:), 'LineWidth', lw);
hold on
plot(time_line, vi_grad_Ki_Gp_m, 'Color',colorblind(2,:), 'LineWidth', lw);
plot(time_line, vi_Div_Pit_Gp_m, 'Color',colorblind(6,:), 'LineWidth', lw);
plot(time_line, Ki_Div_vi_Gp_m, 'Color',colorblind(7,:), 'LineWidth', lw);
plot(time_line, -qi_ni_vi_E_Gp_m, 'Color',colorblind(5,:), 'LineWidth', lw);
plot(time_line, miuiXi_i, 'Color',colorblind(8,:), 'LineWidth', lw);
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
ylabel('Kinetic energy rate terms','Interpreter','latex','FontSize',fs)
xlabel('$t \omega_{pi}$','Interpreter','latex','FontSize',fs)
legend('$d\varepsilon^{k}_{i}/dt$','$u_{i}\cdot \nabla \varepsilon^{k}_{i}$','$u_{i} \cdot \nabla \cdot \overline{P}_{i}$',...
    '$\varepsilon^{k}_{i} \nabla \cdot u_{i}$','$-q_{i}n_{i}u_{i}\cdot E$','$m_{i}u_{i}\Xi_{i}^{1}$','Interpreter','latex','FontSize',fs)
xlim([time_line(1), time_line(end)])
hold off

f903 = figure(903);
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(time_line, dtUe_in_Gp_m, 'Color',colorblind(1,:), 'LineWidth', lw);
hold on
plot(time_line, ve_grad_Ue_Gp_m, 'Color',colorblind(2,:), 'LineWidth', lw);
plot(time_line, dkQiik_05_e_Gp_m, 'Color',colorblind(6,:), 'LineWidth', lw);
plot(time_line, Pet_grad_ve_Gp_m, 'Color',colorblind(7,:), 'LineWidth', lw);
plot(time_line, Pet_grad_ve_ii_Gp_m, 'Color',colorblind(5,:), 'LineWidth', lw);
plot(time_line, Pet_grad_ve_ij_Gp_m, 'Color',colorblind(8,:), 'LineWidth', lw);
plot(time_line, Ue_Div_ve_Gp_m, 'Color',colorblind(9,:), 'LineWidth', lw);
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
ylabel('Kinetic energy rate terms','Interpreter','latex','FontSize',fs)
xlabel('$t \omega_{pi}$','Interpreter','latex','FontSize',fs)
%legend('$dU_{e}/dt$','$u_{e}\cdot \nabla U_{e}$','$\nabla \cdot Q_{e}$',...
 %   '$P_{e} : \nabla u_{e}$','$(P_{e} : \nabla u_{e})_{ii}$',...
  %  '$(P_{e} : \nabla u_{e})_{ij}$','$U_{e}\ \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs)
legend('$d\varepsilon^{th}_{e}/dt$','$u_{e}\cdot \nabla \varepsilon^{th}_{e}$','$\nabla \cdot h_{e}$',...
    '$\overline{P}_{e} : \nabla u_{e}$','$p\theta_{e}$',...
    '$PiD_{e}$','$\varepsilon^{th}_{e}\ \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs)
xlim([time_line(1), time_line(end)])
hold off
h2=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(time_line, dtUi_in_Gp_m, 'Color',colorblind(1,:), 'LineWidth', lw);
hold on
plot(time_line, vi_grad_Ui_Gp_m, 'Color',colorblind(2,:), 'LineWidth', lw);
plot(time_line, dkQiik_05_i_Gp_m, 'Color',colorblind(6,:), 'LineWidth', lw);
plot(time_line, Pit_grad_vi_Gp_m, 'Color',colorblind(7,:), 'LineWidth', lw);
plot(time_line, Pit_grad_vi_ii_Gp_m, 'Color',colorblind(5,:), 'LineWidth', lw);
plot(time_line, Pit_grad_vi_ij_Gp_m, 'Color',colorblind(8,:), 'LineWidth', lw);
plot(time_line, Ui_Div_vi_Gp_m, 'Color',colorblind(9,:), 'LineWidth', lw);
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
ylabel('Kinetic energy rate terms','Interpreter','latex','FontSize',fs)
xlabel('$t \omega_{pi}$','Interpreter','latex','FontSize',fs)
legend('$d\varepsilon^{th}_{i}/dt$','$\varepsilon^{th}_{i}\cdot \nabla \varepsilon^{th}_{i}$','$\nabla \cdot h_{i}$',...
    '$\overline{P}_{i} : \nabla u_{i}$','$p\theta_{i}$',...
    '$PiD_{i}$','$\varepsilon^{th}_{i}\ \nabla \cdot u_{i} $','Interpreter','latex','FontSize',fs)
xlim([time_line(1), time_line(end)])
hold off
%}
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f902,strcat('Kinetic_energy_time_Xi.png'));
saveas(f903,strcat('Internal_energy_time_Xi.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f1131,strcat('Blines_ve_vi_B'+ string(N_steps) +'.png'));
%cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%clf(f1131)
%--------------------------------------------------------------------------


%}
%--------------------------------------------------------------------------

%cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
%scatter
%view(3);

%f78=figure(78)
%scatter(ccc(:,1) , ccc(:,2))
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
sv=7;
%sv=5;
%--------------------------------------------------------------------------

%PLOTS!
%--------------------------------------------------------------------------
%Plots on the slide. The velocities ions, electrons, current and magnetic field components
%--------------------------------------------------------------------------

%Vi Ve J B
f1=figure(1);
for i=1:12
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(4,3,i,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=vi_Gp{1,i};
    titl = 'vi';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=ve_Gp{1,i2};
    titl = 've';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=J_Gp{1,i3};
    titl = 'J';
    s=i3;
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=B_Gp{1,i4};
    if (i4==2)
        dum_paux=B_Gp{1,i4};
        dum_p = dum_paux ;%+ 0.1;
    end
    titl = 'B';
    s=i4;    
end
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p)))); 
set(hc,'edgecolor','none')
hcb=colorbar;
if(i==11)
    colormap(h1,hot);
    %caxis([-lim_yp 0])
else
    colormap(h1,BWR);
    caxis([-lim_yp lim_yp])
end
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
%hcb.Title.String='';
%set(get(hcb,'Title'),'String',titl + string(s) ,'FontSize',14)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = '$p / d_{i}$';
    ax.YTick = [2 4 6 8 10];
end
if (9<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];    
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', colorFVBl);
%set(h2,'AutoScale','on', 'AutoScaleFactor', 2, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
%}

%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f1,strcat('vi_ve_J_B_vev_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

%end



%Remove this to see the rest. This is important!!
%--------------------------------------------------------------------------
%sv=3;


clf(f77)
clf(f78)
clf(f13)
%}

%--------------------------------------------------------------------------
