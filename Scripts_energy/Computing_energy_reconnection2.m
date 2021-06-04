%--------------------------------------------------------------------------
%This script is to compute the energy budget
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data'; %This is important because the xdmf files are in that directory
path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

% useful parameters used in the simulation
%--------------------------------------------------------------------------
kb = 1.; mu0 =1.; mi_over_me_ =100.; vA_over_c_ = 0.1; B0_ = vA_over_c_;
mi = 1.; me = 1. / mi_over_me_;
omega_pi=1;
Omega_e=vA_over_c_*sqrt(mi_over_me_)*omega_pi;
Omega_i=Omega_e / mi_over_me_;
TOmega_e=1/Omega_e; 
TOmega_i=1/Omega_i;
dt=0.06;
%--------------------------------------------------------------------------
S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
i=19;
%i=2;

%for i=1:14
%--------------------------------------------------------------------------
disp(strcat('Computing step ...',S(i).name))
%--------------------------------------------------------------------------
%N_steps = 2000;
%time_t = dt*N_steps; %(In terms of 1/omega_pi)
str = string(S(i).name);
newStr = extractBetween(str,"pfd.","_p000000.h5");
N_steps = str2double(newStr);
time_t = dt*N_steps;
%--------------------------------------------------------------------------

fileID =  S(i).name; %change the 1 per i
%h5disp(fileID);
info = h5info(fileID);
T_1st = info.Groups(1).Name;
E_1st = info.Groups(17).Name;
B_1st = info.Groups(18).Name;
J_1st = info.Groups(19).Name;
V_1st = info.Groups(25).Name;
n_1st = info.Groups(23).Name;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
nx1=166;    nx2=333;
ny1=232;    ny2=400;
nz1=1200;   nz2=1367;

px2 = 168; py2 = 168; pz2 = 168;
start = [nx1 ny1 nz1];
count = [px2 py2 pz2];
%------------------------------------------------------------------

resx=0.06; %0.06081
resy=0.06;
resz=0.06; %0.06152

Lx=px2*resx; k0x=1/Lx;
Ly=px2*resy; k0y=1/Ly;
Lz=px2*resz; k0z=1/Lz;
k00=sqrt(k0x*k0x + k0y*k0y + k0z*k0z);

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
nec{1,1} = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'),start,count); %ne=ne(nx1:nx2,ny1:ny2,nz1:nz2);
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
%vec_a = [-3.3333 0 1]; %is this correct
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
ni=ni{1,1};
ne=ne{1,1};
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

%For the partial time derivative <----------- This is where I was trying to
%get a different times
%------------------------------------------------------------------
Ki_1 = Ki;
Ke_1 = Ke;
Ui_1 = Ui;
Ue_1 = Ue;
%clearvars -except Ki_0 Ke_0 Ui_0 Ue_0 

% This requeires to store first Ks_0 and Us_0
partialKi = (Ki_1 - Ki_0)/6;
partialKe = (Ke_1 - Ke_0)/6;
partialUi = (Ui_1 - Ui_0)/6;
partialUe = (Ue_1 - Ue_0)/6;
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
% Pit_grad_vi
Pit_grad_vi = Pij_i{1,1}.*grad_vi{1,1} + Pij_i{1,2}.*grad_vi{1,2} + Pij_i{1,3}.*grad_vi{1,3} +...
              Pij_i{1,2}.*grad_vi{2,1} + Pij_i{2,2}.*grad_vi{2,2} + Pij_i{2,3}.*grad_vi{2,3} +...
              Pij_i{1,3}.*grad_vi{3,1} + Pij_i{2,3}.*grad_vi{3,2} + Pij_i{3,3}.*grad_vi{3,3};

% Pet_grad_ve
Pet_grad_ve = Pij_e{1,1}.*grad_ve{1,1} + Pij_e{1,2}.*grad_ve{1,2} + Pij_e{1,3}.*grad_ve{1,3} +...
              Pij_e{1,2}.*grad_ve{2,1} + Pij_e{2,2}.*grad_ve{2,2} + Pij_e{2,3}.*grad_ve{2,3} +...
              Pij_e{1,3}.*grad_ve{3,1} + Pij_e{2,3}.*grad_ve{3,2} + Pij_e{3,3}.*grad_ve{3,3};

% vi_Div_Pit
vi_Div_Pit = vi{1,1}.*Div_Pi_t{1,1} + vi{1,2}.*Div_Pi_t{1,2} + vi{1,3}.*Div_Pi_t{1,3}; 

% ve_Div_Pet
ve_Div_Pet = ve{1,1}.*Div_Pe_t{1,1} + ve{1,2}.*Div_Pe_t{1,2} + ve{1,3}.*Div_Pe_t{1,3}; 

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

%HEre I am changing stuff... The following step is to plot in the same
%slices the energy terms and analise them

%--------------------------------------------------------------------------


% This part is to make the plots to see the data
%--------------------------------------------------------------------------
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

%
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

partialKic{1,1}=partialKi;
Gp=partialKic;
[partialKi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
partialKec{1,1}=partialKe;
Gp=partialKec;
[partialKe_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

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

partialUic{1,1}=partialUi;
Gp=partialUic;
[partialUi_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP
partialUec{1,1}=partialUe;
Gp=partialUec;
[partialUe_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi); clearvars GP

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



% 1D cuts to see how the terms releates each other (keep just the energy terms)
%{
%--------------------------------------------------------------------------
ypc_i=127;

cut_1d=ve_Gp{1,3}';
ve_2_1d=cut_1d(:,ypc_i); 

% Kinetic
cut_1d=partialKe_Gp{1,1}';
partialKe_1d = cut_1d(:,ypc_i);
cut_1d=ve_grad_Ke_Gp{1,1}';
ve_grad_Ke_1d = cut_1d(:,ypc_i);
cut_1d=ve_Div_Pet_Gp{1,1}';
ve_Div_Pet_1d = cut_1d(:,ypc_i);
cut_1d=Ke_Div_ve_Gp{1,1}';
Ke_Div_ve_1d = cut_1d(:,ypc_i);
cut_1d=qe_ne_ve_E_Gp{1,1}';
qe_ne_ve_E_Gp_1d = cut_1d(:,ypc_i);

f213 = figure(213);
%subplot(1,2,1)
plot(xll,partialKe_1d + ve_grad_Ke_1d,'k')
hold on
plot(xll,ve_Div_Pet_1d,'-+r')
plot(xll,Ke_Div_ve_1d,'b')
plot(xll,-qe_ne_ve_E_Gp_1d,'m')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
legend('dKe/dt','vedivPe','Kedivve','-qEe')
%plot(xll,ve_2_1d)
hold off
%subplot(1,2,2)
%dum_p=ve_Gp{1,2};
%hc = pcolor(XLL,YLL,dum_p);
%set(hc,'edgecolor','none')

f214 = figure(214);
%subplot(1,2,1)
plot(xll,ve_2_1d)


% Internal
cut_1d=partialUe_Gp{1,1}';
partialUe_1d = cut_1d(:,ypc_i);
cut_1d=ve_grad_Ue_Gp{1,1}';
ve_grad_Ue_1d = cut_1d(:,ypc_i);
%}
%--------------------------------------------------------------------------



%PLOTS!
%--------------------------------------------------------------------------
%Plots on the slide. The velocities ions, electrons, current and magnetic field components
%--------------------------------------------------------------------------
%colorFVBl=[0.5 0.5 0.5 0.5];
colorFVBl=[0 0 0 0.5];
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

%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f1,strcat('vi_ve_J_B_vev_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

%clearvars ni ne
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

sv=2;
f77=figure(77);
hlines=streamline(aaa);
hold on
%plot(bbb(:,1) , bbb(:,2), '*k')
%plot(ccc(:,1) , ccc(:,2), '*r')
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', 5, 'Color', 'b');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
xlabel('$$r / d_{i}$$','Interpreter','latex','FontSize',20);
ylabel('$$p / d_{i}$$','Interpreter','latex','FontSize',20);
title('$$v_{e} \ vectors \ and \ Blines$$','Interpreter','latex','FontSize',20)
set(gca,'XScale','lin','YScale','lin','FontSize',18)
hold off

%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f77,strcat('ve_vec_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
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
    caxis([-lim_yp 0])
elseif(i==2)
    caxis([-lim_yp lim_yp])
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
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
set(h2,'AutoScale','on', 'AutoScaleFactor', 6, 'Color', 'k');
%colorbar('off')
end
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)



%{
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f13,strcat('Ba_Ja_vec_ele_slideRF2_per_2'+ string(N_steps) +'.png'));
%saveas(f1,strcat('vi_ve_J_B_vev_no_Bl_slideRF2_per'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%}

% For ni and ne
%--------------------------------------------------------------------------
f2=figure(2);
for i=1:2
h1=subaxis(1,2,i,'SV',0.035,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i==1)
    i1=i;
    dum_p=ni_Gp{1,1};
    titl = 'ni';
    s=i1;
elseif (i==2)
    i2=i-3;
    dum_p=ne_Gp{1,1};
    titl = 'ne';
    s=i2;
end
%dum_p = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
if (i==1 || i==2 )
    caxis([0 lim_yp])
    colormap(h1,jet); 
end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl ,0,[0.5 0]})
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
if (3>i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
contour(XLL,YLL,dum_p);
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------

% For partial Ke and partial Ue
%--------------------------------------------------------------------------
f27=figure(27);
for i=1:4
h1=subaxis(2,2,i,'SV',0.035,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i==1)
    dum_p=partialKe_Gp{1,1};
    titl = 'dKe/dt';
elseif (i==2)
    dum_p=partialUe_Gp{1,1};
    titl = 'dUe/dt';
elseif (i==3)
    dum_p=ve_grad_Ke_Gp{1,1};
    titl = 'vegradKe';
elseif (i==4)
    dum_p=ve_grad_Ue_Gp{1,1};
    titl = 'vegradUe';
end
%dum_p = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
%if (i==1 || i==2 )
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
%end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl ,0,[0.5 0]})
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
if (3>i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
contour(XLL,YLL,dum_p);
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------


% For Ex, Ey, Ez, and the poynting terms
%--------------------------------------------------------------------------
f3=figure(3);
for i=1:6
h1=subaxis(2,3,i,'SV',0.035,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i==1)
    i1=i;
    dum_p=E_Gp{1,1};
    titl = 'E_{1}';
elseif (i==2)
    dum_p=E_Gp{1,2};
    titl = 'E_{2}';
elseif (i==3)
    dum_p=E_Gp{1,3};
    titl = 'E_{3}';
elseif (i==4)
    dum_p=JdotE_Gp{1,1};
    titl = '\bf{J} \cdot \bf{E}';
elseif (i==5)
    dum_p=Div_Poyn_v_Gp{1,1};
    titl = '\nabla \cdot \bf{S}';
elseif (i==6)
    dum_p=dt_E_em_Gp{1,1};
    %dum_p=Div_Poyn_v_Gp{1,1} - JdotE_Gp{1,1};
    titl = 'dE_{em}/dt';   
end
%dum_p = sign(dum_p).*log10(abs(dum_p));
dum_p = medfilt2(dum_p);
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
elseif (i==4 )
    caxis([-0.005 0.005])
    colormap(h1,BWR); 
elseif (i==5 || i==6)
    caxis([-0.02 0.02])
    colormap(h1,BWR); 
end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Title.String='';
hcb.Title.Position = [-120 258];
hcb.Ruler.Exponent = -2;
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl ,0,[0.5 0]})
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
if (3<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
contour(XLL,YLL,dum_p);
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------

% Poynting components
%--------------------------------------------------------------------------
f31=figure(31);
for i=1:3
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(1,3,i,'SV',0.01,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i==1)
    i1=i;
    dum_p=Poyn_v_Gp{1,1};
    titl = 'S_{r}';
    s=i1;
elseif (i==2)
    i2=i-3;
    dum_p=Poyn_v_Gp{1,2};
    titl = 'S_{a}';
    s=i2;
elseif (i==3)
    i3=i-6;
    dum_p=Poyn_v_Gp{1,3};
    titl = 'S_{p}';
    s=i3;
end
hc = pcolor(XLL,YLL,dum_p);
lim_yp=0.8*max(max(max(abs(dum_p))));
if(i==1) 
    colormap(h1,BWR);
    caxis([-lim_yp lim_yp])
elseif(i>1)
    caxis([-lim_yp lim_yp])
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
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1.5, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 2, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
dCp=25;
%--------------------------------------------------------------------------
% This is the plot for the energy terms in the kinetic equation
%--------------------------------------------------------------------------
GPijf{1,1}=-Ki_Div_vi_Gp{1,1}; GPijf{1,2}=-vi_Div_Pit_Gp{1,1}; GPijf{1,3}=qi_ni_vi_E_Gp{1,1};
GPijf{2,1}=-Ke_Div_ve_Gp{1,1}; GPijf{2,2}=-ve_Div_Pet_Gp{1,1}; GPijf{2,3}=qe_ne_ve_E_Gp{1,1};

f4=figure(4);
for i=1:6
h1=subaxis(2,3,i,'SV',0.04,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=GPijf{1,i};
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=GPijf{2,i2};
    s=i2;
end
if (i==1)
    titl = '- K_{i} \nabla \cdot \bf{v}_{i}';
elseif(i==2)
    titl = '- \bf{v}_{i} \cdot \nabla \cdot \bf{P}_{i}';
elseif(i==3)
    titl = 'q_{i}n_{i}\bf{v}_{i} \cdot \bf{E}';
elseif(i==4)
    titl = '- K_{e} \nabla \cdot \bf{v}^{e}';
elseif(i==5)
    titl = '- \bf{v}_{e} \cdot \nabla \cdot \bf{P}_{e}';
elseif(i==6)
    titl = 'q_{e}n_{e}\bf{v}_{e} \cdot \bf{E}';
end 
%dum_p = sign(dum_p).*log10(abs(dum_p)); % this is to get the log scale with negative values
%hc = pcolor(XLL,YLL,dum_p);
hc = contourf(XLL,YLL,dum_p);
lim_yp=0.7*max(max(max(abs(dum_p))));
%lim_yp=1*max(max(max(abs(dum_p))));
%set(hc,'edgecolor','none')
if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
elseif (i==4 || i==5 || i==6)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
end
hcb=colorbar;
hcb.Title.Position = [-120 258];
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Ruler.Exponent = -2;
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[1 0]})
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
if (3<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2, \ sg(y)log(y)$$','Interpreter','latex','FontSize',20)
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
clearvars GPijf
%--------------------------------------------------------------------------


% This is the plot for the energy terms in the potential equation
%--------------------------------------------------------------------------
GPijf{1,1}=Ui_Div_vi_Gp{1,1}; GPijf{1,2}=Pit_grad_vi_Gp{1,1}; GPijf{1,3}=dkQiik_05_i_Gp{1,1};
GPijf{2,1}=Ue_Div_ve_Gp{1,1}; GPijf{2,2}=Pet_grad_ve_Gp{1,1}; GPijf{2,3}=dkQiik_05_e_Gp{1,1};

f5=figure(5);
for i=1:6
h1=subaxis(2,3,i,'SV',0.04,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=GPijf{1,i};
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=GPijf{2,i2};
    s=i2;
end
if (i==1)
    %titl = 'iUdivV';
    titl = 'U_{i} \nabla \cdot \bf{v}_{i}';
elseif(i==2)
    titl = '\bf{P}_{i} : \nabla \bf{v}_{i}';
elseif(i==3)
    titl = '0.5 \nabla \cdot \bf{Q}_{i} ';
elseif(i==4)
    titl = 'U_{e} \nabla \cdot \bf{v}_{e}';
elseif(i==5)
    titl = '\bf{P}_{e} : \nabla \bf{v}_{e}';
elseif(i==6)
    titl = '0.5 \nabla \cdot \bf{Q}_{e}';
end 
dum_p = sign(dum_p).*log10(abs(dum_p)); % this is to get the log scale with negative values
%hc = pcolor(XLL,YLL,dum_p);
hc = contourf(XLL,YLL,dum_p);
%lim_yp=0.8*max(max(max(abs(dum_p))));
lim_yp=1*max(max(max(abs(dum_p))));
%set(hc,'edgecolor','none')
if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,jet(dCp)); 
elseif (i==4 || i==5 || i==6)
    caxis([-lim_yp lim_yp])
    colormap(h1,jet(dCp)); 
end
hcb=colorbar;
hcb.Label.FontSize = 14;
%hcb.Title.String=''; hcb.Label.String = 'Elevation (ft in 1000s) $\frac{1}{2}$';
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
%hcb.Ruler.Exponent = -2;
hcb.Title.Position = [-120 258];
set(get(hcb,'Title'),'String',titl,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[1 0]})
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
if (3<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2, \ sg(y)log(y)$$','Interpreter','latex','FontSize',20)
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
clearvars GPijf
%--------------------------------------------------------------------------




% For the total time derivatives of the kinetic and potential energy 
%--------------------------------------------------------------------------
f6=figure(6);
for i=1:4
h1=subaxis(2,2,i,'SV',0.04,'SH',0.004,'MR',0.02,'ML',0.06,'PL',0.005,'PR',0.006);
if (i==1)
    i1=i;
    dum_p=dtKi_Gp{1,1};
    titl = 'dK_{i}/dt';
elseif (i==2)
    dum_p=dtKe_Gp{1,1};
    titl = 'dK_{e}/dt';
elseif (i==3)
    dum_p=dtUi_Gp{1,1};
    titl = 'dU_{i}/dt';
elseif (i==4)
    dum_p=dtUe_Gp{1,1};
    titl = 'dU_{e}/dt';
end
%dum_p = medfilt2(dum_p);
dum_p = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p);
%lim_yp=0.5*max(max(max(abs(dum_p))));
lim_yp=1*max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
%if (i==1 || i==2 )
    caxis([-lim_yp lim_yp])
    colormap(h1,jet(dCp)); 
%end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Title.Position = [-120 205];
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
%hcb.Ruler.Exponent = -2;
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl ,0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.TickLabelInterpreter = 'latex';ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.FontSize=14;
ax.XTick = [];
ax.YTick = [];
if (mod(i,2)==1)
    ax.YLabel.String = '$p / d_{i}$';
    ax.YTick = [2 4 6 8 10];
end
if (2<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
contour(XLL,YLL,dum_p);
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2 \ sg(y)log(y)$$','Interpreter','latex','FontSize',20)
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------


%Kinetic and potential energy of the electrons
%--------------------------------------------------------------------------
GPijf{1,1}= - Ke_Div_ve_Gp{1,1}; GPijf{1,2}= - ve_Div_Pet_Gp{1,1}; GPijf{1,3}= qe_ne_ve_E_Gp{1,1};
GPijf{2,1}= - Ue_Div_ve_Gp{1,1}; GPijf{2,2}= - Pet_grad_ve_Gp{1,1}; GPijf{2,3}= -dkQiik_05_e_Gp{1,1};
f7=figure(7);
for i=1:6
h1=subaxis(2,3,i,'SV',0.04,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=GPijf{1,i};
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=GPijf{2,i2};
    s=i2;
end
if (i==1)
    titl = '- K_{e} \nabla \cdot \bf{v}_{e}';
elseif(i==2)
    titl = '- \bf{v}_{e} \cdot \nabla \cdot \bf{P}_{e}';
elseif(i==3)
    titl = 'q_{e}n_{e}\bf{v}_{e} \cdot \bf{E}';
elseif(i==4)
    titl = '- U_{e} \nabla \cdot \bf{v}_{e}';
elseif(i==5)
    titl = '- \bf{P}_{e} : \nabla \bf{v}_{e}';
elseif(i==6)
    titl = '- 0.5 \nabla \cdot \bf{Q}_{e}';
end 
%dum_p = sign(dum_p).*log10(abs(dum_p)); % this is to get the log scale with negative values
%hc = contourf(XLL,YLL,dum_p);
hc = pcolor(XLL,YLL,dum_p);
set(hc,'edgecolor','none')
lim_yp=0.1*max(max(max(abs(dum_p))));
if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
elseif (i==4 || i==5 || i==6)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
end
hcb=colorbar;
hcb.Title.Position = [-120 258];
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Ruler.Exponent = -2;
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[1 0]})
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
if (3<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
%set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2, \ sg(y)log(y)$$','Interpreter','latex','FontSize',20)
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
clearvars GPijf

%Kinetic and potential energy of the ions
%--------------------------------------------------------------------------
GPijf{1,1}= - Ki_Div_vi_Gp{1,1}; GPijf{1,2}= - vi_Div_Pit_Gp{1,1}; GPijf{1,3}= qi_ni_vi_E_Gp{1,1};
GPijf{2,1}= - Ui_Div_vi_Gp{1,1}; GPijf{2,2}= - Pit_grad_vi_Gp{1,1}; GPijf{2,3}= -dkQiik_05_i_Gp{1,1};
f8=figure(8);
for i=1:6
h1=subaxis(2,3,i,'SV',0.04,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=GPijf{1,i};
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=GPijf{2,i2};
    s=i2;
end
if (i==1)
    titl = '- K_{i} \nabla \cdot \bf{v}_{i}';
elseif(i==2)
    titl = '- \bf{v}_{i} \cdot \nabla \cdot \bf{P}_{i}';
elseif(i==3)
    titl = 'q_{i}n_{i}\bf{v}_{i} \cdot \bf{E}';
elseif(i==4)
    titl = '- U_{i} \nabla \cdot \bf{v}_{i}';
elseif(i==5)
    titl = '- \bf{P}_{i} : \nabla \bf{v}_{i}';
elseif(i==6)
    titl = '- 0.5 \nabla \cdot \bf{Q}_{i}';
end 
%dum_p = sign(dum_p).*log10(abs(dum_p)); % this is to get the log scale with negative values
%hc = contourf(XLL,YLL,dum_p);
hc = pcolor(XLL,YLL,dum_p);
set(hc,'edgecolor','none')
lim_yp=0.5*max(max(max(abs(dum_p))));
if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
elseif (i==4 || i==5 || i==6)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
end
hcb=colorbar;
hcb.Title.Position = [-120 258];
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Ruler.Exponent = -2;
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[1 0]})
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
if (3<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2, \ sg(y)log(y)$$','Interpreter','latex','FontSize',20)
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
clearvars GPijf


%--------------------------------------------------------------------------

%These are the better plots but for any reason the streamlines and vectors do not appear 
%--------------------------------------------------------------------------
%Vi Ve J B
%--------------------------------------------------------------------------
%{
f11=figure(11);
for i=1:12
h1=subaxis(4,3,i,'SV',0.002,'SH',0.004,'MR',0.002,'ML',0.04,'PL',0.002,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=vi_Gp{1,i};
    titl = 'vi';
    s=i1;
    tic=0.07;  
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=ve_Gp{1,i2};
    titl = 've';
    s=i2;
    if(i==4 || i==5)
    tic=0.25; 
    elseif(i==6)
    tic=1;   
    end
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=J_Gp{1,i3};
    titl = 'J';
    s=i3;
    if(i==7 || i==8)
    tic=0.2;  
    elseif(i==9)
    tic=0.8;  
    end
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=B_Gp{1,i4};
    %if (i4==3)
    %    dum_paux=B{1,i4};
    %    dum_p = dum_paux;% - 0.1;
    %end
    titl = 'B';
    s=i4;
    if(i==10 || i==11)   
    tic=0.12;  
    elseif(i==12)
    tic=0.07;  
    end
end
%imagesc(xl,yl,dum_p(:,:,slide_n)') % This was from the old set up
imagesc(xl,yl,dum_p(:,:)')
%hc = pcolor(XLL,YLL,dum_p);
%lim_yp=max(max(max(abs(dum_p))));
%caxis([-lim_yp lim_yp])
caxis([-tic tic])  
hcb=colorbar;
colormap(BWR);
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
hcb.Location = 'eastoutside';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
t=get(hcb,'Limits');
 T=[-3*tic/4 -tic/2, -tic/4, 0, tic/4, tic/2 3*tic/4];
set(hcb,'Ticks',T)
TL=arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
set(hcb,'TickLabels',TL)
%pbaspect([1 1 1])
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (9<i)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
hold on
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xl(1) xl(end)]);
ylim([yl(1) yl(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1} \ ; \ xy-plane$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------


%Pij_i
%--------------------------------------------------------------------------
Pij_e_d{1,1} = Pij_i{1,2}; Pij_e_d{1,2} = Pij_i{1,3}; Pij_e_d{1,3} = Pij_i{2,3};
Pij_e_d{2,1} = Pij_i{1,1}; Pij_e_d{2,2} = Pij_i{2,2}; Pij_e_d{2,3} = Pij_i{3,3};
f2=figure(2);
k=1;
for i=1:6
h1=subaxis(2,3,k,'SV',0.002,'SH',0.004,'MR',0.002,'ML',0.04,'PL',0.002,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=Pij_e_d{1,i};
    k=i+1;
elseif (4<=i && i<=6)
    i2=i-3;
    dum_p=Pij_e_d{2,i2};
    k=k+1;
end
if (i==1)
    titl = 'Pi12';
elseif(i==2)
    titl = 'Pi13';
elseif(i==3)
    titl = 'Pi23';
elseif(i==4)
    titl = 'Pi11';
elseif(i==5)
    titl = 'Pi22';
elseif(i==6)
    titl = 'Pi33';
end 
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,slide_n)))));
hcb=colorbar;
if (i<=3)
   colormap(h1,BWR); 
   %caxis([0 lim_yp])
   caxis([-0.005 0.005])
   set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[0 0.05]})
else
   colormap(h1,jet);
   %caxis([-lim_yp lim_yp])
   caxis([0 0.025])
   set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[0.0125 0.05]})
end
hcb.Label.Interpreter = 'latex';
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
hcb.Location = 'northoutside';
axis tight
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (i==1 || i==4)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (3<i)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
hold on
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), uex(1:sv:end,1:sv:end), vey(1:sv:end,1:sv:end), 0);
hlines=streamline(xm,ym,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1.3, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xl(1) xl(end)]);
ylim([yl(1) yl(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}\ ; \ xy-plane$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------

%Pij_e
%--------------------------------------------------------------------------
Pij_e_d{1,1} = Pij_e{1,2}; Pij_e_d{1,2} = Pij_e{1,3}; Pij_e_d{1,3} = Pij_e{2,3};
Pij_e_d{2,1} = Pij_e{1,1}; Pij_e_d{2,2} = Pij_e{2,2}; Pij_e_d{2,3} = Pij_e{3,3};
f3=figure(3);
k=1;
for i=1:6
h1=subaxis(2,3,k,'SV',0.002,'SH',0.004,'MR',0.002,'ML',0.04,'PL',0.002,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=Pij_e_d{1,i};
    k=i+1;
elseif (4<=i && i<=6)
    i2=i-3;
    dum_p=Pij_e_d{2,i2};
    k=k+1;
end
if (i==1)
    titl = 'Pe12';
elseif(i==2)
    titl = 'Pe13';
elseif(i==3)
    titl = 'Pe23';
elseif(i==4)
    titl = 'Pe11';
elseif(i==5)
    titl = 'Pe22';
elseif(i==6)
    titl = 'Pe33';
end 
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,slide_n)))));
hcb=colorbar;
if (i<=3)
   colormap(h1,BWR); 
   %caxis([0 lim_yp])
   caxis([-0.005 0.005])
   set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[0 0.05]})
else
   colormap(h1,jet);
   %caxis([-lim_yp lim_yp])
   caxis([0 0.01])
   set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[0.005 0.05]})
end
hcb.Label.Interpreter = 'latex';
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
hcb.Location = 'northoutside';
axis tight
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (i==1 || i==4)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (3<i)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
hold on
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), uex(1:sv:end,1:sv:end), vey(1:sv:end,1:sv:end), 0);
hlines=streamline(xm,ym,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1.3, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xl(1) xl(end)]);
ylim([yl(1) yl(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}\ ; \ xy-plane$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------

%------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f1,strcat('vi_ve_J_B_vev_Bl_cte_'+ string(N_steps) +'.png'));
saveas(f2,strcat('Pij_i_vev_Bl_cte_'+ string(N_steps) +'.png'));
saveas(f3,strcat('Pij_e_vev_Bl_cte_'+ string(N_steps) +'.png'));
cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
%}
%end
%--------------------------------------------------------------------------




%Lets check that the rotation is happening. This works. 
%-----------------------------------------------------------------------
%{
x=linspace(0,2*pi,100); y=x; z=x;
[X_dum, Y_dum, Z_dum] = meshgrid(x,y,z);
v_dum_x = cos(X_dum) + sin(Y_dum); 
v_dum_y = cos(Y_dum);
v_dum_z = cos(Z_dum);

v_dum_mag=sqrt(v_dum_x.*v_dum_x + v_dum_y.*v_dum_y + v_dum_z.*v_dum_z);

vp_dum_x=v_dum_x.*DCM(1,1) + v_dum_y.*DCM(1,2) + v_dum_z.*DCM(1,3);
vp_dum_y=v_dum_x.*DCM(2,1) + v_dum_y.*DCM(2,2) + v_dum_z.*DCM(2,3);
vp_dum_z=v_dum_x.*DCM(3,1) + v_dum_y.*DCM(3,2) + v_dum_z.*DCM(3,3);

vp_dum_mag=sqrt(vp_dum_x.*vp_dum_x + vp_dum_y.*vp_dum_y + vp_dum_z.*vp_dum_z);


f38=figure(38);
ax(1)=subplot(1,2,1);
pcolor(x, y,v_dum_x(:,:,1)); %Xd
colormap(BWR);
ax(2)=subplot(1,2,2);
pcolor(x, y,vp_dum_x(:,:,1)); %Xd
colormap(BWR);

[slicex, sliceInd1,subX1,subY1,subZ1] = extractSlice(v_dum_x,pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);
[slicey, sliceInd1,subX1,subY1,subZ1] = extractSlice(v_dum_y,pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);
[slicez, sliceInd1,subX1,subY1,subZ1] = extractSlice(v_dum_z,pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);
[slice2, sliceInd2,subX2,subY2,subZ2] = extractSlice(vp_dum_x,pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);
X1 = slicex;
X2 = slice2;

size_X=size(X1);
xll=linspace(1,size_X(2),size_X(2));
yll=linspace(1,size_X(1),size_X(1));
[XLL, YLL] = meshgrid(xll,yll);

ubx = slicex;   vby = slicey; wbz = slicez;
bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

Vev_x = slicex; Vev_y = slicey; Vev_z = slicez; 
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

NN=80; xstart = max(xll)*rand(NN,1); ystart = max(yll)*rand(NN,1); %Do this just once
sv=7;

f39=figure(39);
ax(1)=subplot(1,2,1);
hc1 = pcolor(subX1,subY1,X1); %Xd
hold on
colormap(BWR);
lim_yp=max(max(max(abs(X1))));
caxis([-lim_yp lim_yp])
set(hc1,'edgecolor','none')
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
hold off
ax(2)=subplot(1,2,2);
hc2 = pcolor(subX2,subY2,X2); %Xd %Xd
colormap(BWR);
lim_yp=max(max(max(abs(X2))));
caxis([-lim_yp lim_yp])
set(hc2,'edgecolor','none')
%-----------------------------------------------------------------------

% This is to confirm or check with one slice
%---------------------------------------------------------------------
[slice3, sliceInd3,subX3,subY3,subZ3] = extractSlice(G{4,2},pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);
[slice4, sliceInd4,subX4,subY4,subZ4] = extractSlice(Gp{4,2},pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);
X3 = slice3;
X4 = slice4;
f39=figure(39);
ax(3)=subplot(1,2,1);
hc3 = pcolor(subX3,subY3,X3); %Xd
colormap(BWR);
lim_yp=max(max(max(abs(X3))));
caxis([-lim_yp lim_yp])
set(hc3,'edgecolor','none')
ax(4)=subplot(1,2,2);
hc4 = pcolor(subX4,subY4,X4); %Xd %Xd
colormap(BWR);
lim_yp=max(max(max(abs(X4))));
caxis([-lim_yp lim_yp])
set(hc4,'edgecolor','none')
%}
%--------------------------------------------------------------------------



% Testing Plot the field over the slide single variable
%--------------------------------------------------------------------------
%{
Xd=Div_Poyn_v_Gp{1,1};
%Xd=Ke_Div_ve_Gp{1,1};
%Xd=dt_E_em_Gp{1,1}-Div_Poyn_v_Gp{1,1};
%Xd=dt_E_em_Gp{1,1}+JdotE_Gp{1,1}+Div_Poyn_v_Gp{1,1};
f39=figure(39);
%Yd = sign(Xd).*log10(abs(Xd));
Yd=Xd;
hc = pcolor(XLL,YLL,Yd);
cb=colormap(BWR);
set(gca,'colorscale','lin')
%lim_yp=max(max(max(abs(Yd))));
%caxis([-lim_yp lim_yp])
colorbar
set(hc,'edgecolor','none')
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
%set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
%set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
title('New FRAME')
hold off


Yd=Div_Poyn_v_Gp{1,1};
f39=figure(39);
subplot(1,3,1) 
hc = pcolor(XLL,YLL,Yd);
cb=colormap(BWR);
caxis([-0.0225 0.025])
set(gca,'colorscale','lin')
colorbar
set(hc,'edgecolor','none')
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
title('New FRAME')
hold off
subplot(1,3,2) 
Yd=Div_Poyn_v_Gp{1,1};
Kd = wiener2(Yd,[5 5]); %This appears to be a really good filter
hc = pcolor(XLL,YLL,Kd);
cb=colormap(BWR);
caxis([-0.0225 0.025])
set(gca,'colorscale','lin')
colorbar
set(hc,'edgecolor','none')
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
title('wiener2,[5 5]')
hold off
subplot(1,3,3) 
Yd=Div_Poyn_v_Gp{1,1};
Kd = medfilt2(Yd,[5,5]); %This appears to be a really good filter
hc = pcolor(XLL,YLL,Kd);
cb=colormap(BWR);
caxis([-0.0225 0.025])
set(gca,'colorscale','lin')
colorbar
set(hc,'edgecolor','none')
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
title('medfilt2')
hold off



%--------------------------------------------------------------------------
f77=figure(77)
subplot(2,2,1)
dump=squeeze(dt_E_em(:,87,:));
hc=pcolor(dump);
colorbar
%caxis([-0.05 0.05])
set(hc,'edgecolor','none')
title('dtEem')
subplot(2,2,2)
dump=squeeze(-JdotE(:,87,:));
hc=pcolor(dump);
colorbar
%caxis([-0.05 0.05])
title('-JE')
set(hc,'edgecolor','none')
subplot(2,2,3)
dump=squeeze(- Div_Poyn_v(:,87,:));
hc=pcolor(dump);
colorbar
set(hc,'edgecolor','none')
colormap(BWR);
%caxis([-0.05 0.05])
title('-DivS')
subplot(2,2,4)
dump =dt_E_em + JdotE + Div_Poyn_v;
dump=squeeze(dump(:,87,:));
hc=pcolor(dump);
colorbar
set(hc,'edgecolor','none')
colormap(BWR);
%caxis([-0.05 0.05])
title('zero')


f78=figure(78)
subplot(2,2,1)
dump=dt_E_em_Gp{1,1};
hc=pcolor(dump);
colorbar
%caxis([-0.05 0.05])
set(hc,'edgecolor','none')
title('dtEem')
subplot(2,2,2)
dump=-JdotE_Gp{1,1};
hc=pcolor(dump);
colorbar
%caxis([-0.05 0.05])
title('-JE')
set(hc,'edgecolor','none')
subplot(2,2,3)
dump= - Div_Poyn_v_Gp{1,1};
hc=pcolor(dump);
colorbar
set(hc,'edgecolor','none')
colormap(BWR);
%caxis([-0.05 0.05])
title('-DivS')
subplot(2,2,4)
dump=dt_E_em_Gp{1,1} + JdotE_Gp{1,1} + Div_Poyn_v_Gp{1,1};
%dump=Gdum_Gp{1,1};
hc=pcolor(dump);
colorbar
set(hc,'edgecolor','none')
colormap(BWR);
%caxis([-0.05 0.05])
title('Zero')
%}
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------


%
%--------------------------------------------------------------------------
% This is just for the ion pressure tensor plot in the new reference frame
%--------------------------------------------------------------------------
nc=3; mc=3;
GPijf=cell(nc,mc);
GPijp=Pij_i_Gp;
%--------------------------------------------------------------------------
% Pij_i
GPijf{1,1}=GPijp{1,2}; GPijf{1,2}=GPijp{1,3}; GPijf{1,3}=GPijp{2,3};
% Pii_i
GPijf{2,1}=GPijp{1,1}; GPijf{2,2}=GPijp{2,2}; GPijf{2,3}=GPijp{3,3};
% Pji_i
GPijf{3,1}=GPijp{2,1}; GPijf{3,2}=GPijp{3,1}; GPijf{3,3}=GPijp{3,2};
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This is just for the electron pressure tensor plot in the new reference frame
%--------------------------------------------------------------------------
GPijf_e=cell(nc,mc);
GPijp=Pij_e_Gp;
%--------------------------------------------------------------------------
% Pij_e
GPijf_e{1,1}=GPijp{1,2}; GPijf_e{1,2}=GPijp{1,3}; GPijf_e{1,3}=GPijp{2,3};
% Pii_e
GPijf_e{2,1}=GPijp{1,1}; GPijf_e{2,2}=GPijp{2,2}; GPijf_e{2,3}=GPijp{3,3};
% Pji_e
GPijf_e{3,1}=GPijp{2,1}; GPijf_e{3,2}=GPijp{3,1}; GPijf_e{3,3}=GPijp{3,2};
%--------------------------------------------------------------------------

% Pij for ion in the slide F2
%--------------------------------------------------------------------------
f42=figure(42);
for i=1:6
h1=subaxis(2,3,i,'SV',0.04,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=GPijf{1,i};
    %titl = 'iP_{ij}';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=GPijf{2,i2};
    %titl = 'iP_{ii}';
    s=i2;
end
if (i==1)
    %titl = 'iUdivV';
    titl = 'P_{12}_{i}';
elseif(i==2)
    titl = 'P_{13}_{i}';
elseif(i==3)
    titl = 'P_{23}_{i} ';
elseif(i==4)
    titl = 'P_{11}_{i}';
elseif(i==5)
    titl = 'P_{22}_{i}';
elseif(i==6)
    titl = 'P_{33}_{i}';
end 
%dum_p = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
elseif (i==4 || i==5 || i==6)
    caxis([0 lim_yp])
    colormap(h1,jet); 
end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Title.Position = [-120 258];
hcb.Ruler.Exponent = -2;
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
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
if (3<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f42,strcat('Pij_i_slideRF2_per'+ string(N_steps) +'.png'));
%saveas(f42,strcat('Pij_i_slideRF2_per_no_Bl'+ string(N_steps) +'.png'));
%cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%--------------------------------------------------------------------------


% Pij for the electrons in the slide F2
%--------------------------------------------------------------------------
f43=figure(43);
for i=1:6
h1=subaxis(2,3,i,'SV',0.04,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=GPijf_e{1,i};
    %titl = 'eP_{ij}';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=GPijf_e{2,i2};
    %titl = 'eP_{ii}';
    s=i2;
end
if (i==1)
    %titl = 'iUdivV';
    titl = 'P_{12}_{e}';
elseif(i==2)
    titl = 'P_{13}_{e}';
elseif(i==3)
    titl = 'P_{23}_{e} ';
elseif(i==4)
    titl = 'P_{11}_{e}';
elseif(i==5)
    titl = 'P_{22}_{e}';
elseif(i==6)
    titl = 'P_{33}_{e}';
end 
hc = pcolor(XLL,YLL,dum_p);
lim_yp=0.8*max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,BWR); 
elseif (i==4 || i==5 || i==6)
    caxis([0 lim_yp])
    colormap(h1,jet); 
end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Title.Position = [-120 258];
hcb.Ruler.Exponent = -2;
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%hcb.Title.String='';
%set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
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
if (3<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
%h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
%set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------

% dk_Qijk_i_Gp
%--------------------------------------------------------------------------
f44=figure(44);
for i=1:9
h1=subaxis(3,3,i,'SV',0.035,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=dk_Qijk_i_Gp{1,i};
    titl = 'idkQijk';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=dk_Qijk_i_Gp{2,i2};
    titl = 'idkQijk';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=dk_Qijk_i_Gp{3,i3};
    titl = 'idkQijk';
    s=i3;
end
if (i==1)
    %titl = 'iUdivV';
    titl = 'dQ_{11,i}';
elseif(i==2)
    titl = 'dQ_{12,i}';
elseif(i==3)
    titl = 'dQ_{13,i} ';
elseif(i==4)
    titl = 'dQ_{21,i}';
elseif(i==5)
    titl = 'dQ_{22,i}';
elseif(i==6)
    titl = 'dQ_{23,i}';
elseif(i==7)
    titl = 'dQ_{31,i}';
elseif(i==8)
    titl = 'dQ_{32,i}';
elseif(i==9)
    titl = 'dQ_{33,i}';
end 
%dum_p = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
%if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,jet(dCp)); 
%elseif (i==4 || i==5 || i==6)
   % caxis([-lim_yp lim_yp])
   % colormap(h1,BWR); 
%end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Title.Position = [-100 185];
hcb.Ruler.Exponent = -2;
%hcb.Title.String='';
set(get(hcb,'Title'),'String',titl ,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
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
if (6<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2, \ sg(y)log(y)$$','Interpreter','latex','FontSize',20)
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------


% dk_Qijk_e_Gp
%--------------------------------------------------------------------------
f45=figure(45);
for i=1:9
h1=subaxis(3,3,i,'SV',0.035,'SH',0.004,'MR',0.02,'ML',0.05,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=dk_Qijk_e_Gp{1,i};
    titl = 'edkQijk';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=dk_Qijk_e_Gp{2,i2};
    titl = 'edkQijk';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=dk_Qijk_e_Gp{3,i3};
    titl = 'edkQijk';
    s=i3;
end

if (i==1)
    %titl = 'iUdivV';
    titl = 'dQ_{11,e}';
elseif(i==2)
    titl = 'dQ_{12,e}';
elseif(i==3)
    titl = 'dQ_{13,e} ';
elseif(i==4)
    titl = 'dQ_{21,e}';
elseif(i==5)
    titl = 'dQ_{22,e}';
elseif(i==6)
    titl = 'dQ_{23,e}';
elseif(i==7)
    titl = 'dQ_{31,e}';
elseif(i==8)
    titl = 'dQ_{32,e}';
elseif(i==9)
    titl = 'dQ_{33,e}';
end 
%dum_p = sign(dum_p).*log10(abs(dum_p));
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
set(hc,'edgecolor','none')
%if (i==1 || i==2 || i==3)
    caxis([-lim_yp lim_yp])
    colormap(h1,jet(dCp)); 
%elseif (i==4 || i==5 || i==6)
   % caxis([-lim_yp lim_yp])
   % colormap(h1,BWR); 
%end
hcb=colorbar;
hcb.Label.FontSize = 14;
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Title.Position = [-100 185];
hcb.Ruler.Exponent = -2;
%hcb.Title.String=titl + string(s);
set(get(hcb,'Title'),'String',titl,'FontSize',14)
%set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
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
if (6<i)
    ax.XLabel.String = '$r / d_{i}$';
    ax.XTick = [2 4 6 8];
end
hold on
%plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', colorFVBl);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xll(1) xll(end)]);
ylim([yll(1) yll(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
%sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2, \ sg(y)log(y)$$','Interpreter','latex','FontSize',20)
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}, \ Slide RF2$$','Interpreter','latex','FontSize',20)

%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
%{
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f40,strcat('vi_ve_J_B_vev_Bl_slideRF1_per'+ string(N_steps) +'.png'));
%saveas(f41,strcat('vi_ve_J_B_vev_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%saveas(f41,strcat('vi_ve_J_B_vev_no_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%saveas(f11,strcat('vi_ve_J_B_vev_slideRF2_no_bl_per'+ string(N_steps) +'.png'));
saveas(f1,strcat('vi_ve_J_B_vev_Bl_slideRF2_per'+ string(N_steps) +'.png'));
saveas(f2,strcat('ni_ne_slideRF2_bl_per'+ string(N_steps) +'.png'));
saveas(f3,strcat('E_poyntin_slideRF2_bl_pe'+ string(N_steps) +'.png'));
saveas(f4,strcat('Kinetic_slideRF2_bl_pe_log_'+ string(N_steps) +'.png'));
saveas(f5,strcat('Potential_slideRF2_bl_pe_log_'+ string(N_steps) +'.png'));
saveas(f6,strcat('Total derivatives_slideRF2_bl_pe_log_'+ string(N_steps) +'.png'));
saveas(f42,strcat('Pij_i_slideRF2_bl_per'+ string(N_steps) +'.png'));
saveas(f43,strcat('Pij_e_slideRF2_bl_per'+ string(N_steps) +'.png'));
saveas(f44,strcat('dkQijk_i_slideRF2_bl_per_'+ string(N_steps) +'.png'));
saveas(f45,strcat('dkQijk_e_slideRF2_bl_per_'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%--------------------------------------------------------------------------
%}


%--------------------------------------------------------------------------

% For the magnetic field lines in 3D
%--------------------------------------------------------------------------
%{
u=B{1,1};
v=B{1,2};
w=B{1,3}-0.1;
[xm3,ym3, zm3] = meshgrid(xl,yl,zl);
%[sx,sy,sz] = meshgrid(10:2:20,10:2:24,72:2:82);
dx=xl(2)-xl(1); dy=yl(2)-yl(1); dz=zl(2)-zl(1); 
[sx,sy,sz] = meshgrid(xl(1):8*dx:xl(end),yl(1):8*dy:yl(end),zl(83));
%[sx,sy,sz] = meshgrid(xl,yl,72);
XYZ2=stream3(xm3,ym3,zm3,u,v,w,sx,sy,sz);
dum_j=J{1,3};

f68=figure(68);
hc=pcolor(xl,yl,dum_j(:,:,slide_n));
%hc.ZData = zl(slide_n) + hc.ZData; 
hc.ZData = 77 + hc.ZData; 
set(hc,'edgecolor','none')
hold on
hlines=streamline(XYZ2);
view(3);
hold off

dum_e=E{1,3};
f68=figure(68);
hc=pcolor(xl,yl,dum_e(:,:,slide_n));
hc.ZData = zl(slide_n) + hc.ZData; 
set(hc,'edgecolor','none')
hold on
hlines=streamline(XYZ2);
view(3);
hold off
%--------------------------------------------------------------------------


%Now, Calculating the squshing factor
%--------------------------------------------------------------------------
%the indices run along y direction first and then along x direction
xyzinit = XYZ2(1,:);
xyzend = XYZ2(end,:);

xyzinit = reshape(xyzinit,[21,21]); %it is like this (x1,:,z)
xyzend = reshape(xyzend,[21,21]);


Xx = xyzinit; Xx=cell2mat(Xx); 

Xx = xyzinit(1); Xx=cell2mat(Xx); 
xyzinit = XYZ2{22}(1,:)

xyzinit = XYZ2{144}(1,:)
xyzend = XYZ2{144}(end,:)

colorstring = 'kbgry';
for i = 1:5
  set(hlines(i), 'Color', colorstring(i))
end
%}
%--------------------------------------------------------------------------


% 2d histograms
%--------------------------------------------------------------------------
%{
%--------------------------------------------------------------------------
Nx=zeros(100,168);
nbins=100;
edgesz = linspace(zl(1),zl(end),168);


% For the z component of the magnetic field
f64=figure(64);
ax(1)=subplot(2,2,1);
dum_p=B{1,3}-0.1;
maxd=max(max(max(abs(dum_p))));
edges = linspace(-maxd,maxd,101);
for slide_n=1:168
dmh=dum_p(:,:,slide_n);
[Nx(:,slide_n),edgesx] = histcounts(dmh,edges);
end
[EZ, EX]=meshgrid(edgesz(1:168),edgesx(1:100)); 
hpc=pcolor(EZ, EX, Nx);
set(hpc,'edgecolor','none')
% For the z component of the ELECTRIC field
ax(2)=subplot(2,2,2);
dum_p=E{1,3};
maxd=max(max(max(abs(dum_p))));
edges = linspace(-maxd,maxd,101);
for slide_n=1:168
dmh=dum_p(:,:,slide_n);
[Nx(:,slide_n),edgesx] = histcounts(dmh,edges);
end
[EZ, EX]=meshgrid(edgesz(1:168),edgesx(1:100)); 
hpc=pcolor(EZ, EX, Nx);
set(hpc,'edgecolor','none')
% For the electron cross helicity
ax(3)=subplot(2,2,3);
dum_p=ve{1,1}.*B{1,1} + ve{1,2}.*B{1,2} + ve{1,3}.*B{1,3};
maxd=max(max(max(abs(dum_p))));
edges = linspace(-maxd,maxd,101);
for slide_n=1:168
dmh=dum_p(:,:,slide_n);
[Nx(:,slide_n),edgesx] = histcounts(dmh,edges);
end
[EZ, EX]=meshgrid(edgesz(1:168),edgesx(1:100)); 
hpc=pcolor(EZ, EX, Nx);
set(hpc,'edgecolor','none')

% This does the 1d plots
plot(edgesx(2:end),Nx)
%--------------------------------------------------------------------------

%}
%--------------------------------------------------------------------------


% Other Plots energy, tempereature, etc%%
%--------------------------------------------------------------------------
%{

%Lets plot the electromagnetic energy and its contributions
%------------------------------------------------------------------
f4=figure(4);
for i=1:3
h1=subaxis(2,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
if (i<=3)
    i1=i;
    dum_p=Poyn_v{1,i};
    titl = 'S';
    s=i1;
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
caxis([-lim_yp lim_yp])
hcb=colorbar;
colormap(BWR);
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[3 0]})
%pbaspect([1 1 1])
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (9<i)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%-------------------------------------------------------------------------
h1=subaxis(2,3,4,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
dum_p=E_em;
titl = 'Eem';
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
hcb=colorbar;
colormap(h1,jet);
caxis([0 lim_yp])
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[3 0]})
%pbaspect([1 1 1])
axis tight
ax = gca;
set(ax,'YDir','normal')
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
%--------------------------------------------------------------------------
h1=subaxis(2,3,5,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
dum_p=Div_Poyn_v;
titl = 'DivS';
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
hcb=colorbar;
colormap(h1,BWR);
caxis([-lim_yp lim_yp])
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[3 0]})
%pbaspect([1 1 1])
axis tight
ax = gca;
set(ax,'YDir','normal')
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
%------------------------------------------------------------------
h1=subaxis(2,3,6,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
dum_p=dt_E_em;
titl = 'dtEem';
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
hcb=colorbar;
colormap(h1,BWR);
caxis([-lim_yp lim_yp])
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl,0,[3 0]})
%pbaspect([1 1 1])
axis tight
ax = gca;
set(ax,'YDir','normal')
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
%------------------------------------------------------------------

%Tij_i
%------------------------------------------------------------------
f5=figure(5);
k=1;
for i=1:9
h1=subaxis(2,3,k,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
if (i<=3)
    i1=i;
    dum_p=Tij_i{1,i};
    titl = 'Ti1';
    s=i1;
    k=i+1;
elseif (i==4)
    continue
elseif (4<i && i<=6)
    i2=i-3;
    dum_p=Tij_i{2,i2};
    titl = 'Ti2';
    s=i2;
    k=k+1;
elseif (6<i && i<=8)
    continue
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=Tij_i{3,i3};
    titl = 'Ti3';
    s=i3;
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
hcb=colorbar;
if (i==1 || i==5 || i==9)
   colormap(h1,jet); 
   caxis([0 lim_yp])
else
   colormap(h1,BWR);
   caxis([-lim_yp lim_yp])
end
hcb.Label.Interpreter = 'latex';
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[5 0]})
%pbaspect([1 1 1])
axis tight
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (i==1 || i==5)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (3<k)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%------------------------------------------------------------------


%Tij_e
%------------------------------------------------------------------
f6=figure(6);
k=1;
for i=1:9
h1=subaxis(2,3,k,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
if (i<=3)
    i1=i;
    dum_p=Tij_e{1,i};
    titl = 'Te1';
    s=i1;
    k=i+1;
elseif (i==4)
    continue
elseif (4<i && i<=6)
    i2=i-3;
    dum_p=Tij_e{2,i2};
    titl = 'Te2';
    s=i2;
    k=k+1;
elseif (6<i && i<=8)
    continue
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=Tij_e{3,i3};
    titl = 'Te3';
    s=i3;
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
hcb=colorbar;
if (i==1 || i==5 || i==9)
   colormap(h1,jet); 
   caxis([0 lim_yp])
else
   colormap(h1,BWR);
   caxis([-lim_yp lim_yp])
end
hcb.Label.Interpreter = 'latex';
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[5 0]})
%pbaspect([1 1 1])
axis tight
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (i==1 || i==5)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (3<k)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%------------------------------------------------------------------


%dk_Qij_i    
%------------------------------------------------------------------    
f7=figure(7);
k=1;
for i=1:9
h1=subaxis(2,3,k,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
if (i<=3)
    i1=i;
    dum_p=dk_Qijk_i{1,i};
    titl = 'dkQi1';
    s=i1;
    k=i+1;
elseif (i==4)
    continue
elseif (4<i && i<=6)
    i2=i-3;
    dum_p=dk_Qijk_i{2,i2};
    titl = 'dkQi2';
    s=i2;
    k=k+1;
elseif (6<i && i<=8)
    continue
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=dk_Qijk_i{3,i3};
    titl = 'dkQi3';
    s=i3;
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
hcb=colorbar;
if (i==1 || i==5 || i==9)
   colormap(h1,jet); 
   caxis([0 lim_yp])
else
   colormap(h1,BWR);
   caxis([-lim_yp lim_yp])
end
hcb.Label.Interpreter = 'latex';
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[5 0]})
%pbaspect([1 1 1])
axis tight
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (i==1 || i==5)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (3<k)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%------------------------------------------------------------------

%dk_Qij_e    
%------------------------------------------------------------------    
f8=figure(8);
k=1;
for i=1:9
h1=subaxis(2,3,k,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
if (i<=3)
    i1=i;
    dum_p=dk_Qijk_e{1,i};
    titl = 'dkQe1';
    s=i1;
    k=i+1;
elseif (i==4)
    continue
elseif (4<i && i<=6)
    i2=i-3;
    dum_p=dk_Qijk_e{2,i2};
    titl = 'dkQe2';
    s=i2;
    k=k+1;
elseif (6<i && i<=8)
    continue
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=dk_Qijk_e{3,i3};
    titl = 'dkQe3';
    s=i3;
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
hcb=colorbar;
if (i==1 || i==5 || i==9)
   colormap(h1,jet); 
   caxis([0 lim_yp])
else
   colormap(h1,BWR);
   caxis([-lim_yp lim_yp])
end
hcb.Label.Interpreter = 'latex';
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[5 0]})
%pbaspect([1 1 1])
axis tight
ax = gca;
set(ax,'YDir','normal')
ax.XTick = [];
ax.YTick = [];
if (i==1 || i==5)
    ax.YLabel.String = 'y / d_{i}';
    ax.YTick = [16 20 24];
end
if (3<k)
    ax.XLabel.String = 'x / d_{i}';
    ax.XTick = [10 15 20];
end
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%------------------------------------------------------------------



%------------------------------------------------------------------

hold on
sv=10;
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), u(1:sv:end,1:sv:end), v(1:sv:end,1:sv:end), 0);
%h2=quiver3(xmm(1:sv:endn,1:sv:endn,1:sv:endn),ymm(1:sv:endn,1:sv:endn,1:sv:endn),...
 %   zmm(1:sv:endn,1:sv:endn,1:sv:endn),dp_x(1:sv:endn,1:sv:endn,1:sv:endn),dp_y(1:sv:endn,1:sv:endn,1:sv:endn),dp_z(1:sv:endn,1:sv:endn,1:sv:endn),0);
hlines=streamline(xm,ym,u,v,xstart,ystart);
set(h2,'AutoScale','on', 'AutoScaleFactor', 1.5, 'Color', 'k');
set(hlines, 'Color', 'k');
[M,c] = contour(xl,yl,w,[mean(w,'all') 0.5*mean(w,'all')]); 
c.LineWidth = 1;
c.LineColor = 'k';
hold off

%-----------------------------------------------------------------------
f9=figure(9)
dum_p=J{1,3};
%[M,c] = contour(xl,yl,dum_p(:,:,slide_n)','ShowText','on');
[M,c] = contour(xl,yl,dum_p(:,:,slide_n)',[0.001 0.0015]);
c.LineWidth = 2;
colormap(BWR);
hold on
quiver3(xmm(1:sv:endn,1:sv:endn,1:sv:endn),ymm(1:sv:endn,1:sv:endn,1:sv:endn),...
    zmm(1:sv:endn,1:sv:endn,1:sv:endn),dp_x(1:sv:endn,1:sv:endn,1:sv:endn),dp_y(1:sv:endn,1:sv:endn,1:sv:endn),dp_z(1:sv:endn,1:sv:endn,1:sv:endn))
hold off

[xm,ym] = meshgrid(xl,yl);
dum_p_x=B{1,1}; dum_p_y=B{1,2}; dum_p_z=B{1,3}; dum_p_jz=J{1,3};
u=dum_p_x(:,:,slide_n)'; v=dum_p_y(:,:,slide_n)'; w=dum_p_z(:,:,slide_n)';
wj=dum_p_jz(:,:,slide_n)';


du_x=J{1,3}; du_y=J{1,3}; du_z=B{1,3}; du_jz=J{1,3};
du_u=du_x(:,:,slide_n)'; du_v=du_y(:,:,slide_n)'; du_w=du_z(:,:,slide_n)';
%rng('default')
%s = rng
f13=figure(13);
%quiver(xl,yl,u,v)
%N=100; xstart = max(xl)*rand(N,1); ystart = max(yl)*rand(N,1); 
%streamline(xm,ym,u,v,xstart,ystart)
streamline(xm,ym,du_u,du_v,xstart,ystart)
hold on
[M,c] = contour(xl,yl,wj);
%[M,c] = contour(xl,yl,w);
hold off

f11=figure(11);
quiver3(xl,yl,xl,u,v,w)
%N=100; xstart = max(xl)*rand(N,1); ystart = max(yl)*rand(N,1); 
%streamline(xm,ym,u,v,xstart,ystart)
%quiver3(X,Y,Z,U,V,W)
hold on
[M,c] = contour(xl,yl,w);
hold off


hold on
sv=10;
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), u(1:sv:end,1:sv:end), v(1:sv:end,1:sv:end), 0);
hlines=streamline(xm,ym,u,v,xstart,ystart);
set(h2,'AutoScale','on', 'AutoScaleFactor', 1.5)
set(hlines, 'Color', 'k');
hold off

endn=168;
xstart = max(xl)*rand(N,1); ystart = max(yl)*rand(N,1); zstart = max(yl)*rand(N,1); 
[sx,sy,sz] = meshgrid(xstart,ystart,zstart);

[xmm, ymm, zmm] = meshgrid(xl,yl,xl);
f44=figure(44)
dp_x=B{1,1}; dp_y=B{1,2}; dp_z=B{1,3}; 
quiver3(xmm(1:sv:endn,1:sv:endn,1:sv:endn),ymm(1:sv:endn,1:sv:endn,1:sv:endn),...
    zmm(1:sv:endn,1:sv:endn,1:sv:endn),dp_x(1:sv:endn,1:sv:endn,1:sv:endn),dp_y(1:sv:endn,1:sv:endn,1:sv:endn),dp_z(1:sv:endn,1:sv:endn,1:sv:endn))
axis equal

%--------------------------------------------------------------------------

%Save the plots
%--------------------------------------------------------------------------
saveas(f1,strcat('vi_ve_J_B','.png'));
saveas(f2,strcat('Pij_i','.png'));
saveas(f3,strcat('Pij_e','.png'));
saveas(f4,strcat('EM_energy','.png'));
saveas(f5,strcat('Tij_i','.png'));
saveas(f6,strcat('Tij_e','.png'));
%--------------------------------------------------------------------------

%}
%--------------------------------------------------------------------------


% Practice plots
%--------------------------------------------------------------------------
%{
%--------------------------------------------------------------------------
%Non presenting plots
%------------------------------------------------------------------
f5=figure(5);
title('$\omega$ vs $k$','Interpreter','latex')
for i=1:12
h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'padding', 0.02, 'margin', 0.02);
imagesc(magic(2*i))
ax=gca;
end
%------------------------------------------------------------------

%------------------------------------------------------------------
x=linspace(0,2*pi,100);
h=(2*pi/100);
y=sin(x);
dy1=diff(y)/h;
dy2=gradient(y)/h;
f3=figure(3);
plot(x(1:end-1),y(1:end-1),'k')
hold on
plot(x(1:end-1),dy1,'r')
plot(x(1:end-1),dy2(1:end-1),'b')
legend('y','dy1','dy2')
hold off

f2=figure(2);
histogram(ni,100)
hold on
histogram(ne,100)
%plot(x(1:end-1),dy2(1:end-1),'b')
legend('ni','ne')
hold off

%------------------------------------------------------------------

%------------------------------------------------------------------

 f3=figure(3);
 %ax(1)=subplot(2,1,1);
 dum_p=vi{1,1};
 h=pcolor(dum_p(:,:,1));
 hold on
 set(h, 'EdgeColor', 'none');
 lim = caxis;  
 
 caxis([-23.5 36])
 cmap = [pink(23.5*2); gray(36*2)];
 colormap(cmap)
 
 %caxis([Brange 1]); 
 %colormap(colors_p3)
 set(gca,'ColorScale','lin','FontSize',18)
 colorbar
 %xlim([-1.2 1.2]); ylim([0 1.2]);
 set(gca,'XScale','lin','YScale','lin','FontSize',18)
 xlabel('$$v_{e\|}$$','Interpreter','latex','FontSize',20)
 ylabel('$$v_{e\perp}$$','Interpreter','latex','FontSize',20)
 title('$$\Delta VDF_{e} / VDF_{e,max0}$$','Interpreter','latex','FontSize',20)
 %[M,c1]=contour(X,Y,Z,v); c1.LineWidth = 3;
 %xline(vthe_c_par_,'--k','LineWidth',1);
 %xline(-vthe_c_par_,'--k','LineWidth',1);
 %[M,c2]=contour(X,Y,Z,v2); c2.LineWidth = 3;
 grid off
 %pbaspect([1 0.5 0.5])
 axis square
 hold off


 % Create sample data:
f3=figure(4);
y_plot=dum_p(:,:,4);
imagesc(y_plot);
lim_yp=max(max(max(abs(dum_p(:,:,:)))));
caxis([-lim_yp lim_yp])
colorbar
%colorMap = [redColorMap; blueColorMap; zeros(1, 256)]';
colormap(BWR);
set(gca,'YDir','normal') % this is the one that reverse the yaxis


 % Create sample data:
correlations = peaks(300);
minValue = min(correlations(:));
maxValue = max(correlations(:));
% Scale the data to between -1 and +1.
correlations = (correlations-minValue) * 2 / (maxValue - minValue) - 1;
% Display - will use some weird color map to start with.
imagesc(correlations);
colorbar
% Create colormap that is green for negative, red for positive,
% and a chunk inthe middle that is black.
%greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
% Apply the colormap.
colormap(colorMap);
 
load clown
clims = [10 60];
imagesc(flipud(dum_p(:,:,1)),clims)
colormap(jet)
set(gca,'YDir','normal')
 


cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
f1=figure(1)
bwr1(1)=subplot(1,2,1);
imagesc(peaks(150), [0 10])
bwr1(2)=subplot(1,2,2);
imagesc(peaks(150), [0 10])
colorbar(bwr1(1));colormap(bwr1(1),hot)
colorbar(bwr1(2));colormap(bwr1(2),BWR)

figure;
imagesc(peaks(150));
colormap(BWR(20)), colorbar

%}
%--------------------------------------------------------------------------

% Extra stuff
%--------------------------------------------------------------------------
%{
f37=figure(37);
surf(subX,subY,subZ,slice,'FaceColor','texturemap','EdgeColor','none');
colormap(BWR);
caxis([-1 1])
%[Bob,x,y,z] = obliqueslice(G,point,normal); This on;y works with the
%package
%------------------------------------------------------------------------

[xmm, ymm, zmm] = meshgrid(xl,yl,xl);
dp_x=B{1,1}; dp_y=B{1,2}; dp_z=B{1,3}; 
endn=168;
f34=figure(34);
% Define 3D data
[x,y,z] = meshgrid(linspace(12,yl(end),168));
v = x.^2 + y.^2 + z.^2;
% Define the slice plane
[xi, yi] = meshgrid(linspace(12,yl(end),200));
zi = yi + 3 ;%(xi+yi)./2; Note: yl with zi-xi
% Slice it
hs=slice(x,y,z,J{1,3},xi,yi,zi);
set(hs,'edgecolor','none')
colormap(BWR);
caxis([-1 1])


f35=figure(35);
% Define 3D data
dd=J{1,3};
hs=pcolor(dd(:,:,slide_n)'); %This is the same (zi=yi)
set(hs,'edgecolor','none')
colormap(BWR);
caxis([-1 1])

hold on 
h2=quiver3(xmm(1:sv:endn,1:sv:endn,1:sv:endn),ymm(1:sv:endn,1:sv:endn,1:sv:endn),...
    zmm(1:sv:endn,1:sv:endn,1:sv:endn),dp_x(1:sv:endn,1:sv:endn,1:sv:endn),dp_y(1:sv:endn,1:sv:endn,1:sv:endn),dp_z(1:sv:endn,1:sv:endn,1:sv:endn),0);
hold off

f34=figure(34);
% Define 3D data
[x,y,z] = meshgrid(linspace(xl(1),xl(end),168));
v = x.^2 + y.^2 + z.^2;
% Define the slice plane
[xi, yi] = meshgrid(linspace(xl(1),xl(end),200));
%zi =  (zl(slide_n)-zl(1)) + 0.001.*xi;
zi =  10 + 0.001.*xi;
% Slice it
hs=slice(x,y,z,vi{1,3},xi,yi,zi);
set(hs,'edgecolor','none')


[xm3,ym3, zm3] = meshgrid(xl,yl,zl);
[xi, yi] = meshgrid(xl,yl);
zi = xi;
%hs=slice(x,y,z,vi{1,3},xi,yi,zi);
hs=slice(xm3,ym3,zm3,vi{1,3},xi,yi,zi);
set(hs,'edgecolor','none')

%-------------------------------------------------------------------------
%i-th xy slice:
G=dd;
G_yx=G(:,:,slide_n); % Y-by-X array
%i-th xz slice:
G_xz=permute(G(slide_n,:,:),[2 3 1]); % X-by-Z array
%i-th yz slice:
G_yz=permute(G(:,slide_n,:),[1 3 2]); % Y-by-Z array

f36=figure(36);
% Define 3D data
hxy=pcolor(G_yz); %This is the same (zi=yi)
set(hxy,'edgecolor','none')
colormap(BWR);
caxis([-1 1])
hold on
hxz=pcolor(G_xz); %This is the same (zi=yi)
set(hxz,'edgecolor','none')
hyz=pcolor(G_yz); %This is the same (zi=yi)
set(hyz,'edgecolor','none')
hold off
%-----------------------------------------------------------

%This calculates the rotation matrix between two vectors 
%--------------------------------------------------------------------------
% two random 3D vectors
p0 = randi(10,3,1);
p1 = randi(10,3,1);
% calculate cross and dot products
C = cross(p0, p1) ; 
D = dot(p0, p1) ;
NP0 = norm(p0) ; % used for scaling
if ~all(C==0) % check for colinearity    
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
    R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
else
    R = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
end
% R is the rotation matrix from p0 to p1, so that (except for round-off errors) ...
R * p0      % ... equals p1 
inv(R) * p1 % ... equals p0
%--------------------------------------------------------------------------
%}
%--------------------------------------------------------------------------
%}
