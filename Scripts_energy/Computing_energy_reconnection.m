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

%for i=2:30
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
ni = h5read(fileID,strcat(n_1st,'/n_i/p0/3d'),start,count); %ni=ni(nx1:nx2,ny1:ny2,nz1:nz2);
ne = h5read(fileID,strcat(n_1st,'/n_e/p0/3d'),start,count); %ne=ne(nx1:nx2,ny1:ny2,nz1:nz2);
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

%vectors as 1 arrow 3 columns

%mi=1;
mime=100;
%me=mi/mime;
qi=1;
qe=-1;

% This is to unlock. These are the lines to calculate the energy terms
%{  

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
dt_E_em = JdotE - Div_Poyn_v;  

% Time change of the electromagnetic energy 
%------------------------------------------------------------------

%}

%HEre I am changing stuff...


% This part is to make the plots to see the data
%------------------------------------------------------------------

cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';

xl=linspace(nx1,nx1+px2,px2)*resx;
yl=linspace(ny1,ny1+py2,py2)*resy;
zl=linspace(nz1,nz1+pz2,pz2)*resz;

%Pick the slide number along z-axes to get the plots
slide_n =83;

%----------------------------------------
[xm,ym] = meshgrid(xl,yl);

[xm3,ym3, zm3] = meshgrid(xl,yl,zl);

% For the magnetic field lines and quivers considering the projection on
% the plane. This works only for the xy plane.
dum_p_x=B{1,1}; dum_p_y=B{1,2}; dum_p_z=B{1,3};
ub=dum_p_x(:,:,slide_n)'; vb=dum_p_y(:,:,slide_n)'; wb=dum_p_z(:,:,slide_n)';
bm=sqrt(ub.*ub + vb.*vb + wb.*wb);
ubx=ub.*sqrt(1-((wb./bm).*(wb./bm)));
vby=vb.*sqrt(1-((wb./bm).*(wb./bm)));

% For the electron velocity lines and quivers
dum_p_x=ve{1,1}; dum_p_y=ve{1,2}; dum_p_z=ve{1,3}; 
duex=dum_p_x(:,:,slide_n)'; duey=dum_p_y(:,:,slide_n)'; duez=dum_p_z(:,:,slide_n)';
vem=sqrt(duex.*duex + duey.*duey + duez.*duez);
uex=duex.*sqrt(1-((duez./vem).*(duez./vem)));
vey=duey.*sqrt(1-((duez./vem).*(duez./vem)));

% For the Jz
dum_p_jz=J{1,3};
wj=dum_p_jz(:,:,slide_n)';
%----------------------------------------

%----------------------------------------
%hold on
rng('default')
s = rng;
NN=80; xstart = max(xl)*rand(NN,1); ystart = max(yl)*rand(NN,1); %Do this just once
sv=7;
%streamline(xm,ym,u,v,xstart,ystart)
%[M,c] = contour(xl,yl,w)
%hold off
%----------------------------------------

% old plots
%--------------------------------------------------------------------------
%{
% Plot of Vi, Ve, J, B
%------------------------------------------------------------------
f1=figure(1);
for i=1:12
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(4,3,i,'SV',0,'SH',0,'MR',0,'ML',0.04,'PL',0.002,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=vi{1,i};
    titl = 'vi';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=ve{1,i2};
    titl = 've';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=J{1,i3};
    titl = 'J';
    s=i3;
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=B{1,i4};
    if (i4==3)
        dum_paux=B{1,i4};
        dum_p = dum_paux - 0.1;
    end
    titl = 'B';
    s=i4;    
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,slide_n)))));
caxis([-lim_yp lim_yp])
hcb=colorbar;
colormap(BWR);
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
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
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), uex(1:sv:end,1:sv:end), vey(1:sv:end,1:sv:end), 0);
hlines=streamline(xm,ym,ubx,vby,xstart,ystart);
set(hlines, 'Color', [0.5 0.5 0.5]);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xl(1) xl(end)]);
ylim([yl(1) yl(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1} \ ; \ xy-plane$$','Interpreter','latex','FontSize',20)
%------------------------------------------------------------------



%------------------------------------------------------------------
%length = 20;
%red = [1, 0, 0]; white = [1, 1, 1]; blue = [0, 0, 1];
%pink = [255, 192, 203]/255;
%colors_p1 = [linspace(white(1),red(1),length)', linspace(white(2),red(2),length)', linspace(white(3),red(3),length)'];
%colors_p2 = [linspace(white(1),blue(1),length)', linspace(white(2),blue(2),length)', linspace(white(3),blue(3),length)'];
%------------------------------------------------------------------

%Pij_i    
%------------------------------------------------------------------    
f2=figure(2);
k=1;
for i=1:9
%h1=subaxis(2,3,k,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(2,3,k,'SV',0.002,'SH',0.002,'MR',0.002,'ML',0.04,'PL',0.002,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=Pij_i{1,i};
    titl = 'Pi1';
    s=i1;
    k=i+1;
elseif (i==4)
    continue
elseif (4<i && i<=6)
    i2=i-3;
    dum_p=Pij_i{2,i2};
    titl = 'Pi2';
    s=i2;
    k=k+1;
elseif (6<i && i<=8)
    continue
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=Pij_i{3,i3};
    titl = 'Pi3';
    s=i3;
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
%contourf(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,slide_n)))));
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
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
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
hold on
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), uex(1:sv:end,1:sv:end), vey(1:sv:end,1:sv:end), 0);
hlines=streamline(xm,ym,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xl(1) xl(end)]);
ylim([yl(1) yl(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}\ ; \ xy-plane$$','Interpreter','latex','FontSize',20)
%------------------------------------------------------------------

%Pij_e
%------------------------------------------------------------------
f3=figure(3);
k=1;
for i=1:9
%h1=subaxis(2,3,k,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(2,3,k,'SV',0.002,'SH',0.002,'MR',0.002,'ML',0.04,'PL',0.002,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=Pij_e{1,i};
    titl = 'Pe1';
    s=i1;
    k=i+1;
elseif (i==4)
    continue
elseif (4<i && i<=6)
    i2=i-3;
    dum_p=Pij_e{2,i2};
    titl = 'Pe2';
    s=i2;
    k=k+1;
elseif (6<i && i<=8)
    continue
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=Pij_e{3,i3};
    titl = 'Pe3';
    s=i3;
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,slide_n)))));
hcb=colorbar;
if (i==1 || i==5 || i==9)
   colormap(h1,jet); 
   %caxis([0 lim_yp])
   caxis([0 0.01])
else
   colormap(h1,BWR);
   %caxis([-lim_yp lim_yp])
   caxis([-0.005 0.005])
end
hcb.Label.Interpreter = 'latex';
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.6 0]})
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
hold on
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), uex(1:sv:end,1:sv:end), vey(1:sv:end,1:sv:end), 0);
hlines=streamline(xm,ym,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xl(1) xl(end)]);
ylim([yl(1) yl(end)]);
hold off
end
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}\ ; \ xy-plane$$','Interpreter','latex','FontSize',20)
%}
%--------------------------------------------------------------------------

%These are the better plots 
%--------------------------------------------------------------------------

%Vi Ve J B
%--------------------------------------------------------------------------
%
f1=figure(1);
for i=1:12
h1=subaxis(4,3,i,'SV',0.002,'SH',0.004,'MR',0.002,'ML',0.04,'PL',0.002,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=vi{1,i};
    titl = 'vi';
    s=i1;
    tic=0.07;  
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=ve{1,i2};
    titl = 've';
    s=i2;
    if(i==4 || i==5)
    tic=0.25; 
    elseif(i==6)
    tic=1;   
    end
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=J{1,i3};
    titl = 'J';
    s=i3;
    if(i==7 || i==8)
    tic=0.2;  
    elseif(i==9)
    tic=0.8;  
    end
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=B{1,i4};
    if (i4==3)
        dum_paux=B{1,i4};
        dum_p = dum_paux - 0.1;
    end
    titl = 'B';
    s=i4;
    if(i==10 || i==11)   
    tic=0.12;  
    elseif(i==12)
    tic=0.07;  
    end
end
imagesc(xl,yl,dum_p(:,:,slide_n)')
lim_yp=max(max(max(abs(dum_p(:,:,slide_n)))));
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
h2=quiver(xl(1:sv:end,1:sv:end), yl(1:sv:end,1:sv:end), uex(1:sv:end,1:sv:end), vey(1:sv:end,1:sv:end), 0);
hlines=streamline(xm,ym,ubx,vby,xstart,ystart);
set(hlines,'LineWidth',1.3, 'Color', [0.5 0.5 0.5]);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
xlim([xl(1) xl(end)]);
ylim([yl(1) yl(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
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
%}
%------------------------------------------------------------------

%{
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f1,strcat('vi_ve_J_B_vev_Bl_cte_'+ string(N_steps) +'.png'));
saveas(f2,strcat('Pij_i_vev_Bl_cte_'+ string(N_steps) +'.png'));
saveas(f3,strcat('Pij_e_vev_Bl_cte_'+ string(N_steps) +'.png'));
cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';
%}
%end
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% This sligthly works to produce the slice. However, I still need to 
% carefully check the limits or dimensions so they match the simulation ones.
% Also, it is necesary to define the correct plane and see how to plot the quick vec and streamlines 
% in that plane 
%--------------------------------------------------------------------------


% This helps to display several slides of the same quantity 
%-----------------------------------------------------------
%f37=figure(37);
%montage(G','Size',[9 9]);
%colormap(BWR);
%--------------------------------------------------------


% This is to get the quantities in the slide that is related to the correct plane
% or it seems that way.
%--------------------------------------------------------

% Choose the field to get the slide from.
dd_Bz = B{1,3}; dd_Bz0 = dd_Bz - 0.1;
nc = 4; mc = 3;

G=cell(nc,mc);
%--------------------------------------------------------------------------
%
% Magnetic field
G{1,1}=B{1,1}; G{1,2}=B{1,2}; G{1,3}=dd_Bz0;
% Electron velocity
G{2,1}=ve{1,1}; G{2,2}=ve{1,2}; G{2,3}=ve{1,3};
% Ion velocity
G{3,1}=vi{1,1}; G{3,2}=vi{1,2}; G{3,3}=vi{1,3};
% current
G{4,1}=J{1,1}; G{4,2}=J{1,2}; G{4,3}=J{1,3};
%}
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
pt = [x_pi y_pi z_pi];
%pt = [0 0 0];

% This is the normal vector to the plane perpendicular to the current (Vec adjacent)
vec_a = [1 0 -3.3333];
mag_va = sqrt(vec_a(1)*vec_a(1) + vec_a(2)*vec_a(2) + vec_a(3)*vec_a(3));
% This is the normal vector to the plane along the current (vec perpendicular) 
vec_p = [vec_a(1) vec_a(2) -(vec_a(1) + vec_a(2))/vec_a(3)]; 
mag_vp = sqrt(vec_p(1)*vec_p(1) + vec_p(2)*vec_p(2) + vec_p(3)*vec_p(3));
% This is the vector from the right hand rule
vec_rh = cross(vec_a, vec_p);
mag_vrh = sqrt(vec_rh(1)*vec_rh(1) + vec_rh(2)*vec_rh(2) + vec_rh(3)*vec_rh(3));

%--------------------------------------------------------------------------
%To get the slides using the normal vector along the current filament
vec=vec_a;
%To get the slides using the normal vector perpendiculr direction (z'')
%vec=vec_p;
%--------------------------------------------------------------------------

%The unitary vectors of the reference system 1
ui=[1 0 0];
uj=[0 1 0];
uk=[0 0 1];

% The unitary vectors of the new reference frame 2
uip = vec_rh/mag_vrh;
ujp = vec_a/mag_va;
ukp = vec_p/mag_vp;

% dummy proof
%vec=[0 0 1];
%uip = [-1 0 0];
%ujp = [0 -1 0];
%ukp = [0 0 1];

%Lets define the DCM
DCM = [dot(uip,ui) dot(uip,uj) dot(uip,uk);...
       dot(ujp,ui) dot(ujp,uj) dot(ujp,uk);...
       dot(ukp,ui) dot(ukp,uj) dot(ukp,uk)];

% Now lets compute the quantities in the new reference frame
%--------------------------------------------------------------------------
Gp = cell(nc,mc);
% Magnetic field in the prime reference frame
Gp{1,1}=G{1,1}.*DCM(1,1) + G{1,2}.*DCM(1,2) + G{1,3}.*DCM(1,3);
Gp{1,2}=G{1,1}.*DCM(2,1) + G{1,2}.*DCM(2,2) + G{1,3}.*DCM(2,3);
Gp{1,3}=G{1,1}.*DCM(3,1) + G{1,2}.*DCM(3,2) + G{1,3}.*DCM(3,3);
% Electron velocity
Gp{2,1}=G{2,1}.*DCM(1,1) + G{2,2}.*DCM(1,2) + G{2,3}.*DCM(1,3);
Gp{2,2}=G{2,1}.*DCM(2,1) + G{2,2}.*DCM(2,2) + G{2,3}.*DCM(2,3);
Gp{2,3}=G{2,1}.*DCM(3,1) + G{2,2}.*DCM(3,2) + G{2,3}.*DCM(3,3);
% Ion velocity
Gp{3,1}=G{3,1}.*DCM(1,1) + G{3,2}.*DCM(1,2) + G{3,3}.*DCM(1,3);
Gp{3,2}=G{3,1}.*DCM(2,1) + G{3,2}.*DCM(2,2) + G{3,3}.*DCM(2,3);
Gp{3,3}=G{3,1}.*DCM(3,1) + G{3,2}.*DCM(3,2) + G{3,3}.*DCM(3,3);
% Current
Gp{4,1}=G{4,1}.*DCM(1,1) + G{4,2}.*DCM(1,2) + G{4,3}.*DCM(1,3);
Gp{4,2}=G{4,1}.*DCM(2,1) + G{4,2}.*DCM(2,2) + G{4,3}.*DCM(2,3);
Gp{4,3}=G{4,1}.*DCM(3,1) + G{4,2}.*DCM(3,2) + G{4,3}.*DCM(3,3);
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



%--------------------------------------------------------------------------


% calculate the quantities on the slides
%-----------------------------------------------------------------------
slice_G = cell(nc,mc);
for i=1:nc
for j=1:mc
[slice, sliceInd,subX,subY,subZ] = extractSlice(G{i,j},pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);

% This removes the nans keeping the correct values
X = slice;
X(isnan(X(:,1)),:) = []; X=X';
X(isnan(X(:,1)),:) = []; X=X';

slice_G{i,j} = X;
end
end
%-----------------------------------------------------------------------

%--------------------------------------------------------------------------
slice_Gp = cell(nc,mc);
for i=1:nc
for j=1:mc
[slice, sliceInd,subX,subY,subZ] = extractSlice(Gp{i,j},pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);

% This removes the nans keeping the correct values
X = slice;
X(isnan(X(:,1)),:) = []; X=X';
X(isnan(X(:,1)),:) = []; X=X';

slice_Gp{i,j} = X;
end
end
%--------------------------------------------------------------------------


% built the grid to do the right plots
size_X=size(X);
xll=linspace(1,size_X(2),size_X(2))*0.06;
yll=linspace(1,size_X(1),size_X(1))*0.06;
[XLL, YLL] = meshgrid(xll,yll);


% Plot the field over the slide and the vectors in the original frame. 
%--------------------------------------------------------------------------

% The point where the reconnection is happening is around 
xp_mr = xll(106);
yp_mr = yll(127);
%--------------------------------------------------------------------------
Xd = slice_G{1,3};
% This are just the components and not the projections for the magnetic
% field
ubx = slice_G{1,1};   vby = slice_G{1,2}; wbz = slice_G{1,3};
bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
%ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
%vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

% This is for the electron velocities
Vev_x = slice_G{2,1}; Vev_y = slice_G{2,2}; Vev_z = slice_G{2,3}; 
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
%Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
%Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

NN=80; xstart = max(xll)*rand(NN,1); ystart = max(yll)*rand(NN,1); %Do this just once
sv=7;

f38=figure(38);
hc = pcolor(XLL,YLL,Xd); %Xd
colormap(BWR);
lim_yp=max(max(max(abs(Xd))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none')
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
title('Original frame')
hold off
%--------------------------------------------------------------------------

%Plots on the slide
%--------------------------------------------------------------------------
f40=figure(40);
for i=1:12
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(4,3,i,'SV',0,'SH',0,'MR',0,'ML',0.04,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=slice_G{3,i};
    titl = 'vi';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=slice_G{2,i2};
    titl = 've';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=slice_G{4,i3};
    titl = 'J';
    s=i3;
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=slice_G{1,i4};
    if (i4==3)
        dum_paux=slice_G{1,i4};
        dum_p = dum_paux;% - 0.1;
    end
    titl = 'B';
    s=i4;    
end
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none')
hcb=colorbar;
colormap(BWR);
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = 'p / d_{i}';
    ax.YTick = [2 4 6 8 10];
end
if (9<i)
    ax.XLabel.String = 'rh / d_{i}';
    ax.XTick = [2 4 6 8];
    
end
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
%xlim([xl(1) xl(end)]);
%ylim([yl(1) yl(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1} \ Slide RF1$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------



% Plot the field over the slide and the vectors in the NEW frame. 
%--------------------------------------------------------------------------
Xd = slice_Gp{1,3};

% This are the projections for the magnetic field
%ubx = slice_Gp{1,1};   vby = slice_Gp{1,2}; wbz = slice_Gp{1,3};
ubx = slice_Gp{1,3};   vby = slice_Gp{1,1}; wbz = slice_Gp{1,2};
bm=sqrt(ubx.*ubx + vby.*vby + wbz.*wbz);
ubx=ubx.*sqrt(1-((wbz./bm).*(wbz./bm)));
vby=vby.*sqrt(1-((wbz./bm).*(wbz./bm)));

% This is for the electron velocities
%Vev_x = slice_Gp{2,1}; Vev_y = slice_Gp{2,2}; Vev_z = slice_Gp{2,3}; 
Vev_x = slice_Gp{2,3}; Vev_y = slice_Gp{2,1}; Vev_z = slice_Gp{2,2}; 
Vevm=sqrt(Vev_x.*Vev_x + Vev_y.*Vev_y + Vev_z.*Vev_z);
Vev_x=Vev_x.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));
Vev_y=Vev_y.*sqrt(1-((Vev_z./Vevm).*(Vev_z./Vevm)));

NN=80; xstart = max(xll)*rand(NN,1); ystart = max(yll)*rand(NN,1); %Do this just once
sv=7;

f39=figure(39);
hc = pcolor(XLL,YLL,Xd);
colormap(BWR);
lim_yp=max(max(max(abs(Xd))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none')
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
title('New FRAME')
hold off
%--------------------------------------------------------------------------


%Plots on the slide. The velocities, current and magnetic field
%--------------------------------------------------------------------------
f41=figure(41);
for i=1:12
%h1=subaxis(4,3,i,'SV',0.04,'SH',0.04,'MR',0.02,'ML',0.05,'PL',0.05,'PR',0.05);
h1=subaxis(4,3,i,'SV',0,'SH',0,'MR',0,'ML',0.04,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=slice_Gp{3,i};
    titl = 'vi';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=slice_Gp{2,i2};
    titl = 've';
    s=i2;
elseif (6<i && i<=9)
    i3=i-6;
    dum_p=slice_Gp{4,i3};
    titl = 'J';
    s=i3;
elseif (9<i && i<=12)
    i4=i-9;
    dum_p=slice_Gp{1,i4};
    if (i4==3)
        dum_paux=slice_Gp{1,i4};
        dum_p = dum_paux;% - 0.1;
    end
    titl = 'B';
    s=i4;    
end
hc = pcolor(XLL,YLL,dum_p);
lim_yp=max(max(max(abs(dum_p))));
caxis([-lim_yp lim_yp])
set(hc,'edgecolor','none')
hcb=colorbar;
colormap(BWR);
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = 'p / d_{i}';
    ax.YTick = [2 4 6 8 10];
end
if (9<i)
    ax.XLabel.String = 'rh / d_{i}';
    ax.XTick = [2 4 6 8];
    
end
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
%set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
%xlim([xl(1) xl(end)]);
%ylim([yl(1) yl(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1} \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f40,strcat('vi_ve_J_B_vev_Bl_slideRF1_per'+ string(N_steps) +'.png'));
%saveas(f41,strcat('vi_ve_J_B_vev_Bl_slideRF2_per'+ string(N_steps) +'.png'));
saveas(f41,strcat('vi_ve_J_B_vev_no_Bl_slideRF2_per'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';





% Now the transformation and the plots for the tensors
%--------------------------------------------------------------------------
%Now the proper transformation of the tensor.
%--------------------------------------------------------------------------
%For ions
GPij=cell(3,3);
%--------------------------------------------------------------------------
% Choose the tensor
GPij{1,1}=Pij_i{1,1}; GPij{1,2}=Pij_i{1,2}; GPij{1,3}=Pij_i{1,3};
GPij{2,1}=Pij_i{1,2}; GPij{2,2}=Pij_i{2,2}; GPij{2,3}=Pij_i{2,3};
GPij{3,1}=Pij_i{1,3}; GPij{3,2}=Pij_i{2,3}; GPij{3,3}=Pij_i{3,3};
%--------------------------------------------------------------------------
% % Transform the tensor. Make this a function
%--------------------------------------------------------------------------
GPijp = cell(3,3);
% DCM*Pij
GPijp{1,1}=GPij{1,1}.*DCM(1,1) + GPij{2,1}.*DCM(1,2) + GPij{3,1}.*DCM(1,3);
GPijp{1,2}=GPij{1,2}.*DCM(1,1) + GPij{2,2}.*DCM(1,2) + GPij{3,2}.*DCM(1,3);
GPijp{1,3}=GPij{1,3}.*DCM(1,1) + GPij{2,3}.*DCM(1,2) + GPij{3,3}.*DCM(1,3);
GPijp{2,1}=GPij{1,1}.*DCM(2,1) + GPij{2,1}.*DCM(2,2) + GPij{3,1}.*DCM(2,3);
GPijp{2,2}=GPij{1,2}.*DCM(2,1) + GPij{2,2}.*DCM(2,2) + GPij{3,2}.*DCM(2,3);
GPijp{2,3}=GPij{1,3}.*DCM(2,1) + GPij{2,3}.*DCM(2,2) + GPij{3,3}.*DCM(2,3);
GPijp{3,1}=GPij{1,1}.*DCM(3,1) + GPij{2,1}.*DCM(3,2) + GPij{3,1}.*DCM(3,3);
GPijp{3,2}=GPij{1,2}.*DCM(3,1) + GPij{2,2}.*DCM(3,2) + GPij{3,2}.*DCM(3,3);
GPijp{3,3}=GPij{1,3}.*DCM(3,1) + GPij{2,3}.*DCM(3,2) + GPij{3,3}.*DCM(3,3);
% DCM*Pij*DCM'
GPijp2{1,1}=GPijp{1,1}.*DCM(1,1) + GPijp{1,2}.*DCM(1,2) + GPijp{1,3}.*DCM(1,3);
GPijp2{1,2}=GPijp{1,1}.*DCM(2,1) + GPijp{1,2}.*DCM(2,2) + GPijp{1,3}.*DCM(2,3);
GPijp2{1,3}=GPijp{1,1}.*DCM(3,1) + GPijp{1,2}.*DCM(3,2) + GPijp{1,3}.*DCM(3,3);
GPijp2{2,1}=GPijp{2,1}.*DCM(1,1) + GPijp{2,2}.*DCM(1,2) + GPijp{2,3}.*DCM(1,3);
GPijp2{2,2}=GPijp{2,1}.*DCM(2,1) + GPijp{2,2}.*DCM(2,2) + GPijp{2,3}.*DCM(2,3);
GPijp2{2,3}=GPijp{2,1}.*DCM(3,1) + GPijp{2,2}.*DCM(3,2) + GPijp{2,3}.*DCM(3,3);
GPijp2{3,1}=GPijp{3,1}.*DCM(1,1) + GPijp{3,2}.*DCM(1,2) + GPijp{3,3}.*DCM(1,3);
GPijp2{3,2}=GPijp{3,1}.*DCM(2,1) + GPijp{3,2}.*DCM(2,2) + GPijp{3,3}.*DCM(2,3);
GPijp2{3,3}=GPijp{3,1}.*DCM(3,1) + GPijp{3,2}.*DCM(3,2) + GPijp{3,3}.*DCM(3,3);
GPij_i=GPijp2;
clearvars GPijp2;
%--------------------------------------------------------------------------

%For electrons
GPij=cell(3,3);
%--------------------------------------------------------------------------
% Choose the tensor
GPij{1,1}=Pij_e{1,1}; GPij{1,2}=Pij_e{1,2}; GPij{1,3}=Pij_e{1,3};
GPij{2,1}=Pij_e{1,2}; GPij{2,2}=Pij_e{2,2}; GPij{2,3}=Pij_e{2,3};
GPij{3,1}=Pij_e{1,3}; GPij{3,2}=Pij_e{2,3}; GPij{3,3}=Pij_e{3,3};
%--------------------------------------------------------------------------
% % Transform the tensor. Make this a function
%--------------------------------------------------------------------------
GPijp = cell(3,3);
% DCM*Pij
GPijp{1,1}=GPij{1,1}.*DCM(1,1) + GPij{2,1}.*DCM(1,2) + GPij{3,1}.*DCM(1,3);
GPijp{1,2}=GPij{1,2}.*DCM(1,1) + GPij{2,2}.*DCM(1,2) + GPij{3,2}.*DCM(1,3);
GPijp{1,3}=GPij{1,3}.*DCM(1,1) + GPij{2,3}.*DCM(1,2) + GPij{3,3}.*DCM(1,3);
GPijp{2,1}=GPij{1,1}.*DCM(2,1) + GPij{2,1}.*DCM(2,2) + GPij{3,1}.*DCM(2,3);
GPijp{2,2}=GPij{1,2}.*DCM(2,1) + GPij{2,2}.*DCM(2,2) + GPij{3,2}.*DCM(2,3);
GPijp{2,3}=GPij{1,3}.*DCM(2,1) + GPij{2,3}.*DCM(2,2) + GPij{3,3}.*DCM(2,3);
GPijp{3,1}=GPij{1,1}.*DCM(3,1) + GPij{2,1}.*DCM(3,2) + GPij{3,1}.*DCM(3,3);
GPijp{3,2}=GPij{1,2}.*DCM(3,1) + GPij{2,2}.*DCM(3,2) + GPij{3,2}.*DCM(3,3);
GPijp{3,3}=GPij{1,3}.*DCM(3,1) + GPij{2,3}.*DCM(3,2) + GPij{3,3}.*DCM(3,3);
% DCM*Pij*DCM'
GPijp2{1,1}=GPijp{1,1}.*DCM(1,1) + GPijp{1,2}.*DCM(1,2) + GPijp{1,3}.*DCM(1,3);
GPijp2{1,2}=GPijp{1,1}.*DCM(2,1) + GPijp{1,2}.*DCM(2,2) + GPijp{1,3}.*DCM(2,3);
GPijp2{1,3}=GPijp{1,1}.*DCM(3,1) + GPijp{1,2}.*DCM(3,2) + GPijp{1,3}.*DCM(3,3);
GPijp2{2,1}=GPijp{2,1}.*DCM(1,1) + GPijp{2,2}.*DCM(1,2) + GPijp{2,3}.*DCM(1,3);
GPijp2{2,2}=GPijp{2,1}.*DCM(2,1) + GPijp{2,2}.*DCM(2,2) + GPijp{2,3}.*DCM(2,3);
GPijp2{2,3}=GPijp{2,1}.*DCM(3,1) + GPijp{2,2}.*DCM(3,2) + GPijp{2,3}.*DCM(3,3);
GPijp2{3,1}=GPijp{3,1}.*DCM(1,1) + GPijp{3,2}.*DCM(1,2) + GPijp{3,3}.*DCM(1,3);
GPijp2{3,2}=GPijp{3,1}.*DCM(2,1) + GPijp{3,2}.*DCM(2,2) + GPijp{3,3}.*DCM(2,3);
GPijp2{3,3}=GPijp{3,1}.*DCM(3,1) + GPijp{3,2}.*DCM(3,2) + GPijp{3,3}.*DCM(3,3);
GPij_e=GPijp2;
clearvars GPijp2;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% This is just for the ion pressure tensor plot in the new reference frame
%--------------------------------------------------------------------------
nc=3; mc=3;
GPijf=cell(nc,mc);
GPijp=GPij_i;
%--------------------------------------------------------------------------
% Pij_i
GPijf{1,1}=GPijp{1,2}; GPijf{1,2}=GPijp{1,3}; GPijf{1,3}=GPijp{2,3};
% Pii_i
GPijf{2,1}=GPijp{1,1}; GPijf{2,2}=GPijp{2,2}; GPijf{2,3}=GPijp{3,3};
% Pji_i
GPijf{3,1}=GPijp{2,1}; GPijf{3,2}=GPijp{3,1}; GPijf{3,3}=GPijp{3,2};
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
slice_GPijp = cell(nc,mc);
for i=1:nc
for j=1:mc
[slice, sliceInd,subX,subY,subZ] = extractSlice(GPijf{i,j},pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);

% This removes the nans keeping the correct values
X = slice;
X(isnan(X(:,1)),:) = []; X=X';
X(isnan(X(:,1)),:) = []; X=X';
slice_GPijp{i,j} = X;
end
end
%--------------------------------------------------------------------------

% Pij for ion in the slide F2
%--------------------------------------------------------------------------
f42=figure(42);
for i=1:6
h1=subaxis(2,3,i,'SV',0,'SH',0,'MR',0,'ML',0.04,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=slice_GPijp{1,i};
    titl = 'iPij';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=slice_GPijp{2,i2};
    titl = 'iPii';
    s=i2;
end
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
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = 'p / d_{i}';
    ax.YTick = [2 4 6 8 10];
end
if (9<i)
    ax.XLabel.String = 'rh / d_{i}';
    ax.XTick = [2 4 6 8];
end
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
%set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
%xlim([xl(1) xl(end)]);
%ylim([yl(1) yl(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1} \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f42,strcat('Pij_i_slideRF2_per'+ string(N_steps) +'.png'));
saveas(f42,strcat('Pij_i_slideRF2_per_no_Bl'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% This is just for the electron pressure tensor plot in the new reference frame
%--------------------------------------------------------------------------
nc=3; mc=3;
GPijf=cell(nc,mc);
GPijp=GPij_e;
%--------------------------------------------------------------------------
% Pij_i
GPijf{1,1}=GPijp{1,2}; GPijf{1,2}=GPijp{1,3}; GPijf{1,3}=GPijp{2,3};
% Pii_i
GPijf{2,1}=GPijp{1,1}; GPijf{2,2}=GPijp{2,2}; GPijf{2,3}=GPijp{3,3};
% Pji_i
GPijf{3,1}=GPijp{2,1}; GPijf{3,2}=GPijp{3,1}; GPijf{3,3}=GPijp{3,2};
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
slice_GPijp = cell(nc,mc);
for i=1:nc
for j=1:mc
[slice, sliceInd,subX,subY,subZ] = extractSlice(GPijf{i,j},pt(1),pt(2),pt(3),vec(1),vec(2),vec(3),radi);

% This removes the nans keeping the correct values
X = slice;
X(isnan(X(:,1)),:) = []; X=X';
X(isnan(X(:,1)),:) = []; X=X';
slice_GPijp{i,j} = X;
end
end
%--------------------------------------------------------------------------

% Pij for the electrons in the slide F2
%--------------------------------------------------------------------------
f43=figure(43);
for i=1:6
h1=subaxis(2,3,i,'SV',0,'SH',0,'MR',0,'ML',0.04,'PL',0.005,'PR',0.006);
if (i<=3)
    i1=i;
    dum_p=slice_GPijp{1,i};
    titl = 'ePij';
    s=i1;
elseif (3<i && i<=6)
    i2=i-3;
    dum_p=slice_GPijp{2,i2};
    titl = 'ePii';
    s=i2;
end
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
hcb.Label.FontSize = 12;
hcb.Label.Interpreter = 'latex';
hcb.Title.String='';
%set(get(hcb,'Title'),'String','','FontSize',10)
set(hcb.XLabel,{'String','Rotation','Position'},{titl + string(s),0,[0.5 0]})
%pbaspect([1 1 1])
ax = gca;
ax.XTick = [];
ax.YTick = [];
if (mod(i,3)==1)
    ax.YLabel.String = 'p / d_{i}';
    ax.YTick = [2 4 6 8 10];
end
if (9<i)
    ax.XLabel.String = 'rh / d_{i}';
    ax.XTick = [2 4 6 8];
end
hold on
plot(xp_mr,yp_mr,'mo','MarkerSize',10)
h2=quiver(xll(1:sv:end,1:sv:end), yll(1:sv:end,1:sv:end), Vev_x(1:sv:end,1:sv:end), Vev_y(1:sv:end,1:sv:end), 0);
%hlines=streamline(XLL,YLL,ubx,vby,xstart,ystart);
%set(hlines,'LineWidth',1.5, 'Color', 'm');
%set(hlines,'LineWidth',1, 'Color', [0.5 0.5 0.5]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 3, 'Color', 'k');
%xlim([xl(1) xl(end)]);
%ylim([yl(1) yl(end)]);
hold off
%set(gca,'XScale','lin','YScale','lin','FontSize',10)
end
%sgtitle('Subplot Grid Title')
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1} \ Slide RF2$$','Interpreter','latex','FontSize',20)
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
%saveas(f40,strcat('vi_ve_J_B_vev_Bl_slideRF1_per'+ string(N_steps) +'.png'));
%saveas(f41,strcat('vi_ve_J_B_vev_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%saveas(f41,strcat('vi_ve_J_B_vev_no_Bl_slideRF2_per'+ string(N_steps) +'.png'));
%saveas(f42,strcat('Pij_i_slideRF2_per'+ string(N_steps) +'.png'));
saveas(f43,strcat('Pij_e_slideRF2_no_bl_per'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------








% for the magnetic field lines in 3D
%--------------------------------------------------------------------------
u=B{1,1};
v=B{1,2};
w=B{1,3};
[xm3,ym3, zm3] = meshgrid(xl,yl,zl);
%[sx,sy,sz] = meshgrid(10:2:20,10:2:24,72:2:82);
dx=xl(2)-xl(1); dy=yl(2)-yl(1); dz=zl(2)-zl(1); 
[sx,sy,sz] = meshgrid(xl(1):8*dx:xl(end),yl(1):8*dy:yl(end),zl(83));
%[sx,sy,sz] = meshgrid(xl,yl,72);
XYZ2=stream3(xm3,ym3,zm3,u,v,w,sx,sy,sz);
dum_j=J{1,3};


f67=figure(67);
hc=pcolor(xl,yl,dum_j(:,:,slide_n));
hc.ZData = zl(slide_n) + hc.ZData; 
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



% 2d histograms 
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