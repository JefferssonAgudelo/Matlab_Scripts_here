%--------------------------------------------------------------------------
%This script is to compute the energy budget
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cd '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data'; %This is important because the xdmf files are in that directory
path = '/Volumes/PSC_DiRAC_DATA/DATACB104/raw_data';

%--------------------------------------------------------------------------
S = dir(fullfile(path,'*.h5'));
N=numel(S); % number of files to use
H=zeros(N);
i=19;

%for i=2:N
%--------------------------------------------------------------------------
disp(strcat('Computing step ...',S(i).name))
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

mi=1;
mime=100;
me=mi/mime;
qi=1;
qe=-1;

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
%vivi
vivi{1,1} = vi{1,1}.*vi{1,1}; vivi{1,2} = vi{1,1}.*vi{1,2}; vivi{1,3} = vi{1,1}.*vi{1,3};  
vivi{2,1} = vi{1,2}.*vi{1,1}; vivi{2,2} = vi{1,2}.*vi{1,2}; vivi{2,3} = vi{1,2}.*vi{1,3};
vivi{3,1} = vi{1,3}.*vi{1,1}; vivi{3,2} = vi{1,3}.*vi{1,2}; vivi{3,3} = vi{1,3}.*vi{1,3};

%veve
veve{1,1} = ve{1,1}.*ve{1,1}; veve{1,2} = ve{1,1}.*ve{1,2}; veve{1,3} = ve{1,1}.*ve{1,3};  
veve{2,1} = ve{1,2}.*ve{1,1}; veve{2,2} = ve{1,2}.*ve{1,2}; veve{2,3} = ve{1,2}.*ve{1,3};
veve{3,1} = ve{1,3}.*ve{1,1}; veve{3,2} = ve{1,3}.*ve{1,2}; veve{3,3} = ve{1,3}.*ve{1,3};

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


% To calculate the grad of the heat flux vector wang et al, 2015
%------------------------------------------------------------------
p_i = (Pij_i{1,1} + Pij_i{2,2} + Pij_i{3,3})./3;  
p_e = (Pij_e{1,1} + Pij_e{2,2} + Pij_e{3,3})./3;  
%------------------------------------------------------------------


% From the definition in the code I think Tij is actually Pij as
% int( fvivj )dv^3 and not PPij=int f(vi - ui)(vj - uj) dv^3 which is
% PPij = Pij - nmuiuj
%------------------------------------------------------------------
PPij_i{1,1} = Pij_i{1,1} - ni.*mi.*vi{1,1}.*vi{1,1}; PPij_i{1,2} = Pij_i{1,2} - ni.*mi.*vi{1,1}.*vi{1,2}; PPij_i{1,3} = Pij_i{1,3} - ni.*mi.*vi{1,1}.*vi{1,3}; 
PPij_i{2,2} = Pij_i{2,2} - ni.*mi.*vi{1,2}.*vi{1,2}; PPij_i{2,3} = Pij_i{2,3} - ni.*mi.*vi{1,2}.*vi{1,3}; 
PPij_i{3,3} = Pij_i{3,3} - ni.*mi.*vi{1,3}.*vi{1,3}; 

PPij_e{1,1} = Pij_e{1,1} - ne.*me.*ve{1,1}.*ve{1,1}; PPij_e{1,2} = Pij_e{1,2} - ne.*me.*ve{1,1}.*ve{1,2}; PPij_e{1,3} = Pij_e{1,3} - ne.*me.*ve{1,1}.*ve{1,3}; 
PPij_e{2,2} = Pij_e{2,2} - ne.*me.*ve{1,2}.*ve{1,2}; PPij_e{2,3} = Pij_e{2,3} - ne.*me.*ve{1,2}.*ve{1,3}; 
PPij_e{3,3} = Pij_e{3,3} - ne.*me.*ve{1,3}.*ve{1,3}; 


PPij_i=cell(3,3);
PPij_e=cell(3,3);
for i=1:3
    for j=1:3
        a=isempty(Pij_i{i,j});
        if(a ==1)
            continue
        end
        PPij_i{i,j} = Pij_i{i,j} - ni.*mi.*vi{1,i}.*vi{1,j};
        PPij_e{i,j} = Pij_e{i,j} - ne.*me.*ve{1,i}.*ve{1,j};
    end
end


%------------------------------------------------------------------

% Calculating the temperatures
%------------------------------------------------------------------
TTij_i=cell(3,3);
TTij_e=cell(3,3);
for i=1:3
    for j=1:3
        if(isempty(PPij_i{i,j}) ==1)
            continue
        end
        TTij_i{i,j} =PPij_i{i,j}./(3.*ni);
        TTij_e{i,j} =PPij_e{i,j}./(3.*ne);
    end
end



%------------------------------------------------------------------

%HEre I am changing stuff...

%Lets compute the averages to see what is what
%------------------------------------------------------------------
PPij_i_bar=zeros(3);
for i=1:3
    for j=1:3
        PPij_i_bar(i,j)=mean(PPij_i{i,j}, 'all');
    end 
end
PPij_e_bar=zeros(3);
for i=1:3
    for j=1:3
        PPij_e_bar(i,j)=mean(PPij_e{i,j}, 'all');
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
for i=1:3
    for j=1:3
        Pij_i_bar(i,j)=mean(Pij_i{i,j}, 'all');
    end 
end
Pij_e_bar=zeros(3);
for i=1:3
    for j=1:3
        Pij_e_bar(i,j)=mean(Pij_e{i,j}, 'all');
    end 
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


