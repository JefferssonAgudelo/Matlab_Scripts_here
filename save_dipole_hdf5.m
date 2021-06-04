%--------------------------------------------------------------------------
% This is to create the simple theoretical case in the .h5 format
%--------------------------------------------------------------------------
Lx=2; Ly=2; Lz=8;
nx = 100; ny = 100; nz = 400;

xl=linspace(-Lx,Lx,nx);
yl=linspace(-Ly,Ly,ny);
zl=linspace(-Lz,Lz,nz);
%rl=sqrt(xl.^2 + yl.^2 + zl.^2); This doesn't work

[XL, YL, ZL]=meshgrid(xl,yl,zl);
RL=sqrt(XL.^2 + YL.^2 + ZL.^2);
M=1*ones(size(XL));

% Dipole magnetic field cartesian
Bx=3.*M.*XL.*ZL./(RL.^5);
By=3.*M.*YL.*ZL./(RL.^5);
Bz=M.*(3.* (ZL.^2) - RL.^2)./(RL.^5);

divB = divergence(XL,YL,ZL,Bx,By,Bz); 

% For the streamlines
fr=1;
NN=200;  %Do this just once

ax=xl(1)+3.5; bx=xl(end)-3.5;
ay=yl(1)+3.5; by=yl(end)-3.5;
az=zl(1)+7.5; bz=zl(end)-7.5;

xstart = ax + (bx-ax)*rand(NN,1);
ystart = ay + (by-ay)*rand(NN,1);
zstart = az + (bz-az)*rand(NN,1);


%dx=xl(2)-xl(1); dy=yl(2)-yl(1); dz=zl(2)-zl(1); 
%[sx,sy,sz] = meshgrid(xl(1):8*dx:xl(end),yl(1):8*dy:yl(end),zl(83));
%XYZ2=stream3(XL,YL,ZL,Bx,By,Bz,sx,sy,sz);
XYZ2=stream3(XL,YL,ZL,Bx,By,Bz,xstart,ystart,zstart);

f341=figure(341);
histogram(divB)


f342=figure(342);
hline=streamline(XYZ2);
view(3);

f343=figure(343);
scatter3(xstart,ystart,zstart); 
view(3);

% This works
%--------------------------------------------------------------------------
% Create the data file
if isfile('dipole.h5'), delete geometry.h5; end
h5create('dipole.h5','/g1/Bx',size(Bx))
h5create('dipole.h5','/g1/By',size(By))
h5create('dipole.h5','/g1/Bz',size(Bz))

% write the data file
h5write('dipole.h5', '/g1/Bx', Bx)
h5write('dipole.h5', '/g1/By', By)
h5write('dipole.h5', '/g1/Bz', Bz)

h5disp('dipole.h5')
Bx2=h5read('dipole.h5','/g1/Bx');
%--------------------------------------------------------------------------







%--------------------------------------------------------------------------
% Define dimension of the 3D-array
Lx = 1.0;  nx = 121; 
Ly = 1.0;  ny = 111;
Lz = 1.0;  nz = 101;
% Build axis of the mesh
dx = Lx/(nx-1);  x = 0:dx:Lx;
dy = Ly/(ny-1);  y = 0:dy:Ly;
dz = Lz/(nz-1);  z = 0:dz:Lz;
% Build mesh rectangular grid
[X,Y,Z] = meshgrid(y,x,z);
% Initial condition of pressure field
pressure = 1.0*exp(-( (X-0.5).^2+(Y-0.5).^2+(Z-0.5).^2 ) ./ (0.31)^2 );
% Create geometry file
if isfile('geometry.h5'), delete geometry.h5; end
% x-axis
h5create('geometry.h5','/x_nodes',nx);
h5write('geometry.h5','/x_nodes',x);
h5disp('geometry.h5');
% y-axis
h5create('geometry.h5','/y_nodes',ny);
h5write('geometry.h5','/y_nodes',y);
h5disp('geometry.h5');
% z-axis
h5create('geometry.h5','/z_nodes',nz);
h5write('geometry.h5','/z_nodes',z);
h5disp('geometry.h5');
% Create data file
if isfile('fields_it0.h5'), delete fields_it0.h5; end
% Pressure field
h5create('fields_it0.h5','/p',[nx,ny,nz]);
h5write('fields_it0.h5','/p',pressure);
h5disp('fields_it0.h5');
% Numbers precision
precision = 8; % bit for real double 
time = 0.0; % assuming t=0
it = 0; % assuming this is the initial condition
% write associated XDMF file
fileID = fopen('openWithParaView.xmf','w');
fprintf(fileID,'<?xml version="1.0" ?>\n');
fprintf(fileID,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fileID,'<Xdmf Version="2.0">\n');
fprintf(fileID,' <Domain>\n');
fprintf(fileID,'  <Grid Name="mesh" GridType="Uniform">\n');
fprintf(fileID,'   <Topology TopologyType="3DRectMesh" NumberOfElements="%d %d %d"/>\n',nz,ny,nx);
fprintf(fileID,'   <Geometry GeometryType="VXVYVZ">\n');
fprintf(fileID,'    <DataItem Name="coordx" Dimensions="%g" NumberType="Float" Precision="%g" Format="HDF">\n',nx,precision);
fprintf(fileID,'     geometry.h5:/x_nodes\n');
fprintf(fileID,'    </DataItem>\n');
fprintf(fileID,'    <DataItem Name="coordy" Dimensions="%g" NumberType="Float" Precision="%g" Format="HDF">\n',ny,precision);
fprintf(fileID,'     geometry.h5:/y_nodes\n');
fprintf(fileID,'    </DataItem>\n');
fprintf(fileID,'    <DataItem Name="coordz" Dimensions="%g" NumberType="Float" Precision="%g" Format="HDF">\n',nz,precision);
fprintf(fileID,'     geometry.h5:/z_nodes\n');
fprintf(fileID,'    </DataItem>\n');
fprintf(fileID,'   </Geometry>\n');
fprintf(fileID,'   <Time TimeType="Single" Value="%g"/>\n',time);
fprintf(fileID,'   <Attribute Name="p" AttributeType="Scalar" Center="Node">\n');
fprintf(fileID,'    <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="%d" Format="HDF">\n',nz,ny,nx,precision);
fprintf(fileID,'     fields_it%g.h5:/p\n',it);
fprintf(fileID,'    </DataItem>\n');
fprintf(fileID,'   </Attribute>\n');
fprintf(fileID,'  </Grid>\n');
fprintf(fileID,' </Domain>\n');
fprintf(fileID,'</Xdmf>\n');
fclose(fileID);

%-------------------------------------------------------------------------

mydata = rand(100,100,400);


f88=figure(88)
pcolor(Bz(:,:,4))


pt=[50,50,100];
rl=sqrt((xl-pt(1)).*(xl-pt(1)) + (yl-pt(2)).*(yl-pt(2)) + (zl-pt(3)).*(zl-pt(3)));
Rhat=[XL./RL, YL./RL, ZL./RL];
ui=[1,0,0]; uj=[0,1,0]; uk=[0,0,1];
