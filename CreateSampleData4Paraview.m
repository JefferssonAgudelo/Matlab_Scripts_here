%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Visualize a Matlab 3D-arrays in Paraview/VisIt using HDF5 & XDMF
%
%               Coded by Manuel Diaz, ENSMA, 2020.07.28.
%                   Copyright (c) 2020, Manuel Diaz.
%                           All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks: 
%    1. This snippet exemplifies how to export a 3D array into ParaView as
%    a 3DRectMesh. This is a classical approach to export large volumes of
%    data preferred in High Perfomance Computations (HPC).  
%             ( like reading/writting [10,000]^3 arrays )
%    2. This example outputs 3 files: the geometry.h5, the field.h5 and
%    the solution.xmf. The 3 files are required to be in the same folder.
%    But for opening the array in ParaView/Visit, only the *.xmf file is
%    initially required. ParaView/VisIt will then load the two h5-files.
%    3. When exporting data from Matlab/Fortran observe that the data is
%    transposed. This is because Matlab (like Fortran) store data in memory
%    using a column-major order, but HDF5 and ParaView are C-based
%    applications that store data using a row-major order. More details can
%    be found at:  
%        https://support.hdfgroup.org/HDF5/doc1.6/UG/12_Dataspaces.html
%    4. XDMF stands for eXtensible Data Model and Format. Is a standardized
%    approach to exchange scientific data between HPC codes. For further
%    details, I refer to: 
%        http://www.xdmf.org/
%    5. Matlab supports natively HDF5 format, actually *.mat files use the
%    hdf5 technology on the background. So using h5-files are highly
%    encouraged to read/write/share data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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