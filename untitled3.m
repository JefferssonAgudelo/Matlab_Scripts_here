% Example 1
load mri
D = squeeze(D);
%DD=1./D;
vtkwrite('mri.vtk', 'structured_points', 'mri1', D)
%vtkwrite('mri.vtk', 'structured_points', 'mri2', DD)% This doen't work

load wind
[cu,cv,cw] = curl(x, y, z, u, v, w);
div = divergence(x, y, z, u, v, w);
vtkwrite('wind.vtk', 'structured_grid', x, y, z, ...
'vectors', 'vector_field', u, v, w, 'vectors', 'vorticity', cu, cv, cw, 'scalars', 'divergence', div);



