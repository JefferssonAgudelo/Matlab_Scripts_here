function [slice_Gp, sliceInd,subX,subY,subZ] = Vars_on_Slice(Gp,pt,vec,radi)
% This function is to calculate the slides. It requeires the quantity either scalar, vector or tensor
% The inputs are the quantity: Gp, the point: pt, the normal vector: vec and
% the radii:radi

%--------------------------------------------------------------------------
% calculate the quantities on the slides for scalars
%--------------------------------------------------------------------------
if (size(Gp)==[1 1])
nc=1; mc=1;

%Gp=nic; i=1; j=1;

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
end
%--------------------------------------------------------------------------
% calculate the quantities on the slides for vectors
%--------------------------------------------------------------------------
if (size(Gp)==[1 3])
nc=1; mc=3;
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
end
%--------------------------------------------------------------------------
% calculate the quantities on the slides for tensors
%--------------------------------------------------------------------------
if (size(Gp)==[3 3])
nc=3; mc=3;
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
end
%--------------------------------------------------------------------------

end