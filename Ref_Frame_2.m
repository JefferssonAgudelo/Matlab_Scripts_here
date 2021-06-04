function out = Ref_Frame_2(uip,ujp,ukp,B)
% This is to calculate the reference frame transformation of the different
% elements. The output is just the quantity in the new reference frame.

% uip, ujp, ukp are the unitary vectors in for the new reference frame.
% B as a cell is either the vector to be transformed or the tensor 
% The unitary vectors of the new reference frame 2
% uip = vec_rh/mag_vrh; ujp = vec_a/mag_va; ukp = vec_p/mag_vp;

%--------------------------------------------------------------------------
%The unitary vectors of the reference system 1
ui=[1 0 0]; uj=[0 1 0]; uk=[0 0 1];

%Lets define the DCM
DCM = [dot(uip,ui) dot(uip,uj) dot(uip,uk);...
       dot(ujp,ui) dot(ujp,uj) dot(ujp,uk);...
       dot(ukp,ui) dot(ukp,uj) dot(ukp,uk)];
   
   % Now lets compute the quantities in the new reference frame
%--------------------------------------------------------------------------
% In case of imput as a vector
if (size(B)==[1 3])
G=B;
Gp = cell(1,3);
Gp{1,1}=G{1,1}.*DCM(1,1) + G{1,2}.*DCM(1,2) + G{1,3}.*DCM(1,3);
Gp{1,2}=G{1,1}.*DCM(2,1) + G{1,2}.*DCM(2,2) + G{1,3}.*DCM(2,3);
Gp{1,3}=G{1,1}.*DCM(3,1) + G{1,2}.*DCM(3,2) + G{1,3}.*DCM(3,3);
out = Gp;  
clearvars GP;

% In case of imput as a tensor
elseif((size(B)==[3 3]))    
%--------------------------------------------------------------------------
GPij=B;
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
out=GPijp2;
clearvars GPijp2;
end
%--------------------------------------------------------------------------
end


