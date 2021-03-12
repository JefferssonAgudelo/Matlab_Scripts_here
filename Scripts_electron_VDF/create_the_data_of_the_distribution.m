% This is the information with the other file, the one that doesn't have
% the data below 30 ev

load('JeffA.mat')
f44=figure(44);
pcolor(VparL315212,VperL315212,log10(DistL315212))
colorbar

load('JeffA151220190401.txt')
vpar1d=JeffA151220190401(:,2);
vper1d=JeffA151220190401(:,1);
hper1d=JeffA151220190401(:,3);
% To go from the 1D to 2D knowing the number of pitch angles 
vpar2D=reshape(vpar1d,19,12);
vper2D=reshape(vper1d,19,12);
hper2D=reshape(hper1d,19,12);

f55=figure(55);
pcolor(vpar2D,vper2D,log10(hper2D))
colorbar

max(hper1d)
%This in hear generates the data point that correspond to the distribution
d_funct_3=hper1d*1e14;
Drealtotal2=[];

for j=1:length(d_funct_3)
    dtotal_j2=zeros(fix(d_funct_3(j,1)),1);
    for i=1:fix(d_funct_3(j,1))
        dtotal_j2(i,2)=vper2(j,1); 
        dtotal_j2(i,1)=vpar2(j,1); %par per
        %p_per=vper2(j,1);%q_par=vpar2(j,1);
        %save('per_par.mat','p_per','q_par','-append')
    end
    dtotal_j2=dtotal_j;
    Drealtotal2=vertcat(Drealtotal2, dtotal_j2);    
end

% Now pick the random Numbers
n_points=length(Drealtotal2);
%x = randn(n_points,1); % This returns a uniform distribution
x=Drealtotal2;
arrayA2=zeros(n_points,2);
for i = 1:n_points
    xi=randi(length(x));
    arrayA2(i,1) = x(xi,1);
    arrayA2(i,2) = x(xi,2);
end

[hist_xhat2,Xedges2,Yedges2] = histcounts2(arrayA2(:,1),arrayA2(:,2),[19 12]);


[qxa2,qya2] = meshgrid(linspace(min(min(Xedges2)),max(max(Xedges2)),19),linspace(min(min(Yedges2)),max(max(Yedges2)),12));
f3=figure(3);
pcolor(qxa2,qya2,log10(hist_xhat2'))
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Users/jeffersson_agudelo/Downloads/Electron_project'
load('DV201904011512.txt')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vper=DV201904011512(:,1);
vpar=DV201904011512(:,2);
d_funct_2=DV201904011512(:,3);
f88=figure(88);
histogram(d_funct_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpar2D_D=reshape(vpar,6,71);
vper2D_D=reshape(vper,6,71);
d_funct_2D=reshape(d_funct_2,6,71);
f77=figure(77);
pcolor(vper2D_D,vpar2D_D,log10(d_funct_2D))
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Drealtotal3=Drealtotal;

% With the uninterpolated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This in hear generates the data point that correspond to the distribution
%d_funct_L=1e21*d_funct_2;
%dlog3=log(d_funct_L);
d_funct_3=1e17*d_funct_2;
Drealtotal=[];
%x='v_per'; y='v_par'; %save per_par.mat x y

% This generates the data tha represents the distribution
for j=1:length(d_funct_3)
    dtotal_j=zeros(fix(d_funct_3(j,1)),1);
    for i=1:fix(d_funct_3(j,1))
        dtotal_j(i,2)=vper(j,1); 
        dtotal_j(i,1)=vpar(j,1); %par per
        %p_per=vper2(j,1);%q_par=vpar2(j,1);
        %save('per_par.mat','p_per','q_par','-append')
    end
    dtotal_j=dtotal_j;
    Drealtotal=vertcat(Drealtotal, dtotal_j);    
end

Drealtotal3=Drealtotal;
% Now pick the pairs of random Numbers
n_points=length(1*Drealtotal3);
%x = randn(n_points,1); % This returns a uniform distribution
x=Drealtotal3; %3
arrayA=zeros(n_points,2);
for i = 1:n_points
    xi=randi(length(x));
    arrayA(i,1) = x(xi,1);
    arrayA(i,2) = x(xi,2);
end

f288=figure(288);
h = histogram2(arrayA(:,1),arrayA(:,2));



% This generates the 2D counts (This might not be working properly)
[hist_xhat,Xedges,Yedges] = histcounts2(arrayA(:,1),arrayA(:,2),50);
% Make some fake data
xb = double(arrayA(:,1));
yb = double(arrayA(:,2));
%zb = hist_xhat(:); This might not be quite right to do
% Put data onto a grid
[qxb,qyb] = meshgrid(linspace(min(xb),max(xb),50),linspace(min(yb),max(yb),50));
%Fb = TriScatteredInterp(xb,yb,zb); This doens't work for the 2D output of the histcounts2 
qzb = Fb(qxb,qyb);
% Make contour plot
f66=figure(66);
pcolor(qxb,qyb,log10(hist_xhat))
colorbar


% With the interpolated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xa=DV201904011512;
% Make some fake data
xa = double(Xa(:,1));
ya = double(Xa(:,2));
za = double(Xa(:,3));
% Put data onto a grid
[qxa,qya] = meshgrid(linspace(min(xa),max(xa),40),linspace(min(ya),max(ya),40));
Fa = TriScatteredInterp(xa,ya,za);
qza = Fa(qxa,qya);

% Make contour plot
f11=figure(11);
pcolor(qxa,qya,log(qza))
colorbar

d_funct_3I=1e16*qza;
Drealtotal=[];
%x='v_per'; y='v_par'; %save per_par.mat x y

% this is not doing anything yet. 
for j=1:length(d_funct_3I)
    dtotal_j=zeros(fix(d_funct_3I(j,1)),1);
    for i=1:fix(d_funct_3(j,1))
        dtotal_j(i,2)=vper(j,1); 
        dtotal_j(i,1)=vpar(j,1);
    end
    dtotal_j=dtotal_j;
    Drealtotal=vertcat(Drealtotal, dtotal_j);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f1=figure(1);
histogram(d_funct_2,100);
f2=figure(2);
histogram(d_funct_L,100);
f3=figure(3);
histogram(d_funct_3,100);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%using the external function
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';
%Generate_rand_numbers(pxe1, pxe2, num_bins, times_l)
[xhat1, xhat2] = Generate_rand_numbers(Drealtotal(:,1), Drealtotal(:,2), 100, 4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Users/jeffersson_agudelo/Downloads/Electron_project/myplots';
saveas(f1,'log_distributions_par_per.png');
%saveas(f2,'sub_electron_data.png');


  
  
  
% Put data onto a grid
[qvpar2,qvper2] = meshgrid(linspace(min(vpar2),max(vpar2),60),linspace(min(vper2),max(vper2),60));
Fa = TriScatteredInterp(vpar2,vper2,d_funct_2);
%Fa = scatteredInterpolant(vpar2,vper2,d_funct_2); %This works but pcolor complains becasue non real values in the log
qd_funct_2 = Fa(qvpar2,qvper2);
% Make contour plot
f55=figure(55);
pcolor(qvpar2,qvper2,log(qd_funct_2))
colorbar