% To find structures in a 3d array

%
rng('default') % For reproducibility
r = normrnd(0,1);
%

% In 2D this works
f1=figure(1);
c=randn(5,5,5);
pcolor(c(:,:,1))
colorbar

v2=c(:,:,1) > 0.5;
v4=find(c>0.5);
v3=double(v2);
pcolor(v3)
colorbar

%To get the indexes and loose the matrix form
[r1,c1,v1] = ind2sub(size(c),find(c >= 0.5));

% in general more than 3D
[subI{1:ndims(c)}] = ind2sub(size(c), find(c>0.5));
B=squeeze(QQ(:,:,1));

%In 3D is as simple as this

c=randn(30,30,30);
f2=figure(2);
pcolor(c(:,:,10))
colorbar

x = c>=0; %create the matrix with the same dimention that fulfils the condition
A=c.*x; % get the values in that position
f3=figure(3);
%pcolor(A(:,:,10))
pcolor(x(:,:,10))
colorbar

f4=figure(4);
pcolor(A(:,:,10))
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=randn(10,10);
f2=figure(2);
pcolor(c)
colorbar

y=find(c>=0.5);
length(y) % this gives the number of cell fullfillin the condition

x = c>=0.5; %create the matrix with the same dimention that fulfils the condition
A=c.*x; % get the values in that position

xx=0; xmax=[];
yy=0; ymax=[];
rmax=[];

x1=0; y1=0; x2=0; y2=0;
for j=1:length(A(:,1))-1
    for i=1:length(A(1,:))-1
        if x(i,j)==0
            xmax(i,j)=x2;
            ymax(i,j)=y2;
            x1=0; y1=0; x2=0; y2=0;
            continue
        elseif x(i,j)==1 && x(i+1,j)==0 && x(i,j+1)==0           
            x1=1; y1=1;
            x2=x1; y2=y1;
             
            disp('one')
        elseif x(i,j)==1 && x(i+1,j)==1 && x(i,j+1)==0
            x1=1; y1=1; 
            x2=x2+x1; y2=y1;
            
            disp('two')
        elseif x(i,j)==1 && x(i+1,j)==1 && x(i,j+1)==1
            x1=1; y1=1;
            x2=x2+x1; y2=y2+y1;
            
            disp('three')      
            
        elseif x(i,j)==1 && x(i+1,j)==0 && x(i,j+1)==1
            x1=1; y1=1;
            x2=x1; y2=y2+y1;
              
            disp('four')
            
        else 

            disp('strange')
        end 
       
    end
end

f3=figure(3);
pcolor(x) %This way of ploting neglets the borders
colorbar

f4=figure(4);
pcolor(xmax)
colorbar

f5=figure(5);
pcolor(ymax)
colorbar








f6=figure(6);
pcolor(A)
colorbar

%to plot isosurfaces
fv = isosurface(X,Y,Z,V,isovalue)

