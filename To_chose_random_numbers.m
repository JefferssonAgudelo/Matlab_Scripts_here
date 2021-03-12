%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is to generate and pick random numbers fron a given distribution
%function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Produce the probability distribution function
% Create a standard normal distribution object with the mean ? equal to 0
% and the standard deviation ? equal to 1.

mu = 0;
sigma = 1;
pd = makedist('Normal','mu',mu,'sigma',sigma);

% Define the input vector x to contain the values at which to calculate the pdf.
ax = 0;
bx = 1;
dx = 0.0001; %The thickness
x = -ax:dx:bx;
y = pdf(pd,x);

%check the distribution
f1=figure(1);
plot(x,y)

%%%%% Extract the numbers from the array 
y_max = max(y);

%xk=zeros(5,1); yk=zeros(5,1);
A=zeros(length(x),2);
A2=[];
% randomnly chose the xk and its probability yk
%for i=1:(length(x)-1)/4
i=0;
j=0;
while A(length(x),1)==0
    xpos = randi(length(x)); % choose the random position of x
    x1 = ax+(bx-ax)*x(xpos);
    ypos = randi(length(y)); %chose the randonm position of y
    y1=y_max*y(ypos);
    %y1=y_max*ypos;
    
    if y1 < y(xpos) %If y1 < y(x1), keep the x1. Otherwise choose a couple again
        %xk(i) = x(xpos); yk(i) = y(ypos);
        i=i+1;
        A(i,1) = x(xpos);
        A(i,2) = y(ypos);
    else
        j=j+1;
        A2(i,1) = x(xpos);
        A2(i,2) = y(ypos);
    end
end


%check the distribution of the picked random numbers
f2=figure(2);
scatter(A(:,1),A(:,2))

%
x = rand(n_points,1); % This returns a uniform distribution
pi_gn = randn(n_points,1); %unnormalized function 
psi_max=max(pi_gn); %bounding value
pos = randi(length(x)); 
x1 = a+(b-a)*psi(pos);


% This shows that using a uniform distribution to extract values from a normal distribution works 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_points=50000;
a= 0;
b= 1;
x = randn(n_points,1); % This returns a uniform distribution
arrayA=zeros(n_points/2,1);
for i = 1:n_points/2
    arrayA(i) = x(randi(length(x)));
end
f11=figure(11);
histogram(x,100);
f12=figure(12);
histogram(arrayA,100);

r = randi([-10 10],1,1000);
f13=figure(13);
histogram(r,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points=50000;
x = randn(n_points,1); % This returns a normal distribution
a = -0.5; b = 0.5;
xab = a + (b-a)*x; 

x_sum=sum(x);
dim_x=length(x);
xmean=mean(x,'all');
xmx2_mean= mean((x-xmean).^2, 'all');
xstd=sqrt((dim_x/(dim_x-1))* xmx2_mean);
% r = normrnd(mu,sigma,sz1,...,szN)
rng('default') % For reproducibility
s = rng;
r = normrnd(xmean,xstd,[1,10000]);


% Check the normal distribution of the values
f1=figure(1);
x_histo = histogram(x);
f2=figure(2);
xab_histo = histogram(xab);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots 1D

f1=figure(1);
x_histo = histogram(x2);
hold on
r_histo = histogram(r);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate values from a bivariate normal distribution with specified mean
% vector and covariance matrix.

mu = [1 2]; %bivariate 
sigma = [1 0.5; 0.5 2];
R = chol(sigma);
%repmat(A,n) returns an array containing n copies of A in the row and column dimensions.
%The size of B is size(A)*n when A is a matrix.
z = repmat(mu,10,1) + randn(10,2)*R;


x = rand(100);
x_max=max(max(x));
x_min=min(min(x));

dim_x=length(x)*length(x);

xmean=mean(x,'all');
xmx2_mean= mean((x-xmean).^2, 'all');
xstd=sqrt((dim_x/(dim_x-1))* xmx2_mean);

% r = normrnd(mu,sigma,sz1,...,szN)
rng('default') % For reproducibility
s = rng;
r = normrnd(xmean,xstd,[100,100]);

z1=size(zeros([2 20]));
x_hist = hist3(x);
r_hist = hist3(r);

f1=figure(1);
pcolor(x)
colorbar

f2=figure(2);
pcolor(r)
colorbar

f1=figure(1);
p(x)
colorbar

f2=figure(2);
pcolor(r)
colorbar