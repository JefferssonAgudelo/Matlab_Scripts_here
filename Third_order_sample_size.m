%--------------------------------------------------------------------------
% THIS IS TO ESTIMATE THE NUMBER OF SAMPLES REQUIRED TO ACCURETALY ESTIMATE
% THE HIGHGER ORDER MOMENTS
%--------------------------------------------------------------------------

% Define the first distribution
s = rng;
A=1000000;
mu_t=3;
sigma_t=10;
r = normrnd(mu_t,sigma_t,[1,A]);
r2 = normrnd(0,sigma_t,[1,A]);

%For an easy normalised case
%--------------------------------------------------------------------------
a=-5*sigma_t;
b=5*sigma_t;
dx=(b-a)/A;
x = [a:dx:b];
pd = makedist('Normal','mu',mu_t,'sigma',sigma_t);
y = pdf(pd,x);
dx2=(b-a)/500;
x2 = [a:dx2:b];

[Ny,edgesy] = histcounts(y,x2);

sum(Ny.*edgesy(2:end));
%--------------------------------------------------------------------------

% Select the samples 
%--------------------------------------------------------------------------
%R1 = randsample(r,A*0.01);
%R2 = randsample(r,A*0.1);
%R3 = randsample(r,A*0.2);
%R4 = randsample(r,A*0.5);
%R5 = randsample(r,A*1.0);
%--------------------------------------------------------------------------

% Comopute the moments
%--------------------------------------------------------------------------

sig_3_to_mu_3_1=zeros(1,5);
r_mean=zeros(1,5);
r_std=zeros(1,5);
r_third=zeros(1,5);
NN=zeros(1,i);

for i=1:5
    rs=randsample(r, A*(1/(10^i)) );
N=length(rs);
r_m1 = sum(rs.^1)/N;
r_m2 = sum(rs.^2)/N;
r_m3 = sum(rs.^3)/N;
r_m4 = sum(rs.^4)/N;
r_m5 = sum(rs.^5)/N;
r_m6 = sum(rs.^6)/N;
r_mean(1,i)=r_m1;
r_std(1,i)=sqrt(r_m2);
r_third(1,i)=r_m3;
sig_3_to_mu_3_1(1,i) = (N^(-1/2)) * (abs( r_m6/(r_m3^2) - 1)).^(1\2);
NN(1,i)=N;
end

%--------------------------------------------------------------------------


f21=figure(21);
subplot(2,2,1)
plot(NN,r_mean)
subplot(2,2,2)
plot(NN,r_std)
subplot(2,2,3)
plot(NN,r_third)
subplot(2,2,4)
semilogy(NN,sig_3_to_mu_3_1)


plot(edgesy(2:end),Ny,'o')

