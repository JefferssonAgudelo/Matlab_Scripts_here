%--------------------------------------------------------------------------
% THIS IS TO ESTIMATE THE NUMBER OF SAMPLES REQUIRED TO ACCURETALY ESTIMATE
% THE HIGHGER ORDER MOMENTS
%--------------------------------------------------------------------------

% Define the first distribution
s = rng;
A=1000000;
mu_t=3;
sigma_t=2;
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
ii=6;
kk=500;
sig_3_to_mu_3_1=zeros(kk,ii);
r_mean=zeros(kk,ii);
r_std=zeros(kk,ii);
r_third=zeros(kk,ii);
NN=zeros(kk,ii);

for k=1:kk
    for i=1:ii
    rs=randsample(r, A*(1/(10^(6-i))) );
    N=length(rs);
    r_m1 = sum(rs.^1)/N;
    r_m2 = sum((rs-r_m1).^2)/N;
    r_m3 = sum(rs.^3)/N;
    r_m4 = sum(rs.^4)/N;
    r_m5 = sum(rs.^5)/N;
    r_m6 = sum(rs.^6)/N;
    r_mean(k,i)=r_m1;
    r_std(k,i)=sqrt(r_m2);
    r_third(k,i)=r_m3;
    sig_3_to_mu_3_1(k,i) = (N^(-1/2)) * (abs( (r_m6/(r_m3^2)) - 1)).^(1\2);
    NN(k,i)=N;
    end
end
%--------------------------------------------------------------------------
lw=1.5;
f21=figure(21);
subplot(2,2,1)
ax=gca;
plot(NN,mean(r_mean),'-*b','LineWidth', lw);
ylabel('mean'); set(gca,'XScale','log','YScale','lin','FontSize',12);
ax.XTick = [10^1 10^2 10^3 10^4 10^5];
subplot(2,2,2)
ax=gca;
plot(NN,mean(r_std),'-*b','LineWidth', lw);
ylabel('std'); set(gca,'XScale','log','YScale','lin','FontSize',12)
ax.XTick = [10^1 10^2 10^3 10^4 10^5];
subplot(2,2,3)
ax=gca;
plot(NN,mean(r_third),'-*b','LineWidth', lw);
ylabel('Third'); set(gca,'XScale','log','YScale','lin','FontSize',12)
ax.XTick = [10^1 10^2 10^3 10^4 10^5];
subplot(2,2,4)
ax=gca;
plot(NN,mean(sig_3_to_mu_3_1),'-*b','LineWidth', lw);
ylabel('ratio'); set(gca,'XScale','log','YScale','log','FontSize',12)
yline(0.1,'--k','LineWidth',1)
xline(3000,'--k','LineWidth',1);
ax.XTick = [10^1 10^2 10^3 10^4 10^5];

saveas(f21,strcat('Third_order_thin.png'));
%--------------------------------------------------------------------------


plot(edgesy(2:end),Ny,'o')

