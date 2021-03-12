%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These commands from the code indeed produce a normal (lognormal) distribution 

% ran1 = random() / ((float) RAND_MAX + 1);
% ran2 = random() / ((float) RAND_MAX + 1);

% double pxi = npt->p[0] +
%     sqrtf(-2.f*npt->T[0]/npt->m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);


% The cos() function actually is necesary to make the distribution normal
% without it the distribution is shifted to the left


%Thus, if the random variable X is log-normally distributed, then Y = ln(X) has a normal distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0=1*ones(100000,1);
p02=1*randn(100000,1);
T=1*ones(100000,1);
m=ones(100000,1);
beta=ones(100000,1);
k=5;
rng(k);
ran1=rand(100000,1); %randn doesnt allow to use pd.
ran2=rand(100000,1);
p0cos=p0.*cos(0.5);

%gaussian around a constant value
pxi_1 = p0+sqrt(-2*(T./m).*(beta.*beta).*log(1-ran1)).* cos(2*pi*ran2);
%p02 is a normal distributed variable
pxi_2 = p02+sqrt(-2*(T./m).*(beta.*beta).*log(1-ran1)).* cos(2*pi*ran2);

pxi_3 = p0+sqrt(-2*(T./m).*(beta.*beta).*log(1-ran1));

pxi_4 = p0+cos(2*pi*ran2);

%y_2=log(pxi_2); %normal log
f888=figure(888);
subplot 221, 
pd_1 = fitdist(pxi_1,'Normal');
hispxi_1=histogram(pxi_1);
ylabel('P1','Interpreter','latex','FontSize',18)
xlabel('$$|Bins|$$','Interpreter','latex','FontSize',18)
subplot 222, 
pd_2 = fitdist(pxi_2,'Normal');
hispxi_2=histogram(pxi_2);
ylabel('P2','Interpreter','latex','FontSize',18)
xlabel('$$|Bins|$$','Interpreter','latex','FontSize',18)
subplot 223, 
pd_3 = fitdist(pxi_3,'Normal');
hispxi_3=histogram(pxi_3);
ylabel('P3','Interpreter','latex','FontSize',18)
xlabel('$$|Bins|$$','Interpreter','latex','FontSize',18)
 subplot 224, 
pd_4 = fitdist(pxi_4,'Normal');
hispxi_4=histogram(pxi_4);
ylabel('P4','Interpreter','latex','FontSize',18)
xlabel('$$|Bins|$$','Interpreter','latex','FontSize',18)
 



 pd_2 = fitdist(pxi_2,'Normal');
hispxi_2=histogram(pxi_2);


f778=figure(778);

histfit(pxi_2);


y_2 = y(~isnan(y_2));
f2=figure(2);
scatter(pxi_2,y_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This distribution produces a Normal distribution
   %    mu = 0.000313782   [-0.00588749, 0.00651505]
   %    sigma =     1.00052   [0.996158, 1.00493]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
