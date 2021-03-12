P2D=spectravi.P2D;
kper=spectravi.kper;
kpar=spectravi.kpar;

kpermax=max(max(kper));
kparmax=max(max(kpar));

kper_new=linspace(0,kpermax, 500);
kpar_new=linspace(0,kparmax, 500);
k0_new=sqrt(kper_new.*kper_new + kpar_new.*kpar_new);

kper2D = (kper).* ones(size(P2D));
ko=sqrt(kper.*kper + kpar.*kpar);


i=2;

P2Dk = P2D((k0_new(i) < kper) & ( < k0_new(i+1)));

f1=figure(1)
pcolor(P2D)