f1=figure(1);
hold on
set(gca,'FontSize',18)
plot(histoEparppc50(:,1),histoEparppc50(:,2)/rms(histoEparppc50(:,2),'all'),'r')
plot(histoEparppc100(:,1),histoEparppc100(:,2)/rms(histoEparppc100(:,2),'all'), 'g')
plot(histoEparppc150(:,1),histoEparppc150(:,2)/rms(histoEparppc150(:,2),'all'),'b')
plot(histoEparppc200(:,1),histoEparppc200(:,2)/rms(histoEparppc200(:,2),'all'),'c')
plot(histoEparppc250(:,1),histoEparppc250(:,2)/rms(histoEparppc250(:,2),'all'),'m')
plot(histoEparppc300(:,1),histoEparppc300(:,2)/rms(histoEparppc300(:,2),'all'),'k')
plot(histoEparsubset(:,1),histoEparsubset(:,2)/rms(histoEparsubset(:,2),'all'),'*k')
legend({'$ppc_{50}$','$ppc_{100}$','$ppc_{150}$','$ppc_{200}$','$ppc_{250}$','$ppc_{300}$','$CB104_{100}$'},'Interpreter','latex')
    %xlabel('$k_\perp d_i$','Interpreter','latex')
    ylabel('$E_{par}$','Interpreter','latex')
    title('Normalised to rms')
hold off

m=200*200*200;
mm=400*400*800;

f2=figure(2);
hold on
set(gca,'FontSize',18)
plot(histoEparppc50(:,1),histoEparppc50(:,2)/m,'r')
plot(histoEparppc100(:,1),histoEparppc100(:,2)/m, 'g')
plot(histoEparppc150(:,1),histoEparppc150(:,2)/m,'b')
plot(histoEparppc200(:,1),histoEparppc200(:,2)/m,'c')
plot(histoEparppc250(:,1),histoEparppc250(:,2)/m,'m')
plot(histoEparppc300(:,1),histoEparppc300(:,2)/m,'k')
plot(histoEparsubset(:,1),histoEparsubset(:,2)/mm,'*k')
legend({'$ppc_{50}$','$ppc_{100}$','$ppc_{150}$','$ppc_{200}$','$ppc_{250}$','$ppc_{300}$','$CB104_{100}$'},'Interpreter','latex')
    %xlabel('$k_\perp d_i$','Interpreter','latex')
    ylabel('$E_{par}$','Interpreter','latex')
    title('Normalised to number of cells')
hold off

 cd '/Volumes/PSC_DiRAC_DATA/TEST2020_1/changing_ppc/Images';
    % Save the plots
    %-------------------------------------------------------------------------
    saveas(f1,'Epar_histogram.png');
    saveas(f2,'Epar_histogram_subset.png');