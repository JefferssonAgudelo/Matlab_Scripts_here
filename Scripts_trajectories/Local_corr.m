function [t,dXdY]=Local_corr(t2,Xx,Yx,dd,window)
 % time/distance, var1, var2, increment, window  
 %Define the arrays
 %Ex: t2=d; Xx=vim/Vai; Yx=hm/B0;
 %Define the increment/window
 %dd=10; window=10; % 10 is the optimun value that do not shift the plot
 %Filter data
 b = (1/window)*ones(1,window);
 a = 1;
 Xx = filter(b,a,Xx);
 Yx = filter(b,a,Yx);  
 %Calculate the increments (derivatives) 
 Xx2=Xx(1+dd:length(Xx)); Xx1=Xx(1:length(Xx)-dd);
 dXx=(Xx2 - Xx1)/dd;
 Yx2=Yx(1+dd:length(Yx)); Yx1=Yx(1:length(Yx)-dd);
 dYx=(Yx2 - Yx1)/dd;
 %outputs
 dXdY=dXx.*dYx;
 %dXdY=dXdY/max(abs(dXdY));
 t=t2(1:length(t2)-dd);
 end