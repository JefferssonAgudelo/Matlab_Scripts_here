function [t,covXxYx] = Local_corr2(t2,Xs,Ys,bin1,window)
 % time/distance, var1, var2, increment, window  
 %Define the arrays
 %Ex: t2=d; Xx=vim/Vai; Yx=hm/B0;
 %Define the increment/window
 
 %Xs=vim/Vai; Ys=hm/B0; bin1=10;
 %window=10; % 10 is the optimun value that do not shift the plot
 %Filter data
 b = (1/window)*ones(1,window);
 a = 1;
 Xs = filter(b,a,Xs);
 Ys = filter(b,a,Ys);  
  
 dta =bin1; dt0=1; lenXs=length(Xs); Nf=fix((lenXs-1)/dta);
 covXxYx=zeros(Nf-1,1);
 %t=zeros(lenXs-1,1);
 t=zeros(Nf-1,1);
 ted=zeros(Nf,1);
 for i=1:Nf
     Xx = Xs(dt0:dta); 
     Yx = Ys(dt0:dta);
     tx = median(t2(dt0:dta));
     ted_1=t2(dta);
     covXxYx_i=sum(((Xx-mean(Xx)).*(Yx-mean(Yx))))./((length(Xx)-1)*std(Xx)*std(Yx));
     covXxYx(i)=covXxYx_i; 
     t(i) = tx;
     ted(i+1)= ted_1;
     dt0 = dta; 
     dta = bin1*(i+1);
     covXxYx(covXxYx==0)=nan;
 end
     
 %outputs
end
