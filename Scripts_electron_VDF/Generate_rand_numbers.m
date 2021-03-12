function [xhat1, xhat2] = Generate_rand_numbers(pxe1, pxe2, num_bins, times_l)
% Jeff 19/05/2020 
% Generate two 1D vectores of data with a PDF for two given data distribution following the rejection method
%. Note that the larger
% Data_Num is, the more Out_PDF will resemble to CHist
% Input
% pxe1, pxe2: Distribution functions
% num_bins: number of bins to be used more 
% times_l: how many times the size of the given PDF data
% Output
% xhat1, xhat2 generated distribution functions


%num_bins=100;
[hist_pxe1,edges_pxe1] = histcounts(pxe1,num_bins); % create the histogram counts
[hist_pxe2,edges_pxe2] = histcounts(pxe2,num_bins);

hist_pxe1_max=max(hist_pxe1); % Define the maximum upper boundary
hist_pxe2_max=max(hist_pxe2);

a1=edges_pxe1(1); % define the edges
b1=edges_pxe1(num_bins+1);
a2=edges_pxe2(1); % define the edges
b2=edges_pxe2(num_bins+1);

%times_l=3; %Chose how big the size of the new data
xhat1=zeros(times_l*length(pxe1),1); % allocate the array to keep the new distribution
xhat2=zeros(times_l*length(pxe2),1);

%first distribution
count=0;
count2=0;
while (xhat1(times_l*length(pxe1),1)==0)% && count~=100000000) 
    xhati=a1 +(b1-a1)*rand; %calculate the random x
    yhati=rand*hist_pxe1_max; % calculate the random distribution yhat    
    for i=1:num_bins % search the position
        if edges_pxe1(1,i) <= xhati && xhati <edges_pxe1(1,i+1)  
            foo=hist_pxe1(1,i); % this is pi_tilde(xhat)
        else
            continue
        end
        if yhati < foo
            count2=count2+1;
            xhat1(count2,1)=xhati; % Save the value of the distribution
            %yhat(count2,1)=yhati;
        else
            continue
        end
    end
        count=count+1;
end
   
%Second distribution
count=0;
count2=0;
while (xhat2(times_l*length(pxe2),1)==0)% && count~=100000000) 
    xhati=a2 +(b2-a2)*rand; %calculate the random x
    yhati=rand*hist_pxe2_max; % calculate the random distribution yhat    
    for i=1:num_bins % search the position
        if edges_pxe2(1,i) <= xhati && xhati <edges_pxe2(1,i+1)  
            foo=hist_pxe2(1,i); % this is pi_tilde(xhat)
        else
            continue
        end
        if yhati < foo
            count2=count2+1;
            xhat2(count2,1)=xhati; % Save the value of the distribution
            %yhat(count2,1)=yhati;
        else
            continue
        end
    end
        count=count+1;
end

[hist_xhat1,edges_xhat1] = histcounts(xhat1,num_bins); % create the histogram counts of the new variable
[hist_xhat2,edges_xhat2] = histcounts(xhat2,num_bins);

% make plots to check the results and see how it works

f1=figure(1);
subplot 221, bar(hist_pxe1); 
ylabel('Log vpar','Interpreter','latex','FontSize',18)
 xlabel('$$|v_{par}|$$','Interpreter','latex','FontSize',18)
subplot 222, bar(hist_pxe2); 
ylabel('Log vper','Interpreter','latex','FontSize',18)
 xlabel('$$|v_{per}|$$','Interpreter','latex','FontSize',18)
subplot 223, bar(hist_xhat1);
ylabel('Extended vpar','Interpreter','latex','FontSize',18)
 xlabel('$$|v_{vpar}|$$','Interpreter','latex','FontSize',18)
subplot 224, bar(hist_xhat2);
ylabel('Extended vper','Interpreter','latex','FontSize',18)
 xlabel('$$|v_{vper}|$$','Interpreter','latex','FontSize',18)

end

