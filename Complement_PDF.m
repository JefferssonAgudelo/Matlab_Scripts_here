function [Out_Data, Out_PDF, CHist] = Complement_PDF(Hist, Data_Num, p)
% Generate a 1D vector of data with a PDF specified as the complementary PDF of input historgram. Note that the larger
% Data_Num is, the more Out_PDF will resemble to CHist
% Input
% Hist: PDF/Histogram of data
% Data_Num: Desired number of data to be generated
% p: Precision given by number of digits after 0
% Output
% Out_Data: Generated data as per the complementary PDF
% Out_PDF: The complementary PDF as per Out_Data
% CHist: The complementary PDF as per Hist
% Example
% Hist = [1, 6, 7, 100, 0, 0, 0, 2, 3, 5];
% Data_Number = 100000;
% p = 3


Hist = Hist/sum(Hist);
CHist = 1- Hist;
CHist = CHist/sum(CHist);
CDF_CHist = cumsum(CHist);
CDF_CHist = double(int32(CDF_CHist*10^p))/10^p;

Out_Data = zeros(1, Data_Num);
Out_PDF = zeros(1, length(CDF_CHist));
for i = 1:Data_Num
    % Generate a uniformly distributed variable
    x = double(int32(rand*10^p))/10^p;
    % Inversely index CDF
    Out_Data(i) = Inverse_CDF(x, CDF_CHist);
    temp = floor(Out_Data(i) * length(CDF_CHist));
    Out_PDF(temp) = Out_PDF(temp) + 1;
end

figure;
subplot 221, bar(Hist); 
subplot 222, bar(CHist); 
subplot 223, plot(CDF_CHist);
subplot 224, bar(Out_PDF);

end

function [y] = Inverse_CDF(x, CDF_CHist)

CDF_CHist_Ext = [0, CDF_CHist];
y = 1;
for ind = 1:length(CDF_CHist)
    if (x >= CDF_CHist_Ext(ind)) && (x < CDF_CHist_Ext(ind+1))
        y = ind/length(CDF_CHist);
        break;
    end
end


end