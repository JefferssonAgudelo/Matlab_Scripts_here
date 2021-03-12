clear all
close all
clc

%% Define how many rows/columns you want as well as the percent spacing you
%% wish to have between each subplot

rows = 3;
cols = 16;
space = 0;

%% Creates a matrix of handles to your custom-spaced subplots.
ax = getCustomAxesPos(rows,cols,space);

%% Put what you want to plot here.
for i = 1:rows,
     for j = 1:cols,
         
         x = 1:3;
         y = rand(1,length(x));
         
         plot(ax(i,j),x,y);      % Plots (x,y) in the (i,j)th subplot
         set(ax(i,j),'xtick',[],'ytick',[]);
     end
end