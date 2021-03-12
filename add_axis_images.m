

img = imread('/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images/Bxy0_1.png');
img2 = imread('/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images/Bxy24_1.png');
img3 = imread('/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images/Bz24_1.png');
img4 = imread('/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images/J24_1.png');
img5 = imread('/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images/J243_1.png');
img6 = imread('/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images/T_811_1.png');
img7 = imread('/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images/T_10_1.png');


f1=figure(1);
image(img)    % display Earth image
axis image
ax = gca;
properties(ax)
ax.FontSize=18;
ax.TickLength = [0.04, 0.04];
ax.LineWidth = 1.5;
%set(gca,'FontSize',15)
oldTick_x = get(ax,'XTick');
oldTick_y = get(ax,'YTick'); %get the current tick points of the y axis
newTickStr_x = cellstr(num2str(floor(oldTick_x'/24)));
newTickStr_y = cellstr(num2str(floor(oldTick_y'/24))); %create a cell array of strings
set(ax,'XTickLabel',newTickStr_x)
set(ax,'YTickLabel',newTickStr_y)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',20)
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',20)


f2=figure(2);
image(img2)    % display Earth image
axis image
ax = gca;
properties(ax)
ax.FontSize=18;
ax.TickLength = [0.04, 0.04];
ax.LineWidth = 1.5;
%set(gca,'FontSize',15)
oldTick_x = get(ax,'XTick');
oldTick_y = get(ax,'YTick'); %get the current tick points of the y axis
newTickStr_x = cellstr(num2str(floor(oldTick_x'/24)));
newTickStr_y = cellstr(num2str(floor(oldTick_y'/24))); %create a cell array of strings
set(ax,'XTickLabel',newTickStr_x)
set(ax,'YTickLabel',newTickStr_y)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',20)
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',20)

f3=figure(3);
image(img3)    % display Earth image
axis image
ax = gca;
properties(ax)
ax.FontSize=18;
ax.TickLength = [0.04, 0.04];
ax.LineWidth = 1.5;
%set(gca,'FontSize',15)
oldTick_x = get(ax,'XTick');
oldTick_y = get(ax,'YTick'); %get the current tick points of the y axis
newTickStr_x = cellstr(num2str(floor(oldTick_x'/24)));
newTickStr_y = cellstr(num2str(floor(oldTick_y'/24))); %create a cell array of strings
set(ax,'XTickLabel',newTickStr_x)
set(ax,'YTickLabel',newTickStr_y)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',20)
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',20)


f4=figure(4);
image(img4)    % display Earth image
axis image
ax = gca;
properties(ax)
ax.FontSize=18;
ax.TickLength = [0.04, 0.04];
ax.LineWidth = 1.5;
%set(gca,'FontSize',15)
oldTick_x = get(ax,'XTick');
oldTick_y = get(ax,'YTick'); %get the current tick points of the y axis
newTickStr_x = cellstr(num2str(floor(oldTick_x'/24)));
newTickStr_y = cellstr(num2str(floor(oldTick_y'/24))); %create a cell array of strings
set(ax,'XTickLabel',newTickStr_x)
set(ax,'YTickLabel',newTickStr_y)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',20)
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',20)


f5=figure(5);
image(img5)    % display Earth image
axis image
ax = gca;
properties(ax)
ax.FontSize=18;
ax.TickLength = [0.01, 0.01];
ax.LineWidth = 1.3;
%set(gca,'FontSize',15)
oldTick_x = get(ax,'XTick');
oldTick_y = get(ax,'YTick'); %get the current tick points of the y axis
newTickStr_x = cellstr(num2str(floor(oldTick_x'/22)));
newTickStr_y = cellstr(num2str(floor(oldTick_y'/10))); %create a cell array of strings
set(ax,'XTickLabel',newTickStr_x)
set(ax,'YTickLabel',newTickStr_y)
xlabel('$z/d_{i}$','Interpreter','latex','FontSize',20)
ylabel('$x/d_{i}$','Interpreter','latex','FontSize',20)


f6=figure(6);
image(img6)    % display Earth image
axis image
ax = gca;
properties(ax)
ax.FontSize=18;
ax.TickLength = [0.04, 0.04];
ax.LineWidth = 1.5;
%set(gca,'FontSize',15)
oldTick_x = get(ax,'XTick');
oldTick_y = get(ax,'YTick'); %get the current tick points of the y axis
newTickStr_x = cellstr(num2str(floor(oldTick_x'/24)));
newTickStr_y = cellstr(num2str(floor(oldTick_y'/24))); %create a cell array of strings
set(ax,'XTickLabel',newTickStr_x)
set(ax,'YTickLabel',newTickStr_y)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',20)
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',20)

f7=figure(7);
image(img7)    % display Earth image
axis image
ax = gca;
properties(ax)
ax.FontSize=18;
ax.TickLength = [0.04, 0.04];
ax.LineWidth = 1.5;
%set(gca,'FontSize',15)
oldTick_x = get(ax,'XTick');
oldTick_y = get(ax,'YTick'); %get the current tick points of the y axis
newTickStr_x = cellstr(num2str(floor(oldTick_x'/24)));
newTickStr_y = cellstr(num2str(floor(oldTick_y'/24))); %create a cell array of strings
set(ax,'XTickLabel',newTickStr_x)
set(ax,'YTickLabel',newTickStr_y)
xlabel('$x/d_{i}$','Interpreter','latex','FontSize',20)
ylabel('$y/d_{i}$','Interpreter','latex','FontSize',20)


%----------------------------------------------------------------------
    cd '/Volumes/PSC_DiRAC_DATA/Analysis_CB104_1/tmp_images/modify_images';
    % Save the plots
%-------------------------------------------------------------------------

  saveas(f1,'Bxy_2.png');
  saveas(f2,'Bxy24_2.png');
  saveas(f3,'Bz24_2.png');
  saveas(f4,'J24_2.png');
  saveas(f5,'J243_2.png');
  saveas(f6,'T_811_2.png');
  saveas(f7,'T_10_2.png');
  
  saveas(f1,'Bxy_2','epsc');
  saveas(f2,'Bxy24_2','epsc');
  saveas(f3,'Bz24_2','epsc');
  saveas(f4,'J24_2','epsc');
  saveas(f5,'J243_2','epsc');
  saveas(f6,'T_811_2','epsc');
  saveas(f7,'T_10_2','epsc');
  'epsc'