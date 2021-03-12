%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is to make the video along the box

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis tight manual 
set(gca,'nextplot','replacechildren'); 

v = VideoWriter('peaks_2.avi');
open(v);

for k = 1:2015 
   contourf(squeeze(B_zero(:,:,k)))
   hold on
   scatter(191.599143698154,394.816519489787)
   hold off
   frame = getframe(gcf);
   writeVideo(v,frame);
end

close(v);