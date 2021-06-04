cd '/Volumes/PSC_DiRAC_DATA/DATACB103_1/Images_things/Bfield_gif'
% load the images
 images    = cell(59,1);
 %n_lides=59;
 n_lides=168;
 for i=1:n_lides
     %images{i} = imread(sprintf('Bfield_%d.png',i));
     f1=figure(1);
     y_plot=dum_p(:,:,i);
    imagesc(y_plot);
    lim_yp=max(max(abs(y_plot)));
    caxis([-lim_yp lim_yp])
    colorbar
    %colorMap = [redColorMap; blueColorMap; zeros(1, 256)]';
    colormap(BWR);
    set(gca,'YDir','normal')
 end
 % create the video writer with 30 fps
 writerObj = VideoWriter('video_1.avi');
 writerObj.FrameRate = n_lides;
 % open the video writer
 open(writerObj);
 % write the frames to the video
    for u=1:n_lides    
       % convert the image to a frame
       %frame = im2frame(images{u});
       frame = im2frame(dum_p(:,:,i));
       %for v=1:30
       writeVideo(writerObj, frame);
       %end
   end
   % close the writer object
   close(writerObj);
   implay('video_1.avi');