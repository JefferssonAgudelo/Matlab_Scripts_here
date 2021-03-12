cd '/Volumes/PSC_DiRAC_DATA/DATACB103_1/Images_things/Bfield_gif'
% load the images
 images    = cell(59,1);
 for i=1:59
     images{i} = imread(sprintf('Bfield_%d.png',i));
 end
 % create the video writer with 30 fps
 writerObj = VideoWriter('Bfield.avi');
 writerObj.FrameRate = 59;
 % open the video writer
 open(writerObj);
 % write the frames to the video
    for u=1:59    
       % convert the image to a frame
       frame = im2frame(images{u});
       %for v=1:30
       writeVideo(writerObj, frame);
       %end
   end
   % close the writer object
   close(writerObj);
   implay('Bfield.avi');