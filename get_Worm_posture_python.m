function [ spline, head, tail ] = get_Worm_posture_python( worm, bg, x_offset, y_offset )
%get_Worm_posture_python this is just a wrapper that calls the python
%function called getposture_ELIAS_WIP.get_posture();

% I ravel X to a row vector, and unravel with Numpy
np_worm = py.numpy.reshape(worm(:)',size(worm),'F'); 
%this can now be handed
np_worm = py.numpy.asarray(np_worm,'float32');

% I ravel X to a row vector, and unravel with Numpy
np_bg = py.numpy.reshape(bg(:)',size(bg),'F'); 
%this can now be handed
np_bg = py.numpy.asarray(np_bg,'float32');

posturestuff = py.getposture_ELIAS_WIP.get_posture(np_worm,np_bg);

midline_y = int32(posturestuff{1})'+1;%+y_offset; %add 1 to all to change from python to matlab indexing
midline_x = int32(posturestuff{2})'+1;%+x_offset;
spline = [midline_x midline_y];
head_1 = int32(posturestuff{3})+1;%+y_offset;
head_2 = int32(posturestuff{4})+1;%+x_offset;
head = [head_2 head_1];
tail_1 = int32(posturestuff{5})+1;%+y_offset;
tail_2 = int32(posturestuff{6})+1;%+x_offset;
tail = [tail_2 tail_1];


end

