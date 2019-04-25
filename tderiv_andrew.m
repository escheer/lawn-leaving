function derivative = tderiv_andrew(vector,window,fps)
%TDERIV Takes the time derivative of VECTOR using dt = WINDOW. WINDOW
%should be an odd number.  If not, the program will automatically add one
%to it.If no WINDOW is defined, it is assumed WINDOW = 1. If frames per
%second (FPS) is defined, then the output is in dX/dt (AU/sec).  If not,
%output is dX/dframe.

if nargin < 2
    window = 5;
end
if nargin < 3
    fps = 10;
end

derivative = NaN(size(vector));

if mod(window,2)==0
    window = window+1;
end

stepsize = floor(window/2);

for m=stepsize+1:length(vector)-stepsize
    derivative(m) = (vector(m+stepsize) - vector(m-stepsize))/window;
end

derivative = derivative*fps;

end
