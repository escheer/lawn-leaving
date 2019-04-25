function [m,b] = fit_line(x,y)
% [m,b] = fit_line(x,y)

if(nargin==0)
    disp('[m,b] = fit_line(x,y)')
    return
end

if(length(x)<=3)
    m = (y(end)-y(1))/(x(end)-x(1));
    b = y(1)-m*x(1);
    return;
end

    s = warning('off','all');

    x(isnan(x)) = [];
    y(isnan(y)) = [];
    
    X = [ones(size(x))'  x'];
    
    p = X\y';
    b  = p(1); m  = p(2);

    clear('p'); clear('X');
    
    warning(s);
    
return;
end
