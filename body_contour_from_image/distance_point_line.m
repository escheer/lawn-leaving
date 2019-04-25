function d = distance_point_line(x,y, m, b)
% d = distance_point_line(x,y, m, b)
% shortest distance from a point x,y and a line y=mx+b

d = abs(y - m*x - b)/sqrt(m^2 +1);

return;
end
