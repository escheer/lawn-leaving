function d = point_to_line(pt, v1, v2)
v1 = repmat(v1,size(pt,1),1);
v2 = repmat(v2,size(pt,1),1);
a = v1 - v2;
b = pt - v2;
c_pdt = cross(a,b,2);
c_pdt_norm = sqrt(sum(abs(c_pdt).^2,2));
a_norm = sqrt(sum(abs(a).^2,2));
d = c_pdt_norm ./ a_norm;
end
