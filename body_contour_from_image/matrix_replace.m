function x = matrix_replace(y, relationship, a,  b)
% x = matrix_replace(y, relationship, a,  b)
% replaces the elements relationship a in matrix y w/ b
% >> q = [1 2 3; 4 3 6; 3 8 9]
% q =
%      1     2     3
%      4     3     6
%      3     8     9
%      
% >> matrix_replace(q,'==',3,-1)
% ans =
%      1     2    -1
%      4    -1     6
%     -1     8     9
%     
% >> matrix_replace(q,'<',3,-1)
% ans =
%     -1    -1     3
%      4     3     6
%      3     8     9    

if(nargin==0)
    disp('x = matrix_replace(y, relationship, a,  b)')
    return
end

x = y;

if(isnan(a))
    st = sprintf('x(isnan(x)) = %f;',b);
else
    if(isinf(a))
        st = sprintf('x(isinf(x)) = NaN;');
    else
        if(strcmp(relationship,'=') || strcmp(relationship,'==')) % for replacing values of a w/ b need to deal w/ floating point stuff
            st = sprintf('x(find(abs(x-%f)<1e-6)) = %f;',a, b);
        else
            st = sprintf('x(x %s %f) = %f;',relationship,a,b);
        end
    end
end

eval(st);

return;
end

