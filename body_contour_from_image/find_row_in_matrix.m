function row_idx = find_row_in_matrix(a, target_row)

num_row = size(a,1);
num_col = size(a,2);

row_idx = [];

for(i=1:num_row)
    found_flag=1;
    for(j=1:num_col)
       if(target_row(j) ~= a(i,j))
           found_flag=0;
           break;
       end
    end
    if(found_flag==1)
        row_idx = [row_idx i];
    end
end

return;
end
