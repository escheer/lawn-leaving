function [dist, path_idx] = branch_length(neighbor_matrix, body_contour_coords,branch_coords, branch_idx, tip_idx)
% dist = branch_length(neighbor_matrix, body_contour_coords,branch_coords, branch_idx,tip_idx)
% given neighbor_matrix (neighbor_matrix(i,j)==1 means i,j are neigbors)
% return number of points between branch_idx,tip_idx ;
% if any point between branch_idx and tip_idx is a branch (num_neighbors>2), return 1e10

dist=0;
path_idx = [];

current_idx = tip_idx;
next_idx = find(neighbor_matrix(current_idx,:)>0);
if(length(next_idx)>2) % a branch
    if(isempty(find(next_idx==branch_idx)))
        if(~isempty(find_row_in_matrix(branch_coords, body_contour_coords(next_idx,:))))
            dist=1e10;
            return;
        end
    else
        next_idx(find(next_idx==branch_idx)) = [];
    end
end


prev_idx = tip_idx;
path_idx = [];
while(next_idx~=branch_idx)
    dist = dist+1;
    path_idx = [path_idx current_idx];
    next_idx = find(neighbor_matrix(current_idx,:)>0);
    
    if(length(next_idx)>2) % a branch
        if(isempty(find(next_idx==branch_idx)))
            if(~isempty(find_row_in_matrix(branch_coords, body_contour_coords(next_idx,:))))
                dist=1e10;
                return;
            end
        else
            next_idx(find(next_idx==branch_idx)) = [];
        end
    end
    
    next_idx(find(next_idx==prev_idx)) = [];
    
    if(length(next_idx)>1)
        return;
    end
    
    prev_idx = current_idx;
    current_idx = next_idx;
end

return;
end
