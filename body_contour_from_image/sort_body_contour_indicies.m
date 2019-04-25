function sorted_indicies = sort_body_contour_indicies(body_contour, neighbor_matrix, contour_ends_idx)

sorted_indicies = [];
if(isempty(contour_ends_idx))
    return;
end

length_contour = length(body_contour.x);

if(contour_ends_idx(1) > 0)
    % pick one of the body ends and find its neighbor
    sorted_indicies(1) = contour_ends_idx(1);
    neighbor_matrix(find(neighbor_matrix==sorted_indicies(1))) = 0;
    
    idx_prev = neighbor_matrix(sorted_indicies(1), find(neighbor_matrix(sorted_indicies(1),:)>0));
    
    if(~isempty(idx_prev))
        sorted_indicies(2) = idx_prev(1);
        neighbor_matrix(find(neighbor_matrix==sorted_indicies(2))) = 0;
        
        k=3;
        while(k<=length_contour)
            idx = neighbor_matrix(sorted_indicies(k-1), find(neighbor_matrix(sorted_indicies(k-1),:)>0));
            
            if(~isempty(idx))
                for(p=1:length(idx))
                    sorted_indicies(k) = idx(p);
                    neighbor_matrix(find(neighbor_matrix==sorted_indicies(k))) = 0;
                    k=k+1;
                end
            else
                idx = neighbor_matrix(sorted_indicies(k-2), find(neighbor_matrix(sorted_indicies(k-2),:)>0));
                if(~isempty(idx))
                    for(p=1:length(idx))
                        sorted_indicies(k) = idx(p);
                        neighbor_matrix(find(neighbor_matrix==sorted_indicies(k))) = 0;
                        k=k+1;
                    end
                else
                    break
                end
            end
        end
        
    end
    
end

% if(length(contour_ends_idx)>1)
% if(contour_ends_idx(2) > 0)
%     sorted_indicies(length_contour) = contour_ends_idx(2);
%     idx_prev = neighbor_matrix(sorted_indicies(length_contour), find(neighbor_matrix(sorted_indicies(length_contour),:)>0));
%
%     sorted_indicies(length_contour-1) = idx_prev;
%
%     k=length_contour-2;
%     while(k>=1)
%         idx = neighbor_matrix(sorted_indicies(k+1), find(neighbor_matrix(sorted_indicies(k+1),:)>0 & neighbor_matrix(sorted_indicies(k+1),:)~=sorted_indicies(k+2)));
%
%         if(~isempty(idx))
%             idx = idx(1);
%             sorted_indicies(k) = idx;
%             idx_prev = idx;
%         else
%             break
%         end
%         k=k-1;
%     end
% end
% end


return;
end

