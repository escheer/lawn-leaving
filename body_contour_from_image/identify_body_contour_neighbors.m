function neighbor_matrix = identify_body_contour_neighbors(body_contour_coords)

len_contour = max(size(body_contour_coords));

neighbor_matrix=[];

if(len_contour > 2)
    
    neighbor_matrix(len_contour,len_contour) = 0;
    
    for(i=1:len_contour)
        for(j=i+1:len_contour)
            if( abs(body_contour_coords(i,1) - body_contour_coords(j,1)) <= 1 && ...
                    abs(body_contour_coords(i,2) - body_contour_coords(j,2)) <= 1 )
                neighbor_matrix(i,j) = j;
                neighbor_matrix(j,i) = i;
            end
        end
    end
    
end

return;
end
