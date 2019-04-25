function body_contour_struct = body_contour_from_image(Image, bound_box_corner)


% creates the thin line trace of the body from an image bounding_box_corner

body_contour_struct.x=[];
body_contour_struct.y=[];
% body_contour_struct.midbody=0;
% body_contour_struct.head=0;
% body_contour_struct.tail=0;
% body_contour_struct.neck=0;
% body_contour_struct.lumbar=0;
body_contour_struct.curvature = [];

if(isempty(Image))
    return;
end

if(nargin<2)
    bound_box_corner = [0 0];
end

perimeter_image = [];
skeleton = [];


if(sum(sum(Image))>100)
    % find the smooth perimeter
    perimeter_image = bwperim(bwmorph(bwmorph(bwmorph(bwmorph(imfill(Image,'holes'),'spur'),'majority'),'clean'),'remove'));
    %ELIAS EDIT 09/09/17
%     perimeter_image = bwperim(Image);
    if(sum(sum(perimeter_image))==0)
        return;
    end
    % and fill it in for the working local_image
    skeleton = bwmorph(bwmorph(  bwmorph(imfill(perimeter_image,'holes'),'majority')   ,'thin',Inf),'clean');
else
    skeleton = bwmorph(bwmorph(Image,'spur'),'thin',Inf);
end
if(sum(sum(skeleton))==0)
    return;
end



[body_contour_coords(:,2), body_contour_coords(:,1)] = find(skeleton==1);
neighbor_matrix = identify_body_contour_neighbors(body_contour_coords);

endpoints = bwmorph(skeleton, 'endpoints');
[end_coords(:,2), end_coords(:,1)] = find(endpoints==1);
if(size(end_coords,1)>2) % too many endpoints ... find the likely real ones
    
    branchpoints = bwmorph(skeleton, 'branchpoints');
    r = []; c = [];
    [r, c] = find(branchpoints==1);
    
    if(~isempty(r))
        branch_coords(:,2) = r;
        branch_coords(:,1) = c;
    else
        return;
    end
    
    % %close all
    % figure(2);
    % imshow(skeleton);
    % hold on
    % plot(end_coords(:,1), end_coords(:,2),'.r');
    % plot(branch_coords(:,1), branch_coords(:,2),'.g');
    % hold off
    % pause
    % %     x=10;
    
    % too many ends, so delete the shortest branch from each branchpoint
    ctr=1;
    prev_num_endpoints = size(end_coords,1);
    while(size(end_coords,1)>2 && ctr<=3*size(branch_coords,1))
        ctr = ctr+1;
        
        deleted_branch_ctr=0;
        for(b=1:size(branch_coords,1))
            mindist = 1e11;
            path_idx = [];
            
            b1 = find_row_in_matrix(body_contour_coords, branch_coords(b,:));
            
            for(e=1:size(end_coords,1))
                
                e1 = find_row_in_matrix(body_contour_coords, end_coords(e,:));
                
                [dist, p_idx] = branch_length(neighbor_matrix,  body_contour_coords,branch_coords, b1,e1);
                if(dist<mindist)
                    mindist = dist;
                    path_idx = p_idx;
                end
                
                
            end
            % wipe out pixels for the short branch affiliated w/ branchpoints(b)
            for(p=1:length(path_idx))
                deleted_branch_ctr = deleted_branch_ctr+1;
                skeleton(body_contour_coords(path_idx(p),2), body_contour_coords(path_idx(p),1)) = 0;
            end
        end
        
        body_contour_coords=[];
        [body_contour_coords(:,2), body_contour_coords(:,1)] = find(skeleton==1);
        neighbor_matrix = identify_body_contour_neighbors(body_contour_coords);
        end_coords =[];
        endpoints = bwmorph(skeleton, 'endpoints');
        [end_coords(:,2), end_coords(:,1)] = find(endpoints==1);
        branch_coords = [];
        branchpoints = bwmorph(skeleton, 'branchpoints');
        [branch_coords(:,2), branch_coords(:,1)] = find(branchpoints==1);
        
        % should have fewer endpoints now since we have deleted branch(s)
        % if we have more endpoints than we should, a gap was created
        if(deleted_branch_ctr>0)
            if(size(end_coords,1) > (prev_num_endpoints-deleted_branch_ctr))
                for(i=1:size(end_coords,1)-1)
                    for(j=i+1:size(end_coords,1))
                        d = sqrt((end_coords(i,1) - end_coords(j,1))^2 + (end_coords(i,2) - end_coords(j,2))^2);
                        if(d <= 3)
                            body_contour_coords = [body_contour_coords; round([(end_coords(i,1) + end_coords(j,1))/2 (end_coords(i,2) + end_coords(j,2))/2])];
                            skeleton(body_contour_coords(end,2), body_contour_coords(end,1))=1;
                        end
                    end
                end
                
                body_contour_coords=[];
                [body_contour_coords(:,2), body_contour_coords(:,1)] = find(skeleton==1);
                neighbor_matrix = identify_body_contour_neighbors(body_contour_coords);
                end_coords =[];
                endpoints = bwmorph(skeleton, 'endpoints');
                [end_coords(:,2), end_coords(:,1)] = find(endpoints==1);
                branch_coords = [];
                branchpoints = bwmorph(skeleton, 'branchpoints');
                [branch_coords(:,2), branch_coords(:,1)] = find(branchpoints==1);
                
            end
        end
        
        
        prev_num_endpoints = size(end_coords,1);
        
        %     x=x+1; figure(x);
        %     imshow(skeleton);
        %     hold on
        %     plot(body_contour_coords(:,1), body_contour_coords(:,2),'ob');
        %     plot(end_coords(:,1), end_coords(:,2),'.r');
        %     plot(branch_coords(:,1), branch_coords(:,2),'.g');
        %     hold off
        
    end
    
    % figure(1);
    % subplot(1,2,1);
    % imshow(Image);
    % hold on;
    % plot(body_contour_coords(:,1), body_contour_coords(:,2),'.r')
    % hold off
    
    % end_coords
    
    ctr = 0;
    if(size(end_coords,1)>2 && ctr<=2*size(end_coords,1)) % still too many ends ... are they seperated by one pixel?
        ctr=ctr+1;
        
        for(i=1:size(end_coords,1)-1)
            for(j=i+1:size(end_coords,1))
                d = sqrt((end_coords(i,1) - end_coords(j,1))^2 + (end_coords(i,2) - end_coords(j,2))^2);
                if(d <= 3)
                    body_contour_coords = [body_contour_coords; round([(end_coords(i,1) + end_coords(j,1))/2 (end_coords(i,2) + end_coords(j,2))/2])];
                    skeleton(body_contour_coords(end,2), body_contour_coords(end,1))=1;
                end
            end
        end
        
        neighbor_matrix = identify_body_contour_neighbors(body_contour_coords);
        end_coords =[];
        endpoints = bwmorph(skeleton, 'endpoints');
        [end_coords(:,2), end_coords(:,1)] = find(endpoints==1);
        branch_coords = [];
        branchpoints = bwmorph(skeleton, 'branchpoints');
        [branch_coords(:,2), branch_coords(:,1)] = find(branchpoints==1);
    end
    
end

bodyend_indicies = [];
for(i=1:size(end_coords,1))
    bodyend_indicies = [bodyend_indicies find_row_in_matrix(body_contour_coords, end_coords(i,:))];
end

body_contour_struct.x = body_contour_coords(:,1) + floor(bound_box_corner(1)) + 0.5;
body_contour_struct.y = body_contour_coords(:,2) + floor(bound_box_corner(2)) + 0.5;

if(~isempty(neighbor_matrix))
    
    if(~isempty(bodyend_indicies))
        sorted_indicies = sort_body_contour_indicies(body_contour_struct, neighbor_matrix, bodyend_indicies);
        
        body_contour_struct.x = body_contour_struct.x(sorted_indicies);
        body_contour_struct.y = body_contour_struct.y(sorted_indicies);
        
        % include two more points ... the point on the perimeter closest
        % to the tangent from each bodyend
        if(~isempty(perimeter_image) && length(body_contour_struct.x)>3)
            [perimeter_coords(:,2), perimeter_coords(:,1)] = find(perimeter_image==1);
            firstcoord = [body_contour_struct.x(1) body_contour_struct.y(1)];
            lastcoord = [body_contour_struct.x(end) body_contour_struct.y(end)];
            headcoord = [];
            tailcoord = [];
            [mhead,bhead] = fit_line(body_contour_struct.x(1:3), body_contour_struct.y(1:3));
            [mtail,btail] = fit_line(body_contour_struct.x(end-2:end), body_contour_struct.y(end-2:end));
            minheaddist = 1e10; mintaildist = 1e10;
            
            for(j=1:size(perimeter_coords,1))
                
                d = distance_point_line(perimeter_coords(j,1), perimeter_coords(j,2), mhead, bhead);
                if(d<=1)
                    d = ((firstcoord(1) - perimeter_coords(j,1) )^2 + (firstcoord(2) - perimeter_coords(j,2) )^2);
                    if(d<minheaddist)
                        headcoord = perimeter_coords(j,:);
                        minheaddist = d;
                    end
                end
                
                d = distance_point_line(perimeter_coords(j,1), perimeter_coords(j,2), mtail, btail);
                if(d<=1)
                    d = ((lastcoord(1) - perimeter_coords(j,1) )^2 + (lastcoord(2) - perimeter_coords(j,2) )^2);
                    if(d<mintaildist)
                        tailcoord = perimeter_coords(j,:);
                        mintaildist = d;
                    end
                end
            end
            
            if(~isempty(headcoord))
                body_contour_struct.x = [headcoord(1);  body_contour_struct.x];
                body_contour_struct.y = [headcoord(2);  body_contour_struct.y];
            end
            if(~isempty(tailcoord))
                body_contour_struct.x = [body_contour_struct.x; tailcoord(1)];
                body_contour_struct.y = [body_contour_struct.y; tailcoord(2)];
            end
        end
        
        clear('sorted_indicies');
    end
    clear('bodyend_indicies');
end

% interpolate between contour points
if(~isempty(body_contour_struct.x))
    x1 = [];
    y1 = [];
    for(j=2:length(body_contour_struct.x))
        x1 = [x1 body_contour_struct.x(j-1) (body_contour_struct.x(j-1) + body_contour_struct.x(j))/2 ];
        y1 = [y1 body_contour_struct.y(j-1) (body_contour_struct.y(j-1) + body_contour_struct.y(j))/2 ];
    end
    x1 = [x1 body_contour_struct.x(end)];
    y1 = [y1 body_contour_struct.y(end)];
    body_contour_struct.x = x1;
    body_contour_struct.y = y1;
end
% hard code to get 52 body contour points
[body_contour_struct.curvature, body_contour_struct.x, body_contour_struct.y]  = relative_curvature_angles_from_coordinates(body_contour_struct.x, body_contour_struct.y, 50);

% body_contour_struct.midbody = round(length(body_contour_struct.x)/2);

clear('body_contour');
clear('body_contour_coords');
clear('neighbor_matrix');

% figure(2);
% subplot(1,2,2);
% %imshow(bwmorph(bwmorph(bwmorph(imfill(Image,'holes'),'spur'),'majority'),'remove'));
% imshow(Image);
% %imshow(perimeter_image);
% hold on;
% plot(body_contour_struct.x, body_contour_struct.y,'.b')
% plot(body_contour_struct.x(1), body_contour_struct.y(1),'or')
% plot(body_contour_struct.x(end), body_contour_struct.y(end),'og')
% hold off

return;
end