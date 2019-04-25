function midline = getmidline( wormcrop, framenum )
%getmidline.m This function computes the midline from a cropped binarized
%image of the worm.
pad_img = padarray(wormcrop, [2 2], 0); %pad image to make sure you get the whole outline
% BW = edge(pad_img,'canny'); %get outline points by edge detection
pad_img = imdilate(pad_img,strel('disk',1));
pad_img = imfill(pad_img,'holes');
BW = bwperim(pad_img);
[edge_y, edge_x] = find(BW); %get these points
edge_x = edge_x-2; edge_y = edge_y-2; %remove padding
D = squareform(pdist([edge_x edge_y])); %get pairwise distance matrix between all points

order = 1; %random choice for seed
sorted_points = [edge_x(order(end)) edge_y(order(end)) order(end)];%third entry is the index
next_pt = [];
reached_end = 0;
counter = 1;

startplot = 0;
if framenum >= 529
    disp('debug');
    startplot = 1;
end

while ~reached_end
    %     disp(counter);
    
    %     if counter == 216 && framenum == 45
    %         disp('debug');
    %     end
    %
    if counter==1
        curr_pt = sorted_points(1,:);
    else
        curr_pt = next_pt;
    end
    [sortedD,sortind] = sort(D(curr_pt(3),:)); %sort other points' distance from the reference point (ascending)
    sortind = sortind(2:end); %don't include self
    sortedD = sortedD(2:end);
    while ismember(sortind(1),order) %find closest ind not already in order (and if the distance to this ind is less than back to the start)
        sortind = sortind(2:end); %pop that element
        sortedD = sortedD(2:end);
        if length(sortind)==1 %|| sortedD(1) > 10 %D(order(1),sortind(1))
            reached_end = 1;
            break;
        end
    end
    if ~reached_end %even though we have already reached the end, don't add the point in
        order = [order; sortind(1)];
        next_pt = [edge_x(sortind(1)) edge_y(sortind(1)) sortind(1)];
        sorted_points = [sorted_points; next_pt];
        counter = counter+1;
    end
end
%make two circularizations of this outline so you don't run into boundary
%errors with curvature peak finding
start1 = floor(length(sorted_points)/4); start2 = floor(length(sorted_points)*(1/2));
sorted1 = [sorted_points(start1:end,:);(sorted_points(1:start1-1,:))];
sorted2 = [sorted_points(start2:end,:);(sorted_points(1:start2-1,:))];
%apply polynomial filtering to this outline shape
smth1 = sgolayfilt(sorted1,3,19);
smth2 = sgolayfilt(sorted2,3,19);

strike1 = 0;
strike2 = 0;

Vertices=smth1;
Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
k_1=LineCurvature2D(Vertices,Lines); %gets the radius of the tangent circle at every point
kfilt_1 = sgolayfilt(k_1,3,7); %do a little more filtering
warning('off','signal:findpeaks:largeMinPeakHeight');
[pks_up_1,locs_up_1,w_up_1,p_up_1] = findpeaks(kfilt_1,'MinPeakHeight',0.2,'MinPeakDistance',83,'SortStr','descend','NPeaks',2);
[pks_dn_1,locs_dn_1,w_dn_1,p_dn_1] = findpeaks(-1*kfilt_1,'MinPeakHeight',0.2,'MinPeakDistance',83,'SortStr','descend','NPeaks',2);
locs_1 = [];
% get rid of peaks that are too far apart (there is no maxpeakdistance
% argument to findpeaks)
if length(locs_up_1)==2
    if abs(locs_up_1(2)-locs_up_1(1))>140
        locs_up_1 = [];
    end
end
if length(locs_dn_1)==2
    if abs(locs_dn_1(2)-locs_dn_1(1))>140
        locs_dn_1 = [];
    end
end
if length(locs_up_1)==2 && length(locs_dn_1)~=2
    locs_1 = locs_up_1;
elseif length(locs_dn_1)==2 && length(locs_up_1)~=2
    locs_1 = locs_dn_1;
else
    warning(['SMOOTH_1: ' num2str(length(locs_up_1)) 'up peaks found; ' num2str(length(locs_dn_1)) 'dn peaks found']);
    strike1 = 1;
end
if startplot
    figure(1); clf; hold on;
    plot(kfilt_1);
    for i = 1:length(locs_1) 
        scatter(locs_1(i),kfilt_1(locs_1(i)),10,'r+');
    end
    figure(2); clf; hold on;
    N=LineNormals2D(Vertices,Lines);
    plot([Vertices(:,1) Vertices(:,1)+kfilt_1.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+kfilt_1.*N(:,2)]','g');
    plot(smth1(:,1), smth1(:,2),'b');
    for i = 1:length(locs_1)
        scatter(Vertices(locs_1(i),1),Vertices(locs_1(i),2),10,'r*');
    end
    scatter(Vertices(1,1),Vertices(1,2),'bs');
    scatter(Vertices(end,1),Vertices(end,2),'ms');
    axis equal;
end

Vertices=smth2;
Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
kfilt_2=LineCurvature2D(Vertices,Lines); %gets the radius of the tangent circle at every point
kfilt_2 = sgolayfilt(kfilt_2,3,7); %do a little more filtering
warning('off','signal:findpeaks:largeMinPeakHeight');
[pks_up_2,locs_up_2,w_up_2,p_up_2] = findpeaks(kfilt_2,'MinPeakHeight',0.2,'MinPeakDistance',83,'SortStr','descend','NPeaks',2);
[pks_dn_2,locs_dn_2,w_dn_2,p_dn_2] = findpeaks(-1*kfilt_2,'MinPeakHeight',0.2,'MinPeakDistance',83,'SortStr','descend','NPeaks',2);
locs_2 = [];
% get rid of peaks that are too far apart (there is no maxpeakdistance
% argument to findpeaks)
if length(locs_up_2)==2
    if abs(locs_up_2(2)-locs_up_2(1))>140
        locs_up_2 = [];
    end
end
if length(locs_dn_2)==2
    if abs(locs_dn_2(2)-locs_dn_2(1))>140
        locs_dn_2 = [];
    end
end
if length(locs_up_2)==2 && length(locs_dn_2)~=2
    locs_2 = locs_up_2;
elseif length(locs_dn_2)==2 && length(locs_up_2)~=2
    locs_2 = locs_dn_2;
else
    warning(['SMOOTH_2: ' num2str(length(locs_up_2)) 'up peaks found; ' num2str(length(locs_dn_2)) 'dn peaks found']);
    strike2 = 1;
end

if startplot
    figure(3); clf; hold on;
    plot(kfilt_2);
    for i = 1:length(locs_2)
        scatter(locs_2(i),kfilt_2(locs_2(i)),10,'r+');
    end
    figure(4); clf; hold on;
    N=LineNormals2D(Vertices,Lines);
    plot([Vertices(:,1) Vertices(:,1)+kfilt_2.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+kfilt_2.*N(:,2)]','g');
    plot(smth2(:,1), smth2(:,2),'b');
    for i = 1:length(locs_2)
        scatter(Vertices(locs_2(i),1),Vertices(locs_2(i),2),10,'r*');
    end
    scatter(Vertices(1,1),Vertices(1,2),'bs');
    scatter(Vertices(end,1),Vertices(end,2),'ms');
    axis equal;
%     pause();
end

if strike1 && strike2
    warning('-----------------NEITHER OPENING COULD FIND ENPOINTS!------------------');
end





% [pks,locs,w,p] = findpeaks(abs(kfilt),'MinPeakHeight',0.25,'MinPeakWidth',3,'SortStr','descend','NPeaks',2)%,'MinPeakProminence',0.4)'MinPeakHeight',0.3
% warning('off','signal:findpeaks:largeMinPeakHeight');
% [pks_up,locs_up,w_up,p_up] = findpeaks(kfilt,'MinPeakHeight',0.25,'MinPeakWidth',3,'SortStr','descend','NPeaks',2);
% [pks_dn,locs_dn,w_dn,p_dn] = findpeaks(-1*kfilt,'MinPeakHeight',0.25,'MinPeakWidth',3,'SortStr','descend','NPeaks',2);
% locs = [];
% if length(locs_up)==2 && length(locs_dn)~=2
%     locs = locs_up;
% elseif length(locs_dn)==2 && length(locs_up)~=2
%     locs = locs_dn;
% else
%     warning([num2str(length(locs_up)) 'up peaks found; ' num2str(length(locs_dn)) 'dn peaks found']);
% end

% if startplot
%     figure(1); clf; hold on;
%     plot(k);
%     for i = 1:length(locs)
%         scatter(locs(i),k(locs(i)),5,'r+');
%     end
%     %     for i = 1:length(locs_up)
%     %         scatter(locs_up(i),k(locs_up(i)),5,'r+');
%     %     end
%     %     for i = 1:length(locs_dn)
%     %         scatter(locs_dn(i),k(locs_dn(i)),5,'m+');
%     %     end
%     %     figure(2); clf; plot(abs(kfilt)); hold on;
%     %     plt = abs(kfilt);
%     %     for i = 1:length(locs)
%     %         scatter(locs(i),plt(locs(i)),'r+');
%     %     end
%     figure(3); clf; hold on;
%     N=LineNormals2D(Vertices,Lines);
%     plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
%     % plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
%     plot(smooth_outline(:,1), smooth_outline(:,2),'b');
%     
%     for i = 1:length(locs)
%         scatter(Vertices(locs(i),1),Vertices(locs(i),2),7,'r*');
%     end
% %         for i = 1:length(locs_up)
% %             scatter(Vertices(locs_up(i),1),Vertices(locs_up(i),2),5,'r*');
% %         end
% %         for i = 1:length(locs_dn)
% %             scatter(Vertices(locs_dn(i),1),Vertices(locs_dn(i),2),5,'m*');
% %         end
%     scatter(Vertices(1,1),Vertices(1,2),'bs');
%     scatter(Vertices(end,1),Vertices(end,2),'ms');
%     axis equal;
%     pause();
% end


% [~, sort_k_ind] = sort(k,'ascend');
% endpt1 = smooth_outline(sort_k_ind(1),:); endpt2 = smooth_outline(sort_k_ind(2),:);%the lowest 2 curvatures should correspond to head and tail
% %get each side of the animal
% D = squareform(pdist(smooth_outline)); %get pairwise distance matrix between all points
% side1 = endpt1; %seed
% curr_pt = endpt1;
% order1 = sort_k_ind(1);
% while ~isequal(curr_pt,endpt2)
%     curr_ind = order1(end); %the last ones you added in
%     %     curr_pt = side1(end);
%     [~,sortind] = sort(D(curr_ind,:)); %find closest points in order
%     sortind = sortind(2:end);
%     while ismember(sortind(1),order1) %find closest ind not already in order
%         sortind = sortind(2:end); %pop that element
%     end
%     curr_pt = smooth_outline(sortind(1),:);
%     if ~isequal(curr_pt,endpt2)
%         side1 = [side1; curr_pt];
%         order1 = [order1; sortind(1)]; %append a point to side 1
%     end
% end
%
% side2 = endpt2; %seed
% curr_pt = endpt2;
% order2 = sort_k_ind(2);
% counter = 1;
% while length(side1)+length(side2)<length(smooth_outline)%abs(order2(end)-sort_k_ind(1))>1 %as long as you are more than one index away from the start of side 1's index
%     %     disp(counter);
%     %     if counter == 5
%     %         disp('debug');
%     %     end
%
%     curr_ind = order2(end); %the last ones you added in
%     [~,sortind] = sort(D(curr_ind,:)); %find closest points in order
%     sortind = sortind(2:end);
%     while ismember(sortind(1),order2) || ismember(sortind(1),order1) %find closest ind not already in order (or the other side)
%         sortind = sortind(2:end); %pop that element
%     end
%     curr_pt = smooth_outline(sortind(1),:);
%     if ~isequal(curr_pt,endpt1)
%         side2 = [side2; curr_pt];
%         order2 = [order2; sortind(1)]; %append a point to side 2
%     end
%     counter = counter+1;
% end
%
% % %the other side is all of the other points
% % side2 = setdiff(smooth_outline,side1,'rows','stable'); %this is in order but starts in the middle for some reason
%
% figure(1); clf; imshow(BW); hold on;
% scatter(side1(:,1),side1(:,2),'g*');
% scatter(side2(:,1),side2(:,2),'r+');
% scatter(endpt1(1),endpt1(2),10,'b*');
% scatter(endpt2(1),endpt2(2),10,'y*');
% figure(2); plot(k);
% pause();

midline = [];
end

