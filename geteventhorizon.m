function [ ev_ho_x, ev_ho_y ] = geteventhorizon(background,old_ev_ho_x,old_ev_ho_y)
%GETEVENTHORIZON.M this function takes in the background image (not
%blurred) and extracts the event horizon (lawn-leaving boundary)
blurbg = imgaussfilt(background,6);

th = 0.0007;
% th = 0.0005;

BW = edge(blurbg,'sobel',th,'nothinning');
BW2 = bwareaopen(BW,5000); %get rid of schmutz
% BWnobord = imclearborder(BW2, 8);
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
se = strel('disk',3);
BWdil = imdilate(BWnobord, [se90 se0]);
BWdil = imdilate(BWdil, [se90 se0]);
BWdil = imdilate(BWdil,se);
BWdil = imdilate(BWdil,se);
BWdfill = imfill(BWdil, 'holes');
BWerode = imerode(BWdfill,se);
BWerode = imerode(BWerode,se);
BWerode = imerode(BWerode,[se90 se0]);
BWerode = imerode(BWerode,[se90 se0]);
se = strel('disk',7);
BWerode = imerode(BWerode,se);%restore from Gaussian filter
BWfinal = bwareaopen(BWerode,5000); %get rid of any schmutz created by the image erosion

%find the largest connected component, this should be the lawn
CC = bwconncomp(BWfinal);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);

if isempty(idx) && sum(isnan(old_ev_ho_x))==0 && sum(isnan(old_ev_ho_y))==0 %if there is a problem finding a connected component, use the event horizon from the last video (this is a hack for 11/20)
    ev_ho_x = old_ev_ho_x;
    ev_ho_y = old_ev_ho_y;
    return;
elseif isempty(idx) && sum(isnan(old_ev_ho_x))>0 && sum(isnan(old_ev_ho_y))>0
    figure(2); ax2 = axes();
    imshow(background); hold on;
    h = imfreehand(ax2); wait(h);
    pos = h.getPosition();
    h_poly = impoly(gca,pos);
    h_poly.setVerticesDraggable(false);
    ev_ho = getPosition(h_poly);
    plot(ev_ho(:,1),ev_ho(:,2),'LineWidth',5);
    pause();
    close all;
    ev_ho_x = ev_ho(:,1);
    ev_ho_y = ev_ho(:,2);
    return;
end

BWbig = zeros(size(BWfinal));
BWbig(CC.PixelIdxList{idx})=1;


BWbig = bwareaopen(BWbig,5000);
BWbig = imfill(BWbig, 'holes');

boundaries = bwboundaries(BWbig);

if size(boundaries,1)==1
    firstBoundary = boundaries{1};
    x = firstBoundary(:, 2);
    y = firstBoundary(:, 1);

else
    error('there must be a connected component in order to get the event horizon!');
end
polynomialOrder = 3;
windowWidth = 101;
smoothX = sgolayfilt(x, polynomialOrder, windowWidth);
smoothY = sgolayfilt(y, polynomialOrder, windowWidth);
smoothX = smoothX(10:end-10,:);%get rid of weird stuff at the ends
smoothY = smoothY(10:end-10,:);

% figure();
% imshow(background); hold on;
% plot(smoothX,smoothY,'b');
% pause(); close

ev_ho_x = smoothX;
ev_ho_y = smoothY;

% good = input('EVENT HORIZON LOOKS GOOD?');
% if ~good %if not, choose manually
%         clf;
%         ax2 = axes();
%         imshow(background); hold on;
%         ev_ho = imfreehand(ax2); wait(ev_ho);
%         ev_ho_line = getVertices(ev_ho);
%         ev_ho_x = ev_ho_line(:,1);
%         ev_ho_y = ev_ho_line(:,2);
%         plot(ev_ho_line(:,1),ev_ho_line(:,2),'LineWidth',5);
%         pause();
% else
%     ev_ho_x = smoothX;
%     ev_ho_y = smoothY;
% end



% K = convhull(smoothX,smoothY);
% ev_ho_x = smoothX(K);
% ev_ho_y = smoothY(K);

% imshow(background); hold on
% imshow(background); hold on 
% plot(smoothX,smoothY,'b');
% plot(ev_ho_x,ev_ho_y,'g');

end