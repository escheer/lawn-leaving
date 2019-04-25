function pixpermm = getpixpermm( frame, diam_mm )
%GETPIXPERMM.m Extracts the conversion factor between pixels and
%millimeters -- helpful later on for computing the speed of animals, etc.

bw = edge(frame,'sobel',0.035);%was 0.07
bw = bwareaopen(bw,10);
h1 = figure(1); ax1 = axes();
imshow(bw); hold on;
[well_edge_y, well_edge_x] = find(bw);
H = imfreehand;
Position = wait(H);
pos_x = Position(:,1); pos_y = Position(:,2);
in = inpolygon(well_edge_x,well_edge_y,pos_x,pos_y);
inside_x = well_edge_x(in);
inside_y = well_edge_y(in);
K = convhull(inside_x,inside_y);
well_x = inside_x(K); well_y = inside_y(K);
[~,~,Rfit,~] = circfit(well_x,well_y);%Rfit is the fit radius
pixpermm = Rfit*2/diam_mm; %wells are 10 mm - answer should be around 111.6 pix/mm
close all;
end

