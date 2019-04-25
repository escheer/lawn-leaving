function img=triangleFill(img,p1x,p1y,p2x,p2y,p3x,p3y,colr)
%img=triangleFill(img,p1x,p1y,p2x,p2y,p3x,p3y,colr)
%
% rasterize a list of flat-colored triangles and blit it onto an image
% triangles specified as vertex coordinate arrays
%
%Saul Kato, saul@kato.com, 2014-03-31
%adapted from http://www.gdchaos.net/sites/default/files/Primitive2D.cs__0.txt

[xoffset,yoffset]=size(img);
xoffset=int16(floor(xoffset/2));
yoffset=int16(floor(yoffset/2));

for i=1:size(p1x,2)
    
    p1.X=p1x(i);
    p1.Y=p1y(i);

    p2.X=p2x(i);
    p2.Y=p2y(i);

    p3.X=p3x(i);
    p3.Y=p3y(i);

    % Order points by swapping them around
    if (p2.Y > p1.Y)  tp = p1; p1 = p2; p2 = tp; end;
    if (p3.Y > p1.Y)  tp = p1; p1 = p3; p3 = tp; end;
    if (p3.Y > p2.Y)  tp = p2; p2 = p3; p3 = tp; end;

    % Number of steps for each edge
    steps13 = int16(ceil(p1.Y - p3.Y));
    steps12 = int16(ceil(p1.Y - p2.Y));
    steps23 = int16(ceil(p2.Y - p3.Y));

    sx13 = double(p1.X - p3.X) / (double(steps13));
    sx12 = double(p1.X - p2.X) / (double(steps12));
    sx23 = double(p2.X - p3.X) / (double(steps23));

    % Draw lower part of triangle
    x13 = double(p1.X);
    x12 = double(p1.X);
    dx = double(1);   

    for j=0:steps12-1

        % Decrease by step
        x13 = x13 - sx13;           
        x12 = x12 - sx12;
        % Calculate new distance between points of the two edges
        dx=x13 - x12;   
        % Draw horizontal line 
        yl=int16(p1.Y) - j + yoffset;
        if (yl <= 2*yoffset-1) && (yl >=1) %sauls clipping
            if (dx>0) 
                img(max([1 int16(x12) + xoffset]):min([2*xoffset int16(x12)+int16(abs(dx)) + xoffset]),yl)= colr;
            else
                img(max([1 int16(x13) + xoffset]):min([2*xoffset int16(x13)+int16(abs(dx)) + xoffset]),yl)= colr;
            end
        end
    end

    % Draw upper part of triangle
    % x13 should still hold a midpoint value from last loop so we need to 
    % connect this with points of x23 starting at point 2 (p2)
    x23 = double(p2.X);
    for j=0:steps23-1

        % Decrease both by step as we move towards top corner
        x13=x13 - sx13;
        x23=x23 - sx23;
        % Length of horizontal line to draw
        dx = x13 - x23;
        yl= int16(p2.Y) - j + yoffset;

        if (yl <= 2*yoffset-1) && (yl >=1)  %saul's clipping
            if (dx > 0) 
               img(max([1 int16(x23) + xoffset]):min([2*xoffset int16(x23)+int16(abs(dx)) + xoffset]),yl)=colr;
            else
               img(max([1 int16(x13) + xoffset]):min([2*xoffset int16(x13)+int16(abs(dx)) + xoffset]),yl)=colr;
            end
        end
    end

end
