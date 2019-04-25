% funtion caltulates heading and Bearing of XY trajectories to a reference
% point
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2015

function [ bearingangles, dBearing1, curvature,L]= getbearing_d15(X,Y,Xref,Yref,Tracks,T)

hx=X;
hy=Y;
if nargin==6
speed=Tracks(1,T).Speed(~isnan(X));
end
L=length(hx);
beelineX=NaN(L,1);
beelineY=NaN(L,1);
headingX=NaN(L,1);
headingY=NaN(L,1);
nenner=NaN(L,1);
Brecher=NaN(L,1);
cosinus_angle=NaN(L,1);
bearingangles=NaN(L,1);
headingangles=NaN(L,1);
dBearing1=NaN(L,1);

if L>15
    for ii=1:L-15
        beelineX(ii)=(Xref-(hx(ii)));
        if Yref==0;
            beelineY(ii)=0;
        else
        beelineY(ii)=(Yref-(hy(ii)));
        end
        headingX(ii)=(hx(ii+15))-(hx(ii));
        headingY(ii)=(hy(ii+15))-(hy(ii));
        nenner(ii)=((beelineX(ii)*headingX(ii))+(beelineY(ii)*headingY(ii)));
        Brecher(ii)=(sqrt((beelineX(ii)).^2+(beelineY(ii)).^2))*(sqrt((headingX(ii)).^2+(headingY(ii)).^2));
        cosinus_angle(ii)=nenner(ii)/Brecher(ii);
        bearingangles(ii+7,1)=(acos(cosinus_angle(ii)))*(360/(2*pi));
        headingangles(ii+7) = atan2(headingY(ii),headingX(ii));
    
    end
    
    dBearing1=vertcat(diff(bearingangles) , NaN);
    if nargin==6
    Speed1=speed(1:end-3);
    end
    curvature=(diff(headingangles));
    bi=find(curvature>0.9|curvature<-0.9);% account for jumps at 3-3 rad border
    curvature(bi)=NaN;
    
else
    dBearing1=NaN;
    bearing1=NaN;
    curvature=NaN;
    Speed1=NaN;
    
end
