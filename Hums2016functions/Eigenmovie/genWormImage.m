function [imgOut,outlinex,outliney]=genWormImage(p,s,opts)
% [imgOut,outlinex,outliney]=genWormImage(p,s,opts)
% 
% generate synthetic worm renderings from angle data
%
% (c) Saul Kato, saul@kato.com
% Created 2014-3-27

colr=255;

if nargin<3
    opts.width=200;
    opts.height=200;
    opts.interpolationFactor=1;
end

if ~isfield(opts,'interpolationFactor')
    opts.interpolationFactor=1;
end

imgOut=zeros(opts.width,opts.height);

numEdgePoints=s.numSegs*opts.interpolationFactor;


p.angles=1.15*p.angles;

%generate centerline
[ds.sx_nonInterp,ds.sy_nonInterp]=computeJointPositionsFromAngles(p.angles,s.segLen,s.fixSeg,s.numJoints);

bodyPos=(0:s.segLen:s.segLen*(s.numJoints+1))';

%interpolate centerline
numNewJoints=(s.numJoints+1)*opts.interpolationFactor-1;

dshi=spline(bodyPos,[ds.sx_nonInterp'; ds.sy_nonInterp'],linspace(0,bodyPos(end),numNewJoints+2));
ds.sx=dshi(1,:)';
ds.sy=dshi(2,:)';

%interpolate widths
widthVecInterp=interp1([0 0.5+(0:s.numSegs-1) s.numSegs],[0 s.wth 0],linspace(1/(opts.interpolationFactor+1),s.numSegs-1/(opts.interpolationFactor+1),s.numSegs*opts.interpolationFactor),'PCHIP');

%increase fixSeg to make commensurate with new number of segments
fixSegInterp=s.fixSeg*opts.interpolationFactor;

%generate dorsal and ventral edges
[ds.sxD,ds.syD,ds.sxV,ds.syV]=GWIcomputeEdgePositionsFromSegments(ds.sx,ds.sy,fixSegInterp,widthVecInterp,numNewJoints);

%apply rotation matrix
ds.sDRot=(p.rm*[ds.sxD ds.syD  ]')';
ds.sVRot=(p.rm*[ds.sxV ds.syV  ]')';
ds.sRot =(p.rm*[ds.sx  ds.sy ]')';

%merge into single interleaved point list and add head and tail points
pointlistx=zeros(1,2+2*numEdgePoints);
pointlisty=zeros(1,2+2*numEdgePoints);

pointlistx(1)=p.cx+s.scale*ds.sRot(1,1);
pointlistx(2:2:2*numEdgePoints)=p.cx+s.scale*ds.sDRot(:,1)';
pointlistx(3:2:2*numEdgePoints+1)=p.cx+s.scale*ds.sVRot(:,1)';
pointlistx(2*numEdgePoints+2)=p.cx+s.scale*ds.sRot(end,1);

pointlisty(1)=p.cy+s.scale*ds.sRot(1,2);
pointlisty(2:2:2*numEdgePoints)=p.cy+s.scale*ds.sDRot(:,2)';
pointlisty(3:2:2*numEdgePoints+1)=p.cy+s.scale*ds.sVRot(:,2)';
pointlisty(2*numEdgePoints+2)=p.cy+s.scale*ds.sRot(end,2);

pointlistx=opts.scale*pointlistx;
pointlisty=opts.scale*pointlisty;

trilist=zeros(numEdgePoints*2,3);
for i=1:(numEdgePoints*2)
    trilist(i,:)=[i i+1 i+2];
end

%rasterize triangles [currently into global img]
imgOut=triangleFill(imgOut,pointlistx(trilist(:,1)),pointlisty(trilist(:,1)),...
       pointlistx(trilist(:,2)),pointlisty(trilist(:,2)),...
       pointlistx(trilist(:,3)),pointlisty(trilist(:,3)),colr);

%outlines for patch drawing
outlinex=p.cx+s.scale*[ds.sRot(1,1); ds.sDRot(:,1); ds.sRot(end,1); ds.sVRot(end:-1:1,1)];
outliney=p.cy+s.scale*[ds.sRot(1,2); ds.sDRot(:,2); ds.sRot(end,2); ds.sVRot(end:-1:1,2)];

end %main function

function [sx,sy]=computeJointPositionsFromAngles(angles,segLen,fixSeg,numJoints)

     sx=zeros(numJoints+2,1);
     sy=zeros(numJoints+2,1);
          
     sx(fixSeg+1)=0;  sy(fixSeg+1)=0;  %fix segment 5 in space
     sx(fixSeg)=0;  sy(fixSeg)=segLen(fixSeg);
     
     angcum=0;
     
     %compute head segments
     for i=(fixSeg-1):-1:1
         angcum=angcum+angles(i); 
         sx(i)=sx(i+1)+segLen(i)*sin(angcum);
         sy(i)=sy(i+1)+segLen(i)*cos(angcum);
     end
     
     %compute tail segments
     angcum=0;
     for i=(fixSeg+2):numJoints+2
         angcum=angcum+angles(i-2); 
         sx(i)=sx(i-1)+segLen(i-1)*sin(angcum);
         sy(i)=sy(i-1)-segLen(i-1)*cos(angcum);
     end
end

function [sxD,syD,sxV,syV]=GWIcomputeEdgePositionsFromSegments(sx,sy,fixSeg,widthVec,numJoints)

     sxD=zeros(numJoints+1,1);
     syD=zeros(numJoints+1,1);
     sxV=zeros(numJoints+1,1);
     syV=zeros(numJoints+1,1);
     
     segx=zeros(numJoints+1,1);
     segy=zeros(numJoints+1,1);
     
     segPerp=zeros(numJoints+1,2);
     
     fullvec=1:numJoints+1;

     for j=fullvec
         
         segx(j)=sx(j+1)-sx(j);
         segy(j)=sy(j+1)-sy(j);
         
         segPerp(j,:)=[-segy(j) segx(j)];
         segPerp(j,:)=segPerp(j,:)/norm(segPerp(j,:));
         
         segMidptx(j)=(sx(j)+sx(j+1))/2;
         segMidpty(j)=(sy(j)+sy(j+1))/2;
         
         sxD(j)=segMidptx(j) - widthVec(j)*segPerp(j,1);
         syD(j)=segMidpty(j) - widthVec(j)*segPerp(j,2);
         sxV(j)=segMidptx(j) + widthVec(j)*segPerp(j,1);
         syV(j)=segMidpty(j) + widthVec(j)*segPerp(j,2);
         
     end
     
      sxD(fixSeg)=-widthVec(fixSeg);
      sxV(fixSeg)=+widthVec(fixSeg);
      syD(fixSeg)=segMidpty(fixSeg);
      syV(fixSeg)=segMidpty(fixSeg);
end

