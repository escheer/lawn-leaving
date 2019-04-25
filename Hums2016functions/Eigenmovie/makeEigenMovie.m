function makeEigenMovie(dataStructOrMatfile,movieFileName,opts)
%makeEigenMovie(dataStructOrMatfile,movieFileName,opts)
%
% Make an EigenMovie
%
% struct for input dataStructOrMatfile (shortened here to data):
% data.numSegs    (scalar)
% data.anglesFull (numFrames,numSegs)
% data.anglesCBW  (numFrames,numSegs)
% data.anglesRBW  (numFrames,numSegs)
% data.wormSegmentVector (numFrames,2,numSegs)
% data.wormMoviePath (numFrames,2)
% data.movieRange (integer x 2)
% data.movieStartTime (scalar)
% data.movieFPS (scalar)
% data.arrowTimes (scalar column array)
% data.palette       (either a string, or a Nx3 array) (redblue default)
% data.paletteRanges (2x3 array)  (each row is a range)
%
% movieFileName (string)
% 
% opts.outputDir (pwd)
% opts.imageSize (101)
% opts.figureSize ([1000 1000])
% opts.figureBackgroundColor ([1 1 1])
% opts.title (none)
% opts.palette (overrides data.palette)
% opts.outputStillPDFs (true)
%    
% how to call:
% makeEigenMovie(inputData,movieFileName,opts)
% 
% (c) Saul Kato, saul@kato.com
% Created 2014-3-27

% code dependencies:
%  computeRotationMatrix.m
%  fastsmooth.m
%  genWormImage.m
%   -> triangeFill.m
%  mtv.m
%  proplot.m
%  textul.m
%  textur.m
%
% third party code dependencies:
%  export_fig.m
%  VideoPlayer (VideoUtils 1.2.4 package)
%  subtightplot.m
%
% to output single frame PDFs, ghostscript needs to be installed (export_fig)


%parse first argument
if isstruct(dataStructOrMatfile)
    data=dataStructOrMatfile;
elseif ischar(dataStructOrMatfile)
    data=load(dataStructOrMatfile);
else
    disp('makeEigenMovie> First argument must be a .mat filename or a struct. quitting.');
    return;
end


if nargin<3
    opts=[];
end

if nargin<2
    movieFileName=[];
end

if ~isfield('opts','outputDir') || isempty(opts.outputDir)
    opts.outputDir=pwd;
end

if ~isfield('opts','imageSize') || isempty(opts.imageSize)
    opts.imageSize=101;
end

if ~isfield('opts','smoothOrientationFlag') || isempty(opts.smoothOrientationFlag)
    opts.smoothOrientationFlag=false;
end

if ~isfield('opts','title')
    opts.title=[];
end

if ~isfield('opts','figureSize')
    opts.figureSize=[1000 1000];
end

if ~isfield('opts','palette')  || isempty(opts.palette)
    if isfield('data','palette')
        opts.palette=data.palette;
    else
        opts.palette='redblue';
    end
end

if ~isfield('opts','figureBackgroundColor') || isempty(opts.figureBackgroundColor)
    opts.figureBackgroundColor=[1 1 1];
end

if ~isfield('opts','flipHorizontalMovie') || isempty(opts.flipHorizontalMovie)
    opts.flipHorizontalMovie=false;
end

if ~isfield('opts','outputMovieFrameRate') || isempty(opts.outputMovieFrameRate)
    opts.outputMovieFrameRate=10;
end

if ~isfield('opts','outputStillPDFs') || isempty(opts.outputStillPDFs)
    opts.outputStillPDFs=true;
end

if ~isfield('data','movieFPS')
    opts.movieFPS=10;
end

if ~isfield('data','movieFPS')
    opts.movieFPS=10;
end


%debug flags
debugDrawLines=false;


%preliminaries
halfmovie=floor(opts.imageSize-1)/2;
numFrames=size(data.anglesFull,1);
timeVec=mtv(1:numFrames,1/opts.movieFPS)+data.movieStartTime;
numJoints=size(data.anglesFull,2); %should be one less then data.numSegs
fixseg=ceil(data.numSegs/2);  %which body segment is fixed in the center of the images


%figure setup
figure('Position',[0 0 opts.figureSize(1) opts.figureSize(2)],'Color',opts.figureBackgroundColor);
whitebg(opts.figureBackgroundColor);
gap=[.03 .03];
marg_h=[0.1 .00];
marg_w=[0.04 0.02];


%s struct: static worm model properties
SWs.numJoints=numJoints;
SWs.scale=4; %3.8;
SWs.numSegs=numJoints+1;
SWs.segLen=2*ones(1,SWs.numSegs)/data.numSegs*9;
SWs.wth=[2 3 3.4 3.4 4 4 3.5 3.0 2.0 1.5 0]/6; %10 segments

SWs.wth=interp1(1:length(SWs.wth),SWs.wth,linspace(1,length(SWs.wth),SWs.numSegs),'linear');

SWs.fixSeg=fixseg;

%m struct: genWormImage properties
SWm.imageActual=zeros(opts.imageSize,opts.imageSize,'uint8');
SWm.scale=2;  %controls image resampling factor
SWm.width=opts.imageSize*2;
SWm.height=opts.imageSize*2;

SWm.xoffset=floor(SWm.width/2);
SWm.yoffset=floor(SWm.height/2);
SWm.interpolationFactor=2;  %controls segment uprez factor

%p struct: dynamic worm model properties
SWp.cx=0;
SWp.cy=0;

%fix and smooth angles

angs0=data.anglesFull;
angs1=data.anglesCBW;
angs2=data.anglesRBW;
    
if opts.smoothOrientationFlag && size(data.anglesFull,1)>2

    for j=1:size(angs0,2)
        anglesFull_sm(:,j)=fastsmooth(angs0(:,j),3,3,1);
        anglesCBW_sm(:,j)=fastsmooth(angs1(:,j),3,3,1);
        anglesRBW_sm(:,j)=fastsmooth(angs2(:,j),3,3,1);
    end

else
    
    anglesFull_sm=angs0;
    anglesCBW_sm=angs1;
    anglesRBW_sm=angs2;
    
end

%Real Worm Movie
if ~isempty(movieFileName)
    
    %Create Movie Reader object
    vp = VideoPlayer(movieFileName, 'Verbose', false, 'ShowTime', false,'InitialFrame',data.movieRange(1));
    %vp = VideoReader(movieFileName);
    
    subtightplot(5,4,[1 5],gap,marg_h,marg_w);
    handles.im(1)=image(zeros(opts.imageSize,opts.imageSize));
    axis image;axis off;
    set(gca,'YDir','Reverse');
    title('Actual Worm');
    handles.frameCountText=textul('1');
    handles.timeStampText=textur([num2str(data.movieStartTime) ' s']);
    
    if debugDrawLines
        %plot midpt to nose line
        handles.imline=line([0 10],[0 10],'Color','b');
    end
   
end

%Eigenworm images
subtightplot(5,4,[2 6],gap,marg_h,marg_w);
handles.wormpatch(1)=patch(10*[1 2 3 1 2],10*[2 4 2 4 5],'k');
set(handles.wormpatch(1),'FaceColor',[ 0 0 0]);
xlim([-50 50]);
ylim([-50 50]);
axis square;axis off;
set(gca,'YDir','normal');
title('Undulation Mode');

subtightplot(5,4,[3 7],gap,marg_h,marg_w);
handles.wormpatch(2)=patch(10*[1 2 3 1 2],10*[2 4 2 4 5],'k');
set(handles.wormpatch(2),'FaceColor',[ 0 0 0]);
xlim([-50 50]);
ylim([-50 50]);
axis square;axis off;
set(gca,'YDir','normal');
title('Turning Mode');

subtightplot(5,4,[4 8],gap,marg_h,marg_w);
handles.wormpatch(3)=patch(10*[1 2 3 1 2],10*[2 4 2 4 5],'k');
set(handles.wormpatch(3),'FaceColor',[0 0 0]);
xlim([-50 50]);
ylim([-50 50]);
axis square;
axis off;
set(gca,'YDir','normal');
title('Full Reconstruction');
hold on;

%KYMOGRAPHS
subtightplot(5,4,9:12,gap,marg_h,marg_w);
handles.im(5)=imagesc(timeVec,1:numJoints,data.anglesFull',data.paletteRanges(1,:));
colormap(opts.palette);
colorbar;
set(gca,'YTick',2:2:numJoints); 
set(gca,'XTick',timeVec(1):30:timeVec(end));
set(gca,'XTickLabel',[]);
set(gca,'XMinorTick','on');
ylim([0 numJoints]+0.5);
xlim([timeVec(1) timeVec(end)]);
ylabel('segment angle')
handles.line(1)=line([timeVec(1) timeVec(1)],[0 numJoints+1],'Color','r','LineWidth',3);
hold on;
plot([0 0],[0 numJoints]+0.5,'k--');
title('Body Posture');
proplot;

subtightplot(5,4,13:16,gap,marg_h,marg_w);
handles.im(6)=imagesc(timeVec,1:numJoints, data.anglesCBW',data.paletteRanges(2,:));
colorbar;
set(gca,'YTick',2:2:numJoints); 
set(gca,'XTick',timeVec(1):30:timeVec(end));
set(gca,'XTickLabel',[]);
set(gca,'XMinorTick','on');
ylim([0 numJoints]+0.5);
xlim([timeVec(1) timeVec(end)]);
ylabel('segment angle')
handles.line(2)=line([timeVec(1) timeVec(1)],[0 numJoints+1],'Color','r','LineWidth',3);
hold on;
plot([0 0],[0 numJoints]+0.5,'k--');
title('Undulation Mode');
proplot;

subtightplot(5,4,17:20,gap,marg_h,marg_w);
handles.im(7)=imagesc(timeVec,1:numJoints, data.anglesRBW',data.paletteRanges(3,:));
colorbar;
set(gca,'YTick',2:2:numJoints); 
set(gca,'XTick',10*floor(timeVec(1)/10):10:timeVec(end));
set(gca,'XMinorTick','on');
ylim([0 numJoints]+0.5);
xlim([timeVec(1) timeVec(end)]);
ylabel('segment angle')
handles.line(3)=line([timeVec(1) timeVec(1)],[0 numJoints+1],'Color','r','LineWidth',3);
hold on;
plot([0 0],[0 numJoints]+0.5,'k--');
title('Turning Mode');
xlabel('time (s)');
proplot;

%Annotation Arrows
for i=1:length(data.arrowTimes)
    arrowFrames(i)=(data.arrowTimes(i)-data.movieStartTime)*data.movieFPS;
    annotation('arrow',.89*([arrowFrames(i) arrowFrames(i)])/(data.movieRange(2)-data.movieRange(1))+.04,[.66 .63]);
end

%Figure Title
if ~isempty(opts.title)
    annotation('textbox','Position',[.1 .97 .4 .02],'EdgeColor','none','String',strrep(opts.title,'_','\_'));
end

%Compute corrected nose angle of real worm movie
for i=1:numFrames
    
        svx=squeeze(data.wormSegmentVector(i,1,:));
        svy=squeeze(data.wormSegmentVector(i,2,:));

        segLen=sqrt(svx.^2+svy.^2);

        ws.bodyPos=cumsum(segLen);
        ws.bodyPos=[0; ws.bodyPos];

        ws.sx=cumsum(svx);
        ws.sy=cumsum(svy);

        ws.sx=[0; ws.sx];
        ws.sy=-[0; ws.sy];
         
        ws.sx=ws.sx-mean(ws.sx);
        ws.sy=ws.sy-mean(ws.sy);

        videoNoseAngle(i)=atan2( -ws.sx(1), ws.sy(1));
        
        ws.sxLast=ws.sx;
        ws.syLast=ws.sy;

end

%Compute smoothed rotation angles
for i=1:numFrames
    
        SWs2=SWs;
       %SWp.e=anglesDecompTr(ds.tr_rel).loadings(i,:);    
        SWp.angles=anglesCBW_sm(i,:)'; 
        SWp2=SWp;
        [cs.sx,cs.sy]=computeJointPositionsFromAngles(SWp2.angles,SWs2.segLen,SWs2.fixSeg,SWs2.numJoints);
        globalAngleE12(i)=(atan2(-(cs.sx(1)-mean(cs.sx)),cs.sy(1)-mean(cs.sy) )+videoNoseAngle(i));
      
        SWp.angles=anglesRBW_sm(i,:)';    
        SWp2=SWp;
        [cs.sx,cs.sy]=computeJointPositionsFromAngles(SWp2.angles,SWs2.segLen,SWs2.fixSeg,SWs2.numJoints);
        globalAngleE39(i)=(atan2(-(cs.sx(1)-mean(cs.sx)),cs.sy(1)-mean(cs.sy) )+videoNoseAngle(i));
        
        SWp.angles=anglesFull_sm(i,:)';    
        SWp2=SWp;
        [cs.sx,cs.sy]=computeJointPositionsFromAngles(SWp2.angles,SWs2.segLen,SWs2.fixSeg,SWs2.numJoints);  
        globalAngleE19(i)=(atan2(-(cs.sx(1)-mean(cs.sx)),cs.sy(1)-mean(cs.sy) )+videoNoseAngle(i));
        
end

gAsmooth12=globalAngleE12;
gAsmooth19=globalAngleE19;
gAsmooth39=globalAngleE39;

subFrameRange=[1 numFrames];

%Set up Output Movie
if ischar(dataStructOrMatfile)
    movieOutName=dataStructOrMatfile(1:end-4);
else
    movieOutName='EigenMovie';
end
videoOutObj=VideoWriter(movieOutName,'MPEG-4');
videoOutObj.FrameRate=opts.outputMovieFrameRate;
videoOutObj.Quality=10;
open(videoOutObj);
disp('makeEigenMovie> video file opened for writing.');

%%%%%%%%%%%%%%%%%%
%Make Output Movie
for i=subFrameRange(1):subFrameRange(end)
    
        %video update
        
        if ~isempty(movieFileName)
            
            VidFrame=vp.getFrameAtNum(i+data.movieRange(1)-1);
             
            cptrx = uint16(max([1 round(data.wormMoviePath(i,1))-halfmovie]):min([round(data.wormMoviePath(i,1)+halfmovie) vp.Width]));
            cptry = uint16(max([1 round(data.wormMoviePath(i,2))-halfmovie]):min([round(data.wormMoviePath(i,2)+halfmovie) vp.Height]));
            m_oneframe = VidFrame(cptry,cptrx,:);
            
            if opts.flipHorizontalMovie
                m_oneframe=m_oneframe(:,end:-1:1,:);
                
            end
            set(handles.im(1),'CData',m_oneframe);

            %plot centroid to nose
            svx=squeeze(data.wormSegmentVector(i,1,1:end));
            svy=squeeze(data.wormSegmentVector(i,2,1:end));

            segLen=sqrt(svx.^2+svy.^2);

            ws.bodyPos=cumsum(segLen);
            ws.bodyPos=[0; ws.bodyPos];

            ws.sx=cumsum(svx);
            ws.sy=cumsum(svy);

            ws.sx=[0; ws.sx];
            ws.sy=-[0; ws.sy];

            ws.sx=ws.sx-mean(ws.sx);
            ws.sy=ws.sy-mean(ws.sy);

            set(handles.frameCountText,'String',num2str(i));

            set(handles.timeStampText,'String',[num2str(timeVec(i-subFrameRange(1)+1)) ' s']);

            if debugDrawLines
                 set(handles.imline,'XData',51-[0 50*sin(videoNoseAngle(i))]);
                 set(handles.imline,'YData',51-[0 50*cos(videoNoseAngle(i))]);
            end

            ws.sxLast=ws.sx;
            ws.syLast=ws.sy;


        end
        
        %red lines
        for j=1:3
           set(handles.line(j),'XData',[timeVec(i-subFrameRange(1)+1) timeVec(i-subFrameRange(1)+1)]);
        end
        
        %synthworms
        SWs2=SWs;
   
        SWp.angles=data.anglesCBW(i,:); 
        SWp2=SWp;
        SWp2.angles=SWp.angles';  
        SWp2.rm=computeRotationMatrix(gAsmooth12(i));
        [~,outlinex,outliney]=genWormImage(SWp2,SWs2,SWm);
        set(handles.wormpatch(1),'XData',-outlinex);
        set(handles.wormpatch(1),'YData',outliney);

        SWp.angles=data.anglesRBW(i,:);    
        SWp2=SWp;
        SWp2.angles=SWp.angles';
        SWp2.rm=computeRotationMatrix(gAsmooth39(i));
        [~,outlinex,outliney]=genWormImage(SWp2,SWs2,SWm);
        set(handles.wormpatch(2),'XData',-outlinex);
        set(handles.wormpatch(2),'YData',outliney);
 
        SWp.angles=data.anglesFull(i,:);    
        SWp2=SWp;
        SWp2.angles=SWp.angles';   
        SWp2.rm=computeRotationMatrix(gAsmooth19(i));
        [~,outlinex,outliney]=genWormImage(SWp2,SWs2,SWm);  
        set(handles.wormpatch(3),'XData',-outlinex);
        set(handles.wormpatch(3),'YData',outliney);

        drawnow;     
        
        framestruct=im2frame(png_cdata(gcf),jet(256));
        writeVideo(videoOutObj,framestruct.cdata); %VideoUtils call
        
        if mod(i,100)==1 disp(['makeEigenMovie> ' num2str(i) ' frames completed.']); end
               
        if opts.outputStillPDFs && ismember(i-1,arrowFrames)
             export_fig([movieOutName '-frame' num2str(i,'%4.4d') '.pdf'],'-nocrop','-painters');
        end
end

close(videoOutObj);

end %MAIN

function cdata = png_cdata(hfig)
    % Get CDATA from hardcopy using opengl
    % Need to have PaperPositionMode be auto 
    orig_mode = get(hfig, 'PaperPositionMode');
    set(hfig, 'PaperPositionMode', 'auto');
    warning('off','MATLAB:hardcopy:DeprecatedHardcopyFunction')
    cdata = hardcopy(hfig, '-Dopengl', '-r0');
    warning('on','MATLAB:hardcopy:DeprecatedHardcopyFunction')
    % Restore figure to original state
    set(hfig, 'PaperPositionMode', orig_mode);
end

function [sx,sy]=computeJointPositionsFromAngles(angles,segLen,fixSeg,numJoints)

     %n+1 segments defined by n internal body angles determining n+1 body
     %points.
     sx=zeros(numJoints+2,1);
     sy=zeros(numJoints+2,1);
          
     sx(fixSeg+1)=0;  sy(fixSeg+1)=0;  %fix a middle segment in space
     sx(fixSeg)=0;  sy(fixSeg)=segLen(fixSeg);

     %compute head segments
     angcum=0;     
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
