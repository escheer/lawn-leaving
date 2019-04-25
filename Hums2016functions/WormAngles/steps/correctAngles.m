function WormGeomTracksC = correctAngles(anglesFileOrFolder,tracks,numIP,~)

%Divides skeleton into segments of equal arc length.
%Revision: Corrects WormTraces (skeletons) by fitting a smoothed spline with same number of points;
%Can be interpolated or downsampled with spline to adjust angle number across frames and worms
%
%saves out a file called XXX-corrected.mat where XXX is the input file name
%containing a corrected WormGeomTracksC. 
%
%Will run on a single .mat angles file or a folder of angles.mat files.
%If specified with no arguments, it will run on the current directory.
%
%runs on .mat files specified by filestring.
%
%examples:
%correctAngles     runs on current directory
%correctAngles('/.../')  runs on current directory
%correctAngles('/.../xxxanglesV8.mat') runs on a file
%
%if a second argument is provided, it will run only on one or more tracks
%for debugging purposes.
%correctAngles([],83);

%(c) Saul Kato, saul@kato.com
%Created 12202013
%Revised 2015 by Ingrid Hums, ingrid.hums@gmail.com
%
% third party code dependencies: smoothn.m


filestring = '*anglesV8.mat';

if nargin<3 || isempty(numIP)
    numIP = 26; % defines how many points the interpolated skeleton should have
end   

%%% how much to smooth the wormtrace with smoothn relative to trace length;
%%% the lower the value, the stronger the smoothing;
%%% can be any positive value
smoothpart = 8; % 8              
                    
if nargin<1 || isempty(anglesFileOrFolder)
    anglesFileOrFolder=pwd;
end

if isdir(anglesFileOrFolder)
    
    disp(['A directory was provided.  correctAngles will try to run on every ' filestring ' file in the folder:']);
    files = dir([anglesFileOrFolder filesep filestring]);
    dirPrefix = [anglesFileOrFolder filesep];
else
    [dirPrefix, files.name, ext] = fileparts(anglesFileOrFolder);
    
    if isempty(dirPrefix) 
        dirPrefix='.'; 
    end
    
    files.name = [files.name ext];
    dirPrefix = [dirPrefix filesep];
end

for f = 1:length(files)
    tic
    clear WormGeomTracks WormGeomTracksC
    WormGeomTracks = {NaN};
    if exist([files(f).name(1:end-5) '8-corrected.mat'])
        display (['already exists...' files(f).name(1:end-5) '8-corrected.mat'])

    else     
    
        disp(['...reparameterizing: ' files(f).name]);
        load([dirPrefix files(f).name]);  %creates WormGeomTracks

        WormGeomTracksC = struct('WormAnglesNewRAD', repmat({NaN(1,numIP-2)},1,length(WormGeomTracks)),...
            'WormSegmentVectorsNew', repmat({NaN(1,2,numIP-1)},1,length(WormGeomTracks)),'WormPointsNew', repmat({{}},1,length(WormGeomTracks)) );


        if ~exist('WormGeomTracks','var')
            disp(['No WormGeomTracks in ' files(f).name ', skipping.']);
        else
        
% %%%--reparameterize skeleton--%%%

            if nargin<2
                tracks = 1:length(WormGeomTracks);
            end
              
            for tr = tracks

                if mod(tr,10)==0
                    display(tr)
                end

                wormlength = NaN(size(WormGeomTracks(tr).WormTraces,2),1);

                WormGeomTracksC(tr).WormSegmentVectorsNew=[];
                WormGeomTracksC(tr).WormAnglesNewRAD=[];

                for fr = 1:size(WormGeomTracks(tr).WormTraces,2);                

                     if ~isempty(WormGeomTracks(1,tr).WormTraces{1,fr}) 

                        X = WormGeomTracks(1,tr).WormTraces{1,fr}(:,1);
                        Y = WormGeomTracks(1,tr).WormTraces{1,fr}(:,2);

                         % fit a smooth spline (cardioid) to the worm trace
                        smoothfactor = length(X)/smoothpart;
                        Z = smoothn({X,Y},smoothfactor); % 'robust' gives slightly smoother spline but chops off more from ends

                        smoothwormX = Z{1};
                        smoothwormY = Z{2};

                        WormDot2X = smoothwormX;
                        WormDot2Y = smoothwormY;   

                        svx = diff(WormDot2X);
                        svy = diff(WormDot2Y);

                        segLen = sqrt(svx.^2+svy.^2);
                        wormlength(fr) = sum(segLen);

                        bodyPos = [0; cumsum(segLen)];

                        zeroSegs = abs(svx)+abs(svy)==0;

                        if all(~isnan(svx)) && all(~zeroSegs)

                            wshi = spline(bodyPos,[WormDot2X'; WormDot2Y'],linspace(0,bodyPos(end),numIP));        

                            WormGeomTracksC(tr).WormPointsNew{fr} = single(wshi');

                            ws_hi.sx=wshi(1,:);
                            ws_hi.sy=wshi(2,:);

                        else
                            ws_hi.sx=zeros(1,numIP);
                            ws_hi.sy=zeros(1,numIP);
                            WormGeomTracksC(tr).WormPointsNew{fr} = {};
                        end

                    else
                            ws_hi.sx=zeros(1,numIP);
                            ws_hi.sy=zeros(1,numIP);
                            WormGeomTracksC(tr).WormPointsNew{fr} = {};

                    end 

                        WormJointPositionsCorr(fr,1,:) = ws_hi.sx;
                        WormJointPositionsCorr(fr,2,:) = ws_hi.sy;


                        WormGeomTracksC(tr).WormSegmentVectorsNew(fr,1,:) = WormJointPositionsCorr(fr,1,2:end)-WormJointPositionsCorr(fr,1,1:end-1);
                        WormGeomTracksC(tr).WormSegmentVectorsNew(fr,2,:) = WormJointPositionsCorr(fr,2,2:end)-WormJointPositionsCorr(fr,2,1:end-1);


                end %fr

                    %compute new joint Angles & convert angles to RAD

                        for fr = 1:size(WormGeomTracks(tr).WormTraces,2)         
                            WormGeomTracksC(tr).WormAnglesNewRAD(fr,:) = asin(computeSinAnglesFromJointPositions(WormJointPositionsCorr(fr,1,:),WormJointPositionsCorr(fr,2,:)));
                        end


                %%% kick out frames where worm skeleton is too short = artifacts
                wormlength = wormlength/nanmean(wormlength);
                tooshort = find(wormlength<=0.85); %0.885 

                WormGeomTracksC(tr).WormSegmentVectorsNew(tooshort,:,:) = NaN;
                WormGeomTracksC(tr).WormAnglesNewRAD(tooshort,:) = NaN;

                for i = 1:length(tooshort)
                    WormGeomTracksC(tr).WormPointsNew{tooshort(i)} = {};
                end

            end % end tracks loop
        
            WormGeomTracksC(1).Status = 'corrected';
            save([dirPrefix filesep files(f).name(1:end-5) '8-smoothn' num2str(smoothpart)  '-corrected.mat'],'WormGeomTracksC','-v6');
        end
    end
   
    toc

end %f

disp('correctAngles done.');

end %main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sinAngles=computeSinAnglesFromJointPositions(xpos,ypos)

    % calculates angles in unit SINE from worm points
    % e.g. if there are 11 worm points, there are 10 body segments, and there are 9 internal angles

    numSegs = size(xpos,3)-1; 
    sinAngles=zeros(numSegs-1,1); 
    xdel=zeros(numSegs,1); 
    ydel=zeros(numSegs,1);

    % segment vectors calculated from worm points
    for i=1:(numSegs)
        xdel(i)=xpos(i+1)-xpos(i);
        ydel(i)=ypos(i+1)-ypos(i);
    end

    for i=1:(numSegs-1)
        sinAngles(i)=-(xdel(i)*ydel(i+1)-xdel(i+1)*ydel(i))/norm([xdel(i) ydel(i)])/norm([xdel(i+1) ydel(i+1)]);
    end

end





    