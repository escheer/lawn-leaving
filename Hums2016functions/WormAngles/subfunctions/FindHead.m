function  [headX, headY,First,Last,CurrentTrack,flips,cc]=FindHead(CurrentTrack,CurrentFrame,Skeleton,prevLast, prevFirst ,Xlast,Ylast, HXvect,flips,cc,breakflag,flipvalue)
%finds head position based on direction of movement and reversal state:
%FindSkeleton assigns first point of skeleton to NE (in
%movie plot upper left corner
%
%(c) Julia Riedl
%Created 2015

cc=cc+1;

%First check if the current frame is a reversal or omega assigned assigned frame:
rFlag=[];

if ~isempty(CurrentTrack.polishedReversals)
RI = (CurrentTrack.polishedReversals(:,4)==0);
CurrentTrack.polishedReversals(RI,:)=NaN;
    for rr=1:size(CurrentTrack.polishedReversals,1)
        ReversalState =sum( CurrentFrame >= CurrentTrack.polishedReversals(rr,1)  & CurrentFrame <= CurrentTrack.polishedReversals(rr,2)) ;
        if ReversalState==1
            rFlag=rr;
            break
        end
    end
    
%if reversal overlaps with end of omega, delete it, probably wrongly assigned

if ~isempty(CurrentTrack.OmegaTrans)

    for rr=1:size(CurrentTrack.OmegaTrans,1)
       OmegaState =sum( CurrentFrame >= CurrentTrack.OmegaTrans(rr,2)-10  & CurrentFrame <= CurrentTrack.OmegaTrans(rr,2)) ;
    end
else
 OmegaState = 0;
end

if OmegaState==1 & ReversalState==1
    CurrentTrack.polishedReversals(rFlag,1:4)=NaN;
    ReversalState=0;
end

else
    
    ReversalState = 0;
end


if ~isempty(Skeleton{1,CurrentFrame})
    First=0;
    Last=0;

    
%%first get the direction of locomotion:

CurrentDirection = round(CurrentTrack.Direction(CurrentFrame));

    if abs(CurrentDirection) <= 90
        DirectionTagNS ='N';
    else
        DirectionTagNS ='S';
    end

    if CurrentDirection >= 0
        DirectionTagEW ='E';
    else
        DirectionTagEW ='W';
    end

Direction=DirectionTagEW;



%then check whether first element in smoothed worm trace is aligned with
%direction of locomotion

WormMajorAxisVecX = Skeleton{1,CurrentFrame}(end,1) - Skeleton{1,CurrentFrame}(1,1) ;

WormMajorAxisVecY = Skeleton{1,CurrentFrame}(end,2) - Skeleton{1,CurrentFrame}(1,2) ;

    if WormMajorAxisVecX >= 0
        HeadTailTagEW = 'W';
        % always west, first point is always on left = west
    else
        HeadTailTagEW = 'E';
        % never east
    end

    if WormMajorAxisVecY >= 0
        HeadTailTagNS = 'N';
    else

        HeadTailTagNS = 'S';
    end

H_T_Direction= HeadTailTagEW;


%Matlab assigns first point in skeleton to W, if worm goes to east(
%Direction=E) while H-T Direction is W, we have to swap and take the last
%point of the skeleton as head.
%This is true only during forward movement, during reversals we don't swap,
%because naturally here the head-tail direction direction inverse to the
%heading direction.
%also if speed goes below reversal threshold keep head assigned as
%previously as worm is pausing or jittering

    if CurrentFrame>1 & CurrentFrame<length(Skeleton)
        Speed=mean(CurrentTrack.Speed(CurrentFrame-1:CurrentFrame+1));
    else
        Speed=mean(CurrentTrack.Speed(CurrentFrame:CurrentFrame));
    end

    if strcmp(HeadTailTagEW, DirectionTagEW)==0 & ReversalState==0
    % same as: strcmp(DirectionTagEW,'E')==1 & strcmp(HeadTailTagEW,'W')==1

        headX=Skeleton{1,CurrentFrame}(end,1);
        headY=Skeleton{1,CurrentFrame}(end,2);
        Last=1;
        First=0;

    elseif strcmp(DirectionTagEW,'W')==1 & strcmp(HeadTailTagEW,'W')==1 & ReversalState==1

        headX=Skeleton{1,CurrentFrame}(end,1);
        headY=Skeleton{1,CurrentFrame}(end,2);
        Last=1;
        First=0;


    else

        headX=Skeleton{1,CurrentFrame}(1,1);
        headY=Skeleton{1,CurrentFrame}(1,2);
        First=1;
        Last=0;

    end



    if CurrentFrame>5
        if Speed<0.055 & CurrentFrame>1 & sum(isnan(HXvect(CurrentFrame-5:CurrentFrame-1)))<3%& ReversalState==0

            if prevLast==1;
                headX=Skeleton{1,CurrentFrame}(end,1);
                headY=Skeleton{1,CurrentFrame}(end,2);
                Last=1;
                First=0;

            elseif prevLast==0;
                headX=Skeleton{1,CurrentFrame}(1,1);
                headY=Skeleton{1,CurrentFrame}(1,2);
                First=1;
                Last=0;

            end

        end

    
        if breakflag<1

    % finally if the distance between new and previous head position is too
    % large, also keep old one, or take the one of the skeleton which is
    % closer (accounts for heading angles close to 0 or 180, were the skeleton
    % flips also the heading direction doesn't).

    d=sqrt(((headX-Xlast).^2)+((headY-Ylast).^2));

     %check for distance of previous head position to alternative ending of
     %skeleton, but only if there was no longer stretch of NaN (e.g. caused by an omega) just before
    Xalt=[];
    Yalt=[];

            if First==1 & Last==0 & sum(isnan(HXvect(CurrentFrame-5:CurrentFrame-1)))<3 %3 for 15fps
                    Xalt=Skeleton{1,CurrentFrame}(end,1);
                    Yalt=Skeleton{1,CurrentFrame}(end,2);
                    pFlag=2;

            elseif First==0 & Last==1 & sum(isnan(HXvect(CurrentFrame-5:CurrentFrame-1)))<3 %3
                    Xalt=Skeleton{1,CurrentFrame}(1,1);
                    Yalt=Skeleton{1,CurrentFrame}(1,2);
                    pFlag=1;
%             else
%                 pFlag=0;
            end


    d1=sqrt(((Xalt-Xlast).^2)+((Yalt-Ylast).^2));

            if d1<d & pFlag==1

                headX=Skeleton{1,CurrentFrame}(1,1);
                headY=Skeleton{1,CurrentFrame}(1,2);
                First=1;
                Last=0;
                flips=flips+1;

            elseif d1<d & pFlag==2

                headX=Skeleton{1,CurrentFrame}(end,1);
                headY=Skeleton{1,CurrentFrame}(end,2);
                Last=1;
                First=0;
                flips=flips+1;
            else
                flips=0;
            end



     %if previous was NaN, take last position of skeleton    
            if CurrentFrame>5
               if isnan(d)& prevLast==1 & sum(isnan(HXvect(CurrentFrame-5:CurrentFrame-1)))<4 %4 for 15fps
                   headX=Skeleton{1,CurrentFrame}(end,1);
                   headY=Skeleton{1,CurrentFrame}(end,2);
                   First=0;
                   Last=1;

               elseif isnan(d)& prevLast==0 & sum(isnan(HXvect(CurrentFrame-5:CurrentFrame-1)))<4 %4
                   headX=Skeleton{1,CurrentFrame}(1,1);
                   headY=Skeleton{1,CurrentFrame}(1,2);
                   First=1;
                   Last=0;

               end
            end
        end
    end
 
 % go back to start of long flip episode to correct:
 if flips > flipvalue
cc=cc-flips;

 end

    

else
    
    headX=NaN;
    headY=NaN;
    Last=prevLast;
    First=prevFirst;
    Direction='.';
    H_T_Direction='.';
    ReversalState=0;
end



