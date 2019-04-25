function [LineSkeleton, RecursionDepth] = Skel2Line(SkelImg,RecursionDepth,CurrentRecursionLimit)
 % SkelImg:                  Input Image that contains a skeleton
 % RecursionDepth:           How many times did the function called itself, start with 0
 % CurrentRecursionLimit:    Max number of recursions
 %
 %In some instances no trimming occurs
 %
 %(c) Manuel Zimmer, manuel.zimmer@imp.ac.at
 %
 % Requires MatLAb Image Processing Toolbox

 
    LineSkeleton = SkelImg;

    Endpoints = bwmorph(SkelImg, 'endpoints');

    Branches = bwmorph(SkelImg,'branchpoints');
    
    %StartingBranches = Branches;
    
    NumTotalBranches= sum(sum(Branches));

    if NumTotalBranches > 0

        while sum(sum(Branches))==NumTotalBranches
            
            OldEndpoints = bwmorph(SkelImg, 'endpoints');
            
            OldBranches = bwmorph(SkelImg,'branchpoints');


            Endpoints(OldEndpoints)=1;

            SkelImg(Endpoints) = 0;

            Branches = bwmorph(SkelImg,'branchpoints');

        end

    BranchSegments = bwlabel(Endpoints | OldBranches);

    BranchSegmentsStats = regionprops(BranchSegments,'Area');

    [~, SmallestSegmentIdx] = max([BranchSegmentsStats.Area]);

    LineSkeleton(BranchSegments==SmallestSegmentIdx)=0;
    
    LineSkeleton(OldBranches) = 1;
    
    LineSkeleton = bwmorph(LineSkeleton, 'skel', 'inf'); % necessary for removing corners, which are detected as branchpoints




if RecursionDepth < CurrentRecursionLimit %recursive function call when RecursionLimit is not reached (there should be a more standard way to do that
    
    RecursionDepth=RecursionDepth +1;

    
    [LineSkeleton, RecursionDepth] = Skel2Line(LineSkeleton,RecursionDepth,CurrentRecursionLimit);
    
%     disp('Current recursion Depth = ');
%     disp(RecursionDepth);
    
else
    
    disp('RecursionLimit was reached !');
    
end

    end

end