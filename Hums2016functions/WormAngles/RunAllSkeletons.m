%%% runs through all 4 codes to get from WormImages to angle time series
%%% per track
% third party code dependencies: smoothn


 clear

%%% (1)
FindAllWormSkeletonsFromTracks;
%%% Calculates the skeletons from binary worm images from tracks_als files.

%%% (2)
correctAngles;
%%% Smoothes skeleton, divides skeleton into segments of equal arc length
%%% and calculates angles in rad.
%%% Set the number of wormpoints = numIP (number of angles + 2) (line 34);
%%% Optionally, edit smoothpart (line 40), if you want to smooth the
%%% skeleton less or more.

%%% (3)
HeadAssign;
%%% Finds head coordinates of worm

%%%(4)
skel_order;
%%% Orders skeleton from head to tail
%%% and fixes left-right flips. 
%%% Optionally, ajdust thresholds ("threshes" in line 7) in this script for correcting left-right
%%% flips.