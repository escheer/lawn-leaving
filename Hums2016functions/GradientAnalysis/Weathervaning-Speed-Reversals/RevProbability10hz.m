% Needs Tracks_als files as input
%calculates bearing dO2/dt and rev onset timepoints for each run and stores them in revinfo file
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2016


clear
warning off
% global rois
% cd('mat')
rate=10;
rate=round((rate*10/3)/10);
home=cd;
nd=(cd);
d= strfind(cd, '\');
name=nd(d(2):d(end));

files=dir('*JR_als.mat');
nd=(cd);
d= strfind(cd, '\');
name=nd(d(2):d(end));


reversals={NaN};

% create gradient matrix
gradient=zeros(1060,2160);
gl=21:-17/(2100):4;
low=4;
high=21;

EC=1;

for i=1:size(gradient,1)
    gradient(i,1:30)=high;
    gradient(i,31:2131)=gl;
    gradient(i,2132:end)=low;
end

% cd ..\
for F=1:length(files)
    
    F
    clearvars -except F exp  home name files bearing dBearing speed_ctr...
        HeadSens HeadSpeed HeadSpeedNet reversals rate gradient sensPath
    cc=0;
    fname=(files(F,1).name)
    load (files(F).name);
    Reversals={NaN};
    dCdT={NaN};
    
    pixToMM=0.0155;
    
    
    %% ---analyze Tracks----
    
    for T= 1:length(Tracks)
        % only tracks which start after establishment of the gradient
        % (ca frame 500)
        if  Tracks(1,T).Frames(end)>3000 && Tracks(1,T).Frames(1)>1000
            
            TX=  (Tracks(1,T).SmoothX);
            TY=  (Tracks(1,T).SmoothY);
            
            
            %----exclude those which occur close to border 
            %(needs to be adapted to reolution of behavior movie if needed):----
            bi1= find( TX<600 | TX>1800);
            bi2= find( TY<50 | TY>980);
            bi= setdiff(bi1, bi2);
            bi=[bi,bi2];
            TX(bi)=NaN;
            TY(bi)=NaN;
            
            %kill runs which have unusual high angular speed
            %and don't cover 1 small worm travel distance (20 px):
            
            meanAngVelocity=(nanmean(abs(Tracks(1,T).AngSpeed)));
            if meanAngVelocity>30
                continue
            else
                d=NaN(1,1);
                ci=1;
                for iii=1:3:length(TY)
                    d(ci)=sqrt(((TX(1)-TX(iii)).^2)+((TY(1)-TY(iii)).^2));
                    ci=ci+1;
                end
                if max(d)<20
                    disp('no displacement')
                    continue
                end
            end
            
            %get O2 concentration at this position of path:
            C_ind=sub2ind(size(gradient),TY(1:rate:end),TX(1:rate:end));
            nan_idx=isnan(C_ind);
            C_ind=C_ind(nan_idx~=1);
            
            %get O2 concentration at this position of path and reinsert Nans:
            sensoryPathNaN=(gradient(round(C_ind)));
            dCdT=NaN(1,length(nan_idx));
            dCdT(nan_idx~=1)=sensoryPathNaN;
            sensP{T}=dCdT;
            
            % reversal vector:
            Revs=zeros(1,ceil(length(TX)/rate));
            Rev=(Tracks(1,T).polishedReversals);
            
            if ~isempty(Rev) & max(d)>=20;% only 1 sec long reversals
                gi=find(Rev(:,4)>0);
                Rev=Rev(gi,:);
                %                 gi=find(Rev(:,2)-Rev(:,1)>(rate*2));
                %                 Rev=Rev(gi);
                Rev=Rev(:,1);
                Rev=round(Rev/rate);
                Revs(Rev)=1;
                
            end
            Reversals{T}=Revs;
            
            
            
            
            %--- get all bearing and heading angles for each left over run relative to
            %%  preferred 02 isocline:------
            if   max(d)>=20
                
                x=TX(1:rate:end);
                y=TY(1:rate:end);
                s=Tracks(1,T).Speed;
                s=s(1:rate:end);
                L=length(x);
                bearingangles=NaN(1,1);
                headingangles=NaN(1,1);
                beelineX=NaN;
                
                
                if F<1
                    CO2=2500;
                else
                    CO2=10;
                end
                
                for ii=1:L-1
                    beelineX=(CO2-(x(ii)));
                    beelineY=0;%(roi(2,2)+100-(x(ii)));  % alternatively bearing towards gas inlet
                    headingX=(x(ii+1))-(x(ii));
                    headingY=(y(ii+1))-(y(ii));
                    nenner=((beelineX*headingX)+(beelineY*headingY));
                    Brecher=(sqrt((beelineX).^2+(beelineY).^2))*(sqrt((headingX).^2+(headingY).^2));
                    cosinus_angle=nenner/Brecher;
                    bearingangles(ii)=(acos(cosinus_angle))*(360/(2*pi));
                    headingangles(ii) = atan2(headingY,headingX);
                end
                
                cc=cc+1;
                
                Ctr_Bearing{T}=[bearingangles NaN ];
                Ctr_dBearing{T}=[diff(bearingangles) NaN NaN ];
                Ctr_Speed{T}=s;
                
            end 
            
        end
        
    end %end Tracks loop
    %%
    
    bearing{F}=Ctr_Bearing;
    dBearing{F}=Ctr_dBearing;
    speed_ctr{F}=Ctr_Speed;
    reversals{F}=Reversals;
    sensPath{F}=sensP;
    
end
%%

home=cd;
nd=(cd);
d= strfind(cd, '\');
name=nd(d(end)+1:end);

display('...save')

save(['revinfo_10hz_JR' fname(1:end-42)  name] , 'bearing' ,'dBearing' ,'speed_ctr','reversals','rate','sensPath');



