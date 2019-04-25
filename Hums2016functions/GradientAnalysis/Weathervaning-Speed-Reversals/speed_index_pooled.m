% use on runinfo files
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2016

clear
home=cd;

%plot?
plotting=1;
SI=NaN(55,2);
if plotting==1
fig=figure;
end
CC=1;

for GT =1:2 % if you have 2 folders with control and gradient data
    tic
    if GT==1
        cd('control')
        
    elseif GT==2
        cd('gradient')
    end
files =dir('*runinfo*');

%parameters:
first=1;
last=179;
binsize=89;

% for each experiment: put all run data into one vector:
mB= NaN(2,2);
mT= NaN(2,2);
mS= NaN(2,2);
c=0;

for batch=1:length(files)
    
    load(files(batch).name);
    disp(files(batch).name);
   

for F=1:length(bearing)
    
    if ~isempty(bearing{F})
    
    c=c+1;
    bearingAll=[];
    dBearingAll=[];
    SpeedAll=[];

    
    for i=1:size(bearing{F},1)
        if ~isempty(dBearing{F}{i,1})
            for rr=1:size(bearing{F}(i,:),2)
                bv=bearing{F}{i,rr};
                tr=dBearing{F}{i,rr};
                sp=speed_ctr{F}{i,rr};
               
                bearingAll=cat(2,bearingAll,bv);
                dBearingAll=cat(2,dBearingAll,tr);
                SpeedAll=cat(2,SpeedAll,sp);

            end
        end
    end
    
    
    bv=round(bearingAll*10)/10;
    tr=round(dBearingAll*10)/10;
    sp=round(SpeedAll*100)/100;

    %kill extremely high turning rates (short reversals etc)
    
    bi=find(tr>19 | tr<-19);
    bi2= find(bv>last| bv<first);
    bi=[bi,bi2];
    tr(bi)=NaN;
    bv(bi)=NaN;
    sp(bi)=NaN;
    
if ~isempty(bv)
    
    [X,v]=hist(bv,10);
    X1(c,:)=X./sum(X);
    cc=1;
    %bin -->for what?
    
    for i=first:binsize:last
        bin_idx= bv>i & bv<=i+binsize;
        pb1=bv(bin_idx);
        pt=tr(bin_idx);
        st=sp(bin_idx);

        mB(c,cc)=nanmean(pb1);
        mT(c,cc)=nanmean(pt);
        mS(c,cc)=nanmean(st);

        cc=cc+1;
    end
    
end
    end
end

end % end files loop

save mS mS
SI(1:length(mS),CC)=-1+(mS(:,2)./mS(:,1));
nd=(cd);
d= strfind(cd, '\');
name=nd(d(end):end);
CC=CC+1;
cd(home)

end

legend('control','gradient ')

if ~exist ('plots', 'dir')
mkdir('plots')
end
%% boxplots:
if plotting==1
figure(fig)

boxplot(-1+SI,'whisker',0.9529,'labels',{'nlp-12 ctrl','nlp-12 gradient'})

ylabel('relative speed change (mm/min)')
cd('plots')
% saveas(gca, 'weathervaning_pooled_olddata.fig')
 saveas(gca, 'speed modulation boxplot.fig')


cd(home)

end


%% (2) bearing  binned by turning rate


