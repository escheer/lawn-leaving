%use on runinfo files
%plots binned turning bias for bearing bin 40-140
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2016
%
%third party code dependencies: errorb



clear
warning off;

%plot?
plot=1;
files =dir('*runinfo*');

%parameters:
first=40;
last=140;
binsize=100;


% for each experiment: put all run data into one vector:
mP= NaN(35,8);
mT= NaN(35,8);
mS= NaN(35,8);
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
    
    bi=find(tr>19 | tr<-19 );
    bi2= find(bv>last| bv<first);
    bi=[bi,bi2];
    tr(bi)=NaN;
    bv(bi)=NaN;
    sp(bi)=NaN;
    tnorm=tr./(sp*67);
    bi=tnorm>10 | tnorm<-10;
    tnorm(bi)=NaN;
     tnorm=tnorm/0.8552;  % 0.8552= conversion from deg/pix-->rad/mm 1/(67*(360/(2*pi)))
    
if ~isempty(bv)
    
    [X,v]=hist(bv,10);
    X1(c,:)=X./sum(X);
    cc=1;
    %bin for bearing value
    
    for i=first :binsize:(last-binsize)
        bin_idx= bv>i & bv<=i+binsize;
        pb1=bv(bin_idx);
        pt=tnorm(bin_idx);
        st=sp(bin_idx);

        mP(c,cc)=nanmean(pb1);
        mT(c,cc)=nanmean(pt);
        mS(c,cc)=nanmean(st);

        cc=cc+1;
    end
    
end
    end
end

end % end files loop



%% plot:(1) turning rate binned by bearing
nd=(cd);
d= strfind(cd, '\');
name=nd(d(end-1):end);
if plot==1

figure
sem=nanstd(mT*3,1)/sqrt(c);
hold on
bar(nanmean(mT*3,1));
errorb(nanmean(mT*3,1),sem)
set(gca,'XTick',[1:length(mP)])
set(gca,'XTickLabel',round(nanmedian(mP,1)));
xlabel('bearing')
ylabel('turning rate (rad/mm)')
ylim([-0.2 0.10])
title ([name ': speed normalized'])


end

save turningbias mT




