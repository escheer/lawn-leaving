% function [ forward, reverse, sorted_posture_angle, sorted_spline ] = getforwardreverse_phase( posture_angle, head, tail, spline )
%GETFORWARDREVERSE_PHASE.M This function looks at the posture angles of
%animals from one frame to the next to see if the phase is advancing toward
%the head or towards the tail

spline = allTracks_slim(1).spline;
posture_angle = allTracks_slim(1).posture_angle;
head = allTracks_slim(1).head;
tail = allTracks_slim(1).tail;

sorted_spline = spline;
sorted_posture_angle = posture_angle;
phase = NaN(size(posture_angle,2),1);

for i = 1:size(posture_angle,2) %loop over frames in the track
    disp(i);
%     if i == 223
%         disp('debug');
%     end
    h = head(i,:);
    t = tail(i,:);
    sp = spline{i};
    pa = posture_angle{i};
    if sum(isnan(pa))==length(pa) || sum(sum(isnan(sp)))==numel(sp) || sum(isnan(h))==2 || sum(isnan(t))==2 %if missing data, skip loop
        continue;
    end
    % first rearrange the posture_angle cell array based on the predetermined
    % head and tail position
    pt1 = sp(1,:); pt2 = sp(end,:);
    d = pdist2([pt1; pt2],[h; t]);
    [~,closer] = min(d);
    if isequal(closer,[2 1]) %this means the spline is sorted tail to head and must be reversed
        disp('FLIP!');
        sp = flipud(sp);
        pa = fliplr(pa);
        sorted_spline{i} = sp;
        sorted_posture_angle{i} = pa;
    end
    % fit a single sine wave of the form: f(x) =  a1*sin(b1*x+c1), where c1 is
    % the phase, which will be used to figure out which way the body wave is
    % propagating. this requires the curve fitting toolbox but there are other
    % ways to do implement this without this toolbox later on if need be.
    f = fit( (1:length(pa))',pa', 'sin1');
    phase(i) = f.c1;
    
    figure(1); clf; hold on;
    plot(1:length(pa),pa); 
    t = 1:length(pa);
    a = f.a1*sin(f.b1*t+f.c1);
    plot(t,a);
    [xaf,yaf] = ds2nfu([t(3)+5 t(3)],[a(3)+0.5 a(3)]); 
    annotation('textarrow',xaf,yaf,'String',['phase: ' num2str(phase(i))]);
    pause();

end



% end

