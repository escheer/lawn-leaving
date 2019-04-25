function timevec=mtv(trace,dt)
%timevec=mtv(trace,dt)
%
% convenience function
% make time vector from t=0, dt steps, to trace length)
%
%(c) Saul Kato 10/3/12, saul@kato.com

if nargin<2
    dt=0.05;
end

timevec=dt*(0:length(trace)-1);