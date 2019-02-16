function newts = intdefine(varargin)
%INTDEFINE interpolate a timeseries to a user-defined timestep.
%   NEWTS = intdefine(T1,V1,Step,Cut) interpolates a timeseries to a 
%   user-defined timestep and user-defined time-gap where
%   interpolated data will be calculated when the difference in time 
%   is greater than this value.  Cut & Step should be expressed in julian
%   date format (e.g., 1 second = 1/24/60/60 = 0.0001574).
%
%   NEWTS is a 2-column matrix, containing the new timestep in the
%   first column, and the Y1 interpolation in the second.
%
%  modified from ALIGNTS.m(Copyright 2010 C. Grant Law)

t1 = varargin{1};
v1 = varargin{2};
Step = varargin{3};
Cut = varargin{4};
tStart = varargin{5};
tEnd = varargin{6};

diff1 = t1(2:end)-t1(1:end-1); diff1(diff1==0)=[];

xi = (tStart:Step:tEnd);

newts = zeros(length(xi),2);
newts(:,1) = xi;
newts(:,2) = interp1(t1,v1,xi,'linear');
%newts(:,2) = interp1_nan(t1,v1,xi,'linear');

% Identify and clean out data gaps from new timeseries
gaps1 = find(diff1 > Cut);

for i = 1:length(gaps1)
    bad = newts(:,1)>t1(gaps1(i))&newts(:,1)<t1(gaps1(i)+1);
    newts(bad,2) = nan;
end
