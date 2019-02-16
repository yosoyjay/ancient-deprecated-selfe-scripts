function [int_kflow, int_kwind, int_vO2] = int_time_NEM(kflow, kwind, vO2)
% [k_flow, k_wind, v02] = int_time_NEM(k_flow, k_wind, vO2)
%
% This function integrates the instantanous k_flow, k_wind, vO2 output by 
% air_water_diffusion.m to hourly values as defined in Needoba et. al 2012.
%
% Assume:
%	Model output is every 900 seconds.  Lots of things depend on that.
%	
% Input:
%	kflow - Kflow values output by air_water_diffusion.m (day,ts,np)
%	kwind - Kwind values output by air_water_diffusion.m (day,ts,np)
%	vO2 - vO2 values output by air_water_diffusion.m (day,ts,np)
%
% Output:
%	int_kflow - Kflow values integrated over an hour for each node (day, node, hours) (cm/h)
%	int_kwind - Kwind values integrated over an hour for each node (cm/h)???
%	vO2 - vO2 values integrated over an hour for each node (cm/h)???
%
% lopezj -3/18/2012
%
% Consts used for dimensioning and time step assumption
nDays = size(kflow,1);
nTS = size(kflow,2)/4;      % Assume 96 output time steps per day		
nNodes = size(kflow,3);
oTS = 900; 					% Model output time step in seconds

% Allocate output arrays
int_kflow = nan(nDays,nTS,nNodes);
int_kwind = nan(nDays,nTS,nNodes);
int_vO2 = nan(nDays,nTS,nNodes);

% Loop over days and hours 
for day = 1:nDays
%tic
	for hour = 1:nTS
		% Take instantanous values from model output time step [kflow cm/s]
		% Multiply by 900 seconds to integrate over output time step
		% Sum over the 4 time steps to make an hour [kflow cm/h]
		% Result: int_flow(days, hours, np)
		int_kflow(day,hour,:) = kflow(day,hour*4,:).*oTS     +...		% ts = 4 - 60 mins
						        kflow(day,(hour*4)-1,:).*oTS +...		% ts = 3 - 45 mins
								kflow(day,(hour*4)-2,:).*oTS +...		% ts = 2 - 30 mins
								kflow(day,(hour*4)-3,:).*oTS;			% ts = 1 - 15 mins
		int_wind(day,hour,:) = kwind(day,hour*4,:).*oTS     +...
						       kwind(day,(hour*4)-1,:).*oTS +...
							   kwind(day,(hour*4)-2,:).*oTS +...
							   kwind(day,(hour*4)-3,:).*oTS;
		int_vO2(day,hour,:) = vO2(day,hour*4,:).*oTS     +...
						      vO2(day,(hour*4)-1,:).*oTS +...
							  vO2(day,(hour*4)-2,:).*oTS +...
							  vO2(day,(hour*4)-3,:).*oTS;
	end
%toc - ~0.2 seconds per day.
end
