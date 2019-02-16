function [sa] = model_data_sa(model, obsInt)
% [sa] = model_data_sa(model, obsInt)
% Calculates skill assessment measures for model data
%
% Input:
%	model - Model data in the from returned by model_data_comp
%	obsInt - Observation data in the struc returned by model_data_comp
% Ouput:
%	sa - Calculated skill assessment 
%

% Addpath for comp_stuct
addpath = '/path/to/directory'

% Check if all stations are in both model and obsInt
[fs1, fs2, er] = comp_struct(model,obsInt)
skipStat = [];
if ~isempty(fs1)
	skipStat = fs1;
end
if ~isempty(fs2)
	skipStat = append(fs1 with fs2)
end

% Now if the station or depth is missing or present
% in skipStat (rename) then you don't do sa for that 
% and fill the array in with a bunch of NaNs.


% Use the model or obs data with least fields.
% Hack to continue when one is missing some data.
names = fieldnames(model);
modStatNames = fieldnames(model.(names{1}));
obsStatNames = fieldnames(obsInt);
if length(modStatNames) <= length(obsStatNames)
	stats = modStatNames;
else 
	stats = obsStatNames;
end

% Loop over model runs and stations in structs and perform skill assessment
for run = 1:length(names)
	for station = 1:length(stats)
		modDepths = fieldnames(model.(names{run}).(stats{station}));
		obsDepths = fieldnames(obsInt.(stats{station}));
		for depth = 1:size(obsDepths)
			temp = sa_combine(obsInt.(stats{station}).(obsDepths{depth})(:,2), ...
							  model.(names{run}).(stats{station}).(modDepths{depth})(:,2));
			sa.(names{run}).(stats{station}).(modDepths{depth}) = temp;
		end
	end
end
