function [average_sa] = calc_average_sa(sa)
%
% Calculates average skill assessment after it has already
% been calculated using model_data_sa.m
%
% Input:
%	sa - Skill asssessment output from model_data_sa.m
%
% Output:
%	average_sa - Average skill assesssment for each run
%

% Deconstrct the runName.stationName.depthName(1,7) struct
% and calculate average skill assessment for every run
rN = fieldnames(sa);
for i = 1:length(rN);
	count = 1;
	sN = fieldnames(sa.(rN{i}));
	average = [0 0 0 0 0 0 0];
	for j = 1:length(sN)
		dN = fieldnames(sa.(rN{i}).(sN{j}));
		for k = 1:length(dN)
			average = average + sa.(rN{i}).(sN{j}).(dN{k});
			count = count + 1;
		end
	end
	average_sa.(rN{i}) = average ./ count;
end
