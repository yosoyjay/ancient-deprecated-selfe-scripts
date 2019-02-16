function [orderedCasts] = order_casts(ctdData)
% Just rearranges casts so they are in order by longitude

% Don't destory original data
ctdCasts = ctdData;

for i = 1:size(ctdCasts,2)
	max = -9999;
	for j = 1:98
		if ctdCasts(j).long > max
			idx = j;
			max = ctdCasts(j).long;
		end
	end
	orderedCasts(i) = ctdCasts(idx);
	ctdCasts(idx).long = -9999;
end
