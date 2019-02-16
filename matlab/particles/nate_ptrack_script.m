% relshape=[46 30 -124 12; 46 30 -124 18; 46 30 -124 24];
% reltime=[2009 5 14 12 0; 2009 5 14 13 0; 2009 5 14 14 0];
relshape=RELEASE_PTS;
reltime=START_TIME;


dims = [-124.5, 46.2; -124.0, 46.8];

names = {'46_30_124_12', '46_30_124_18', '46_30_124_24'}

for i=2:size(relshape,1)
    for j=1:size(reltime,1)
        ttle = sprintf('part_track_%s_t%di%d', names{i}, reltime(j, 4), reltime(j, 5));
        parttrack('dev',relshape(i, :),reltime(j, :),dims,ttle);
    end
end
     