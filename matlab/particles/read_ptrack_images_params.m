function [params] = read_ptrack_images_params(fname)

fid = fopen(fname, 'r');
args = textscan(fid, '%s %s', 'Delimiter', '=');
fclose(fid);

keys = args{1};
vals = args{2};
for i=1:size(keys,1)
    params.(keys{i,1}) = vals{i,1};
end