function [trackid, lon, lat, dnum] = read_drifter(fname)

fid = fopen(fname, 'r');

astr = fgets(fid);
trackid = sscanf(astr, '%s, %*d-%*d-%*d %*d:%*d:%*d-%*d');
vals = sscanf(astr,'%*s%d-%d-%d%d:%d:%d%d');
dnumbase = datenum(vals(1), vals(2), vals(3));

lon=[];
lat=[];
dnum=[];

fgets(fid);
line = fgets(fid);
while (line>=0)
    vals = sscanf(line, '%d, %f, %f');
    dnum(end+1) = dnumbase+(vals(1)-86400)/86400;
    lon(end+1) = vals(2);
    lat(end+1) = vals(3);
    line=fgets(fid);
end
