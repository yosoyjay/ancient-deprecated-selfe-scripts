function [regions, fg] = read_regions(regPath, grPath)
% [regions] = read_regions(regPath, grPath)
%
% This function reads in region information from regPath and
% returns the definition in the regions struture.
%
% Input:
%	regPath - Path to region.dat file definition
%	grPath - Path to hgrid.gr3
%
% Output:
%	regions - Structure with regions information
%		.name - Region names
%		.points - Nodes
%		.area - Area of the region (m^2)
%		.elems - Boolean flags indicating whether element is in region
%	fg - Grid data
%
% lopezj - 3/18/2012
%
% hgrid2fg file
addpath '/home/workspace/users/lopezj/scripts/matlab/general_grid';

% Load in hgrid
fg = hgrid2fg(grPath);

% Deal with region file
fid = fopen(regPath, 'r');
fgets(fid);
nRegions = fscanf(fid,'%i\n',1);
for region = 1:nRegions
	info = fscanf(fid,'%i\n',1);
	fscanf(fid,'%i\n',1);
	regions(region).name = fscanf(fid,'%s\n',1);
	regions(region).points = fscanf(fid,'%f %f\n',[2,info])';

	% Identify nodes defined by polygon defined from region file
	% If all the nodes that define an element in fg.e are in the region, then that element is in region
	inRegion = inpolygon(fg.x,fg.y,regions(region).points(:,1),regions(region).points(:,2));
	% Boolean flag identifying elements in region
	elems = sum(inRegion(fg.e)') == 3;   
	% Calculate the area in this region
	regions(region).area = sum(fg.ar(elems));
	% Send bool of elements in region
	regions(region).elems = elems; 
end


