function nodes = get_node_cords(node_file, grid)
%
%  This script extracts the coordinates of the nodes numbered in the input 
%  file node_file from xmgredit6.  It returns a [n,2] dimensional array with 
%  the x,y coordinates of the nodes listed in the file node_file.  
%
%  The grid number is only used if there is not an hgrid.gr3 file in the cwd.
%
%  Input:
%  node_file - File output from xmgredit defining a list of nodes
%  grid      - Grid number.  [15, 22, 26, or 27]
%
%  Output:
%  nodes - The x,y coordinates of the nodes listed in the node_file
%
%  lopezj - 6/19/2011
%

% Check path and availability of m-elio goodies
	ret_code = exist('gr_readHGrid.m','file');
	if(ret_code ~= 2)
		addpath '/usr/local/cmop/matlab/cmop/m-elio/';
	end

% Check path and availability of hgrid.gr3 
	ret_code = exist('hgrid.gr3','file');
	if(ret_code ~= 2)
		if(grid == 15)
			addpath '/home/workspace/project/lopezj/scripts/selfe_setup/db15/';
		elseif(grid == 22)
			addpath '/home/workspace/project/lopezj/scripts/selfe_setup/db22/';
		elseif(grid == 26)
			addpath '/home/workspace/project/lopezj/scripts/selfe_setup/db26/';
		elseif(grid == 27)
			addpath '/home/workspace/project/lopezj/scripts/selfe_setup/db26c/';
		else
			print 'Unable to identify grid. Please choose 15, 22, 26, or 27'
			exit
		end
	end

% Get the x,y positions of the nodes in the node_file 
	node_nums = importdata(node_file, ' ', 3);
	node_nums = node_nums.data(:,1);
	h_grid = gr_readHGrid('hgrid.gr3');
	nodes = zeros(size(node_nums,1),2);
	for i = 1:size(node_nums)
		nodes(i,1:2) = h_grid.nodes(node_nums(i,1),2:3);
	end
end
