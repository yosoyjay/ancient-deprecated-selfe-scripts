% DRAWELEMS draw 2D/3D FEM element configuration 
%
% DRAWELEMS draws element boundries given a valid grid structure.  
%
%  INPUT : fem_grid_struct - (from LOADGRID, see FEM_GRID_STRUCT) 
%               dim - either '2D' or '3D'
%           
% OUTPUT : hel - handle to the element object.
%
%   CALL : hel=drawelems(fem_grid_struct,dim);
%
% Written by: Brian O. Blanton
% Summer 1997
%                  
function hel=drawelems(fem_grid_struct,dim) 

% DEFINE ERROR STRINGS
err1=['Not enough input arguments; need a fem_grid_struct'];

% check arguments
if nargin ==0 
   error(err1);
end  

if ~is_valid_struct(fem_grid_struct)
   error('    Argument to DRAWELEMS must be a valid fem_grid_struct.')
end

if nargin==1
    dim='2D';
end

% Extract grid fields from fem_grid_struct
%
elems=fem_grid_struct.e;
   % COPY FIRST COLUMN TO LAST TO CLOSE ELEMENTS
   %
   elems=elems(:,[1 2 3 1]);
x=fem_grid_struct.x;
y=fem_grid_struct.y;
if strcmp(dim,'3D')
    z=fem_grid_struct.z;
end

elems=elems';
[m,n]=size(elems);
xt=x(elems);
yt=y(elems);
if strcmp(dim,'3D')
    zt=z(elems);
end
if n~=1 
   if m>n
      xt=reshape(xt,n,m);
      yt=reshape(yt,n,m);
      if strcmp(dim,'3D')
          zt=reshape(zt,n,m);
      end
   else
      xt=reshape(xt,m,n);
      yt=reshape(yt,m,n);
      if strcmp(dim,'3D')
          zt=reshape(zt,m,n);
      end
   end
   xt=[xt
       NaN*ones(size(1:length(xt)))];
   yt=[yt
       NaN*ones(size(1:length(yt)))];
    if strcmp(dim,'3D')
        zt=[zt
           NaN*ones(size(1:length(zt)))];
    end
end
xt=xt(:);
yt=yt(:);
if strcmp(dim,'3D')
    zt=zt(:);
end
%
% DRAW GRID
%
if strcmp(dim,'3D')
    hel=line(xt,yt,zt,'LineWidth',1,'LineStyle','-','Color',[.6 .6 .6]);
    %hel=line(xt,yt,zt,'LineWidth',0.5,'LineStyle','-','Color','k');
else
    hel=line(xt,yt,'LineWidth',0.5,'LineStyle','-','Color',[.6 .6 .6]);
    %hel=line(xt,yt,'LineWidth',0.5,'LineStyle','-','Color','k');
end
set(hel,'Tag','elements');
 
%
%        Brian O. Blanton
%        Curr. in Marine Sciences
%        15-1A Venable Hall
%        CB# 3300
%        Uni. of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        919-962-4466
%        blanton@cuda.chem.unc.edu
%
%        Summer 1997
%
