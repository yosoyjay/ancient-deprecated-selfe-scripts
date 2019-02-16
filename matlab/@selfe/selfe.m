classdef selfe < dynamicprops
    % obj=selfe(DirIn) - A netcdf style data class to assist reading 
    % SELFE output data from "DirIn", for matlab 2008a (and later
    % versions?).  
    %
    % Index reference is strictly in this order:
    % time, elem, vert, u/v
    % e.g. selfeobj.hvel(2,1,1,1) would access u velocity for the 2nd timestep, 1st
    % element and the bottom layer.
    %
    % Note the following important information:
    % hgrid.gr3, vgrid.in and param.in files must be present in "DirIn", as well
    % output data.  
    %
    % The directory must be "clean" with only combined output files,
    % and so uncombined files should be moved to a seperate directory
    % within or outside of the "DirIn" directory.
    %
    % Note that filenames and date format (UK format) are hardwired in
    % "/selfeutility/initSELFEvar.m" as:
    %
    %     result.param_file='param.in';
    %     result.grid_file='hgrid.gr3';
    %     result.vgrid_file='vgrid.in';
    %     result.dateFormat='dd/mm/yyyy HH:MM'
    %
    % start date needs to be at the start of the param file, 1st two lines of the 
    % param file may look like this for e.g.
    % ! Note require start time at start of param file for selfe data object...
    % ! 20/10/2010 10:00
    %
    % WARNING: as file access occurs in the background of this function it may be slow to respond, 
    % particularly for a large number of files, so be careful accessing the time dimension with a ":"
    % e.g. selfeobj.hvel(:,:,:,:) may be slow to respond or run out of memory if there are a large 
    % number of files.
    % 
    % Example of use using provided example files, trisurf plot of surface
    % elevations at time 1:
    % %load SELFE object from Example directory provided
    % dirin=which('selfe');
    % [dirin file]=fileparts(dirin);
    % cd(dirin)
    % slf=selfe('..\Example\');
    % %plot elevations
    % figure; 
    % trisurf(slf.elem,slf.x,slf.y,slf.elev(1,:)'); 
    % view([0 90]);
    % shading interp;
    %
    % ---------------------------------------------------------------------
    % 
    %  2 September 2011, version 1 Beta for testing
    %
    %  Ben Knight, Cawthron Institute, 2011 (Ben.Knight@cawthron.org.nz)
    
    properties
        time
        matTime
        x
        y
        relz
        depth
        elem
        datinfo
    end %properties
    
    methods
        %% constructor
        function obj=selfe(DirIn)
            if nargin<1;
                return;
            end
            dat=initSELFE(DirIn); 
            
            %get critical grid information
            obj.time=dat.time(:,1); %time in seconds
            obj.matTime=dat.datetime(:); %matlab time
            obj.x=dat.hgrid.x;  %node x
            obj.y=dat.hgrid.y;  %node y
            obj.depth=dat.hgrid.depth; %node z
            %get relative depths.. not exact! 
            
            obj.relz=sb_computeZlevels(obj.depth, 0*obj.depth, dat.vgrid)./repmat(obj.depth,[1 dat.vgrid.nvrt]);
            obj.elem=dat.hgrid.elem(:,3:5);
            %clear unnessecary data
            rmfield(dat,{'datetime','hgrid'});
            %keep useful data
            obj.datinfo=dat;
            %clear unnessecary data
            clear dat;
            
            %add file specific properties
%             count=1;
            for vars=obj.datinfo.varNames'
                var=char(vars);
                prop=obj.addprop(var);
                eval(['obj.' var '=selfevar(''' var ''');']); %create selfevar object
%                 obj.findprop(var);
%                  eval(['prop.GetMethod = @get_datafromfile;']);
%                 count=count+1;
            end
        end  %constructor
        
        %% display function
        function disp(obj)
            for fld=fields(obj)';
                fld2=char(fld);
                sz=size(eval(['obj.' fld2]));
                szstr='';
                for i=1:length(sz)
                    szstr=[szstr 'x' num2str(sz(i))];
                end
                szstr(1)=[];
                type=class(eval(['obj.' fld2]));
                switch type
                case 'selfevar'
                    disp(eval(['obj.' fld2]));
                otherwise      
                    disp([char(fld) ': [' szstr ' ' type ']'])
                end %switch
            end
        end
        
        %% subsref for "varibles" with associated files
        function value = subsref(obj, S)
            %get data from oject
            theVar = S(1).subs;
            if ~isstr(theVar)
                error(['SELFEobj: string variable required'])
            end
            if length(S)>1
                eval(['value = subsref(obj.' theVar ', S(2:end));']);
            else
                eval(['value = obj.' theVar ';']);
            end %if length(S)>1
         end
        
        %% access datinfo data
        function value = get.datinfo(obj)
            value=obj.datinfo;
        end
%         
%         function value = get_datafromfile(obj)
% %             S=subsref(obj,substruct('.','blah'));
%             value=1;
%         end
    end
    
    
end %classdef
