classdef selfevar < selfe
    % obj=selfe(varIn) - Data subclass to assist reading SELFE output files from "DirIn" 
    %
    % 2 September 2011, version 1 Beta
    % Ben Knight, Cawthron Institute, 2011
    properties
        name
        file
        dims
    end %properties
    
    methods
        %% constructor
        function obj=selfevar(varIn)
            parent=evalin('caller','obj;');
            %get file type and max dimensions
            whatFiletype=strcmp(varIn,parent.datinfo.varNames);
            theFiletype=parent.datinfo.filetype{whatFiletype};
            
            obj.name=varIn;
            obj.file=[varIn theFiletype];
            %get dims, time,xy,z,u/v
            switch theFiletype
                case '.61' %time,xy
                    obj.dims=[length(parent.matTime) length(parent.x)];
                case '.62' %time,xy,u/v
                    obj.dims=[length(parent.matTime) length(parent.x) 2];
                case '.63' %time,xy,z
                    obj.dims=[length(parent.matTime) length(parent.x) size(parent.relz,2)];
                case '.64' %time,xy,z,u/v
                    obj.dims=[length(parent.matTime) length(parent.x) size(parent.relz,2) 2];
            end

        end  %constructor
        
        %% display function
        function disp(obj)
            sz=obj.dims;
            szstr='';
            for i=1:length(sz)
                szstr=[szstr 'x' num2str(sz(i))];
            end
            szstr(1)=[];
            type='selfevar';
            disp([char(obj.name) ': [' szstr ' ' type ']'])
        end
        
        %% subsref for "varibles" with associated files
        function value = subsref(obj, S)
            parent=evalin('caller','obj;');
            
            theIdx = S(1).subs;  %time,el,z
            %get file type and check number of input dims
            whatFiletype=strcmp(obj.name,parent.datinfo.varNames);
            theFiletype=parent.datinfo.filetype{whatFiletype};

            %check input index and flag error if wrong
            ndim=length(obj.dims); 
            if length(theIdx)>ndim|length(theIdx)<ndim
                error('Wrong number of dims for variable')
            end

            %check dimensions match size of matrix?
            for i=1:ndim
                if isstr(theIdx{i})
                    if theIdx{i}==':'
                        theIdx{i}=1:obj.dims(i);
                        S(1).subs=theIdx;
                    else
                        error('selfevar: Unrecognised index character')
                    end
                end
                if max(theIdx{i})>obj.dims(i)
                    error('selfevar: Index exceeds dimensions')
                end
            end

            %get file(s) for requested time period
            if length(theIdx)==1
                timeIdx=theIdx{1}==parent.datinfo.time(:,4);
            elseif length(theIdx)>=2
                timeIdx=find(parent.datinfo.time(:,4)>=theIdx{1}(1)& ...
                            parent.datinfo.time(:,4)<=theIdx{1}(end));
            end
            whatFileNum=parent.datinfo.time(timeIdx,2);
            whatFileNums=unique(whatFileNum);

            %predim data output
            switch ndim 
                case 2 %time,el
                    datOut2=zeros(length(S(1).subs{1}),length(S(1).subs{2}));
                case 3 %time,el,z
                    datOut2=zeros(length(S(1).subs{1}),length(S(1).subs{2}),length(S(1).subs{3}));
                case 4 %time,el,z,u/v
                    datOut2=zeros(length(S(1).subs{1}),length(S(1).subs{2}),length(S(1).subs{3}),length(S(1).subs{4}));
            end

            %buildup mtx if multiple files
            count=0;
            whatTimeIdx=parent.datinfo.time(timeIdx,3);

            for i=1:length(whatFileNums)
                theFile=[num2str(parent.datinfo.filenums(whatFileNum(i))) '_' obj.file];
                fileTimeIdx=whatTimeIdx(whatFileNum==whatFileNums(i));

                %open the file
                hdr=sz_readHeader([parent.datinfo.SELFE_Dir theFile]);
                [datOut,ts] = sz_readTimeStep(hdr,fileTimeIdx);
                S2=S(1);
                %update for index of file, as used to reference 'cut-down' version of matrix
                S2.subs{1}=fileTimeIdx'-min(fileTimeIdx)+1; 
                idx=count+S2.subs{1};
                 switch theFiletype
                case '.61' %time,el
                    if length(size(datOut))==3, 
                        datOut=permute(datOut,[3 1 2]);
                    else
                        datOut=datOut';
                    end
                    datOut2(idx,:)=subsref(datOut,S2);
                case '.62' %time,el,z
                    datOut=permute(datOut,[3 1 2]);
                    datOut2(idx,:,:)=subsref(datOut,S2);  
                case '.63' %time,el,z
                    datOut=reshape(datOut,[hdr.vgrid.nLevels, hdr.hgrid.np size(datOut,2) size(datOut,3)]);
                    datOut=squeeze(permute(datOut,[4 2 1 3]));
                    datOut2(idx,:,:)=subsref(datOut,S2);  
                case '.64' %time,el,z,u/v
                    datOut=reshape(datOut,[hdr.vgrid.nLevels, hdr.hgrid.np size(datOut,2) size(datOut,3)]);
                    datOut=permute(datOut,[4 2 1 3]);
                    datOut2(idx,:,:,:)=subsref(datOut,S2);
                end
                count=count+length(fileTimeIdx);
            end

            value=datOut2;

%                     %could calc real zlevels here?... but slow!
%                     %e.g. obj.zlevels=sb_computeZlevels(obj.depth, ts{1}.eta,dat.vgrid);
%                     clear ts
%                     
%                     switch theFiletype
%                     case '.61'
%                         value=map_sz2hts(hdr,data);
%                     case '.62'
%                         value=map_sz2hts(hdr,data);
%                     case '.63'
%                         value=map_sz2hts(hdr,data);
%                     case '.64'
%                         value=map_sz2hts(hdr,data);
%                     end
           
        fclose('all')    
        end

    end
    
    
end %classdef