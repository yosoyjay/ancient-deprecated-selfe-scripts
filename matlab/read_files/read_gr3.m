function file14=read_gr3(file)

[grid]=gr_readHGrid(file);
file14.AGRID=grid.name;
file14.NE=grid.ne;
file14.NP=grid.np;
file14.JN=grid.nn;
file14.X=grid.x;
file14.Y=grid.y;
file14.DP=grid.depth;
file14.JE=grid.elem(:,1);
file14.NHY=grid.elem(:,2);
file14.NM=grid.elem(:,3:5);
[x,remain]=strtok(grid.eofLines{1},{'=','!'});
file14.NOPE=str2num(x);
[x,remain]=strtok(grid.eofLines{2},{'=','!'});
file14.NETA=str2num(x);

l=3;
for i=1:file14.NOPE
    [x,remain]=strtok(grid.eofLines{l},{'=','!'});
    file14.NVDLL(i)=str2num(x);
    for j=1:file14.NVDLL(i)
        file14.NBDV(i).node(j)=str2num(grid.eofLines{l+j});
    end
    l=l+str2num(x)+1;
end
    

[x,remain]=strtok(grid.eofLines{l},{'=','!'});
file14.NBOU=str2num(x);

[x,remain]=strtok(grid.eofLines{l+1},{'=','!'});
file14.NVEL=str2num(x);

l=l+1;
for i=1:file14.NBOU
    l=l+1;
    [x,remain]=strtok(grid.eofLines{l},{'=','!'});
    x=str2num(x);
    file14.land(i).NVELL=x(1);
    file14.land(i).IBTYPE=x(2);
    for j=1:x(1)
        file14.land(i).NBVV(j)=str2num(grid.eofLines{l+j});
    end
    l=l+j;
end