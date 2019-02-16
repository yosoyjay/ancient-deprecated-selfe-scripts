function path = partpath(pthfile)

[fid,error] = fopen(pthfile);

if fid<0
    disp(error)
end

[a,count] = fscanf(fid,'%s',1);
[a,count] = fscanf(fid,'%i',1);
step = 0;

while count>0
    step = step+1;
    [a,count] = fscanf(fid,'%f %i',2);
    if count==0
        path.x = x;
        path.y = y;
        path.z = z;
        path.time = time;
        return
    end
    time(step) = a(1);
    pnum = a(2);
    for i = 1:pnum
        [data,count] = fscanf(fid,'%i %f %f %f',4);
        if count~=4
            break
        end
        x(i,step) = data(2);
        y(i,step) = data(3);
        z(i,step) = data(4);
    end
end
