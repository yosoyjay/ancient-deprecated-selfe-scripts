function bnd = input_complex_bnd(file)

f =fopen(file,'r')
fgets(f);
bnd = [];
line = fgets(f);
nbnd = sscanf(line,'%f');
for i= 1:nbnd
    line = fgets(f);
    tst = sscanf(line,'%f',1);
    sz=size(tst);
    if (sz==0)
        line=fgets(f);
    end
    nnodes = sscanf(line,'%f',1);
    d = fscanf(f,'%f %f',[2 nnodes])';
    bnd = [bnd; nan nan; d; d(1,:)];   
end
