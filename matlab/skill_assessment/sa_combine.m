%following Skill Metrics modified from work by Grant Law (2008)
%requires obs & mod to be same length

%apply intdefine as necessary (to model or obs, or both) & remove Nans
%example:
    %Tidei=intdefine(Tide(:,1),Tide(:,2),0.01042,0.0625,xmin,xmax);
    %i=find(~isnan(SalObsS324_15min(:,2)));
    %Obs=SalObsS324_15min(i,:)
    %Mdl=S3S24157(i,:)
%may need to apply 'remove duplicate numbers':
    %diff1 = S13t(2:end)-S13t(1:end-1);
    %S13t(diff1==0)=[];
    %S13v(diff1==0)=[];

function SA = sa_combine(obs,mod)
%define final table for metrics:
SA=zeros([1 7]);
%set obs & mod values:
o = obs;
m = mod;

%remove NaNs from data:
i = find(~isnan(o));
o = o(i);
m = m(i);

%Pearson's Corr Coeff
X(:,1)= o;
X(:,2)= m;
pcorr=corr(X);
pcc=pcorr(1,2);

%Model Bias
Em = mean(m);
Eo = mean(o);
varo = std(o);
if varo~=0;
    mb = (Em-Eo)./varo;
else
    mb = 0;
end

%Correlation Skill
cor = cov(o,m);
cor = cor(1,2);
varo = var(o);
varm = var(m);
varcomp = sqrt(varo*varm);

if varcomp~=0;
    corskill = cor/varcomp;
else
    corskill = 1;
end

%Mean Square Error
mse = mean((m-o).^2);

%Root Mean Square Error
%mse = msecalc(o,m);
rmse = sqrt(mse);

%Normalized Standard Deviation
dev = std(o);
if dev ~= 0;
    nstdev = std(m)./std(o);
else
    nstdev = 1;
end
SA(1)=pcc;
SA(2)=mb;
SA(3)=corskill;
SA(4)=mse;
SA(5)=rmse;
SA(6)=nstdev;
