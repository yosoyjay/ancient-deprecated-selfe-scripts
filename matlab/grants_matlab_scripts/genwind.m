function status = genwind(dt,runtime,u,v)
% DT is the timestep in seconds to be used in the file (must match the
%       timestep value specified in param.in!)
% RUNTIME is the length of the run in seconds
% U is the x component of the windspeed. Use one value for constant winds,
%       or enter an array for time-varying winds
% V is the y component of the windspeed. Use one value for constant winds,
%       or enter an array for time-varying winds

if length(u)==1
    u = [u u];
    v = [v v];
end
t = [dt:(runtime-dt)./(length(u)-1):runtime];
tn = dt:dt:runtime;
un = interp1(t,u,tn);
vn = interp1(t,v,tn);
fid = fopen('wind.th','w');
fprintf(fid,'%f %f\n',[un;vn]);
fclose(fid);
status = exist('wind.th','file');