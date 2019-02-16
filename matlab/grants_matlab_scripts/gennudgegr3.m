function [nudge_fac,status] = gennudgegr3(rate,method)
% 

fg = hgrid2fg('hgrid.gr3');
fid1 = fopen('t_nudge.gr3','w');
fid2 = fopen('s_nudge.gr3','w');
fprintf(fid1,'%s\n',['Nudging rates using ' method ' method']);
fprintf(fid2,'%s\n',['Nudging rates using ' method ' method']);
fprintf(fid1,'%i %i\n',[length(fg.e) length(fg.x)]);
fprintf(fid2,'%i %i\n',[length(fg.e) length(fg.x)]);

switch method
    case 'edge'
        edgex = [];
        edgey = [];
        for ob = 1:length(fg.ob)
            edgex = [edgex;fg.x(fg.ob(ob).ind)]; 
            edgey = [edgey;fg.y(fg.ob(ob).ind)]; 
        end
        dist = zeros(length(fg.x),1);
        nudge_fac = dist;
        for i = 1:length(fg.x)
            d = sqrt((edgex-fg.x(i)).^2 + (edgey-fg.y(i)).^2);
            dist(i) = min(d);
        end
        %nudge_fac(dist<=50000) = (1 - dist(dist<=50000)./50000) .* rate;% Linear
        %nudge_fac = exp(-(dist)./10000) .* rate;% Exponential
        nu_length = 50000;% In meters
        nudge_fac = 1./(1 + exp((dist-(nu_length./2))./(nu_length/10))) .* rate;% Sigmoidal
        block = [1:length(fg.x);fg.x';fg.y';nudge_fac'];
        fprintf(fid1,'%i %f %f %f\n',block);
        fprintf(fid2,'%i %f %f %f\n',block);
    otherwise
        error('Didn''t recognize that method...');
end

block = [1:length(fg.e);ones(1,length(fg.e)).*3;fg.e'];
fprintf(fid1,'%i %i %i %i %i\n',block);
fprintf(fid2,'%i %i %i %i %i\n',block);
fclose(fid1);
fclose(fid2);
status(1) = exist('t_nudge.gr3','file');
status(2) = exist('s_nudge.gr3','file');