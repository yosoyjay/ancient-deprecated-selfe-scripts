function [tot_cross,tot_along,tot_u,tot_v,trDist] = resid_across_transect(h,bpfile,dud,dvd)
% [tot_cross,tot_along] = resid_across_transect(h,bpfile,dud,dvd)

[ob]= ob_ini_fromTrasect(h, bpfile);
badnodes = isnan(ob.xy.w(:,1));
%last = find(~isnan(ob.xy.w(:,1)), 1, 'last' );
trDist = sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2+(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2);
%trLen = [cumsum(trDist)];
sina = repmat(diff(ob.xy.y)./trDist, 1 , size(dud,2));
cosa = repmat(-diff(ob.xy.x)./trDist, 1 , size(dud,2));
trDist = repmat(trDist, 1 , size(dud,2));
%trLen = repmat(trLen, 1 , size(dud,2));

u=ob.xy.H*double(dud);
v=ob.xy.H*double(dvd);
u(badnodes,:) = nan;
v(badnodes,:) = nan;
tot_u = 0.5 * (u(1:end-1,:) + u(2:end,:)).*trDist;
tot_v = 0.5 * (v(1:end-1,:) + v(2:end,:)).*trDist;
tot_cross = tot_u .* sina + tot_v .* cosa;
tot_along = tot_u .* cosa + tot_v .* -sina;
