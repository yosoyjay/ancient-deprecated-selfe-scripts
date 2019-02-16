function zlev = sigmacalc(hc,thb,thf,nlev,sz,z)
%   LEVEL_DATA = SIGMACALC(h_c,theta_b,theta_f,num_levs,sz_trans,max_depth),
%   uses the standard SELFE vgrid.in parameters to build a representative
%   vertical discretization scheme useful in setting up runs. The output of
%   this function can be easily visualized using PLOT.
%  
%   H_C is the isobath beyond which vertical S-level constrictions, defined
%         by theta_f, are imposed
%   THETA_B defines the depth of the sub-surface vertical grid
%         constriction center. A value of 1.0 places the subsurface
%         constriction at the SZ-transition depth. Values from 1.0 to 0.0
%         move the subsurface constriction from the SZ-transition depth
%         towards the surface -- for example, 0.5 puts the constriction
%         half-way between the SZ-transition depth and the surface. A
%         THETA_B value of 0.0 essentially removes the subsurface
%         constriction entirely
%   THETA_F is the constriction factor. commonly used values range from 1
%         to 10
%   NUM_LEVS is the number of vertical levels you wish to use (please note:
%         an additional level, representative of the seafloor, will be
%         added to the output)
%   SZ_TRANS is the depth at which the S-levels give way to Z-levels
%   Z is the maximum depth of your horizontal grid
%
%   EXAMPLE:
%      lev = sigmacalc(30,.7,10,100,100,200);
%      plot(lev);
%

%   C. Grant Law 04-12-2011

if z<0
    z = -z;
end
if sz<0
    sz = -sz;
end
if hc<0
    hc = -hc;
end

sigma = -(0:1/(nlev-1):1);
z = -z .* (1./(1+exp(-6:.1:6))-1); zlev = z;
zlev(zlev>sz) = sz; 
[sigma,zlev] = meshgrid(sigma,zlev); sigma=sigma'; zlev=zlev';

cs=(1-thb).*sinh(thf.*sigma)./sinh(thf) + thb.*(tanh(thf.*(sigma+0.5))-tanh(thf.*0.5))./2./tanh(thf.*0.5);

spec = zlev<=hc;
zlev(spec) = sigma(spec).*zlev(spec);
zlev(~spec) = hc.*sigma(~spec) + (zlev(~spec)-hc).*cs(~spec);
zlev = [zlev;-z]';