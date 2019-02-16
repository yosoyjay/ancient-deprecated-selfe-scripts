function [lon,lat]=convm2ll(x,y,RefLon,RefLat);
%CONVM2LL convert cartesian coordinates (meters) to long/lat (deg)
% CONVM2LL converts from cartesian coordinates in meters to
% long/lat coordinates by inverting the conversion formulas 
% used in the routine CONVLL2M.  This is the inverse of
% converting long/lat coordinates to cartesian coordinates in
% meters. The x,y  coordinates MUST have been computed using
% CONVLL2M, and the  reference long/lat pair MUST be the same
% as those used in the  "forward" conversion (defaults to Boston, MA).
% Any translation in the x,y coordinates that might have been
% done MUST be "undone" prior to conversion back to long/lat
% coordinates.  The inversion for latitude is done with 5 
% Newton-Raphson iterations with an initial guess of the
% reference latitude.  The longitude is trivial.
%
%  Inputs: x,y - cartesian coordinates in meters
%          RefLon,RefLat - reference long/lat values
%
% Outputs: lon,lat - long/lat coordinates.
%
% Call as:  [lon,lat]=convm2ll(x,y,RefLon,RefLat);

% Implement default lat & lon

if nargin == 3
    switch 1
        case strcmp(RefLon,'Atlantic')
            method = 'blant';
            RefLon = 0;
            RefLat = 35.3256;
            R=6367500;  % Approx earth radius in meters
        case strcmp(RefLon,'Pacific')
            method = 'cmop';
            RefLon = -122.555;
            RefLat = 43.585;
            R = 6378206.4;  % Approx earth radius in meters
    end
end
if nargin == 2
    method = 'cmop';
    RefLon = -122.555;% -122.555
    RefLat = 43.585;% 43.585
    R = 6378206.4;  % Approx earth radius in meters
end

if nargin < 2
   error('Too few arguments to CONVM2LL')
end

if nargout ~=2
   error('CONVM2LL requires two return arguments.')
end

% Check sizes of lon,lat
if ~all(size(x)==size(y))
   error('Sizes of lon and lat are NOT equal')
end

switch 1
    
    case strcmp(method,'blan')
        deg2rad=pi/180;
        fac=R*cos(RefLat*deg2rad);

        % Invert for lon (easy)
        lon=x/fac/deg2rad;
        phin=RefLat.*ones(size(x)).*deg2rad;
        CC=exp(y/fac);
        
        % 5 Newton-Raphson iterations, assume convergence!!
        for i=1:5
            phinm1=phin;
            numer=CC - (1+sin(phinm1))./(cos(phinm1));
            denom= -1- (1+sin(phinm1)).*tan(phinm1)./cos(phinm1);
            phin=phinm1-numer./denom;
        end

        lat=phin/deg2rad;
        
    case strcmp(method,'cmop')
        
        n1 = size(x);
        % Memory pre-allocation
        %
        lat=zeros(size(x));
        lon=zeros(size(x));
        
        %sa = 6378137.000000 ; sb = 6356752.314245; % For WGS-84
      	sa = 6378137.000000 ; sb = 6356752.3141; % For WGS-80
  
       	e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
    	e2cuadrada = e2 ^ 2;
     	c = ( sa ^ 2 ) / sb;

     	lat =  y ./ ( 6366197.724 * 0.9996 );                                    
     	v = ( c ./ ( ( 1 + ( e2cuadrada .* ( cos(lat) ) .^ 2 ) ) ) .^ 0.5 ) .* 0.9996;
    	a = (x-500000) ./ v;
     	a1 = sin( 2 .* lat );
       	a2 = a1 .* ( cos(lat) ) .^ 2;
        j2 = lat + ( a1 ./ 2 );
     	j4 = ( ( 3 .* j2 ) + a2 ) ./ 4;
    	j6 = ( ( 5 .* j4 ) + ( a2 .* ( cos(lat) ) .^ 2) ) ./ 3;
      	alfa = ( 3 ./ 4 ) .* e2cuadrada;
     	beta = ( 5 ./ 3 ) .* alfa .^ 2;
     	gama = ( 35 ./ 27 ) .* alfa .^ 3;
     	Bm = 0.9996 .* c .* ( lat - alfa .* j2 + beta .* j4 - gama .* j6 );
      	b = ( y - Bm ) ./ v;
       	Epsi = ( ( e2cuadrada .* a.^ 2 ) ./ 2 ) .* ( cos(lat) ).^ 2;
      	Eps = a .* ( 1 - ( Epsi ./ 3 ) );
       	nab = ( b .* ( 1 - Epsi ) ) + lat;
      	senoheps = ( exp(Eps) - exp(-Eps) ) ./ 2;
       	Delt = atan(senoheps ./ (cos(nab) ) );
       	TaO = atan(cos(Delt) .* tan(nab));
        lon = (Delt .*(180 ./ pi ) ) + RefLon;
     	lat = ( lat + ( 1 + e2cuadrada.* (cos(lat).^ 2) - ( 3 ./ 2 ) .* e2cuadrada .* sin(lat) .* cos(lat) .* ( TaO - lat ) ) .* ( TaO - lat ) ) .* (180 ./ pi) + RefLat;
   
end