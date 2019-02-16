function bldate=isdate(datestr)
%
% use datenum to determine if 'datestr' is a valid datenum date
%
% Note that it will classify any string that can be processed by "datenum", so 
% may 'misclassify' some strings - so will require sanity checking of input
% data.
%
% e.g. isdate('1/1/2010') = true 
%       isdate('blah') = false
%       isdate('!1-12') = true  - will flag this string as a date
%
%ben knight 2011
%
bldate=true;
try
    tst=datenum(datestr);
catch
    bldate=false;
end
 
