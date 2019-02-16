function [d,eta] = cacl_mole_diff(temp)
% [d] = cacl_mole_diff(temp)
%
% Calculates the molecular diffusion of oxygen using the formula described
% in Needoba et. al 2012.
%
% Input:
%	temp - Temperature in celcius.  ROW_MAJOR data
%
% Output:
%	d - Molecular diffusion of oxygen
%

% Constants
x = 2.26; % Association factor of water
M = 18;   % Molar weight of water
V = 25.6; % Molar volume at boiling point

% Convert Celcius to Kelvin
temp = temp + 273.15;

% Calculate the dynamic viscosity of water (\eta)
exponent = (247.8./(temp-140));
eta = (2.414*10^-5)*10.^exponent;

% eta needs to be converted from Pascal to centipoise?
% eta should be about 1 for 20 C
eta = eta.*1000;

% Calculate the molecular diffusion of oxygen
numerator = 7.4*10^-8*(x*M)^(1/2).*temp;
d = numerator./(eta.*V^(0.6));
