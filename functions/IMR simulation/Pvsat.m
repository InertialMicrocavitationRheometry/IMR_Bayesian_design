function [ Pv ] = Pvsat( T )
%Calculates the saturated vapor pressure using temperature 
Pv = 1.17e11*exp(-5200./(T)); 