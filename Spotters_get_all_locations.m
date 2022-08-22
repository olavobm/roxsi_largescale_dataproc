%% Script to get reference location of all Spotters (including
% Smart moorings). Here, I "reference location" means what will be
% used as the location of pressure sensors on Smart Moorings. This
% is NOT the average location because the Spotters does not spend
% time homogenously along its watch circle. Instead a least-squares
% fit is used to calculate the "best circle" and the coordinate at
% its center is taken.
%
% This script uses data in location.csv files, which is generated
% by the Sofar's python script that parses the data for a Spotter
% (note this DOES NOT (!) include the pressure data from Smart Mooring).

clear
close all


%%



