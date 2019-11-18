%% START_ASPIRE

clear all, close all

hdisk = 'E:\';
% hdisk = '/media/fedelet1/Tommaso USZ 1/'

% brainstorm_path = '\\172.16.118.134\TFedele\Tommaso\Matlab Toolboxes\brainstorm_190129\brainstorm3';
brainstorm_path = [hdisk 'Valerii\brainstorm_190129\brainstorm3'];
addpath(genpath(brainstorm_path))

fieldtrip_path = [hdisk '\_USZserver\Matlab Toolbox\FIELDTRIP\fieldtrip-20191028'];
addpath( (fieldtrip_path))
ft_defaults

addpath(genpath([ 'E:\_USZserver\hfEEG\hfEEG software analysis\hfEEGAquistionRoutines']));

workdir = [  'C:\Users\Tom\Documents\GitHub\ASPIRE'];
cd(workdir)