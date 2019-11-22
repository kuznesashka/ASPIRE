%Add path and all subfolders 
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

maindir_data = 'E:\Valerii\EPILEPSY\'

%% Paths -- subject info
% !!! File separator for current platform

paths.anat = [maindir_data filesep, ...
                'MEG_Tommaso' filesep 'anat' filesep]; % anatomy from Brainstorm
paths.data = [maindir_data filesep, ...
                'MEG_Tommaso' filesep 'data' filesep]; % data from Brainstorm
paths.root = 'E:\Valerii\45_cases\'; %45 cases folder

paths.subj_name = 'B1C2'; %when two or more recordings 'B4Z2\Rec_01'
paths.case = 'B1C2';
paths.fname = 'B1C2_ii_run1_raw_tsss_mc_art_corr'; %artifact corrected file name
paths.sh_fname = 'B1C2_ii_run1_raw_tsss_mc'; %short file name


 
 

 