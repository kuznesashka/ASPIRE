% https://natmeg.se/MEEG_course2018/dipole_fitting.html

paths.brainstorm = '/Users/valery/MEG/brainstorm3';
paths.fieldtrip = '/Users/valery/MEG/fieldtrip-20191008';
addpath(genpath(paths.brainstorm));
addpath((paths.fieldtrip))
addpath(([paths.fieldtrip filesep 'plotting']));
addpath(([paths.fieldtrip filesep 'utilities' filesep 'private']));
warning('off', 'MATLAB:MKDIR:DirectoryExists');


%DataFile = sInput.FileName;
%DataMat = in_bst_data(DataFile);
Data = load('/Users/valery/MEG/brainstorm_db/MEG_Tommaso/data/B1C2/B1C2_ii_run1_raw_tsss_mc_art_corr/data_block001.mat');

% Load channel file
% ChannelMat = in_bst_channel(sInput.ChannelFile);
channels = load('/Users/valery/MEG/brainstorm_db/MEG_Tommaso/data/B1C2/B1C2_ii_run1_raw_tsss_mc_art_corr/channel_vectorview306_acc1.mat');
% Get selected sensors
% iChannels = channel_find(ChannelMat.Channel, SensorTypes);
channel_type = 'grad';
if strcmp(channel_type, 'grad') == 1
    grad_idx = setdiff(1:306, 3:3:306);
    channel_idx = grad_idx(Data.ChannelFlag(grad_idx)~=-1);
elseif strcmp(channel_type, 'mag') == 1
    magn_idx = 3:3:306;
    channel_idx = magn_idx(Data.ChannelFlag(magn_idx)~=-1);
end

% Load head model
% HeadModelMat = in_bst_headmodel(sHeadModel.FileName);
G3 = load('/Users/valery/MEG/brainstorm_db/MEG_Tommaso/data/B1C2/B1C2_ii_run1_raw_tsss_mc_art_corr/headmodel_surf_os_meg.mat');

% Convert head model
ftHeadmodel = out_fieldtrip_headmodel(G3, channels, channel_idx);

% Convert data file
ftData = out_fieldtrip_data(Data, channels, channel_idx, 0);
% Generate rough grid for first estimation
%? GridOptions = bst_get('GridOptions_dipfit');
%? GridLoc = bst_sourcegrid(GridOptions, G3.SurfaceFile);
timelockeds.time   = 1:1000;
timelockeds.avg    = ftData.trial{1}(:,1:1000);
timelockeds.label  = ftData.label;
timelockeds.dimord = 'chan_time';
timelockeds        = copyfields(ftData, timelockeds, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});

DipoleModel = 'regional';
NumDipoles = 1;
TimeWindow = [45 65]; % samples (ms)
SymmetryConstraint = [];
% Prepare FieldTrip cfg structure
cfg = [];
cfg.channel     = {channels.Channel(channel_idx).Name};
cfg.headmodel   = ftHeadmodel;
cfg.latency     = TimeWindow;
cfg.numdipoles  = NumDipoles;
cfg.model       = DipoleModel;
cfg.nonlinear   = 'yes';
% cfg.grid.pos    = GridLoc;
% cfg.grid.inside = ones(size(GridLoc,1),1);
% cfg.grid.unit   = 'm';
cfg.grid.resolution = 1;
cfg.grid.unit   = 'cm';
cfg.symmetry    = SymmetryConstraint;
cfg.feedback    = 'textbar';
cfg.gridsearch  = 'yes';
cfg.senstype    = 'MEG';

if exist('fminunc', 'file')
        cfg.dipfit.optimfun = 'fminunc';
else 
        cfg.dipfit.optimfun = 'fminsearch';
end
ftDipole = ft_dipolefitting(cfg, timelockeds);
% dipoles_mag_early{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
nTime = length(ftDipole.time);
switch (DipoleModel)
case 'moving'
            dipPos = cat(1, ftDipole.dip.pos)';
            dipMom = reshape(cat(2, ftDipole.dip.mom), 3, []);
            dipRv  = cat(2, ftDipole.dip.rv);
case 'regional'
            dipPos = repmat(ftDipole.dip.pos, nTime, 1)';
            dipMom = reshape(ftDipole.dip.mom, 3, []);
            dipRv  = ftDipole.dip.rv;
end
Goodness  = 1 - dipRv;
MRI = load('/Users/valery/MEG/brainstorm_db/MEG_Tommaso/anat/B1C2/subjectimage_T1.mat');
Voxels = cs_convert(MRI, 'scs', 'voxel', dipPos);
% voxinds = round(ft_warp_apply(pinv(mri.transform), mnipos));
%% make leadfields eeg and meg

% timelockeds = data

cfg = [];
cfg.elec = timelockeds{1}.elec;
cfg.grid.resolution = 1;
cfg.grid.unit = 'cm';
cfg.senstype = 'meg';
cfg.grad = timelockeds{1}.grad;
cfg.headmodel = headmodel_meg;
cfg.channel = 'megmag';

leadfield_mag = ft_prepare_leadfield(cfg, timelockeds{1});

cfg.channel = 'meggrad';

leadfield_grad = ft_prepare_leadfield(cfg, timelockeds{1});


%% dipole fits

n_events = length(events); % the events to be looped through

% the six cell arrays below are preparation for the fits

dipoles_mag_early =  cell(1, n_events);
dipoles_grad_early = cell(1, n_events);
dipoles_eeg_early =  cell(1, n_events);

dipoles_mag_late =  cell(1, n_events);
dipoles_grad_late = cell(1, n_events);
dipoles_eeg_late =  cell(1, n_events);

early_latency = [0.045 0.065]; % s
late_latency = [0.115 0.155]; % s

for event_index = 1:n_events

    cfg = [];
    cfg.gridsearch = 'yes'; %% search the grid for an optimal starting point
    cfg.grid = leadfield_mag; %% supply the grid/leadfield
    cfg.headmodel = headmodel_meg; %% supply the headmodel
    cfg.dipfit.metric = 'rv'; %% the metric to minimize (the relative residual variance: proportion of variance left unexplained by the dipole model)
    cfg.model = 'regional'; %% we assume that the dipole has a fixed position during the time points in the latency range
    cfg.senstype = 'meg'; %% sensor type
    cfg.channel = 'megmag'; %% which channels to use
    cfg.nonlinear = 'yes'; %% do a non-linear search

    % magnetometer fits
    cfg.latency = early_latency; %% specify the latency
    cfg.numdipoles = 1; %% we only expect contralateral activity
    cfg.symmetry = []; %% empty for single dipole fit
    dipoles_mag_early{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
    cfg.latency = late_latency;
    cfg.numdipoles = 2; %% we expect bilateral activity
    cfg.symmetry = 'x'; %% we expect it to be symmetrical in the x-direction (ear-to-ear)
    dipoles_mag_late{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
 
    % gradiometer fits
    cfg.channel = 'meggrad';
    cfg.grid = leadfield_grad;

    cfg.latency = early_latency;
    cfg.numdipoles = 1;
    cfg.symmetry = [];
    dipoles_grad_early{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
    cfg.latency = late_latency;
    cfg.numdipoles = 2;
    cfg.symmetry = 'x';
    dipoles_grad_late{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
    
end

disp('Done')