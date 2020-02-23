function [IndMax, ValMax, ind_m, spikeind] = dip_fit(spike_ind, Data, G3, channels, channel_idx, t1, t2, t3, MRI, cortex)

Voxels = cs_convert(MRI, 'scs', 'voxel', cortex.Vertices);
IndMax      = zeros(spike_ind)
ValMax      = zeros(spike_ind)
ind_m       = []
spikeind    = zeros(spike_ind)
% Convert head model
ftHeadmodel = out_fieldtrip_headmodel(G3, channels, channel_idx);

% Convert data file
ftData = out_fieldtrip_data(Data, channels, channel_idx, 0);

% Generate rough grid for first estimation
% ??? GridOptions = bst_get('GridOptions_dipfit');
% ??? GridLoc = bst_sourcegrid(GridOptions, G3.SurfaceFile);
DipoleModel = 'moving';
NumDipoles = 1;
TimeWindow = [1]; % samples (ms)
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

for i = 1:length(spike_ind)
    t = spike_ind(i)

    if t1 == 0
        timelockeds.time   = [1];
        timelockeds.avg    = ftData.trial{1}(:, t);
    else
        avr1 = t - (t3-t1);
        avr2 = t - (t3-t2);
        timelockeds.time   = [1];
        timelockeds.avg    = mean(ftData.trial{1}(:, avr1:avr2),2);
    end

    timelockeds.label  = ftData.label;
    timelockeds.dimord = 'chan_time';
    timelockeds        = copyfields(ftData, timelockeds, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});
        

    ftDipole = ft_dipolefitting(cfg, timelockeds);
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
    mri_cord = cs_convert(MRI, 'scs', 'voxel', dipPos);
    %compute Euclidean distances:
    distances =  sum((Voxels - mri_cord) .^ 2, 2);
    %find the smallest distance:
    IndMax(i)   = % find(distances == min(distances));
    ValMax(i)   = Goodness
    spikeind(i) = mean([avr1 avr2])
 
end 



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