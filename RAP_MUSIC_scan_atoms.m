function [Valmax, Indmax, Sources, exp_var_gof] = RAP_MUSIC_scan_atoms(spike, ...
    thresh, G3, channel_idx)
% -------------------------------------------------------
% RAP-MUSIC scan
% -------------------------------------------------------
% INPUTS:
%   spike -- Nch x T MEG data with spike
%   Gain -- Nch x Nsources*3 forward operator
%   G2 -- Nch x Nsources*2 forward operator in tangential plane
%   thresh -- minimal considered correlation value
%
% OUTPUT:
%   Valmax -- subcorr values for all found dipoles higher than threshold
%   IndMax -- location of this sources
%   Sources -- The time series of the source
%   exp_var_gof -- Explained variance of each iteration of RAP MUSIC
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

% 2D forward operator
[G2, ~] = G3toG2(G3, channel_idx);
Gain = G3.Gain(channel_idx,:);

[Ns, Nsrc2] = size(G2);
Nsrc = Nsrc2/2;

Valmax = [];
Indmax = [];
Sources = [];
[U,S,V] = svd(spike);
h = cumsum(diag(S)/sum(diag(S)));
n = find(h>=0.95);
corr = MUSIC_scan(G2, U(:,1:n(1)));
[valmax, indmax] = max(corr);
% indmax = 12053;
[source_ts, gof] = source_reconstruction_atom(spike, G2, indmax, channel_idx);
Sources = [Sources source_ts];
exp_var_gof = gof;

% Error estimation
% A = Gain(:,(indmax*3-2):indmax*3);
% spike_estimated =  (repmat(source_ts,1,3)\(A ))' ;
% error = 1 - norm(spike  - spike_estimated )/norm(spike)

% error = norm(spike/norm(spike) - spike_estimated/norm(spike_estimated))
try
    while valmax > 0.95
        Valmax = [Valmax, valmax];
        Indmax = [Indmax, indmax];
        A = Gain(:,(indmax*3-2):indmax*3);
        P = eye(Ns, Ns)-A*inv(A'*A)*A';
        spike_proj = P*spike;
        G_proj = P*Gain;
        Gain = G_proj;
        
        G2 = zeros(Ns,2*Nsrc);
        range = 1:2;
        for i = 1:Nsrc
            g = [G_proj(:,1+3*(i-1)) G_proj(:,2+3*(i-1)) G_proj(:,3+3*(i-1))];
            [u sv v] = svd(g);
            gt = g*v(:,1:2);
            G2(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
            range = range + 2;
        end
        
        [U,S,V] = svd(spike_proj);
        h = cumsum(diag(S)/sum(diag(S)));
        n = find(h>=0.95);
        corr = MUSIC_scan(G2, U(:,1:n(1)));
        [valmax, indmax] = max(corr);
        [source_ts, gof_n] = source_reconstruction_atom(spike_proj, G2, indmax, channel_idx);
        gof = (1-(exp_var_gof(end)))*gof_n;
        %gof = (1-sum(exp_var_gof))*gof_n;

        exp_var_gof = [exp_var_gof gof];
        Sources = [Sources source_ts];
        
    end
catch
    disp('End of the dipole fitting')
end
end


function [source_ts, gof] = source_reconstruction_atom(Data, G2, dip_ind, channel_idx)
    g = G2(:,(dip_ind*2-1):dip_ind*2); %[sensors x spikes_locationx2]
    [U S ~] = svd(Data); % U:[204 x 204], S:[204xtime]
    h = cumsum(diag(S)/sum(diag(S)));
    n = min(find(h>=0.95));
    [u s v] = svd(U(:,1:n)'*g); % project only main components, SVD them again: V is the oreintation of the source
    g_fixed = g*v(1,:)'; % fixed orientation forward model
    source_ts = Data'*g_fixed; % estimated source signal
    
    % gof 
    source_est = Data'*g;
    linvg =  inv(g'*g)*g';
    Data_est = source_est*linvg;
    % linvg*g %check, it shsource_est*linvgould be I 2 x 2
    gof = norm(Data_est-Data')/norm(Data'); % [0,1] where 1 is perfect and 0 is bad
    % figure, plot([Data_est(:,10)';  Data(10,:)]')
end



