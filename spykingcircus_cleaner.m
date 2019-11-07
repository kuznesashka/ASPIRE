function [spike_ind,spike_clust] = spykingcircus_cleaner(spike_ind,spike_clust)


%% all templates

time_limit = 3;
spike_loc  = sort([spike_ind spike_clust]);
diff_spike_ind = diff(spike_loc(:,1));

debug =0;
if debug
    figure
    time_distr = 1:1000
    distr = hist(diff_spike_ind,time_distr);
    bar(time_distr,distr)
end


ciao  = find(diff_spike_ind< time_limit)+1;
spike_loc(ciao,:)   = [];
spike_ind   = spike_loc(:,1);
spike_clust = spike_loc(:,2);


%% cluster level
% id_clus = unique(spike_clust)';
% 
% 
% time_limit = 10;
% 
% for templ = 1:length(id_clus)
%     
%     spike_templ_loc = spike_ind(spike_clust==id_clus(templ));
%     diff_spike_ind = diff(spike_templ_loc);
%     
% end
% 

    