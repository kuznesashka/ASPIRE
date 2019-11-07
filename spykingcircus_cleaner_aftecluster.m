function cluster = spykingcircus_cleaner_aftecluster(cluster)

time_limit = 30;


%% all templates

for cl = 1:length(cluster)
    
    %%
    matrix = cluster{1, cl};
    
    templ  = unique(matrix(end,:));
    
    matrix_new = [];
    
    for id_clus = 1:length(templ)
        
        %%
        
        GoodOnes = [];
        
        id_spikes         = find(matrix(end,:)==templ(id_clus));
        spike_templ_loc   = matrix(2,id_spikes); % stamps
        spike_gof_loc     = matrix(3,id_spikes); % goodness
        spike_templid_loc = matrix(5,id_spikes); % template
        
        if length(spike_templ_loc)>1
            
            diff_spike_ind = diff(spike_templ_loc);
            
            lowval = find(diff_spike_ind<time_limit);
            
            if numel(lowval)
            
            highval = find(diff_spike_ind>=time_limit);
            
            GoodOnes = [GoodOnes [ spike_templ_loc(:,highval+1)
                spike_gof_loc(:,highval+1)
                spike_templid_loc(:,highval+1)]];
            
            else
                
              GoodOnes = [GoodOnes [ spike_templ_loc 
                                    spike_gof_loc 
                                    spike_templid_loc ]];
            end
            
            [cl id_clus ];
            
            if lowval
                
                guide = (diff_spike_ind<time_limit);                
                Local_train = [];     
                
                for i = 1:length(guide) % 1001000111 
                    
                    if guide(i) == 1
                        ADDED = [ spike_templ_loc(i)   spike_templ_loc(i+1)
                            spike_gof_loc(i)     spike_gof_loc(i+1)
                            spike_templid_loc(i) spike_templid_loc(i+1)];
                        Local_train = [Local_train  ADDED];
                    end
                    
                    if numel(Local_train) % update GoodOnes
                        [value,pos] = max(Local_train(2,:));                        
                        GoodOnes = [GoodOnes Local_train(:,pos)];                        
                    end                    
                    Local_train = [];
                end
            end
            
        else
            
             GoodOnes = [GoodOnes [ spike_templ_loc 
                                    spike_gof_loc 
                                    spike_templid_loc ]];
        end
        
        
        [~,GoodOnes_id ] = unique(GoodOnes(1,:) );
        GoodOnes         = GoodOnes(:,GoodOnes_id);
        matrix_new       = [matrix_new GoodOnes];

    end   %id_clus
    
    [c,ia,ib] = intersect(matrix_new(1,:), matrix(2,:));
    
    cluster{1, cl} =  matrix(:,ib);
    
end %cl

 


    