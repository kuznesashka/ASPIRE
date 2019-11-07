function distr = all2all_prop(cluster,time_w)

time_distr = [-time_w:.01:time_w]*1000;
figure
L = length(cluster);
distr = {}
for pre = 1:L
    stamps_pre = cluster{1,pre}(2,:)'; % column
    for post = 1:L
        stamps_post = cluster{1,post}(2,:); % raw
        diff = repmat(stamps_pre,1,length(stamps_post)) - repmat(stamps_post,length(stamps_pre),1);
        diff  = diff(find(abs(diff)<time_w*1000));
        if pre < post
        subplot(L,L,(pre-1)*L+post)
        distr{pre,post}(:,1) = hist(diff,time_distr);
        bar(time_distr,distr{pre,post})
        xlim([-1 1]*100)
        ylim([0 30])
        end
    end
end
        
