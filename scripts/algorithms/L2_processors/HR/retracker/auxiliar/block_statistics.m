function [std_block,mean_block,mean_std,mean_mean] = block_statistics(input,win_size_detrending)
    N_samples = length(input);
    num_boxes=floor(N_samples/win_size_detrending);
    std_block = zeros(1,num_boxes);
    mean_block = zeros(1,num_boxes);
    for i_box=1:(num_boxes+1)
        init_sample=max([(i_box-1)*win_size_detrending+1,1]);
        last_sample=min([(i_box-1)*win_size_detrending+win_size_detrending,N_samples]);
        std_block(1,i_box)=nanstd(detrend(input(init_sample:last_sample)));
        mean_block(1,i_box)=nanmean((input(init_sample:last_sample)));
    end

    mean_std =nanmean(std_block);
    mean_mean =nanmean(mean_block);

end