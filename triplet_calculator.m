function [ramp_down_triplet_avg,ramp_down_triplet_stdev,ramp_down_triplet_time]...
    = triplet_calculator(ramp_down_avg_start_stop_times,zero_avg_start_stop_times,...
    ramp_down_avg_time_length_of_chunk,ramp_down_avg_gaus_avg,ramp_down_avg_gaus_stdev,...
    zero_avg_gaus_avg,zero_avg_gaus_stdev)
    
    
    ramp_down_triplet_avg = [];
    ramp_down_triplet_stdev = [];
    ramp_down_triplet_time = [];


    for i = 1:length(ramp_down_avg_start_stop_times)



        ramp_down_start_time = ramp_down_avg_start_stop_times(i,1);

        ramp_down_stop_time = ramp_down_avg_start_stop_times(i,2);

        earlier_zero_chunk = find(ramp_down_start_time - zero_avg_start_stop_times(:,2) > 0 &...
            ramp_down_start_time - zero_avg_start_stop_times(:,2) < ramp_down_avg_time_length_of_chunk(i));

        later_zero_chunk = find(zero_avg_start_stop_times(:,1) - ramp_down_stop_time > 0 &...
            zero_avg_start_stop_times(:,1) - ramp_down_stop_time < ramp_down_avg_time_length_of_chunk(i));


        if length(earlier_zero_chunk) == 1 && length(later_zero_chunk) == 1

            triplet_avg = ramp_down_avg_gaus_avg(i) - 0.5*(zero_avg_gaus_avg(earlier_zero_chunk)+...
                zero_avg_gaus_avg(later_zero_chunk));

            triplet_stdev = sqrt(ramp_down_avg_gaus_stdev(i)^2 +zero_avg_gaus_stdev(earlier_zero_chunk)^2 +...
                zero_avg_gaus_stdev(later_zero_chunk)^2);

            triplet_time = ramp_down_start_time + 0.5*ramp_down_avg_time_length_of_chunk(i);

            ramp_down_triplet_time = [ramp_down_triplet_time triplet_time];

            ramp_down_triplet_avg = [ramp_down_triplet_avg triplet_avg];

            ramp_down_triplet_stdev = [ramp_down_triplet_stdev triplet_stdev];

        end

    end