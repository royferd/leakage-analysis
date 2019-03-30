function x_avg_offset = calc_offset(x_avg,x_weight,start_offset,end_offset,offset_length )...

    num_files = length(x_avg(:,1));

    %subtract off the offsets, i.e. when power supply voltage = 0
    x_avg_offset = zeros(num_files,1);
    x_avg_offset = zeros(num_files,1);
    x_avg_offset_wt_sum = zeros(num_files,1);
    x_weight_offset_sum = zeros(num_files,1);
    x_avg_offset_end_wt = zeros(num_files,offset_length);
    x_avg_offset_start_wt = zeros(num_files,offset_length);
    x_avg_offset_start = zeros(num_files,1);
    x_avg_offset_end = zeros(num_files,1);
    x_weight_offset_start = zeros(num_files,1);
    x_weight_offset_end = zeros(num_files,1);

    for i = 1:num_files
        for j = start_offset(i,1):start_offset(i,2)
            x_avg_offset_start_wt(i,j-start_offset(i,1)+1) = x_avg(i,j)*x_weight(i,j);
        end
    end

    for i = 1:num_files
        x_avg_offset_start(i) = sum(x_avg_offset_start_wt(i,:))/sum(x_weight(i,start_offset(i,1):start_offset(i,2)));
        x_weight_offset_start(i) = (std(x_avg(i,start_offset(i,1):start_offset(i,2)),x_weight(i,start_offset(i,1):start_offset(i,2))))^(-2);
    end

    for i = 1:num_files
        for j = end_offset(i,1):end_offset(i,2)
            x_avg_offset_end_wt(i,j-end_offset(i,1)+1) = x_avg(i,j)*x_weight(i,j);
        end
    end

    for i = 1:num_files
        
        x_avg_offset_end(i) = sum(x_avg_offset_end_wt(i,:))/sum(x_weight(i,end_offset(i,1):end_offset(i,2)));
        x_weight_offset_end(i) = (std(x_avg(i,end_offset(i,1):end_offset(i,2)),x_weight(i,end_offset(i,1):end_offset(i,2))))^(-2);

        x_avg_offset_wt_sum(i) = x_avg_offset_start(i)*x_weight_offset_start(i) + x_avg_offset_end(i)*x_weight_offset_end(i);
        x_weight_offset_sum(i) = x_weight_offset_start(i)+x_weight_offset_end(i);

        x_avg_offset(i) = x_avg_offset_wt_sum(i)/x_weight_offset_sum(i);
        x_avg_offset(i) = x_avg_offset(i);

    end

end