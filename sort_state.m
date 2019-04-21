%{
sort_state takes raw leakage current and voltage data inputs and separates
it into data at high negative voltage (ramp up), zero voltage (trash), and
high positive voltage (ramp down). The details of the way data is sorted
depends on the system (e.g. bipoloar MSU or ANL), but in general it looks
for rapid changes in voltage or leakage current.
%}
function [lcm1_avg_trash_raw,lcm1_weight_avg_trash_raw,...
    lcm1_inv_weight_avg_trash_raw,time_trash,num_trash_points,...
    num_trash_chunks,lcm1_avg_ramp_down_raw,...
    lcm1_weight_avg_ramp_down_raw,lcm1_inv_weight_avg_ramp_down_raw,...
    time_ramp_down,num_ramp_down_points,num_down_chunks...
    lcm1_avg_ramp_up_raw,lcm1_weight_avg_ramp_up_raw,...
    lcm1_inv_weight_avg_ramp_up_raw,time_ramp_up,num_ramp_up_points,...
    num_up_chunks,up_chunk_array,...
    down_chunk_array,trash_chunk_array,index_up_chunk_data,index_trash_chunk_data]...
    = sort_state(vmon_avg,vmon_weight,vmon_avg_wt,vmon_weight_raw,...
    start_point,end_point,lcm1_avg_raw,lcm1_weight_raw,time,num_rows,...
    numpoints,sampling_time,EDM_sim,power_supply)

    disp('running sort_state.m');

    %%%%%%%%%%%%%%%%%%%%%%%% ramp test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%08-28-2017 put this into github. implement identical code, both
    %%polarities.

    num_files = length(vmon_avg(:,1));

    vmon_avg_avg = zeros(num_files,1);
    vmon_avg_total_stdev = zeros(num_files,1);
    vmon_ramp_high = zeros(num_files,1);
    vmon_ramp_low = zeros(num_files,1);
    vmon_ramp_high_deviation = zeros(num_files,1);
    vmon_ramp_low_deviation = zeros(num_files,1);
    
    

    % avg around 0. To separate HI/LO values from 0 values, ignore voltages
    % within pm 1 kV of vmon_avg_avg. HI = negative, LO= positive

    
    
    for i = 1:num_files
        
        length_data = length(vmon_avg);
        [vmon_avg_list,sorted_indices] = sort(vmon_avg(i,:));
        vmon_wt_list = vmon_weight(i,sorted_indices);
        vmon_wt_raw_list = vmon_weight_raw(i,sorted_indices);
        vmon_avg_wt_list = vmon_avg_wt(i,sorted_indices);
        
        % look at the bottom 15% and top 15% to find the max and min
        % voltages.
        
        minimum_voltages_index = floor(length_data*.20);
        maximum_voltages_index = floor(length_data*.80);
        
%         vmon_avg_avg(i) = mean(vmon_avg(i,start_point(i):end_point(i)));
%         vmon_avg_total_stdev(i) = std(vmon_avg(i,start_point(i):end_point(i)));
%         vmon_hi_vals = find(vmon_avg(i,:) < vmon_avg_avg(i) - 5);
%         vmon_lo_vals = find(vmon_avg(i,:) > vmon_avg_avg(i) + 5);
        
        % high negative voltage
        vmon_hi_vals = vmon_avg_list(i,1:minimum_voltages_index);
        vmon_hi_weights = vmon_wt_list(i,1:minimum_voltages_index);
        vmon_hi_weights_raw = vmon_wt_raw_list(i,1:minimum_voltages_index);
        vmon_hi_avg_wts = vmon_avg_wt_list(i,1:minimum_voltages_index);
        
        vmon_lo_vals = vmon_avg_list(maximum_voltages_index:end);
        vmon_lo_weights = vmon_wt_list(maximum_voltages_index:end);      
        vmon_lo_weights_raw = vmon_wt_raw_list(maximum_voltages_index:end);
        vmon_lo_avg_wts = vmon_avg_wt_list(maximum_voltages_index:end);
        
%         vmon_ramp_high_deviation(i) = std(vmon_avg(i,vmon_hi_vals),vmon_weight(i,vmon_hi_vals));
%         vmon_ramp_low_deviation(i) = std(vmon_avg(i,vmon_lo_vals),vmon_weight(i,vmon_lo_vals));
%         vmon_ramp_high(i) = sum(vmon_avg_wt(i,vmon_hi_vals))/sum(vmon_weight_raw(i,vmon_hi_vals));
%         vmon_ramp_low(i) = sum(vmon_avg_wt(i,vmon_lo_vals))/sum(vmon_weight_raw(i,vmon_lo_vals));  
        
%         vmon_ramp_high_deviation(i) = std(vmon_hi_vals,vmon_hi_weights);
%         vmon_ramp_low_deviation(i) = std(vmon_lo_vals,vmon_lo_weights);
                
        
        vmon_ramp_high(i) = sum(vmon_hi_avg_wts)/sum(vmon_hi_weights_raw);
        vmon_ramp_low(i) = sum(vmon_lo_avg_wts)/sum(vmon_lo_weights_raw); 

    end

    fprintf('high negative voltage = %f +/- %f kV \n',vmon_ramp_high(i),vmon_ramp_high_deviation(i));
    fprintf('high positive voltage = %f +/- %f kV \n',vmon_ramp_low(i),vmon_ramp_low_deviation(i));
    %ramp_deviation = 0.05;

    time_ramp_up_pass = zeros(num_files,num_rows,1);
    time_ramp_down_pass = zeros(num_files,num_rows,1);

    lcm1_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
    lcm1_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);

    lcm1_weight_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
    lcm1_weight_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);

    time_trash_pass = zeros(num_files,num_rows,1);    
    lcm1_avg_trash_raw_pass = zeros(num_files,num_rows,1);    
    lcm1_weight_avg_trash_raw_pass = zeros(num_files,num_rows,1);



    %Chunk array(file number, chunk number, chunk boundary element): For each 
    %'chunk,' I define a chunk number (first, second, etc) and the elements 
    %that mark the boundaries of each chunk. For example, the second chunk for
    %the first file is given as up_chunk_array(1,1,2). The boundaries of the
    %second chunk would be up_chunk_array(1,2,2) + 1 and up_chunk_array(1,2,3).

    up_chunk_array_pass = zeros(num_files,2,num_rows); 
    down_chunk_array_pass = zeros(num_files,2,num_rows);
    trash_chunk_array_pass = zeros(num_files,2,num_rows);

    up_array_count = zeros(num_files,1);
    down_array_count = zeros(num_files,1);
    trash_array_count = zeros(num_files,1);

    %I start them at one out of necessity and subtract one after the sorting
    num_ramp_up_points = ones(num_files,1);
    num_ramp_down_points = ones(num_files,1);
    num_trash_points = ones(num_files,1);

    count_up_chunks = zeros(num_files,1);
    count_down_chunks = zeros(num_files,1);
    count_trash_chunks = zeros(num_files,1);

    index_up_chunk_data_pass = zeros(num_files,1);
    index_down_chunk_data_pass = zeros(num_files,1);
    index_trash_data_pass = zeros(num_files,1);
    
    % After sorting the up/down/trash data, trim a few more points off the
    % chunks to make double sure that charging data is excluded from the
    % ramp segments.
    extra = 0;

%    if (power_supply == 0) || (power_supply == 1)
    if (power_supply == 0) || (power_supply == 1) || (power_supply == 3) || power_supply == 4
        
        extra = 0;
        for i = 1:num_files
            count_up_chunks(i) = 0;
            count_down_chunks(i) = 0;
            count_trash_chunks(i) = 0;

            for j = start_point(i):numpoints(i)           

%                 if abs(vmon_avg(i,j) - vmon_ramp_high(i)) < ...
%                         3*vmon_ramp_high_deviation(i) && ...
%                         vmon_avg(i,j) < ...
%                         vmon_ramp_high(i) + vmon_ramp_high_deviation(i)
                    
                if abs(vmon_avg(i,j) - vmon_ramp_high(i)) < ...
                        3*vmon_ramp_high_deviation(i) && ...
                        vmon_avg(i,j) < ...
                        vmon_ramp_high(i) + 5*vmon_ramp_high_deviation(i)

                    index_up_chunk_data_pass(i,end+1) = j;
                    num_ramp_up_points(i) = num_ramp_up_points(i) + 1;

                    lcm1_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_avg_raw(i,j);
                    lcm1_weight_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_weight_raw(i,j);

                    time_ramp_up_pass(i,num_ramp_up_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                    up_array_count(i) = up_array_count(i)+1;

                    if (time_ramp_up_pass(i,num_ramp_up_points(i)) - ...
                            time_ramp_up_pass(i,num_ramp_up_points(i)-1) > sampling_time)                

                        count_up_chunks(i) = count_up_chunks(i) + 1;
                        up_chunk_array_pass(i,1, count_up_chunks(i)) = count_up_chunks(i);
                        up_chunk_array_pass(i,2, count_up_chunks(i)) = up_array_count(i) - 1;

                    end

%                 elseif (abs(vmon_avg(i,j) - vmon_ramp_low(i)) < 3*vmon_ramp_low_deviation(i) && ...
%                         vmon_avg(i,j) > vmon_ramp_low(i) - vmon_ramp_low_deviation(i))

                elseif (abs(vmon_avg(i,j) - vmon_ramp_low(i)) < 3*vmon_ramp_low_deviation(i) && ...
                        vmon_avg(i,j) > vmon_ramp_low(i) - 5*vmon_ramp_low_deviation(i))
                    
                    index_down_chunk_data_pass(i,end+1) = j;
                    
                    num_ramp_down_points(i) = num_ramp_down_points(i) + 1;

                    lcm1_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = lcm1_avg_raw(i,j);
                    
                    lcm1_weight_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = lcm1_weight_raw(i,j);

                    time_ramp_down_pass(i,num_ramp_down_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                
                    down_array_count(i) = down_array_count(i) + 1;

                    if (time_ramp_down_pass(i,num_ramp_down_points(i)) - ...
                            time_ramp_down_pass(i,num_ramp_down_points(i)-1) > sampling_time)                

                        count_down_chunks(i) = count_down_chunks(i) + 1;
                        
                        down_chunk_array_pass(i,1, count_down_chunks(i)) = count_down_chunks(i);
                        
                        down_chunk_array_pass(i,2, count_down_chunks(i)) = down_array_count(i)-1;

                    end

                else

                    index_trash_data_pass(i,end+1) = j;
                  
                    num_trash_points(i) = num_trash_points(i) + 1;

                    lcm1_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_avg_raw(i,j);
          
                    lcm1_weight_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_weight_raw(i,j);

                    time_trash_pass(i,num_trash_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
             
                    trash_array_count(i) = trash_array_count(i) + 1;

                    if (time_trash_pass(i,num_trash_points(i)) - ...
                            time_trash_pass(i,num_trash_points(i)-1) > sampling_time)                

                        count_trash_chunks(i) = count_trash_chunks(i) + 1;
                   
                        trash_chunk_array_pass(i,1, count_trash_chunks(i)) = count_trash_chunks(i);
                     
                        trash_chunk_array_pass(i,2, count_trash_chunks(i)) = trash_array_count(i)-1;

                    end

                end

            end

        end

    elseif power_supply == 2
        
        extra = 0;

        for i = 1:num_files

            count_up_chunks(i) = 0;
            
            count_down_chunks(i) = 0;
            
            count_trash_chunks(i) = 0;

            for j = start_point(i):numpoints(i)

                if ((lcm1_avg(i,j) - lcm1_avg_offset(i) < +6*1e3) &&  ...
                        (lcm1_avg(i,j) - lcm1_avg_offset(i) > -6*1e3))

                    num_ramp_up_points(i) = num_ramp_up_points(i) + 1;

                    lcm1_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_avg_raw(i,j);

                    lcm1_weight_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = ...
                        lcm1_weight_raw(i,j);

                    time_ramp_up_pass(i,num_ramp_up_points(i)) = time(i,j); 
                    
                    up_array_count(i) = up_array_count(i)+1;

                    if (time_ramp_up_pass(i,num_ramp_up_points(i)) - ...
                            time_ramp_up_pass(i,num_ramp_up_points(i)-1) > 5*sampling_time)

                        count_up_chunks(i) = count_up_chunks(i) + 1;
                  
                        up_chunk_array_pass(i,1, count_up_chunks(i)) = count_up_chunks(i);
                        
                        up_chunk_array_pass(i,2, count_up_chunks(i)) = up_array_count(i) - 1;

                    end

                else

                    index_trash_data_pass(i,end+1) = j;
                    
                     num_trash_points(i) = num_trash_points(i) + 1;

                    lcm1_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_avg_raw(i,j);
                    
                    lcm1_weight_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_weight_raw(i,j);

                     time_trash_pass(i,num_trash_points(i)) = time(i,j);
                    
                    trash_array_count(i) = trash_array_count(i) + 1;

                    if (time_trash_pass(i,num_trash_points(i)) - ...
                            time_trash_pass(i,num_trash_points(i)-1) > 5*sampling_time)

                        count_trash_chunks(i) = count_trash_chunks(i) + 1;
                        
                        trash_chunk_array_pass(i,1, count_trash_chunks(i)) = count_trash_chunks(i);
                        
                        trash_chunk_array_pass(i,2, count_trash_chunks(i)) = trash_array_count(i)-1;

                    end

                end

            end

        end

    end

    for i = 1:num_files

        num_ramp_up_points(i) = num_ramp_up_points(i) - 1;
        
        num_ramp_down_points(i) = num_ramp_down_points(i) -1;
        
        num_trash_points(i) = num_trash_points(i) -1;

    end    

    num_up_chunks = zeros(num_files,1);

    num_down_chunks = zeros(num_files,1);

    num_trash_chunks = zeros(num_files,1);

    for i = 1:num_files
        %count to third to last trash chunk
        num_trash_chunks(i) = uint8(count_trash_chunks(i) -4);

        %count to second to last down chunk
        num_down_chunks(i) = uint8(count_down_chunks(i) - 2);

        %count up to 2nd to last up chunk
        num_up_chunks(i) = uint8(count_up_chunks(i) - 2);

        if EDM_sim == 1

%            file_number_str = sprintf('%d',i);
%             up_chunks_str = sprintf('%d',num_up_chunks(i));
%             down_chunks_str = sprintf('%d',num_down_chunks);
%             trash_chunks_str = sprintf('%d',num_trash_chunks);
%             mes1 = ['file number: ',file_number_str];
%             mes2 = ['number of up_chunks: ', up_chunks_str];
%             mes3 = ['number of  down chunks: ', down_chunks_str];
%             mes4 = ['number of trash chunks: ', trash_chunks_str];
%             disp(mes1)
%             disp(mes2)
%             disp(mes3)
%             disp(mes4)
%             disp('the number of trash chunks should be the sum of up and down chunks')                                                  

        end

    end
    
    up_chunk_array_pass_pass = zeros(num_files,1); % this triggers an error
    % if more than one file is being read
    
    down_chunk_array = zeros(num_files,2,1);
    
    trash_chunk_array_pass_pass = zeros(num_files,2,1);

    time_ramp_up = zeros(num_files,num_ramp_up_points,1);
%    time_ramp_up = [];    
    
    time_ramp_down = zeros(num_files,num_ramp_down_points,1);
%    time_ramp_down = [];    
    
%    time_trash = zeros(num_files,num_trash_points,1);
    time_trash = [];

    lcm1_avg_ramp_up_raw = zeros(num_files,num_ramp_up_points,1);
%    lcm1_avg_ramp_up_raw = [];

    lcm1_avg_ramp_down_raw = zeros(num_files,num_ramp_down_points,1);
%    lcm1_avg_ramp_down_raw = [];

%    lcm1_avg_trash_raw = zeros(num_files,num_trash_points,1);
    lcm1_avg_trash_raw = [];

    % 3/11/2018 inv_weight is just stdev. I fucked up variable names so I have to use
    %this nomenclature to avoid a lot of work I don't want to do right now.

%    lcm1_weight_avg_ramp_up_raw = zeros(num_files,num_ramp_up_points,1);
    lcm1_weight_avg_ramp_up_raw = [];
    
    lcm1_weight_avg_ramp_down_raw = zeros(num_files,num_ramp_down_points,1);
%    lcm1_weight_avg_ramp_down_raw = [];

%    lcm1_weight_avg_trash_raw = zeros(num_files,num_trash_points,1);
    lcm1_weight_avg_trash_raw = [];    
    
%    lcm1_inv_weight_avg_trash_raw = zeros(num_files,num_trash_points,1);
    lcm1_inv_weight_avg_trash_raw = [];

    index_up_chunk_data = index_up_chunk_data_pass(:,2:end);
    
    index_down_chunk_data = index_down_chunk_data_pass(:,2:end);
    
    index_trash_chunk_data = index_trash_data_pass(:,2:end);

    for i = 1:num_files
        
        for j = 1:num_ramp_up_points(i)

            lcm1_avg_ramp_up_raw(i,j) = lcm1_avg_ramp_up_raw_pass(i,j+1);
            lcm1_weight_avg_ramp_up_raw(i,j) = lcm1_weight_avg_ramp_up_raw_pass(i,j+1);
            lcm1_inv_weight_avg_ramp_up_raw(i,j) = (lcm1_weight_avg_ramp_up_raw(i,j))^(-1/2);
            
            time_ramp_up(i,j) = time_ramp_up_pass(i,j+1);

        end

       

        for j = 1:num_ramp_down_points(i)

            lcm1_avg_ramp_down_raw(i,j) = lcm1_avg_ramp_down_raw_pass(i,j+1);
            lcm1_weight_avg_ramp_down_raw(i,j) = lcm1_weight_avg_ramp_down_raw_pass(i,j+1);
            lcm1_inv_weight_avg_ramp_down_raw(i,j) = (lcm1_weight_avg_ramp_down_raw(i,j))^(-1/2);

            time_ramp_down(i,j) = time_ramp_down_pass(i,j+1);


        end

% 
        for j = 1:num_trash_points(i)

            lcm1_avg_trash_raw(i,j) = lcm1_avg_trash_raw_pass(i,j+1);
      
            lcm1_weight_avg_trash_raw(i,j) = lcm1_weight_avg_trash_raw_pass(i,j+1);
            
            lcm1_inv_weight_avg_trash_raw(i,j) = (lcm1_weight_avg_trash_raw(i,j))^(-1/2);
     
            time_trash(i,j) = time_trash_pass(i,j+1);

        end
        
        %start at 2nd down chunk
        for j =1:num_down_chunks(i) + 1

            down_chunk_array(i,1,j) = down_chunk_array_pass(i,2,j) + extra;
            
            down_chunk_array(i,2,j) = down_chunk_array_pass(i,2,j+1) - extra;
            
            in_this_chunk = down_chunk_array(i,2,j) - down_chunk_array(i,1,j) -1;

        end
        
        %start at 2nd up chunk
        for j =1:num_up_chunks(i) + 1
            
            lookfor_start = up_chunk_array_pass(i,2,j) + extra;
            
            if time_ramp_up(i,lookfor_start+1) > time_ramp_down(i,down_chunk_array(i,1,1)+1)
                
                up_chunk_begin = j;
                
                break;
                
            end
            
        end
        
        
        % now look for last up chunk, which should occur after the last
        % down chunk. 
        for j = up_chunk_begin:num_up_chunks(i) + 1
            
            lookfor_end = down_chunk_array_pass(i,2,j+1) - extra;
            
            fprintf('up chunk number = %d \n',j);
            fprintf('lookfor_end = %d \n',down_chunk_array_pass(i,2,j+1) - extra);
            fprintf('length of time_ramp_up = %d \n',length(time_ramp_up));

            if time_ramp_up(i,lookfor_end) > time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i)))
                
                up_chunk_end = j;
                
                break;
                
            end
            
        end
        
        % If this does not happen, e.g. if the simulation isn't allowed to 
        % finish its last cycle, then reduce the number of down chunks by
        % the appropriate number of cycles.
        if time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i))) > time_ramp_up(i,lookfor_end)
            
            up_chunk_end = num_up_chunks(i)+1;
                
            for j = 1:num_down_chunks(i)-1
                
                if time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i)-j)) < time_ramp_up(i,lookfor_end)
                    
                    num_down_chunks(i) = num_down_chunks(i) - j;
                    
                    break;
                    
                end
                
            end
            
        end
        
        num_up_chunks(i) = up_chunk_end - up_chunk_begin + 1;
        
        for j = up_chunk_begin:num_up_chunks(i)+1

            up_chunk_array_pass_pass(i,1,end+1) = up_chunk_array_pass(i,2,j) + extra;
            
            up_chunk_array_pass_pass(i,2,end) = up_chunk_array_pass(i,2,j+1) - extra;

            in_this_chunk = up_chunk_array_pass_pass(i,2,j) - up_chunk_array_pass_pass(i,1,j) -1;

        end
        

        % first trash chunk starts before first down chunk. Last trash chunk
        % ends after last up chunk.
        for j =1:num_trash_chunks(i) + 1
            
            lookfor_start = trash_chunk_array_pass(i,2,j+1) + extra;
            
            if time_trash(i,lookfor_start+1) > time_ramp_down(i,down_chunk_array(i,1,1)+1) && j > 1
                
                trash_chunk_begin = j-1;
                
                break;
                
            end
            
        end
        
        % now look for last trash chunk, which should occur after the last
        % up chunk. 
        for j = trash_chunk_begin:num_trash_chunks(i) + 1
            
            lookfor_end = trash_chunk_array_pass(i,2,j+1) - extra;
            
            if time_trash(i,lookfor_end) > time_ramp_up(i,up_chunk_array_pass_pass(i,2,num_up_chunks(i)+1))
                
                trash_chunk_end = j;
                
                break;
                
            end
            
        end      
        
        % If this does not happen, e.g. if the simulation isn't allowed to 
        % finish its last cycle, then reduce the number of up and down 
        % chunks by the appropriate number of cycles.
        if time_ramp_up(i,up_chunk_array_pass_pass(i,2,num_up_chunks(i)+1)) > ...
                time_trash(i,lookfor_end)
            
            trash_chunk_end = num_trash_chunks(i)+1;
                
            for j = 1:num_up_chunks(i)-1
                
                if time_ramp_up(i,up_chunk_array_pass_pass(i,2,num_up_chunks(i)+1-j)) < ...
                        time_trash(i,lookfor_end)
                    
                    num_up_chunks(i) = num_up_chunks(i) - j;
                    num_down_chunks(i) = num_down_chunks(i) - j;
                    
                    break;
                    
                end
                
            end
            
        end
        
        num_trash_chunks(i) = trash_chunk_end - trash_chunk_begin + 1;       
        
        
        for j =trash_chunk_begin:num_trash_chunks(i)+1
            
            trash_chunk_array_pass_pass(i,1,end+1) = trash_chunk_array_pass(i,2,j+1) + extra;
            
            trash_chunk_array_pass_pass(i,2,end) = trash_chunk_array_pass(i,2,j+2) + extra;
            
            in_this_chunk = trash_chunk_array_pass_pass(i,2,end) - trash_chunk_array_pass_pass(i,1,end) - 1;

        end
        
        
%         %start at third trash chunk
%         for j =1:num_trash_chunks(i) + 1
%             
%             trash_chunk_array(i,1,j) = trash_chunk_array_pass(i,2,j+1);
%             
%             trash_chunk_array(i,2,j) = trash_chunk_array_pass(i,2,j+2);
%             
%             in_this_chunk = trash_chunk_array(i,2,j) - trash_chunk_array(i,1,j) - 1;
% 
%         end
        
        


%         %start at 2nd down chunk
%         for j =1:num_down_chunks(i) + 1
% 
%             down_chunk_array(i,1,j) = down_chunk_array_pass(i,2,j) + extra;
%             
%             down_chunk_array(i,2,j) = down_chunk_array_pass(i,2,j+1) - extra;
%             
%             in_this_chunk = down_chunk_array(i,2,j) - down_chunk_array(i,1,j) -1;
% 
%         end
        

        
%         %start at 2nd up chunk
%         for j =1:num_up_chunks(i) + 1
% 
%             up_chunk_array(i,1,j) = up_chunk_array_pass(i,2,j) + extra;
%             
%             up_chunk_array(i,2,j) = up_chunk_array_pass(i,2,j+1) - extra;
% 
%             in_this_chunk = up_chunk_array(i,2,j) - up_chunk_array(i,1,j) -1;
% 
%         end
        

        
        
%         for j = 1:num_ramp_up_points(i)
% 
%             lcm1_avg_ramp_up_raw(i,j) = lcm1_avg_ramp_up_raw_pass(i,j+1);
%             lcm1_weight_avg_ramp_up_raw(i,j) = lcm1_weight_avg_ramp_up_raw_pass(i,j+1);
%             lcm1_inv_weight_avg_ramp_up_raw(i,j) = (lcm1_weight_avg_ramp_up_raw(i,j))^(-1/2);
%             
%             time_ramp_up(i,j) = time_ramp_up_pass(i,j+1);
% 
%         end
% 
%        
% 
%         for j = 1:num_ramp_down_points(i)
% 
%             lcm1_avg_ramp_down_raw(i,j) = lcm1_avg_ramp_down_raw_pass(i,j+1);
%             lcm1_weight_avg_ramp_down_raw(i,j) = lcm1_weight_avg_ramp_down_raw_pass(i,j+1);
%             lcm1_inv_weight_avg_ramp_down_raw(i,j) = (lcm1_weight_avg_ramp_down_raw(i,j))^(-1/2);
% 
%             time_ramp_down(i,j) = time_ramp_down_pass(i,j+1);
% 
% 
%         end
% 
% % 
%         for j = 1:num_trash_points(i)
% 
%             lcm1_avg_trash_raw(i,j) = lcm1_avg_trash_raw_pass(i,j+1);
%       
%             lcm1_weight_avg_trash_raw(i,j) = lcm1_weight_avg_trash_raw_pass(i,j+1);
%             
%             lcm1_inv_weight_avg_trash_raw(i,j) = (lcm1_weight_avg_trash_raw(i,j))^(-1/2);
%      
%             time_trash(i,j) = time_trash_pass(i,j+1);
% 
%         end

%         num_ramp_down_points(i) = length(lcm1_avg_ramp_down_raw);
%         num_ramp_up_points(i) = length(lcm1_avg_ramp_up_raw);      
%         lcm1_inv_weight_avg_ramp_up_raw = (lcm1_weight_avg_ramp_up_raw(:,:)).^(-1/2); 
%         lcm1_inv_weight_avg_ramp_down_raw = (lcm1_weight_avg_ramp_down_raw(:,:)).^(-1/2);
       
%        lcm1_inv_weight_avg_trash_raw = (lcm1_weight_avg_trash_raw(:,:)).^(-1/2);
% % %         
%          num_trash_points(i) = length(lcm1_avg_trash_raw);  

    end
    
    up_chunk_array = up_chunk_array_pass_pass(:,:,2:end);
    trash_chunk_array = trash_chunk_array_pass_pass(:,:,2:end);
    
    fprintf('file number: %d \n',i);
    fprintf('number of up chunks (-V) = %d \n', num_up_chunks(i));
    fprintf('number of  down chunks (+V) = %d \n', num_down_chunks(i));
    fprintf('number of trash chunks (0V) = %d \n', num_trash_chunks(i));
    disp('the number of trash chunks should be the sum of up and down chunks') 
    
end