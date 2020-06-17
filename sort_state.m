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

    program = which('sort_state','-all');
    fprintf('running %s \n', program{:});

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
        
        vmon_ramp_high_deviation(i) = 0.05*abs(vmon_ramp_high(i));
        vmon_ramp_low_deviation(i) = 0.05*abs(vmon_ramp_low(i));

    end

    fprintf('high negative voltage = %.1f +/- %.1f kV \n',vmon_ramp_high(i),vmon_ramp_high_deviation(i));
    fprintf('high positive voltage = %.1f +/- %.1f kV \n',vmon_ramp_low(i),vmon_ramp_low_deviation(i));
    %ramp_deviation = 0.05;

%     time_ramp_up_pass = zeros(num_files,num_rows,1);
%     time_ramp_down_pass = zeros(num_files,num_rows,1);
    
    time_ramp_up_pass = zeros(2,num_rows);
    time_ramp_down_pass = zeros(2,num_rows);

    lcm1_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
    lcm1_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);

    lcm1_weight_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
    lcm1_weight_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);

%     time_trash_pass = zeros(num_files,num_rows,1); 
    
    time_trash_pass = zeros(2,num_rows);
    
    % the number of points of points back in time to check to see if a new
    % chunk has been reached
    
    switch_setting_pt_buffer = 3;
    
    lcm1_avg_trash_raw_pass = zeros(num_files,num_rows,1);    
    lcm1_weight_avg_trash_raw_pass = zeros(num_files,num_rows,1);



    %Chunk array(file number, chunk number, chunk boundary element): For each 
    %'chunk,' I define a chunk number (first, second, etc) and the elements 
    %that mark the boundaries of each chunk. For example, the second chunk for
    %the first file is given as up_chunk_array(1,1,2). The boundaries of the
    %second chunk would be up_chunk_array(1,2,2) + 1 and up_chunk_array(1,2,3).

%     up_chunk_array_pass = zeros(num_files,2,num_rows);
    up_chunk_array_pass = zeros(num_files,2,1);
    
    
    down_chunk_array_pass = zeros(num_files,2,num_rows);
    
%     trash_chunk_array_pass = zeros(num_files,2,num_rows);
    
    trash_chunk_array_pass = zeros(num_files,2,1);

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
                        vmon_ramp_high(i) + vmon_ramp_high_deviation(i)

                    index_up_chunk_data_pass(i,end+1) = j;
                    
                    num_ramp_up_points(i) = num_ramp_up_points(i) + 1;

                    lcm1_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_avg_raw(i,j);
                    
                    lcm1_weight_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_weight_raw(i,j);

                    time_ramp_up_pass(1,num_ramp_up_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                    
                    time_ramp_up_pass(2,num_ramp_up_points(i)) = j;
                
                    up_array_count(i) = up_array_count(i)+1;

%                     if (time_ramp_up_pass(i,num_ramp_up_points(i)) - ...
%                             time_ramp_up_pass(i,num_ramp_up_points(i)-1) > sampling_time)
                        
                    if (time_ramp_up_pass(2,num_ramp_up_points(i)) - ...
                            time_ramp_up_pass(2,num_ramp_up_points(i)-1) > switch_setting_pt_buffer) 

                        count_up_chunks(i) = count_up_chunks(i) + 1;
                        
%                         up_chunk_array_pass(i,1, count_up_chunks(i)) = count_up_chunks(i);
%                         up_chunk_array_pass(i,2, count_up_chunks(i)) = up_array_count(i) - 1;
                        
                        up_chunk_array_pass(i,:,count_up_chunks(i)) = ...
                            [count_up_chunks(i) up_array_count(i) - 1];

                    end

%                 elseif (abs(vmon_avg(i,j) - vmon_ramp_low(i)) < 3*vmon_ramp_low_deviation(i) && ...
%                         vmon_avg(i,j) > vmon_ramp_low(i) - vmon_ramp_low_deviation(i))

                elseif (abs(vmon_avg(i,j) - vmon_ramp_low(i)) < 3*vmon_ramp_low_deviation(i) && ...
                        vmon_avg(i,j) > vmon_ramp_low(i) - vmon_ramp_low_deviation(i))
                    
                    index_down_chunk_data_pass(i,end+1) = j;
                    
                    num_ramp_down_points(i) = num_ramp_down_points(i) + 1;

                    lcm1_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = lcm1_avg_raw(i,j);
                    
                    lcm1_weight_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = lcm1_weight_raw(i,j);

                    time_ramp_down_pass(1,num_ramp_down_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                    
                    time_ramp_down_pass(2,num_ramp_down_points(i)) = j;
                
                    down_array_count(i) = down_array_count(i) + 1;

%                     if (time_ramp_down_pass(i,num_ramp_down_points(i)) - ...
%                             time_ramp_down_pass(i,num_ramp_down_points(i)-1) > sampling_time)    
                        
                    if (time_ramp_down_pass(2,num_ramp_down_points(i)) - ...
                            time_ramp_down_pass(2,num_ramp_down_points(i)-1) > switch_setting_pt_buffer)

                        count_down_chunks(i) = count_down_chunks(i) + 1;
                        
                        down_chunk_array_pass(i,1, count_down_chunks(i)) = count_down_chunks(i);
                        
                        down_chunk_array_pass(i,2, count_down_chunks(i)) = down_array_count(i)-1;

                    end

                else

                    index_trash_data_pass(i,end+1) = j;
                  
                    num_trash_points(i) = num_trash_points(i) + 1;

                    lcm1_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_avg_raw(i,j);
          
                    lcm1_weight_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_weight_raw(i,j);

                    time_trash_pass(1,num_trash_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                    
                    time_trash_pass(2,num_trash_points(i)) = j;
             
                    trash_array_count(i) = trash_array_count(i) + 1;

%                     if (time_trash_pass(i,num_trash_points(i)) - ...
%                             time_trash_pass(i,num_trash_points(i)-1) > sampling_time)                

                    if (time_trash_pass(2,num_trash_points(i)) - ...
                            time_trash_pass(2,num_trash_points(i)-1) > switch_setting_pt_buffer) 

                        count_trash_chunks(i) = count_trash_chunks(i) + 1;
                   
%                         trash_chunk_array_pass(i,1, count_trash_chunks(i)) = count_trash_chunks(i);
%                      
%                         trash_chunk_array_pass(i,2, count_trash_chunks(i)) = trash_array_count(i)-1;
                        
                        trash_chunk_array_pass(i,:,count_trash_chunks(i)) = ...
                            [count_trash_chunks(i) trash_array_count(i) - 1];

                    end

                end

            end
            
            % we know the last point in the ramp_up dataset is part of the
            % last chunk.
            up_chunk_array_pass(i,:,count_up_chunks(i)+1) = ...
                [count_up_chunks(i)+1 num_ramp_up_points-1];
            
            trash_chunk_array_pass(i,:,count_trash_chunks(i)+1) = ...
            [count_trash_chunks(i)+1 num_trash_points-1];

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

                    time_ramp_up_pass(1,num_ramp_up_points(i)) = time(i,j); 
                    
                    time_ramp_up_pass(2,num_ramp_up_points(i)) = j;
                    
                    up_array_count(i) = up_array_count(i)+1;

%                     if (time_ramp_up_pass(i,num_ramp_up_points(i)) - ...
%                             time_ramp_up_pass(i,num_ramp_up_points(i)-1) > 5*sampling_time)

                    if (time_ramp_up_pass(2,num_ramp_up_points(i)) - ...
                            time_ramp_up_pass(2,num_ramp_up_points(i)-1) > switch_setting_pt_buffer)
                        
                        count_up_chunks(i) = count_up_chunks(i) + 1;
                  
                        up_chunk_array_pass(i,1, count_up_chunks(i)) = count_up_chunks(i);
                        
                        up_chunk_array_pass(i,2, count_up_chunks(i)) = up_array_count(i) - 1;

                    end

                else
                    

                    index_trash_data_pass(i,end+1) = j;
                    
                    num_trash_points(i) = num_trash_points(i) + 1;

                    lcm1_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_avg_raw(i,j);
                    
                    lcm1_weight_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_weight_raw(i,j);

                    time_trash_pass(1,num_trash_points(i)) = time(i,j);
                     
                    time_trash_pass(2,num_trash_points(i)) = j;
                    
                    trash_array_count(i) = trash_array_count(i) + 1;

%                     if (time_trash_pass(i,num_trash_points(i)) - ...
%                             time_trash_pass(i,num_trash_points(i)-1) > 5*sampling_time)
                        
                    if (time_trash_pass(2,num_trash_points(i)) - ...
                            time_trash_pass(2,num_trash_points(i)-1) > switch_setting_pt_buffer)

                        count_trash_chunks(i) = count_trash_chunks(i) + 1;
                        
%                         trash_chunk_array_pass(i,1, count_trash_chunks(i)) = count_trash_chunks(i);
%                         
%                         trash_chunk_array_pass(i,2, count_trash_chunks(i)) = trash_array_count(i)-1;
                        
                        trash_chunk_array_pass(i,:,count_trash_chunks(i)) = ...
                            [count_trash_chunks(i) trash_array_count(i) - 1];

                    end

                end

            end

        end
        
        % we know the last point in the ramp_up dataset is part of the
        % last chunk.
        up_chunk_array_pass(i,:,count_up_chunks(i)+1) = [count_up_chunks(i)+1 num_ramp_up_points-1];
        
        trash_chunk_array_pass(i,:,count_trash_chunks(i)+1) = ...
            [count_trash_chunks(i)+1 num_trash_points-1];
        
    end

    num_up_chunks = zeros(num_files,1);

    num_down_chunks = zeros(num_files,1);

    num_trash_chunks = zeros(num_files,1);
    
    for i = 1:num_files

        num_ramp_up_points(i) = num_ramp_up_points(i) - 1;
        
        num_ramp_down_points(i) = num_ramp_down_points(i) -1;
        
        num_trash_points(i) = num_trash_points(i) -1;

        %count to third to last trash chunk
%        num_trash_chunks(i) = uint8(count_trash_chunks(i) -4);
        
%         num_trash_chunks(i) = count_trash_chunks(i) -4;
        
        num_trash_chunks(i) = count_trash_chunks(i);
        
%        fprintf('counted %d trash chunks \n',count_trash_chunks(i));

        %count to second to last down chunk
%        num_down_chunks(i) = uint8(count_down_chunks(i) - 2);
        num_down_chunks(i) = count_down_chunks(i) - 2;
%        fprintf('counted %d down chunks \n',count_down_chunks(i));


%         num_up_chunks(i) = count_up_chunks(i) - 2;
        
        num_up_chunks(i) = count_up_chunks(i);



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
            
            time_ramp_up(i,j) = time_ramp_up_pass(1,j+1);

        end

       

        for j = 1:num_ramp_down_points(i)

            lcm1_avg_ramp_down_raw(i,j) = lcm1_avg_ramp_down_raw_pass(i,j+1);
            lcm1_weight_avg_ramp_down_raw(i,j) = lcm1_weight_avg_ramp_down_raw_pass(i,j+1);
            lcm1_inv_weight_avg_ramp_down_raw(i,j) = (lcm1_weight_avg_ramp_down_raw(i,j))^(-1/2);

            time_ramp_down(i,j) = time_ramp_down_pass(1,j+1);


        end

% 
        for j = 1:num_trash_points(i)

            lcm1_avg_trash_raw(i,j) = lcm1_avg_trash_raw_pass(i,j+1);
      
            lcm1_weight_avg_trash_raw(i,j) = lcm1_weight_avg_trash_raw_pass(i,j+1);
            
            lcm1_inv_weight_avg_trash_raw(i,j) = (lcm1_weight_avg_trash_raw(i,j))^(-1/2);
     
            time_trash(i,j) = time_trash_pass(1,j+1);

        end
        
        %start at 2nd down chunk
        for j =1:num_down_chunks(i) + 1

            down_chunk_array(i,1,j) = down_chunk_array_pass(i,2,j) + extra;
            
            down_chunk_array(i,2,j) = down_chunk_array_pass(i,2,j+1) - extra;
            
            in_this_chunk = down_chunk_array(i,2,j) - down_chunk_array(i,1,j) -1;

        end
        
        num_down_chunks(i) = length(down_chunk_array(i,1,:));
        
        up_chunk_begin = 1;
        
        %start at 2nd up chunk
%         for j =1:num_up_chunks(i) + 1
            
        for j =1:length(up_chunk_array_pass(i,2,:))
            
            lookfor_start = up_chunk_array_pass(i,2,j) + extra;
            
            
            
%             fprintf('num_up_chunks = %d \n',j);
%             fprintf('lookfor_start = up_chunk_array_pass(i,2,j)+extra = %d + %d \n',...
%                up_chunk_array_pass(i,2,j),extra);
            
            if time_ramp_up(i,lookfor_start+1) > time_ramp_down(i,down_chunk_array(i,1,1)+1)
                
                up_chunk_begin = j;
%                 fprintf('up_chunk_begin = %d \n',j);
                
                break;
                
            end
            
        end
        
%         time_ramp_down(1,down_chunk_array(1,1,end))
%         up_chunk_array_pass(1,:,1:50)
        
        % now look for last up chunk, which should occur after the last
        % down chunk. 
        
%         length(up_chunk_array_pass)
%         
%         time_ramp_up(1,up_chunk_array_pass(1,2,end-5:end))
%         
%         num_up_chunks
        
%         for j = up_chunk_begin:num_up_chunks(i) + 1
        for j = up_chunk_begin:num_up_chunks(i)
            
            
            if up_chunk_array_pass(i,2,j) + 1 < length(time_ramp_up)
                
                lookfor_end = up_chunk_array_pass(i,2,j) + 1;
                
            else
                
                lookfor_end = length(time_ramp_up);
                
%                 disp('Possible issue finding the last up chunk. Check that chunks are properly identified');
                
            end
             
            up_chunk_end = j+1;
            
%             fprintf('up chunk number = %d \n',j);
%             fprintf('lookfor_end = %d \n',up_chunk_array_pass(i,2,j) +1);
%             fprintf('length of time_ramp_up = %d \n',length(time_ramp_up));

%             time_ramp_up(i,lookfor_end)
            
%             time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i)))

%                 lookfor_end
%                 length(time_ramp_up)

            %if time_ramp_up(i,lookfor_end) > time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i)))
            if time_ramp_up(i,lookfor_end) > time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i)))
      
%                 fprintf('time_ramp_up(lookfor_end) = %f\n',time_ramp_up(lookfor_end));
                
%                 fprintf('time_ramp_down(i,2,num_down_chunks) = %f\n',time_ramp_down(down_chunk_array(i,2,num_down_chunks(i))));
                
%                 up_chunk_end = j;
                
                break;
                
            end
            
        end
        
        % If this does not happen, e.g. if the simulation isn't allowed to 
        % finish its last cycle, then reduce the number of down chunks by
        % the appropriate number of cycles.
        if time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i))) > time_ramp_up(i,lookfor_end)
            
            fprintf('Last down chunk (+V) = %.2f min occurs after last up chunk (-V) = %.2f min \n',...
                time_ramp_down(i,down_chunk_array(i,2,end)),time_ramp_up(i,lookfor_end));
            
            disp('Reducing number of down chunks to correct.');
            
            remove_down_chunks = 0;
            
%             up_chunk_end = num_up_chunks(i)+1;
                
            for j = 1:num_down_chunks(i)-1
                
                remove_down_chunks = remove_down_chunks + 1;
                    
                num_down_chunks(i) = num_down_chunks(i) - j;
                
                if time_ramp_down(i,down_chunk_array(i,2,num_down_chunks(i))) < time_ramp_up(i,lookfor_end)
    
                    break;
                    
                end
                
            end            
            
            down_chunk_array = down_chunk_array(:,:,1:end-remove_down_chunks);
            
            fprintf('Removed %d down chunks. Now, Last down chunk (+V) = %.2f min and last up chunk (-V) = %.2f min\n',...
                remove_down_chunks,time_ramp_down(down_chunk_array(i,2,end)),time_ramp_up(lookfor_end));
            
        end
        
        num_up_chunks(i) = up_chunk_end -1 - up_chunk_begin + 1;
        
%         for j = up_chunk_begin:num_up_chunks(i)+1

%         disp(up_chunk_array_pass(1,2,end-3:end));
            
        for j = up_chunk_begin:up_chunk_end - 1

            up_chunk_array_pass_pass(i,1,end+1) = up_chunk_array_pass(i,2,j) + extra;
            
%             fprintf('up_chunk_array_pass_pass(i,1,end+1) = %d\n',up_chunk_array_pass(i,2,j) + extra);
            
            up_chunk_array_pass_pass(i,2,end) = up_chunk_array_pass(i,2,j+1) - extra;
            
%             fprintf('up_chunk_array_pass_pass(i,2,end) = %d\n',up_chunk_array_pass(i,2,j+1) + extra);

            in_this_chunk = up_chunk_array_pass_pass(i,2,j) - up_chunk_array_pass_pass(i,1,j) -1;

        end
        

        % first trash chunk starts after first down chunk. Last trash chunk
        % ends after last up chunk.
        %for j =1:num_trash_chunks(i) + 1
        for j =1:num_trash_chunks(i) + 1
            
            %lookfor_start = trash_chunk_array_pass(i,2,j+1) + extra;
            lookfor_start = trash_chunk_array_pass(i,2,j) + 1;
            
            %if time_trash(i,lookfor_start+1) > time_ramp_down(i,down_chunk_array(i,1,1)+1) && j > 1
            if time_trash(i,lookfor_start) > time_ramp_down(i,down_chunk_array(i,1,1)+1)
                
                %trash_chunk_begin = j-1;
                trash_chunk_begin = j;
                
%                 fprintf('found first trash chunk! = %d \n',j);
                
                break;
                
            end
            
        end
        
        % now look for last trash chunk, which should occur after the last
        % up chunk. 
%         for j = trash_chunk_begin:num_trash_chunks(i) + 1
        for j = trash_chunk_begin:num_trash_chunks(i)            
%             fprintf('trash chunk cycle = %d\n',j);
            
            %lookfor_end = trash_chunk_array_pass(i,2,j+1) - extra;
            lookfor_end = trash_chunk_array_pass(i,2,j) +1;
            
%             if time_trash(i,lookfor_end) > time_ramp_up(i,up_chunk_array_pass_pass(i,2,num_up_chunks(i)+1))
                
            if time_trash(i,lookfor_end) > time_ramp_up(i,up_chunk_array_pass_pass(i,2,end))
                
%                 trash_chunk_end = j;
                trash_chunk_end = j + 1;
                
%                 fprintf('found last trash chunk! = %d \n',j);
                
                break;                                
                
            end
%             
%             fprintf('looking for last trash chunk, currently on %d, did not reach it yet! \n',j);
            
        end      
        
        % If this does not happen, e.g. if the simulation isn't allowed to 
        % finish its last cycle, then reduce the number of up and down 
        % chunks by the appropriate number of cycles.
        
%         if time_ramp_up(i,up_chunk_array_pass_pass(i,2,num_up_chunks(i)+1)) > ...
%                 time_trash(i,lookfor_end)
            
        if time_ramp_up(i,up_chunk_array_pass_pass(i,2,end)) > ...
                time_trash(i,lookfor_end)
            
            fprintf('Last up chunk (-V) = %d occurs after last trash chunk (0V) = %d. Reducing number of up and down chunks to correct. \n',...
                num_up_chunks,num_trash_chunks);
            
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
        
%         num_trash_chunks(i) = trash_chunk_end - trash_chunk_begin + 1; 
        num_trash_chunks(i) = trash_chunk_end -1 - trash_chunk_begin + 1;
        
%         fprintf('num_trash_chunks = %d - %d + 1 = %d\n',trash_chunk_end,...
%             trash_chunk_begin,num_trash_chunks);

%         for j =trash_chunk_begin:num_trash_chunks(i)+1
%         for j =trash_chunk_begin:trash_chunk_end
        
        for j =trash_chunk_begin:trash_chunk_end-1
            
            
%             trash_chunk_array_pass_pass(i,1,end+1) = trash_chunk_array_pass(i,2,j+1) + extra;
%             
%             trash_chunk_array_pass_pass(i,2,end) = trash_chunk_array_pass(i,2,j+2) + extra;
            
            trash_chunk_array_pass_pass(i,1,end+1) = trash_chunk_array_pass(i,2,j) + extra;
            
            trash_chunk_array_pass_pass(i,2,end) = trash_chunk_array_pass(i,2,j+1) + extra;
            
            in_this_chunk = trash_chunk_array_pass_pass(i,2,end) - trash_chunk_array_pass_pass(i,1,end) - 1;

        end

    end
    
    up_chunk_array = up_chunk_array_pass_pass(:,:,2:end);
    trash_chunk_array = trash_chunk_array_pass_pass(:,:,2:end);
    
%     fprintf('file number: %d \n',i);
    fprintf('number of up (-V) : down (+V) : trash (0V) chunks = %d : %d : %d \n',...
        num_up_chunks(i),num_down_chunks(i),num_trash_chunks(i));
%     fprintf('number of  down chunks (+V) = %d \n', num_down_chunks(i));
%     fprintf('number of trash chunks (0V) = %d \n', num_trash_chunks(i));
%     disp('the number of trash chunks should be the sum of up and down chunks')
    
end