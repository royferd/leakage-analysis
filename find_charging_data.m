%{

find_charging_data() separates data into 1) charging to -V 2) charging to 
+V 3) discharging to 0 from -V and 4) discharging to 0 from +v. It looks at
the behavior of the leakage current sets that include 0 V data plus the 
charging currents

%}
function [lcm1_avg_charge_neg_raw, lcm1_stdev_charge_neg_raw,...
    lcm1_weight_charge_neg_raw, lcm1_charge_neg_time,lcm1_avg_charge_pos_raw,...
    lcm1_stdev_charge_pos_raw,lcm1_weight_charge_pos_raw,...
    lcm1_charge_pos_time,lcm1_avg_discharge_pos_raw,...
    lcm1_stdev_discharge_pos_raw,lcm1_weight_discharge_pos_raw,...
    lcm1_discharge_pos_time,lcm1_avg_discharge_neg_raw,...
    lcm1_stdev_discharge_neg_raw,lcm1_weight_discharge_neg_raw,...
    lcm1_discharge_neg_time,max_length_discharge,max_length_charge,...
    charge_pos_index,charge_neg_index,discharge_pos_index,...
    discharge_neg_index,discharge_index,charge_index,lcm1_avg_zero_raw,...
    lcm1_stdev_zero_raw,lcm1_weight_avg_zero_raw,...
    zero_time,num_zero_chunks,num_zero_points,zero_chunk_array] = ...
    find_charging_data(num_trash_chunks,trash_chunk_array,lcm1_inv_weight_avg_trash,...
    lcm1_avg_trash,lcm1_avg_trash_raw,...
    lcm1_inv_weight_avg_trash_raw,...
    lcm1_weight_avg_trash_raw,time_trash)

    disp('running find_charging_data.m');
    
    % when identifying a charging or discharging segment, include a few
    % extra points in the direction of increasing stability to make double
    % sure we're collecting all transitional data points.
    buffer = 2;

    num_files = length(num_trash_chunks);
    
    lcm1_avg_charge_neg_pass = zeros(num_files,1);
    
    lcm1_stdev_charge_neg_pass = zeros (num_files,1);
    
    lcm1_weight_charge_neg_raw_pass = zeros (num_files,1);
    
    lcm1_avg_charge_neg_raw_pass = zeros(num_files,1);
    
    lcm1_stdev_charge_neg_raw_pass = zeros (num_files,1);
    
    lcm1_charge_neg_pass_time = zeros(num_files,1);

    
    lcm1_avg_charge_pos_pass = zeros(num_files,1);
    
    lcm1_stdev_charge_pos_pass = zeros (num_files,1);
    
    lcm1_weight_charge_pos_raw_pass = zeros (num_files,1);
    
    lcm1_avg_charge_pos_raw_pass = zeros(num_files,1);
    
    lcm1_stdev_charge_pos_raw_pass = zeros (num_files,1);
    
    lcm1_charge_pos_pass_time = zeros(num_files,1);

    
    lcm1_avg_discharge_pos_pass = zeros(num_files,1);
    
    lcm1_stdev_discharge_pos_pass = zeros (num_files,1);
    
    lcm1_avg_discharge_pos_raw_pass = zeros(num_files,1);
    
    lcm1_stdev_discharge_pos_raw_pass = zeros (num_files,1);
    
    lcm1_weight_discharge_pos_raw_pass = zeros (num_files,1);
    
    lcm1_discharge_pos_pass_time = zeros(num_files,1);

    
    lcm1_avg_discharge_neg_pass = zeros(num_files,1);
    
    lcm1_stdev_discharge_neg_pass = zeros (num_files,1);
    
    lcm1_weight_discharge_neg_raw_pass = zeros (num_files,1);
    
    lcm1_avg_discharge_neg_raw_pass = zeros(num_files,1);
    
    lcm1_stdev_discharge_neg_raw_pass = zeros (num_files,1);
    
    lcm1_discharge_neg_pass_time = zeros(num_files,1);
    
    
    lcm1_avg_zero_pass = [];
    
    lcm1_stdev_zero_pass = [];
    
    lcm1_weight_avg_zero_raw_pass = [];
    
    lcm1_avg_zero_raw_pass = [];
    
    lcm1_stdev_zero_raw_pass = [];
    
    lcm1_zero_pass_time = [];
    

% sample +/- steady_state_range to estimate the steady-state leakage
% current
%    steady_state_range = 10;
    steady_state_range = 20;     
    
    for i =1:num_files
        
        for j = 1:num_trash_chunks(i)
            
%            chunk_begin = trash_chunk_array(i,2,j);
            chunk_begin = trash_chunk_array(i,2,j) + 1;
            
            chunk_end = trash_chunk_array(i,2,j+1);
            
            middle = round((chunk_end - chunk_begin)/2,0);
            
            middle_lo = chunk_begin + middle - steady_state_range;
            
            middle_hi = chunk_begin + middle + steady_state_range;
            
            lcm1_avg_trash_steady_avg(i,j) = mean(lcm1_avg_trash(i,middle_lo:middle_hi));
            
            lcm1_avg_trash_steady_stdev(i,j) = ...
                mean(lcm1_inv_weight_avg_trash(i,middle_lo:middle_hi));
            
%             lcm1_avg_trash_steady_stdev(i,j) = ...
%                 std(lcm1_inv_weight_avg_trash(i,middle_lo:middle_hi),...
%                 lcm1_inv_weight_avg_trash(i,middle_lo:middle_hi));
        
        end


        % %track index of discharge currents and charging currents in trash data.
        % %column 1: positive/negative (dis)charge
        % %column 2: end point of previous chunk
        % %column 3: end point of this chunk
        discharge_index_pass = zeros(num_files,num_trash_chunks,3);
        charge_index_pass = zeros(num_files,num_trash_chunks,3);

        %pick out the ramping currents. +V -> 0, 0 -> -V, -V -> 0, 0 -> +V.
        max_length_discharge = 0;

        for j = 1:num_trash_chunks(i)

            %first point might be the start of the ramp, which could be a small
            %stdev. Increment by 1 so we know we're starting on a ramp point
            chunk_begin = trash_chunk_array(i,2,j)+1;

            chunk_end = trash_chunk_array(i,2,j+1);
            middle = round((chunk_end - chunk_begin)/2,0);

            for k = chunk_begin:chunk_begin+ middle

%                 if (lcm1_inv_weight_avg_trash(i,k+1) < ...
%                         1.5 * lcm1_avg_trash_steady_stdev(i,j) && ...
%                         lcm1_inv_weight_avg_trash(i,k+1) > ...
%                         lcm1_inv_weight_avg_trash(i,k+2))
                    
                if (lcm1_inv_weight_avg_trash(i,k+1) < ...
                        1.5 * lcm1_avg_trash_steady_stdev(i,j) && ...
                        abs(lcm1_avg_trash(i,k+1)) > ...
                        abs(lcm1_avg_trash(i,k+2)) && k > chunk_begin)

%                     if (max(lcm1_avg_trash(i,chunk_begin:k)) < ...
%                             lcm1_avg_trash_steady_avg(i,j) && ...
%                             abs(max(lcm1_avg_trash(i,chunk_begin:k))) > ...
%                             3 * abs(lcm1_avg_trash_steady_avg(i,j)))

                    [max_leak_mag,max_leak_index] = max(abs(lcm1_avg_trash(i,chunk_begin:k)));
                    max_leak_value = max_leak_mag(1)*sign(lcm1_avg_trash(i,chunk_begin+max_leak_index(1)-1));
%                     fprintf('max leakage value for chunk %d: %f \n',j,max_leak_value);
                    
                    if max_leak_value > lcm1_avg_trash_steady_avg(i,j)

                        fprintf('%f > %f, this is a discharge ramp from +V \n',...
                            max_leak_value,lcm1_avg_trash_steady_avg(i,j));
                        %discharge from +V
                        discharge_index_pass(i,j,1) = 1;

%                     elseif (max(lcm1_avg_trash(i,chunk_begin:k)) > ...
%                             lcm1_avg_trash_steady_avg(i,j) && ...
%                             abs(max(lcm1_avg_trash(i,chunk_begin:k))) > ...
%                             3 * abs(lcm1_avg_trash_steady_avg(i,j)))
                    
                    elseif max_leak_value < lcm1_avg_trash_steady_avg(i,j)
                        
                        fprintf('%f < %f, this is a discharge ramp from -V \n',...
                            max_leak_value,lcm1_avg_trash_steady_avg(i,j));

                        %discharge from -V
                        discharge_index_pass(i,j,1) = -1;

                    end

                    discharge_index_pass(i,j,2) = chunk_begin;
%                     discharge_index_pass(i,j,3) = k+2;
                    discharge_index_pass(i,j,3) = k+buffer;

                    break            

                end

            end

        % if the condition for the discharging current for a segment is not 
        % met, I simply define the segment to start from chunk_begin and have
        % the same # of points as the previous segment
        if discharge_index_pass(i,j,1) == 0
          
            discharge_index_pass(i,j,1) = - discharge_index_pass(i,j-1,1);            

            discharge_index_pass(i,j,2) = chunk_begin;

            discharge_index_pass(i,j,3) = (chunk_begin + ...
                discharge_index_pass(i,j-1,3) - discharge_index_pass(i,j-1,2));
            
            fprintf('trash chunk # %d: unable to find ramp, assuming discharge from %d voltage. \n',...
                j,sign(discharge_index_pass(i,j,1)) );

        end

        if max_length_discharge < discharge_index_pass(i,j,3) - discharge_index_pass(i,j,2)

            max_length_discharge = discharge_index_pass(i,j,3) - discharge_index_pass(i,j,2);

        end

        end

        max_length_charge = 0;

        for j = 1:num_trash_chunks(i)

            chunk_begin = trash_chunk_array(i,2,j)+1;

            chunk_end = trash_chunk_array(i,2,j+1);

            middle = round((chunk_end - chunk_begin)/2,0);

            for k = chunk_begin+ middle:chunk_end 
                

%                 if (lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k) ...
%                         < 1.5 * lcm1_avg_trash_steady_stdev(i,j) && ...
%                         lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k) ...
%                         > lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k - 1))
                if (lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k) ...
                    < 2.0 * lcm1_avg_trash_steady_stdev(i,j) && ...
                    lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k) ...
                    > lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k - 1))
                    
                    [max_leak_mag,max_leak_index] = max(abs(lcm1_avg_trash(i,chunk_end + chunk_begin + middle - k - 1:chunk_end)));
                    max_leak_value = max_leak_mag(1)*sign(lcm1_avg_trash(i,chunk_end + chunk_begin + middle - k - 1+max_leak_index(1)-1));
                    %fprintf('max leakage value for chunk %d: %f \n',j,max_leak_value);

%                     if (max(lcm1_avg_trash(i,k:chunk_end)) >  lcm1_avg_trash_steady_avg(i,j) && ...
%                             abs(max(lcm1_avg_trash(i,k:chunk_end))) ...
%                             > 3*abs(lcm1_avg_trash_steady_avg(i,j)))
                    if (max_leak_value > lcm1_avg_trash_steady_avg(i,j))

                        % charging to +V
                        charge_index_pass(i,j,1) = 1;

%                     elseif max(lcm1_avg_trash(i,k:chunk_end)) <  lcm1_avg_trash_steady_avg(i,j) && ...
%                             abs(max(lcm1_avg_trash(i,k:chunk_end))) > 3*abs(lcm1_avg_trash_steady_avg(i,j))                
                    elseif max_leak_value <  lcm1_avg_trash_steady_avg(i,j)

                        % charging to -V
                        charge_index_pass(i,j,1) = -1;

                    end

%                     charge_index_pass(i,j,2) = chunk_end + chunk_begin + middle - k -1;
                    charge_index_pass(i,j,2) = chunk_end + chunk_begin + middle - (k + buffer);

                    charge_index_pass(i,j,3) = chunk_end;                

                    break       

                end

            end

        % if the condition for detecting the charging current is not satisfied
        % for a chunk, I simply define the charging segment to be the same #
        % of points as the previous segment
            if charge_index_pass(i,j,1) == 0

                charge_index_pass(i,j,1) = -charge_index_pass(i,j-1,1);

                charge_index_pass(i,j,2) = (chunk_end - (charge_index_pass(i,j-1,3) ...
                    - charge_index_pass(i,j-1,2)));

                charge_index_pass(i,j,3) = chunk_end;
                
                fprintf('trash chunk # %d: unable to find ramp, assuming charge to %d voltage. \n',...
                j,sign(charge_index_pass(i,j,1)) );
            end

            if max_length_charge < charge_index_pass(i,j,3) - charge_index_pass(i,j,2)

                max_length_charge = charge_index_pass(i,j,3) - charge_index_pass(i,j,2);

            end

        end

        %want to make the charge/discharge length long enough so that the last
        %points are almost steady-state
        max_length_charge = max_length_charge+round(0.09*max_length_charge+1,0);

        max_length_discharge = max_length_discharge+round(0.09*max_length_discharge+1,0);

        %make all segments the same length which is determined by the maximum
        %segment size
        for j = 1:num_trash_chunks

            charge_index_pass(i,j,2) = charge_index_pass(i,j,3) - max_length_charge;

            discharge_index_pass(i,j,3) = discharge_index_pass(i,j,2) + max_length_discharge;

        end
    
    end

    % We want all but the first discharge_index and all but the last
    % charge_index.
    %
    % <type>_index_count, where <type> can be charge_pos, charge_neg,
    % discharge_pos, discharge_neg, increments when we iterate thru an
    % index of <type>'s.
    % <type>_index_count_increment is a vector that we will use to check
    % whether the count updated or not.
    
    charge_index = charge_index_pass(:,1:end-1,:);
    discharge_index = discharge_index_pass(:,2:end,:);

    charge_neg_index_count = 0;
    charge_neg_index_count_increment = [0];
    charge_neg_index = zeros(num_files,num_trash_chunks);

    discharge_neg_index_count = 0;
    discharge_neg_index_count_increment = [0];
    discharge_neg_index = zeros(num_files,num_trash_chunks);

    charge_pos_index_count = 0;
    charge_pos_index_count_increment = [0];
    charge_pos_index = zeros(num_files,num_trash_chunks);

    discharge_pos_index_count = 0;
    discharge_pos_index_count_increment = [0];
    discharge_pos_index = zeros(num_files,num_trash_chunks);
    
    num_zero_chunks = 0;
    
    zero_chunk_array_pass = ones (1,2,1);


    %create vectors containing time, leakage avg, and leakage stdev data 
    %for + V -> 0, -V -> 0, 0 -> +V, and 0 -> -V
    
%     'discharge_index: '
%     discharge_index
%     
%     'charge_index: '
%     charge_index
    
    for i = 1:num_files

        for j = 1:num_trash_chunks -1

            if discharge_index(i,j,1) == 1

                discharge_pos_index_count = discharge_pos_index_count+1;
                
%                 discharge_pos_index_count_increment(end+1) = discharge_pos_index_count;
%                 
%                 discharge_neg_index_count_increment(end+1) = discharge_neg_index_count_increment(end);

                discharge_pos_index(i,discharge_pos_index_count+1) = ...
                    discharge_pos_index(i,discharge_pos_index_count) + ...
                    discharge_index(i,j,3) - discharge_index(i,j,2)+1;            

                lcm1_avg_discharge_pos_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = ...
                    lcm1_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));

                lcm1_stdev_discharge_pos_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = ...
                    lcm1_inv_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));

                lcm1_weight_discharge_pos_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = ...
                    lcm1_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));

                lcm1_discharge_pos_pass_time(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = ...
                    time_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));

            elseif discharge_index(i,j,1) == -1

                discharge_neg_index_count = discharge_neg_index_count+1;
                
%                 discharge_neg_index_count_increment(end+1) = discharge_neg_index_count;
%                 
%                 discharge_pos_index_count_increment(end+1) = discharge_pos_index_count_increment(end);

                discharge_neg_index(i,discharge_neg_index_count+1) = discharge_neg_index(i,discharge_neg_index_count) + discharge_index(i,j,3) - discharge_index(i,j,2)+1;

                lcm1_avg_discharge_neg_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = lcm1_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));

                lcm1_stdev_discharge_neg_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) =lcm1_inv_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));

                lcm1_weight_discharge_neg_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) =lcm1_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));

                lcm1_discharge_neg_pass_time(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = time_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));

            end
            
            discharge_pos_index_count_increment(end+1) = discharge_pos_index_count;
                
            discharge_neg_index_count_increment(end+1) = discharge_neg_index_count;
            

            if charge_index(i,j,1) == 1

                charge_pos_index_count = charge_pos_index_count+1; 
                
                %wait for ramp to zero to occur before collecting no
                %voltage data
%                 if discharge_neg_index_count_increment(j+1) > discharge_neg_index_count_increment(j)
%                 
%                     charge_pos_index_count_increment(end+1) = charge_pos_index_count;
%                     
%                 else charge_pos_index_count_increment(end+1) = charge_pos_index_count_increment(end);
%                     
%                 end
                
%                charge_neg_index_count_increment(end+1) = charge_neg_index_count_increment(end);

                lcm1_avg_charge_pos_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = lcm1_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));

                lcm1_weight_charge_pos_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) =lcm1_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));

                lcm1_stdev_charge_pos_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) =lcm1_inv_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));

                lcm1_charge_pos_pass_time(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = time_trash(i,charge_index(i,j,2):charge_index(i,j,3));

                charge_pos_index(i,charge_pos_index_count+1) = charge_pos_index(i,charge_pos_index_count) + charge_index(i,j,3) - charge_index(i,j,2)+1;

            elseif charge_index(i,j,1) == -1
                
                charge_neg_index_count = charge_neg_index_count+1;

                
%                 if discharge_pos_index_count_increment(j+1) > discharge_pos_index_count_increment(j)
%                 
%                     charge_neg_index_count_increment(end+1) = charge_neg_index_count;
%                     
%                 else charge_neg_index_count_increment(end+1) = charge_neg_index_count_increment(end);
% 
%                 end
                    
%                charge_pos_index_count_increment(end+1) = charge_pos_index_count_increment(end);

                lcm1_avg_charge_neg_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = ...
                    lcm1_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));

                lcm1_stdev_charge_neg_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = ...
                    lcm1_inv_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));

                lcm1_weight_charge_neg_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = ...
                    lcm1_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));

                lcm1_charge_neg_pass_time(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = ...
                    time_trash(i,charge_index(i,j,2):charge_index(i,j,3));

                charge_neg_index(i,charge_neg_index_count+1) = charge_neg_index(i,charge_neg_index_count) + ...
                    charge_index(i,j,3) - charge_index(i,j,2)+1;

            end
            
            charge_pos_index_count_increment(end+1) = charge_pos_index_count;
                
            charge_neg_index_count_increment(end+1) = charge_neg_index_count;
            
            
%             if discharge_pos_index_count_increment(j+1) > discharge_pos_index_count_increment(j) && ...
%                     charge_neg_index_count_increment(j+1) > charge_neg_index_count_increment(j) && ...
%                     charge_index(i,j,2) > discharge_index(i,j,3) || ...
%                     discharge_neg_index_count_increment(j+1) > discharge_neg_index_count_increment(j) && ...
%                     charge_pos_index_count_increment(j+1) > charge_pos_index_count_increment(j) && ...
%                     charge_index(i,j,2) > discharge_index(i,j,3)

            if j > 1 && charge_index(i,j,2)-1 > discharge_index(i,j-1,3)+1
                
                num_zero_chunks = num_zero_chunks + 1;
                
                length_of_chunk = (charge_index(i,j,2)-1)-(discharge_index(i,j-1,3)+1);
                
                zero_chunk_array_pass(1,1,end+1) = zero_chunk_array_pass(1,2,end);
                zero_chunk_array_pass(1,2,end) = zero_chunk_array_pass(1,1,end) + length_of_chunk;

                lcm1_avg_zero_raw_pass(end+1:end+1+(charge_index(i,j,2)-1)-(discharge_index(i,j-1,3)+1)) = ...
                    lcm1_avg_trash_raw(i,(discharge_index(i,j-1,3)+1):(charge_index(i,j,2)-1));

                lcm1_stdev_zero_raw_pass(end+1:end+1+(charge_index(i,j,2)-1)-(discharge_index(i,j-1,3)+1)) = ...
                    lcm1_inv_weight_avg_trash_raw(i,(discharge_index(i,j-1,3)+1):(charge_index(i,j,2)-1));

                lcm1_weight_avg_zero_raw_pass(end+1:end+1+(charge_index(i,j,2)-1)-(discharge_index(i,j-1,3)+1)) = ...
                    lcm1_weight_avg_trash_raw(i,(discharge_index(i,j-1,3)+1):(charge_index(i,j,2)-1));

                lcm1_zero_pass_time(end+1:end+1+(charge_index(i,j,2)-1)-(discharge_index(i,j-1,3)+1)) = ...
                    time_trash(i,(discharge_index(i,j-1,3)+1):(charge_index(i,j,2)-1));

            elseif j > 1

                msg = sprintf('at %dth trash chunk, we did not trigger a zero chunk',j)

            end




        end

    end
    
    %get rid of those pesky zeros in front
    
    zero_chunk_array = zero_chunk_array_pass(:,:,2:end);
    
    lcm1_avg_zero_raw = lcm1_avg_zero_raw_pass(2:end);
    
    lcm1_stdev_zero_raw = lcm1_stdev_zero_raw_pass(2:end);
    
    lcm1_weight_avg_zero_raw = lcm1_weight_avg_zero_raw_pass(2:end);
    
    zero_time = lcm1_zero_pass_time(2:end);

    lcm1_avg_charge_neg_raw = lcm1_avg_charge_neg_raw_pass(:,2:end);

    lcm1_stdev_charge_neg_raw = lcm1_stdev_charge_neg_raw_pass(:,2:end);

    lcm1_weight_charge_neg_raw = lcm1_weight_charge_neg_raw_pass(:,2:end);

    lcm1_charge_neg_time = lcm1_charge_neg_pass_time(:,2:end);


    lcm1_avg_charge_pos_raw = lcm1_avg_charge_pos_raw_pass(:,2:end);

    lcm1_stdev_charge_pos_raw = lcm1_stdev_charge_pos_raw_pass(:,2:end);

    lcm1_weight_charge_pos_raw = lcm1_weight_charge_pos_raw_pass(:,2:end);

    lcm1_charge_pos_time = lcm1_charge_pos_pass_time(:,2:end);


    lcm1_avg_discharge_pos_raw = lcm1_avg_discharge_pos_raw_pass(:,2:end);

    lcm1_stdev_discharge_pos_raw = lcm1_stdev_discharge_pos_raw_pass(:,2:end);

    lcm1_weight_discharge_pos_raw = lcm1_weight_discharge_pos_raw_pass(:,2:end);

    lcm1_discharge_pos_time = lcm1_discharge_pos_pass_time(:,2:end);


    lcm1_avg_discharge_neg_raw = lcm1_avg_discharge_neg_raw_pass(:,2:end);

    lcm1_stdev_discharge_neg_raw = lcm1_stdev_discharge_neg_raw_pass(:,2:end);

    lcm1_weight_discharge_neg_raw = lcm1_weight_discharge_neg_raw_pass(:,2:end);

    lcm1_discharge_neg_time = lcm1_discharge_neg_pass_time(:,2:end);
    
    num_zero_points = length(lcm1_avg_zero_raw);
    
    

end