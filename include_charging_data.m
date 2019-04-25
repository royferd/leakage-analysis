%{

include_charging_data takes the charging segments and prepends/appends them
to the existing segments that are already defined for +v and -V and defines
them as new arrays.

%}

function [num_ramp_up_inc_points,num_ramp_down_inc_points,...
    up_chunk_inc_array,down_chunk_inc_array,lcm1_avg_ramp_up_inc_raw,...
    lcm1_weight_ramp_down_inc_raw,lcm1_weight_ramp_up_inc_raw,...
    lcm1_avg_ramp_down_inc_raw,time_ramp_up_inc,time_ramp_down_inc,...
    lcm1_avg_charge_neg_sum_current_raw,lcm1_avg_charge_pos_sum_current_raw,...
    lcm1_avg_discharge_neg_sum_current_raw,lcm1_avg_discharge_pos_sum_current_raw,...
    lcm1_avg_discharge_neg_sum_weight_raw,lcm1_avg_discharge_pos_sum_weight_raw] = ...
    include_charging_data(num_rows,num_up_chunks,num_down_chunks,up_chunk_array,...
    down_chunk_array,lcm1_discharge_pos_time,lcm1_discharge_neg_time,...
    lcm1_weight_discharge_pos_raw,lcm1_weight_discharge_neg_raw,...
    lcm1_avg_discharge_pos_raw,lcm1_avg_discharge_neg_raw,charge_pos_index,...
    discharge_pos_index,charge_neg_index,...
    discharge_neg_index,time_ramp_up,time_ramp_down,...
    lcm1_weight_avg_ramp_up_raw,lcm1_weight_avg_ramp_down_raw,lcm1_avg_ramp_down_raw,...
    lcm1_avg_ramp_up_raw,lcm1_charge_pos_time,...
    lcm1_charge_neg_time,lcm1_weight_charge_pos_raw,lcm1_weight_charge_neg_raw,...
    lcm1_avg_charge_pos_raw,lcm1_avg_charge_neg_raw,...
    max_length_charge,max_length_discharge)

    num_files = length(lcm1_weight_discharge_pos_raw(:,1));

    lcm1_avg_ramp_down_inc_raw_pass = zeros(num_files,1);
    
    lcm1_weight_ramp_down_inc_raw_pass = zeros(num_files,1);
    
    time_ramp_down_inc_pass = zeros(num_files,1);
    
    down_chunk_inc_array_pass = zeros(num_files,2,num_rows);

    lcm1_avg_ramp_up_inc_raw_pass = zeros(num_files,1);
    
    lcm1_weight_ramp_up_inc_raw_pass = zeros(num_files,1);
    
    time_ramp_up_inc_pass = zeros(num_files,1);
    
    up_chunk_inc_array_pass = zeros(num_files,2,num_rows);
    
    
    lcm1_avg_discharge_pos_sum_current_raw = zeros(num_files,max_length_discharge);
    
    lcm1_avg_discharge_pos_sum_weight_raw =  zeros(num_files,max_length_discharge);
    
    lcm1_avg_charge_pos_sum_current_raw = zeros(num_files,max_length_charge);


    
    lcm1_avg_discharge_neg_sum_current_raw = zeros(num_files,max_length_discharge);
    
    lcm1_avg_discharge_neg_sum_weight_raw = zeros(num_files,max_length_discharge);
    
    lcm1_avg_charge_neg_sum_current_raw = zeros(num_files,max_length_charge);

    
    up_chunk_inc_array = zeros(num_files,2,max(num_up_chunks(:)));
    
    down_chunk_inc_array = zeros(num_files,2,max(num_down_chunks(:)));
    
    
    %This builds the inclusive data sets. charging and discharging segments are
    %attached to either end of the appropriate ramp segments.
    
    for i = 1:num_files

%        for j = 1:num_down_chunks(i)
        for j = 1:num_down_chunks(i) - 1 

            down_chunk_inc_array_pass(i,1,j) = length(lcm1_avg_ramp_down_inc_raw_pass(i,:));

            lcm1_avg_ramp_down_inc_raw_pass(i,end+1:end+1+charge_pos_index(i,j+1) - (charge_pos_index(i,j)+1) ) = ...
                lcm1_avg_charge_pos_raw(i,charge_pos_index(i,j)+1:charge_pos_index(i,j+1));

            lcm1_weight_ramp_down_inc_raw_pass(i,end+1:end+1+charge_pos_index(i,j+1) - (charge_pos_index(i,j)+1) ) = ...
                lcm1_weight_charge_pos_raw(i,charge_pos_index(i,j)+1:charge_pos_index(i,j+1));

            time_ramp_down_inc_pass(i,end+1:end+1+charge_pos_index(i,j+1) - (charge_pos_index(i,j)+1) ) = ...
                lcm1_charge_pos_time(i,charge_pos_index(i,j)+1:charge_pos_index(i,j+1));

            lcm1_avg_ramp_down_inc_raw_pass(i,end+1:end+1+down_chunk_array(i,2,j+1) - (down_chunk_array(i,2,j)+1) ) = ...
                lcm1_avg_ramp_down_raw(i,down_chunk_array(i,2,j)+1:down_chunk_array(i,2,j+1));

            lcm1_weight_ramp_down_inc_raw_pass(i,end+1:end+1+down_chunk_array(i,2,j+1) - (down_chunk_array(i,2,j)+1) ) = ...
                lcm1_weight_avg_ramp_down_raw(i,down_chunk_array(i,2,j)+1:down_chunk_array(i,2,j+1));           

            time_ramp_down_inc_pass(i,end+1:end+1+down_chunk_array(i,2,j+1) - (down_chunk_array(i,2,j)+1) ) = ...
                time_ramp_down(i,down_chunk_array(i,2,j)+1:down_chunk_array(i,2,j+1));

            lcm1_avg_ramp_down_inc_raw_pass(i,end+1:end+1+discharge_pos_index(i,j+1) - (discharge_pos_index(i,j)+1) ) = ...
                lcm1_avg_discharge_pos_raw(i,discharge_pos_index(i,j)+1:discharge_pos_index(i,j+1));

            lcm1_weight_ramp_down_inc_raw_pass(i,end+1:end+1+discharge_pos_index(i,j+1) - (discharge_pos_index(i,j)+1) ) = ...
                lcm1_weight_discharge_pos_raw(i,discharge_pos_index(i,j)+1:discharge_pos_index(i,j+1));

            time_ramp_down_inc_pass(i,end+1:end+1+discharge_pos_index(i,j+1) - (discharge_pos_index(i,j)+1) ) = ...
                lcm1_discharge_pos_time(i,discharge_pos_index(i,j)+1:discharge_pos_index(i,j+1));

            down_chunk_inc_array_pass(i,2,j) = length(lcm1_avg_ramp_down_inc_raw_pass(i,:));

%             fprintf('down chunk # = %d \n',num_down_chunks);
%             fprintf('length of lcm1_avg_discharge_pos_sum_current_raw = %d \n',...
%                 length(lcm1_avg_discharge_pos_sum_current_raw));
%             fprintf('length of lcm1_avg_discharge_pos_raw = %d \n',...
%                 length(lcm1_avg_discharge_pos_raw));            
            
            for k = 1:max_length_discharge

                fprintf('down chunk k = %d \n',k);
                fprintf('discharge_pos_index(i,j)+1+k = %d \n',discharge_pos_index(i,j)+1+k);
                
                lcm1_avg_discharge_pos_sum_current_raw(i,k) = ...
                    lcm1_avg_discharge_pos_sum_current_raw(i,k) + ...
                    lcm1_avg_discharge_pos_raw(i,discharge_pos_index(i,j)+1+k);

                lcm1_avg_discharge_pos_sum_weight_raw(i,k) = ...
                    lcm1_avg_discharge_pos_sum_weight_raw(i,k)+ ...
                    lcm1_weight_discharge_pos_raw(i,discharge_pos_index(i,j)+1+k);

            end

            for k = 1:max_length_charge

                lcm1_avg_charge_pos_sum_current_raw(i,k) = ...
                    lcm1_avg_charge_pos_sum_current_raw(i,k)+ lcm1_avg_charge_pos_raw(i,charge_pos_index(i,j)+1+k);               

            end

        end

%        for j = 1:num_up_chunks(i)
        for j = 1:num_up_chunks(i) - 1

            up_chunk_inc_array_pass(i,1,j) = length(lcm1_avg_ramp_up_inc_raw_pass(i,:));          

            lcm1_avg_ramp_up_inc_raw_pass(i,end+1:end+1+charge_neg_index(i,j+1) - (charge_neg_index(i,j)+1) ) = ...
                lcm1_avg_charge_neg_raw(i,charge_neg_index(i,j)+1:charge_neg_index(i,j+1));

            lcm1_weight_ramp_up_inc_raw_pass(i,end+1:end+1+charge_neg_index(i,j+1) - (charge_neg_index(i,j)+1) ) = ...
                lcm1_weight_charge_neg_raw(i,charge_neg_index(i,j)+1:charge_neg_index(i,j+1));

            time_ramp_up_inc_pass(i,end+1:end+1+charge_neg_index(i,j+1) - (charge_neg_index(i,j)+1) ) = ...
                lcm1_charge_neg_time(i,charge_neg_index(i,j)+1:charge_neg_index(i,j+1));

            lcm1_avg_ramp_up_inc_raw_pass(i,end+1:end+1+up_chunk_array(i,2,j+1) - (up_chunk_array(i,2,j)+1) ) = ...
                lcm1_avg_ramp_up_raw(i,up_chunk_array(i,2,j)+1:up_chunk_array(i,2,j+1));

            lcm1_weight_ramp_up_inc_raw_pass(i,end+1:end+1+up_chunk_array(i,2,j+1) - (up_chunk_array(i,2,j)+1) ) = ...
                lcm1_weight_avg_ramp_up_raw(i,up_chunk_array(i,2,j)+1:up_chunk_array(i,2,j+1));           

            time_ramp_up_inc_pass(i,end+1:end+1+up_chunk_array(i,2,j+1) - (up_chunk_array(i,2,j)+1) ) = ...
                time_ramp_up(i,up_chunk_array(i,2,j)+1:up_chunk_array(i,2,j+1));

            lcm1_avg_ramp_up_inc_raw_pass(i,end+1:end+1+discharge_neg_index(i,j+1) - (discharge_neg_index(i,j)+1) ) = ...
                lcm1_avg_discharge_neg_raw(i,discharge_neg_index(i,j)+1:discharge_neg_index(i,j+1));

            lcm1_weight_ramp_up_inc_raw_pass(i,end+1:end+1+discharge_neg_index(i,j+1) - (discharge_neg_index(i,j)+1) ) = ...
                lcm1_weight_discharge_neg_raw(i,discharge_neg_index(i,j)+1:discharge_neg_index(i,j+1));

            time_ramp_up_inc_pass(i,end+1:end+1+discharge_neg_index(i,j+1) - (discharge_neg_index(i,j)+1) ) = ...
                lcm1_discharge_neg_time(i,discharge_neg_index(i,j)+1:discharge_neg_index(i,j+1));

            up_chunk_inc_array_pass(i,2,j) = length(lcm1_avg_ramp_up_inc_raw_pass(i,:));

            for k = 1:max_length_discharge

                lcm1_avg_discharge_neg_sum_current_raw(i,k) = ...
                    (lcm1_avg_discharge_neg_sum_current_raw(i,k) + ...
                    lcm1_avg_discharge_neg_raw(i,discharge_neg_index(i,j)+1+k));

                lcm1_avg_discharge_neg_sum_weight_raw(i,k) = ...
                    (lcm1_avg_discharge_neg_sum_weight_raw(i,k) + ...
                    lcm1_weight_discharge_neg_raw(i,discharge_neg_index(i,j)+1+k));

            end

            for k = 1:max_length_charge
           %     
           %     charge_neg_index(i,j)+1+k;
%            fprintf('k: %d \n',k);
%            fprintf('length of max charge: %d \n',max_length_charge);
%            fprintf('length of lcm1_avg_charge_neg_sum_current_raw: %d \n',...
%                length(lcm1_avg_charge_neg_sum_current_raw));
%            fprintf('length of lcm1_avg_charge_neg_raw: %d \n',...
%                length(lcm1_avg_charge_neg_raw));
%            fprintf('charge_neg_index(i,j)+1+k: %d \n',charge_neg_index(i,j)+1+k);

                lcm1_avg_charge_neg_sum_current_raw(i,k) = ...
                    (lcm1_avg_charge_neg_sum_current_raw(i,k) + ...
                    lcm1_avg_charge_neg_raw(i,charge_neg_index(i,j)+1+k)); 

            end


        end
         
        num_ramp_down_inc_points = length(lcm1_avg_ramp_down_inc_raw_pass(i,:))-1;
      
        for j = 1:num_down_chunks(i)

            down_chunk_inc_array(i,1,j) = down_chunk_inc_array_pass(i,1,j) - 1;
       
            down_chunk_inc_array(i,2,j) = down_chunk_inc_array_pass(i,2,j) - 1;

        end

        num_ramp_up_inc_points = length(lcm1_avg_ramp_up_inc_raw_pass(i,:))-1;
      
        for j = 1:num_up_chunks(i)

            up_chunk_inc_array(i,1,j) = up_chunk_inc_array_pass(i,1,j) - 1;
         
            up_chunk_inc_array(i,2,j) = up_chunk_inc_array_pass(i,2,j) - 1;

        end

    end
    
    lcm1_avg_ramp_down_inc_raw = lcm1_avg_ramp_down_inc_raw_pass(:,2:end);
  
    lcm1_weight_ramp_down_inc_raw = lcm1_weight_ramp_down_inc_raw_pass(:,2:end);

    time_ramp_down_inc = time_ramp_down_inc_pass(:,2:end);
  
    lcm1_avg_ramp_up_inc_raw = lcm1_avg_ramp_up_inc_raw_pass(:,2:end);
 
    lcm1_weight_ramp_up_inc_raw = lcm1_weight_ramp_up_inc_raw_pass(:,2:end);

    time_ramp_up_inc = time_ramp_up_inc_pass(:,2:end);

end