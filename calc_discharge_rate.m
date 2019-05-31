function [dph,dph_std,median_discharge_size,...
    overall_dph,overall_dph_std,overall_median_discharge_size] = ...
    calc_discharge_rate(discharge_times, discharge_vals, ...
    discharge_stdevs,discharges_per_chunk,time_length_of_chunk)

start_time = discharge_times(1,1);

first_discharge_time = discharge_times(1,1);

leftover_time = 0.0;

hour_number = 1;

live_time = 0.0;

    start_time_index = 1;
    
    dph = [];
    overall_dph = 0;
    
    dph_std = [];
    overall_dph_std = 0;
    
    median_discharge_size = [];
    overall_median_discharge_size = 0;
    
    if length(discharge_times(1,:)) < 1
        
        dph = [0];
        overall_dph = 0;
        
        dph_std = [0];
        overall_dph_std = 0;
        
        median_discharge_size = [0];
        overall_median_discharge_size = 0;
        
    elseif length(discharge_times(1,:)) == 1
        
        dph = [1];
        overall_dph = 1;
        
        dph_std = [1];
        overall_dph_std = 1;
        
        median_discharge_size = [discharge_vals(1,1)];
        
        overall_median_discharge_size = discharge_vals(1,1);
        
    elseif length(discharge_times(1,:)) > 1
        
        overall_dph = length(discharge_times(1,start_time_index:end))*60.0/(discharge_times(1,end) - ...
            discharge_times(1,start_time_index));
                
        overall_dph_std = sqrt(overall_dph);
        
        if isempty(find(discharge_vals(1,start_time_index:end),1))
            
            overall_median_discharge_size = 0.0;
            
        else
            
            [row, col, val] = find(discharge_vals(1,start_time_index:end));

            overall_median_discharge_size = median(val);
            
        end
    
        for i =1:length(discharge_times(1,:))

            % time between this discharge and discharge from last dph value
            
            time_delta = discharge_times(1,i) - start_time;
            
%            live_time = live_time + time_length_of_chunk(discharge_times(2,i)); 

            if i==start_time_index
                
                live_time = live_time + ...
                    sum(time_length_of_chunk(start_time_index:discharge_times(2,i)));
            
            elseif discharge_times(2,i) ~= discharge_times(2,i-1)
                    
                    live_time = live_time + ...
                        sum(time_length_of_chunk(discharge_times(2,i-1)+1:discharge_times(2,i)));
                                    
            end
            
%             time_since_last_hour = discharge_times(1,i) - first_discharge_time - 60.0*hour_number;
            
            time_since_last_hour = discharge_times(1,i) - first_discharge_time + ...
                leftover_time - 60.0*hour_number;
            
            % if time_since_last_hour > 0, than over an hour has elapsed
            % since the first discharge in this sequence of discharges, and
            % we compute a new discharge rate
            
            if time_since_last_hour > 0.0 
                
                % if more than an hour passes between discharges, we need
                % to record that zero discharges occurred
                
%                 if time_since_last_hour > 60.0 
                    
                if time_since_last_hour > 60.0 && isempty(find(discharge_vals(1,start_index:i),1))
                    
                    while time_since_last_hour > 60.0
                        
                        time_since_last_hour = time_since_last_hour - 60.0;
                   
                        avg_discharges_this_hour = 0.0;

                        std_this_hour = sqrt(avg_discharges_this_hour);

                        median_discharge_size_this_hour = 0.0;

%                         dph = [dph avg_discharges_this_hour];
% 
%                         dph_std = [dph_std std_this_hour];
% 
%                         median_discharge_size = ...
%                             [median_discharge_size median_discharge_size_this_hour];
                        
                    end
                    
%                 end
                else

                    %discharges per hour

%                     avg_discharges_this_hour = ...
%                         length(discharge_times(1,start_time_index:i))*60.0/live_time;
                    
                    avg_discharges_this_hour = ...
                        length(find(discharge_vals(1,start_time_index:i)))*60.0/live_time;

                    std_this_hour = sqrt(avg_discharges_this_hour);

%                     median_discharge_size_this_hour = ...
%                         median(discharge_vals(1,start_time_index:i));
                    

                    if isempty(find(discharge_vals(1,start_time_index:i),1))
                        
                        median_discharge_size_this_hour = 0.0;
                        
                    else
                        
                        [row, col, val] = find(discharge_vals(1,start_time_index:i));
                       
                        median_discharge_size_this_hour = median(val);
                        
                    end

%                     dph = [dph avg_discharges_this_hour];
% 
%                     dph_std = [dph_std std_this_hour];
% 
%                     median_discharge_size = ...
%                         [median_discharge_size median_discharge_size_this_hour];
                
                end

                % now redefine start_time and start_time_index so that we can
                % look for the avg. # of discharges for the next hour.
                
                % Break loop if we're at the last point                                
                
%                 if i == length(discharge_times(1,:))
%                     
% %                     start_time = discharge_times(1,i);
% 
% %                     start_time_index = i;
%                     
%                     break;
                    
                if i ~= length(discharge_times(1,:))
                    
                    start_time = discharge_times(1,i);

                    start_time_index = i;
                    
                    hour_number = hour_number+1;
                    
                    leftover_time = time_since_last_hour;
                    
                    live_time = 0.0;
                    
                    dph = [dph avg_discharges_this_hour];

                    dph_std = [dph_std std_this_hour];

                    median_discharge_size = ...
                        [median_discharge_size median_discharge_size_this_hour];
                    
                end                                

            end

        end
        
        avg_discharges_this_hour = ...
                length(find(discharge_vals(1,start_time_index:end)))*60/live_time;
            
        std_this_hour = sqrt(avg_discharges_this_hour);

        if isempty(find(discharge_vals(1,start_time_index:end),1))

            median_discharge_size_this_hour = 0.0;

        else
            
            [row, col, val] = find(discharge_vals(1,start_time_index:end));
                       
            median_discharge_size_this_hour = median(val);

        end           
        
        % if we don't get a full last hour, we can extrapolate the last
        % dph
        
        time_delta = discharge_times(1,end) - start_time;
            
%         if time_since_last_hour < 0.0 && start_time_index ~= length(discharge_times(1,:)) && ...
%                 time_delta > 8.0
            
        if time_since_last_hour < 0.0  && time_delta > 15.0
                       
%             avg_discharges_this_hour = ...
%                 length(discharge_times(1,start_time_index:end))*60/live_time;
            
%             avg_discharges_this_hour = ...
%                 length(find(discharge_vals(1,start_time_index:end)))*60/live_time;
%             
%             std_this_hour = sqrt(avg_discharges_this_hour);
%             
%             if isempty(find(discharge_vals(1,start_time_index:end),1))
%                 
%                 median_discharge_size_this_hour = 0.0;
%                 
%             else
%                                             
%                 median_discharge_size_this_hour = ...
%                     median(find(discharge_vals(1,start_time_index:end)));
%                 
%             end    
                
            median_discharge_size = ...
                [median_discharge_size median_discharge_size_this_hour];
            
            dph = [dph avg_discharges_this_hour];
            
            dph_std = [dph_std std_this_hour];
            
%         elseif time_since_last_hour < 0.0 && start_time_index ~= length(discharge_times(1,:)) && ...
%                 time_delta < 8.0
            
        else
          
%             avg_discharges_this_hour = ...
%                 length(discharge_times(1,start_time_index:end))*60.0/live_time;
            
%             avg_discharges_this_hour = ...
%                 length(find(discharge_vals(1,start_time_index:end)))*60.0/live_time;
%          
%             std_this_hour = sqrt(avg_discharges_this_hour);
% 
% %             fprintf('start_time_index = %d \n',start_time_index);
% %             fprintf('length of discharge_stdevs = %d \n',length(discharge_stdevs));
%             
% %             median_discharge_size_this_hour = ...
% %                 median([median_discharge_size(end) discharge_vals(1,start_time_index:end)]);
%             
% 
            if isempty(find([median_discharge_size(end) discharge_vals(1,start_time_index:end)],1))
                
                median_discharge_size_this_hour = 0.0;
                
            else
                
                [row, col, val] = find([median_discharge_size(end) discharge_vals(1,start_time_index:end)]);
                
                median_discharge_size_this_hour = median(val);
                
            end
            
%             dph(end) = (dph(end) + avg_discharges_this_hour)/2.0;
%             
%             dph_std(end) = sqrt((dph_std(end))^2+ (std_this_hour)^2);
%             
%             median_discharge_size(end) = median_discharge_size_this_hour;
            
            dph(end) = (dph(end) + (live_time/60.0)*avg_discharges_this_hour)/(1.0+live_time/60.0);
            
            dph_std(end) = sqrt(dph_std(end));
            
            median_discharge_size(end) = median_discharge_size_this_hour;
    
        end
    
    end