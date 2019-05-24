function [dph,dph_std,median_discharge_size,...
    overall_dph,overall_dph_std,overall_median_discharge_size] = ...
    calc_discharge_rate(discharge_times, discharge_vals, ...
    discharge_stdevs,discharges_per_chunk,time_length_of_chunk)

start_time = discharge_times(1,1);

first_discharge_time = discharge_times(1,1);

leftover_time = 0.0;

hour_number = 1;

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
        
        median_discharge_size = [discharge_stdevs(1,1)];
        overall_median_discharge_size = discharge_stdevs(1,1);
        
    elseif length(discharge_times(1,:)) > 1
        
        overall_dph = length(discharge_times(1,start_time_index:end))*60.0/(discharge_times(1,end) - ...
            discharge_times(1,start_time_index));
                
        overall_dph_std = sqrt(overall_dph);

        overall_median_discharge_size = median(discharge_stdevs(1,start_time_index:end));
    
        for i =1:length(discharge_times(1,:))

            % time between this discharge and discharge from last dph value
            
            time_delta = discharge_times(1,i) - start_time;
            
            time_since_last_hour = discharge_times(1,i) - first_discharge_time - 60.0*hour_number;
            
            time_since_last_hour = discharge_times(1,i) - first_discharge_time + ...
                leftover_time - 60.0*hour_number;

%            if time_delta > 60.0*hour_number
                
            if time_since_last_hour > 0.0       
                
                % if more than an hour passes between discharges, we need
                % to record that zero discharges occurred
                if time_since_last_hour > 60.0
                    
                    while time_since_last_hour > 60.0
                        
                        time_since_last_hour = time_since_last_hour - 60.0;
                   
                        avg_discharges_this_hour = 0.0;

                        std_this_hour = sqrt(avg_discharges_this_hour);

                        median_discharge_size_this_hour = 0.0;

                        dph = [dph avg_discharges_this_hour];

                        dph_std = [dph_std std_this_hour];

                        median_discharge_size = [median_discharge_size median_discharge_size_this_hour];
                        
                    end
                    
                end

                %discharges per hour
                avg_discharges_this_hour = length(discharge_times(1,start_time_index:i))*60.0/time_delta;
                
                std_this_hour = sqrt(avg_discharges_this_hour);
                
                median_discharge_size_this_hour = median(discharge_stdevs(1,start_time_index:i));

                dph = [dph avg_discharges_this_hour];
                
                dph_std = [dph_std std_this_hour];
                
                median_discharge_size = [median_discharge_size median_discharge_size_this_hour];

                % now redefine start_time and start_time_index so that we can
                % look for the avg. # of discharges for the next hour.
                
                % Break loop if we're at the last point                                
                
                if i == length(discharge_times(1,:))
                    
                    start_time = discharge_times(1,i);

                    start_time_index = i;
                    
                    break;
                    
                else
                    
%                     start_time = discharge_times(i+1);
% 
%                     start_time_index = i + 1;
                    
                    start_time = discharge_times(1,i);

                    start_time_index = i;
                    
                    hour_number = hour_number+1;
                    
                    leftover_time = time_since_last_hour;
                                        
                    
                end                                

            end

        end
        
        % if we don't get a full last hour, we can extrapolate the last
        % dph
        %
        time_delta = discharge_times(1,end) - start_time;
        
%         if time_delta < 60.0 && start_time_index ~= length(discharge_times) && ...
%                 time_delta > 8.0
            
        if time_since_last_hour < 0.0 && start_time_index ~= length(discharge_times(1,:)) && ...
                time_delta > 8.0
            
            %avg_discharges_this_hour = length(discharge_times(start_time_index:end))*time_delta/60.0;
            avg_discharges_this_hour = length(discharge_times(1,start_time_index:end))*60/time_delta;
            
            std_this_hour = sqrt(avg_discharges_this_hour);
            
            median_discharge_size_this_hour = median(discharge_stdevs(1,start_time_index:end));
            
            median_discharge_size = [median_discharge_size median_discharge_size_this_hour];
            
            dph = [dph avg_discharges_this_hour];
            
            dph_std = [dph_std std_this_hour];
            
%         elseif time_delta < 60.0 && start_time_index ~= length(discharge_times) && ...
%                 time_delta < 8.0
            
        elseif time_since_last_hour < 0.0 && start_time_index ~= length(discharge_times(1,:)) && ...
                time_delta < 8.0
            
            avg_discharges_this_hour = length(discharge_times(1,start_time_index:end))*60.0/time_delta;
         
            std_this_hour = sqrt(avg_discharges_this_hour);

%             fprintf('start_time_index = %d \n',start_time_index);
%             fprintf('length of discharge_stdevs = %d \n',length(discharge_stdevs));
            
            median_discharge_size_this_hour = ...
                median([median_discharge_size(end) discharge_stdevs(1,start_time_index:end)]);
            
            dph(end) = (dph(end) + avg_discharges_this_hour)/2.0;
            
            dph_std(end) = sqrt((dph_std(end))^2+ (std_this_hour)^2);
            
            median_discharge_size(end) = median_discharge_size_this_hour;
    
        end
    
    end