function [dph,dph_std,median_discharge_size,...
    overall_dph,overall_dph_std,overall_median_discharge_size] = ...
    calc_discharge_rate(discharge_times, discharge_vals, ...
    discharge_stdevs,discharges_per_chunk,time_length_of_chunk,...
    official_sim_start_time,official_sim_end_time,title_string,plotname,save_figs,savepath)

    

    xlabel_string = 'leakage current (pA)';
    ylabel_string = 'frequency';
    %text box assignments
    inside_plot = [0.1375 0.73 0.5 0.19]; % edges: x y width height


    start_time = official_sim_start_time;

    end_time = official_sim_end_time;

    % fprintf('official EDMRS runtime = %.1f hr\n',(official_sim_end_time - official_sim_start_time)/60.0);

    % add an element that indicates to the algorithm when the simulation has
    % officially ended so that we can account for the time between the last
    % discharge and the actual end of the simulation

    discharge_times = [discharge_times [end_time ; length(time_length_of_chunk)] ];
    
    discharge_time_chunk_check_ahead = [discharge_times [end_time; length(time_length_of_chunk)+1] ];

    discharge_vals = [discharge_vals [0 ; length(time_length_of_chunk)] ];

    discharge_stdevs = [discharge_stdevs [ 0 ; length(time_length_of_chunk)] ];

    first_discharge_time = discharge_times(1,1);

    hour_number = 1;

    live_time = 0.0;

    start_time_index = 1;
    
    first_discharge_index = 1;
    
    dph = [];
    overall_dph = 0;
    
    dph_std = [];
    overall_dph_std = 0;
    
    median_discharge_size = [];
    overall_median_discharge_size = 0;
    
    vals_list = [];
    vals_per_hour = [];
    
    live_time_list = [];
    
    
        
    if length(discharge_times(1,:)) < 2
        
        
        
        num_hours = round((end_time - start_time)/60.0);
    
        dph = zeros(1,num_hours);
        overall_dph = zeros(1,num_hours);
        
        dph_std = zeros(1,num_hours);
        overall_dph_std = zeros(1,num_hours);
        
        median_discharge_size = zeros(1,num_hours);
        overall_median_discharge_size = zeros(1,num_hours);
        
    else
        
        overall_dph = ...
            length(discharge_times(1,start_time_index:end))*60.0/(discharge_times(1,end) - ...
            discharge_times(1,start_time_index));
                
        overall_dph_std = sqrt(overall_dph);
        
        
        
        if isempty(find(discharge_vals(1,start_time_index:end),1))
            
            overall_median_discharge_size = 0.0;
            
        else
            
            [row, col, vals] = find(discharge_vals(1,start_time_index:end));

            overall_median_discharge_size = median(vals);
            
        end
        
        
    
        for i =1:length(discharge_times(1,:))
            
            time_since_last_hour = discharge_times(1,i) -official_sim_start_time - 60.0*hour_number; 
            
            % if time_since_last_hour > 0, than over an hour has elapsed
            % since the first discharge in this sequence of discharges, and
            % we compute a new discharge rate
            
            % time between this discharge and discharge from last dph value

            if i==first_discharge_index % if it's the first discharge

                live_time = live_time + ...
                    sum(time_length_of_chunk(discharge_times(2,i)));
%                 fprintf('first discharge, %i\n',discharge_times(2,i))

            elseif discharge_times(2,i) ~= discharge_times(2,i-1)

                live_time = live_time + ...
                    sum(time_length_of_chunk(discharge_times(2,i-1)+1:discharge_times(2,i)));
%                 disp(discharge_times(2,i))

            end
            
%             if time_since_last_hour < 0.0
                
%                 % time between this discharge and discharge from last dph value
% 
%                 if i==first_discharge_index % if it's the first discharge
% 
%                     live_time = live_time + ...
%                         sum(time_length_of_chunk(discharge_times(2,i)));
%                     fprintf('first discharge, %i\n',discharge_times(2,i))
% 
%                 elseif discharge_times(2,i) ~= discharge_times(2,i-1)
% 
%                     live_time = live_time + ...
%                         sum(time_length_of_chunk(discharge_times(2,i-1)+1:discharge_times(2,i)));
%                     disp(discharge_times(2,i))
% 
%                 end
            
%             elseif time_since_last_hour > 0.0 && (discharge_times(2,i) < discharge_times(2,i+1))
            if time_since_last_hour > 0.0 && (discharge_times(2,i) < discharge_time_chunk_check_ahead(2,i+1))
                
                first_discharge_index = i+1;
                
                if length(live_time_list) < length(time_length_of_chunk)  
                    
                    live_time_list = [live_time_list live_time];
                    
                end
                
                % if more than an hour passes between discharges, we need
                % to record that zero discharges occurred
                
                 if time_since_last_hour > 60.0                     
                     
                     hours_between_discharge = floor(time_since_last_hour/60);
                     
                     hour_number = hour_number + hours_between_discharge;
                     
                     start_time = start_time + hours_between_discharge*60.;
                     
                     % check if discharge occured, then if 1++ hours passed
                     % without discharges. In this case that first
                     % discharge needs to be accounted for.
                     
                     if i - 1 > start_time_index
                     
                         [row,col,vals] = find(discharge_vals(1,start_time_index:i-1));
                         
%                          vals_list = [vals_list vals];
%                          
%                          if length(vals_per_hour) < length(time_length_of_chunk)
% 
%                             vals_per_hour = [vals_per_hour length(vals)];
% 
%                          end
                         
%                          avg_discharges_this_hour = length(vals)*60.0/(live_time);
                         avg_discharges_this_hour = ceil(length(vals)*60.0/(live_time));

%                          std_this_hour = sqrt(avg_discharges_this_hour);
                         std_this_hour = ceil(sqrt(length(vals))*60.0/live_time);

                         if isempty(find(discharge_vals(1,start_time_index:i-1),1))

                             median_discharge_size_this_hour = 0.0;

                         else

                             [row, col, vals] = find(discharge_vals(1,start_time_index:i-1));
                             
%                              vals_list = [vals_list vals];
%                          
%                              if length(vals_per_hour) < length(time_length_of_chunk)
%                                  
%                                 vals_per_hour = [vals_per_hour length(vals)];
%                                 
%                              end

                             median_discharge_size_this_hour = median(vals);

                         end
                         
                         
                         
                     else
                         
                         avg_discharges_this_hour = 0.0;
                         
                         std_this_hour = 0.0;
                         
                         median_discharge_size_this_hour = 0.0;
                         
                         % no discharge occured, so vals is empty vector
                         vals = [];
                     
                     end
                     
                     %%%
                     
                     vals_list = [vals_list vals];
                         
                     if length(vals_per_hour) + hours_between_discharge -1 < length(time_length_of_chunk)

                        vals_per_hour = [vals_per_hour zeros(1,hours_between_discharge - 1) length(vals)];
                        
                     elseif length(vals_per_hour) < length(time_length_of_chunk)
                         
                         vals_per_hour = [vals_per_hour length(vals)];
                         
                     elseif length(vals_per_hour) > 0
                         
                         vals_per_hour(end) = vals_per_hour(end) + length(vals);

                     end
                     
                     start_time_index = i;
                     
                     % only build dph if there multiple chunks spanning
                     % hours. A chunk only lasts ~ a minute. Sometimes we
                     % want to look at all the chunks summed together. In
                     % this case we only want one dph value.
                     if length(dph) + hours_between_discharge-1 < length(time_length_of_chunk)  
                     
%                          dph = [dph avg_discharges_this_hour zeros(1,hours_between_discharge-1)];
% 
%                          dph_std = [dph_std std_this_hour zeros(1,hours_between_discharge-1)];
%                          
%                          median_discharge_size = ...
%                              [median_discharge_size median_discharge_size_this_hour 
%                              zeros(1,hours_between_discharge-1)];
                         
                         dph = [dph zeros(1,hours_between_discharge-1) avg_discharges_this_hour];

                         dph_std = [dph_std zeros(1,hours_between_discharge-1) std_this_hour];
                         
                         median_discharge_size = ...
                             [median_discharge_size zeros(1,hours_between_discharge-1) 
                             median_discharge_size_this_hour];
                         
                     elseif length(dph) < length(time_length_of_chunk)
                         
                         dph = [dph avg_discharges_this_hour];

                         dph_std = [dph_std std_this_hour];
                         
                         median_discharge_size = ...
                             [median_discharge_size 
                             median_discharge_size_this_hour];
                         
                     elseif length(vals_per_hour) > 0
                         
                         dph(end) = mean(dph(end) + avg_discharges_this_hour);
                         
                         dph_std(end) = mean(dph_std(end)+std_this_hour);
                         
                         median_discharge_size(end) = mean(median_discharge_size(end)+...
                             median_discharge_size_this_hour);
                     
                     end
                    
                     time_since_last_hour = time_since_last_hour - hours_between_discharge*60.0;                                          

                % if discharges did occur in this hour...                
                else

                    %discharges per hour
                    
                    [row,col,vals] = find(discharge_vals(1,start_time_index:i-1));
                    
                    vals_list = [vals_list vals];
                    
                         
                    if length(vals_per_hour) < length(time_length_of_chunk)

                       vals_per_hour = [vals_per_hour length(vals)];

                    end
                    
%                     avg_discharges_this_hour = ...
%                         length(vals)*60.0/live_time;
                    
                    avg_discharges_this_hour = ...
                        ceil(length(vals)*60.0/live_time);

%                     disp(live_time)
%                     std_this_hour = sqrt(avg_discharges_this_hour);
                    std_this_hour = ceil(sqrt(length(vals))*60.0/live_time);

                    if isempty(find(discharge_vals(1,start_time_index:i),1))
                        
                        median_discharge_size_this_hour = 0.0;
                        
                    else
                        
                        [row, col, vals] = find(discharge_vals(1,start_time_index:i));
                       
                        median_discharge_size_this_hour = median(vals);
                        
                    end
                
                end

                % now redefine start_time and start_time_index so that we can
                % look for the avg. # of discharges for the next hour.
                    
                if i ~= length(discharge_times(1,:))
                    
%                     start_time = discharge_times(1,i) + time_since_last_hour;
                    
                    start_time = official_sim_start_time + 60.*hour_number;

%                     start_time_index = i+1;
                    
                    start_time_index = i;
                    
                    hour_number = hour_number+1;
                    
%                     leftover_time = 0;
                    
                    live_time = 0.0;
                    
                    if length(dph) < length(time_length_of_chunk)  
                    
                        dph = [dph avg_discharges_this_hour];

                        dph_std = [dph_std std_this_hour];

                        median_discharge_size = ...
                            [median_discharge_size median_discharge_size_this_hour];
                    
                    end
                    
                end                                

            end
            

        end                      
        
        if official_sim_end_time > discharge_times(1,end)
                
            time_since_last_hour = time_since_last_hour + (official_sim_end_time - discharge_times(1,end));

        end
        
        % if we don't get a full last hour, we can extrapolate the last
        % dph
        
        %%%
        
        if time_since_last_hour > 60.0                     
                     
             hours_between_discharge = floor(time_since_last_hour/60);

             hour_number = hour_number + hours_between_discharge;

             start_time = start_time + hours_between_discharge*60.;

             % check if discharge occured, then if 1++ hours passed
             % without discharges. In this case that first
             % discharge needs to be accounted for.

             if i - 1 > start_time_index

                 [row,col,vals] = find(discharge_vals(1,start_time_index:i-1));
                 
                 vals_list = [vals_list vals];
                         
                 if length(vals_per_hour) < length(time_length_of_chunk)

                    vals_per_hour = [vals_per_hour length(vals)];

                 end
                 
                 avg_discharges_this_hour = ...
                    ceil(length(vals)*60.0/(live_time));

                
                 std_this_hour = ceil(sqrt(length(vals))*60.0/live_time);

                 if isempty(find(discharge_vals(1,start_time_index:i-1),1))

                     median_discharge_size_this_hour = 0.0;

                 else

                     [row, col, vals] = find(discharge_vals(1,start_time_index:i-1));

                     median_discharge_size_this_hour = median(vals);

                 end

             else

                 avg_discharges_this_hour = 0.0;

                 std_this_hour = 0.0;

                 median_discharge_size_this_hour = 0.0;

             end

             start_time_index = i;

             if length(dph) < length(time_length_of_chunk)  
                 
                 dph = [dph avg_discharges_this_hour zeros(1,hours_between_discharge-1)];

                 dph_std = [dph_std std_this_hour zeros(1,hours_between_discharge-1)];

                 median_discharge_size = ...
                     [median_discharge_size median_discharge_size_this_hour ...
                     zeros(1,hours_between_discharge-1)];

             end
             
             time_since_last_hour = time_since_last_hour - hours_between_discharge*60.0;
             
             

        end
        
        %%%
               
        % if more than a half hour passes, add another point.
        
        if time_since_last_hour >= -30.0

%             fraction_last_hour = abs((60.-abs(time_since_last_hour))/60.);
            
            % scale the live time to what it would have been if we went for
            % a complete last hour
%             live_time = live_time/fraction_last_hour;
            
            if length(live_time_list) < length(time_length_of_chunk)  
            
                live_time_list = [live_time_list live_time];
                
            end
            
            [row,col,vals] = find(discharge_vals(1,start_time_index:end));
            
            vals_list = [vals_list vals];
            
                         
            if length(vals_per_hour) < length(time_length_of_chunk)

                vals_per_hour = [vals_per_hour length(vals)];

            end
            
            avg_discharges_this_hour = ...
                ceil(length(vals)*60.0/live_time);
            
%             disp(live_time)
            
            std_this_hour = ceil(sqrt(length(vals))*60.0/live_time);
            
            if isempty(find(discharge_vals(1,start_time_index:end),1))

                median_discharge_size_this_hour = 0.0;

            else

                [row, col, vals] = find(discharge_vals(1,start_time_index:end));

                median_discharge_size_this_hour = median(vals);

            end     
                           
            if length(dph) < length(time_length_of_chunk)  
                
                median_discharge_size = ...
                    [median_discharge_size median_discharge_size_this_hour];

                dph = [dph avg_discharges_this_hour];

                dph_std = [dph_std std_this_hour];
            
            end

        % if there is less than thirty minutes of unaccounted-for discharge
        % time, just roll it into the last point and re-average that 
        % discharge rate
        else
            
            % check if median_discharge_size is populated
            if isempty(find(median_discharge_size,1))   
                
                if isempty(find(discharge_vals(1,start_time_index:end),1))

                    median_discharge_size_this_hour = 0.0;

                else

%                     [row, col, vals] = find(discharge_vals(1,start_time_index:end));
                    [row, col, last_vals] = find(discharge_vals(1,start_time_index:end));
                    
%                     disp(length(vals))
                    vals = [vals last_vals];
%                     disp(length(vals))

                    median_discharge_size_this_hour = median(vals);
                    
                    vals_list = [vals_list vals];
                         
                    if length(vals_per_hour) < length(time_length_of_chunk)

                       vals_per_hour(end) = length(vals);

                    end

                end
                
            else

                    if isempty(find([median_discharge_size(end) discharge_vals(1,start_time_index:end)],1))

                        median_discharge_size_this_hour = 0.0;

                    else

%                         [row, col, vals] = ...
%                             find([median_discharge_size(end) discharge_vals(1,start_time_index:end)]);
                        
                        [row, col, last_vals] = ...
                            find([median_discharge_size(end) discharge_vals(1,start_time_index:end)]);
                        
%                         disp(length(vals))
                        vals = [vals last_vals];
%                         disp(length(vals))

                        median_discharge_size_this_hour = median(vals);
                        
                        vals_list = [vals_list vals];
                         
                        if length(vals_per_hour) < length(time_length_of_chunk)

                           vals_per_hour(end) = length(vals);

                        end

                    end
                    
            end
            
            previous_live_time = live_time_list(end);            
            
            live_time = live_time + previous_live_time;
            
%             disp(live_time)
            
            avg_discharges_this_hour = ceil(length(vals)*60.0/(live_time));
            
            std_this_hour = ceil(sqrt(length(vals))*60.0/live_time);

            if length(dph) > 0

%                 dph(end) = ceil((dph(end)*previous_live_time/60. + ...
%                     (live_time/60.0)*avg_discharges_this_hour)/((previous_live_time+live_time)/60.0));
% 
%                 dph_std(end) = ceil(sqrt( (dph_std(end)*(previous_live_time/60.0))^2 + ...
%                     (std_this_hour*live_time/60.)^2 )/((previous_live_time+live_time)/60.0));

                median_discharge_size(end) = median_discharge_size_this_hour;
                
                
                
                dph(end) = avg_discharges_this_hour;
                
                dph_std(end) = std_this_hour;
                
            else
                
                if length(dph) < length(time_length_of_chunk)  
                    
%                     dph = [dph ceil(avg_discharges_this_hour/(live_time/60.0))];
% 
%                     dph_std = [ dph_std ceil(sqrt(length(vals))/(live_time/60.0))];
                    
                    dph = [dph avg_discharges_this_hour];

                    dph_std = [dph_std std_this_hour];

                    median_discharge_size = [ median_discharge_size median_discharge_size_this_hour];
                
                end
                
            end
            
        end
        
%         live_time_list
        
 % plot histograms of the discharges. 
 
%         dph
%         
%         vals_per_hour

        

        if save_figs == 1
            
            vals_index = 1;
            
            if length(dph) > 1
                
                save_file_path = fullfile(savepath,sprintf('%s_%d.png',plotname,1));

                if vals_per_hour(1) > 0
                    
                    vals_index = vals_per_hour(1)+1;
                    
                    plot_discharge_hist(vals_list(1:vals_per_hour(1)),...
                                title_string,1,save_file_path,xlabel_string,ylabel_string)
                else
                    
                    plot_discharge_hist([],...
                                title_string,1,save_file_path,xlabel_string,ylabel_string)
                    
                end
                        
                for i = 2:length(dph)

                    save_file_path = fullfile(savepath,sprintf('%s_%d.png',plotname,i));
                    
                    if vals_per_hour(i) > 0
                        
                        plot_discharge_hist(vals_list(vals_index:vals_index+vals_per_hour(i)-1),...
                                    title_string,i,save_file_path,xlabel_string,ylabel_string)

                        vals_index = vals_per_hour(i)+1;
                            
                    else
                        
                        plot_discharge_hist([],...
                                title_string,i,save_file_path,xlabel_string,ylabel_string)
                        
                    end

                end
                
            elseif length(dph) == 1    
                
                save_file_path = fullfile(savepath,sprintf('%s_%d.png',plotname,hour_number));

                plot_discharge_hist(vals_list,...
                            title_string,hour_number,save_file_path,xlabel_string,ylabel_string)
                
            end
                        
        end
    
    end