%{

chunkify takes a set of data with weights and defined sections, or chunks,
and gives an array of weighted and averaged points for each chunk.

%}
function [y_avg_set_avg_chunk, y_avg_set_stdev_chunk,...
    y_avg_set_avg, y_avg_set_stdev,...
    y_avg_set_stdev_avg_chunk, y_avg_set_stdev_stdev_chunk, y_stdev_of_stdev_of_chunk,...
    y_avg_stdev_of_chunk]...
    = chunkify(num_files, num_chunks,...
    num_set_points,chunk_array, y_avg_set_raw,...
    y_weight_avg_set_raw,y_avg_scale)

    % 7/11/2017 can't think of a way to store chunk averaging for more than 
    % 1 data set

    % the first/last few points in each chunk may be intermediate points
    % from a ramp down/up. When chunkifying, trim off these points.
    trim = 2;
    
    %calculate averages
    y_avg_set_raw_avg = zeros(num_files,1);
    y_avg_set_raw_avg_chunk = zeros(num_files,num_chunks,1);
    y_avg_set_raw_stdev = zeros(num_files,1);

    y_avg_set_avg = zeros(num_files,1);
    y_avg_set_avg_chunk = zeros(num_files,num_chunks,1);

    y_avg_set_stdev = zeros(num_files,1);
    y_avg_set_stdev_chunk = zeros(num_files,num_chunks,1);

    y_avg_set_raw_stdev_raw_chunk = zeros(num_files,num_chunks,1);

    y_avg_set_raw_avg_wt = zeros(num_files,num_set_points);

    y_weight_set_raw_stdev_raw_chunk = zeros(num_files,num_chunks);
    y_avg_set_raw_stdev_avg_raw_chunk_wt = zeros(num_files,num_chunks,1);

    y_avg_set_raw_stdev_stdev_raw_chunk = zeros(num_files,1);
    y_avg_set_stdev_stdev_chunk = zeros(num_files,1);

    y_avg_set_raw_stdev_avg_raw_chunk = zeros(num_files,1);
    y_avg_set_stdev_avg_chunk = zeros(num_files,1);
    
    y_stdev_of_stdev_of_chunk = zeros(num_files,1);
    y_stdev_of_stdev_of_chunk_raw = zeros(num_files,1);
    
    y_avg_stdev_of_chunk = zeros(num_files,1);
    y_avg_stdev_of_chunk_raw = zeros(num_files,1);
    
    y_stdev_set_raw = 1./sqrt(y_weight_avg_set_raw);
    
    y_avg_set_raw_trimmed = [];
    
    y_weight_avg_set_raw_trimmed = [];
    
    y_avg_set_raw_avg_wt_trimmed = [];

    for i = 1:num_files
        
% % %         % stdev of the the average of all the chunks
% % %         y_avg_set_raw_stdev(i) = std(y_avg_set_raw(i,:),y_weight_avg_set_raw(i,:));

        for j = 1:num_set_points(i)
            
            % each point avg value * the weight
            y_avg_set_raw_avg_wt(i,j) = y_avg_set_raw(i,j)*y_weight_avg_set_raw(i,j); 
            
        end

% % %         % weighted average of all the points
% % %         y_avg_set_raw_avg(i) = sum(y_avg_set_raw_avg_wt(i,:))/sum(y_weight_avg_set_raw(i,:));
    
    end

    %7/12 change chunk code when I figure out how to handle multiple data sets...
    for i = 1:num_files

        for j = 1:num_chunks(i)

            chunk_begin = chunk_array(i,1,j) +1 + trim;
            
            chunk_end = chunk_array(i,2,j) - trim;
            
%             for k = chunk_begin:chunk_end
                
                y_avg_set_raw_trimmed = ...
                    [y_avg_set_raw_trimmed y_avg_set_raw(chunk_begin:chunk_end)];
                
                y_weight_avg_set_raw_trimmed = ...
                    [y_weight_avg_set_raw_trimmed y_weight_avg_set_raw(chunk_begin:chunk_end)];
                
                % weighted stdev of the average of the chunk
                y_avg_set_raw_stdev_raw_chunk(i,j) = ...
                    ( std(y_avg_set_raw(i,chunk_begin:chunk_end),...
                    y_weight_avg_set_raw(i,chunk_begin:chunk_end)));

                % weight of the chunk
                y_weight_set_raw_stdev_raw_chunk(i,j) = ...
                    (y_avg_set_raw_stdev_raw_chunk(i,j))^(-2);
                
                % stdev of the stdev of the chunk
                y_stdev_of_stdev_of_chunk_raw(i,j) = std(y_stdev_set_raw(chunk_begin:chunk_end));
                
                %avg of the stdev of the chunk
                y_avg_stdev_of_chunk_raw(i,j) = mean(y_stdev_set_raw(chunk_begin:chunk_end));
                
                
%             end
            
            % the weighted stdev of the average chunk times its weight? seems wrong
            y_avg_set_raw_stdev_avg_raw_chunk_wt(i,j) = ...
                y_avg_set_raw_stdev_raw_chunk(i,j)*y_weight_set_raw_stdev_raw_chunk(i,j);
            
            
        end
        
        num_set_points_trimmed = length(y_avg_set_raw_trimmed);
        
        for k = 1:num_set_points_trimmed
            
            % each point avg value * the weight
            y_avg_set_raw_avg_wt_trimmed(k) = ...
                y_avg_set_raw_trimmed(k)*y_weight_avg_set_raw_trimmed(k); 
            
        end
        
        y_avg_set_raw_stdev(i) = ...
            std(y_avg_set_raw_trimmed(:),y_weight_avg_set_raw_trimmed(i,:));

        % weighted average of all the points
%         y_avg_set_raw_avg(i) = ...
%             sum(y_avg_set_raw_avg_wt_trimmed(:))/sum(y_weight_avg_set_raw_trimmed(i,:));
        y_avg_set_raw_avg(i) = ...
            sum(y_avg_set_raw_avg_wt(:))/sum(y_weight_avg_set_raw(i,:));
        
        % stdev of the stdevs of the chunks. Seems correct
        y_avg_set_raw_stdev_stdev_raw_chunk(i) = ...
            1/sqrt(sum(y_weight_set_raw_stdev_raw_chunk(i,:)));

        % the weighted average stdev of all the chunks. this is the thing
        % that seems off.       
%         y_avg_set_raw_stdev_avg_raw_chunk(i) = ...
%             sum(y_avg_set_raw_stdev_avg_raw_chunk_wt(i,...
%             :))/sum(y_weight_set_raw_stdev_raw_chunk(i,:));
        
%         y_avg_set_raw_stdev_avg_raw_chunk(i) = ...
%             mean(y_avg_stdev_of_chunk);
        
                
        y_avg_set_raw_stdev_avg_raw_chunk(i) = ...
            mean(y_avg_stdev_of_chunk_raw);

        last_chunk = 1;
        
        for j = 1:num_chunks(i)
            
%             chunk_begin = chunk_array(i,1,j) +1 + trim;
%             
%             chunk_end = chunk_array(i,2,j) - trim;

            chunk_begin = chunk_array(i,1,j) +1;
            
            chunk_end = chunk_array(i,2,j);
            
            this_chunk = chunk_end - chunk_begin;

            % the weighted average of each chunk
            y_avg_set_raw_avg_chunk(i,j) =( sum(y_avg_set_raw_avg_wt(i,...
                chunk_begin:chunk_end))/sum(y_weight_avg_set_raw(i,chunk_begin:chunk_end)));
            
%              y_avg_set_raw_avg_chunk(i,j) =...
%                  ( sum(y_avg_set_raw_avg_wt_trimmed(last_chunk:this_chunk + (last_chunk-1)))/sum(y_weight_avg_set_raw_trimmed(last_chunk:this_chunk + (last_chunk-1))));
%  
            last_chunk = this_chunk + 1;
            
        end
        
        

        %overall chunk averages and stdevs of avgs
        %y_avg_set_avg(i) = y_avg_set_raw_avg(i)*y_avg_scale;
        y_avg_set_avg(i) = y_avg_set_raw_avg(i)*y_avg_scale;
        
        
        y_avg_set_stdev(i) = y_avg_set_raw_stdev(i)*abs(y_avg_scale);

        %overall chunk averages of stdevs and stdevs of stdevs
        y_avg_set_stdev_avg_chunk(i) = y_avg_set_raw_stdev_avg_raw_chunk(i)*abs(y_avg_scale);
        
        y_avg_set_stdev_stdev_chunk(i) = ...
            y_avg_set_raw_stdev_stdev_raw_chunk(i)*abs(y_avg_scale);

        for j = 1:num_chunks(i)
            
            y_avg_set_avg_chunk(i,j) = y_avg_set_raw_avg_chunk(i,j)*y_avg_scale;

            y_avg_set_stdev_chunk(i,j) = y_avg_set_raw_stdev_raw_chunk(i,j)*abs(y_avg_scale);
            
            y_stdev_of_stdev_of_chunk(i,j) = y_stdev_of_stdev_of_chunk_raw(i,j)*abs(y_avg_scale);
            
            y_avg_stdev_of_chunk(i,j) = y_avg_stdev_of_chunk_raw(i,j)*abs(y_avg_scale);

        end

    end


end