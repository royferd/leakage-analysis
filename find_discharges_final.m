function [discharge_times,discharge_times_cutoff, ...
    discharge_vals, discharge_vals_cutoff,...
    discharge_stdevs,discharge_stdevs_cutoff,...
    discharges_per_chunk,discharges_per_chunk_cutoff,...
    time_length_of_chunk,start_stop_chunk_times,num_not_optimized,...
    gaus_avg_list,gaus_stdev_list] = ...
    find_discharges_final(data_raw,...
    data_weight_raw,time,x_avg_per_chunk,...
    x_stdev_per_chunk,num_chunk,chunk_array,bval,minimize,...
    scale,tail,stdev_discharge_cutoff,title_string,plotname,save_figs,savepath)
    
    %store time and values of the discharges
    xlabel_string = 'leakage current (pA)';
    ylabel_string = 'counts';
    
    discharge_times = [];
    discharge_times_cutoff = [];

    discharge_vals = [];
    discharge_vals_cutoff = [];

    discharge_stdevs = [];
    discharge_stdevs_cutoff = [];

    discharges_per_chunk = [];
    discharges_per_chunk_cutoff = [];
    
    xavg_history = [];
    xstd_history = [];
    xamp_history = [];
    
    gaus_avg_list = [];
    gaus_stdev_list = [];
    
    time_length_of_chunk = zeros(1,num_chunk);
    
    % 2-column vector storing the start and stop times of each chunk.
    start_stop_chunk_times = [];
    
    num_not_optimized = 0;

    % the ramping sorting is imperfect and there is usually a point or two that
    % is an intermediate point that gets grouped in. We don't want to count
    % this as a discharge so we trim the ends of the segments by a few points.
    trim = 5;
    
%     fprintf('---- find_discarges.m: scanning discharges in: %s \n',plotname);

    for i = 1:num_chunk
%     for i = 1:1

        close(gcf);

        % this will tag chunks whose averages and standard deviations
        % cannot be optimized, e.g. if the chunk is too noisy.
        not_optimized = 0;
        
        current_bval = bval;
        
        % we're interested in the number of discharges, the values of the
        % discharges, and the time at which the discharges happened
        count_discharges = 0;
        count_discharges_cutoff = 0;
        
        % first column will store the chunk number that the discharge was
        % found in. The second column will store the discharge value.
        
        discharge_values_pass = zeros(2,1);
        discharge_values_cutoff_pass = zeros(2,1);

        time_of_discharge_pass = zeros(2,1);
        time_of_discharge_cutoff_pass = zeros(2,1);

        stdev_discharge_values_pass = zeros(2,1);
        stdev_discharge_values_cutoff_pass = zeros(2,1);


        in_this_chunk = chunk_array(:,2,i) - chunk_array(:,1,i) - 2*trim;

        xdata = zeros(1,in_this_chunk);

        xdata_stdev_raw = zeros(1,in_this_chunk);

        xdata_stdev = zeros(1,in_this_chunk);

        time_this_chunk = zeros(1,in_this_chunk);
        
        

        for j = (chunk_array(:,1,i)+1)+trim:chunk_array(:,2,i) - trim

            time_this_chunk(j-(chunk_array(:,1,i)+trim)) = ...
                time(j);

            xdata(j-(chunk_array(:,1,i)+trim)) = ...
                data_raw(j)*scale;

            xdata_stdev_raw(j-(chunk_array(:,1,i)+trim)) = ...
                1/sqrt(data_weight_raw(j));

            xdata_stdev(j-(chunk_array(:,1,i)+trim)) = ...
                abs(scale)*xdata_stdev_raw(j-(chunk_array(:,1,i)+trim));
    
        end

        
        min_x = min(xdata);
        
        max_x = max(xdata);                
        
        most_xdata = floor(0.68*length(xdata));
        
        trim_points = ceil((length(xdata)-most_xdata)/2.0);
        
        ordered_xdata = sort(xdata);
        
        if min(xdata) < 0
            
            subset_xdata = ordered_xdata(trim_points:end-trim_points);
            
        else
        
            subset_xdata = ordered_xdata(1:most_xdata);
        
        end
        
        %the minimum number of counts we want in the highest bin
        peak_ht_goal = floor(0.06*length(xdata));
        
        figure0 = figure('visible','off');
        
        h1 = histogram(xdata,'BinWidth',current_bval);

        nbins = h1.NumBins;
        [counts,edges] = histcounts(xdata,nbins);
        
        close(figure0);
        
        for j = 2:11
            
            if max(counts) < peak_ht_goal
                
                current_bval = j*bval;
                
                figure0 = figure('visible','off');
        
                h1 = histogram(xdata,'BinWidth',current_bval);

                nbins = h1.NumBins;
                
                [counts,edges] = histcounts(xdata,nbins);
                
                close(figure0);
                
            else break;
            
            end
            
        end
        
%         if i > 2
        if length(xstd_history) > 2
            
            xstd = median([2.0*std(subset_xdata) xstd_history]);
            
        else
            
            xstd = 2.0*std(subset_xdata);
            
        end
        
        if mode(floor(subset_xdata)) > 0
            
            xavg = mode(current_bval*floor(subset_xdata/current_bval));
            
        else
            
            xavg = mode(current_bval*ceil(subset_xdata/current_bval));
            
        end

        x_hist_data = zeros(1,nbins);
      
        for k = 1:length(edges)-1

            x_hist_data(k) = 0.5*edges(k) + 0.5*edges(k+1);

        end

        [in_range] = find( x_hist_data < xavg + 10*xstd & x_hist_data > xavg - 10*xstd);
        
        opt_counts = [];
        opt_x_hist_data = [];

        for j = 1:length(in_range)

            opt_counts = [opt_counts counts(in_range(j))];

            opt_x_hist_data = [opt_x_hist_data x_hist_data(in_range(j))];

        end

        % view minimization output
      %  options = optimset('Display','iter','FunValCheck','on','TolX',1e-5,'TolFun',1e-5);
        options = optimset('FunValCheck','on','TolX',1e-4,'TolFun',1e-4,'Display','none');

        max_count = find(max(counts));
        
        [in_mid_bin] = find( xdata < xavg + 1.6*current_bval & xdata > xavg - 1.6*current_bval);
        
        guess_amplitude = sqrt(length(in_mid_bin)/3);

        if (minimize == 1)

            % optimize the amplitude of the Gaussian for fixed average and stdev:
        %    fun = @(x)gaus_min_amp(x,opt_x_hist_data,opt_counts,xavg,xstd);
            fun = @(x)gaus_min_amp(x,opt_x_hist_data,opt_counts,xavg,xstd);

            x0 = 0.75*max(counts);

            bestx = fminsearch(fun,x0,options);

            gaus_a = bestx(1);

            gaus_avg = xavg;

            gaus_stdev = xstd;

        elseif (minimize == 2)
            

            % optimize the amplitude and average of the Gaussian for fixed stdev:    
            fun = @(x)gaus_min_amp_avg(x,opt_x_hist_data,opt_counts,xavg,xstd);

            x0 = [max(opt_counts) xavg];

            bestx = fminsearch(fun,x0,options);

            gaus_a = bestx(1);

            gaus_avg = bestx(2);

            gaus_stdev = xstd;
            
        elseif (minimize == 3)

            % optimize the amplitude and average of the Gaussian for fixed stdev:    
            fun = @(x)gaus_min_amp_avg_stdev(x,opt_x_hist_data,opt_counts);
         
            x0 = [sqrt(0.75*max(counts)) xavg xstd];
            
            x0 = [guess_amplitude xavg xstd];

            bestx = fminsearch(fun,x0,options);

            gaus_a = bestx(1)^2;

            gaus_avg = bestx(2);

            gaus_stdev = abs(bestx(3));
            
        end   

        if gaus_stdev > 2.0*xstd || abs(xavg - gaus_avg)/xstd > 2.0
            
            not_optimized = 1;
            
            num_not_optimized = num_not_optimized + 1;
            
            gaus_a = guess_amplitude^2;
            
            gaus_avg = xavg;
            
            gaus_stdev = xstd;
            
        end

        %%%%%  Find discharges with 5 sigma test

        for j = (chunk_array(:,1,i)+1)+trim:chunk_array(:,2,i) - trim

            %one tailed test. we only count values greater than xavg + 5 stdev
            if (tail == 1) 

                if ( xdata(j-(chunk_array(:,1,i)+trim)) > gaus_avg + 5*gaus_stdev)

                    count_discharges = count_discharges + 1;

                    % value of discharge - average leakage current
                    discharge_values_pass(1,end+1) = ...
                        xdata(j-(chunk_array(:,1,i)+trim)) - gaus_avg;

                    discharge_values_pass(2,end) = i;

                    time_of_discharge_pass(1,end+1) = ...
                        time_this_chunk(j-(chunk_array(:,1,i)+trim));

                    time_of_discharge_pass(2,end) = i;

                    stdev_discharge_values_pass(1,end+1) = ...
                        xdata_stdev(j-(chunk_array(:,1,i)+trim)) + gaus_stdev;

                    stdev_discharge_values_pass(2,end) = i;

                end
                
                % cutoff_test
                if ( xdata(j-(chunk_array(:,1,i)+trim)) > gaus_avg + stdev_discharge_cutoff)

                    count_discharges_cutoff = count_discharges_cutoff + 1;

                    % value of discharge - average leakage current
                    discharge_values_cutoff_pass(1,end+1) = ...
                        xdata(j-(chunk_array(:,1,i)+trim)) - gaus_avg;

                    discharge_values_cutoff_pass(2,end) = i;

                    time_of_discharge_cutoff_pass(1,end+1) = ...
                        time_this_chunk(j-(chunk_array(:,1,i)+trim));

                    time_of_discharge_cutoff_pass(2,end) = i;

                    stdev_discharge_values_cutoff_pass(1,end+1) = ...
                        xdata_stdev(j-(chunk_array(:,1,i)+trim)) + gaus_stdev;

                    stdev_discharge_values_cutoff_pass(2,end) = i;
                        
                end

            elseif tail ==2

                if ( xdata(j-(chunk_array(:,1,i)+trim)) < gaus_avg - 5*gaus_stdev || ...
                        xdata(j-(chunk_array(:,1,i)+trim)) > gaus_avg + 5*gaus_stdev)

                    count_discharges = count_discharges + 1;

                    % value of discharge - average leakage current
                    discharge_values_pass(1,end+1) = xdata(j-(chunk_array(:,1,i)+trim));
                    
                    discharge_values_pass(2,end) = i;

                    time_of_discharge_pass(1,end+1) = time_this_chunk(j-(chunk_array(:,1,i)+trim));
                    
                    time_of_discharge_pass(2,end) = i;

                    stdev_discharge_values_pass(1,end+1) = ...
                        xdata_stdev(j-(chunk_array(:,1,i)+trim));
                    
                    stdev_discharge_values_pass(2,end) = i;

                end

            end

        end               

        discharge_values = discharge_values_pass(:,2:end);
        
        time_of_discharge = time_of_discharge_pass(:,2:end);
        
        stdev_discharge_values = stdev_discharge_values_pass(:,2:end);

        discharge_values_cutoff = discharge_values_cutoff_pass(:,2:end);
        
        time_of_discharge_cutoff = time_of_discharge_cutoff_pass(:,2:end);
        
        stdev_discharge_values_cutoff = stdev_discharge_values_cutoff_pass(:,2:end);
        
        %%%%% cutoff test end

        %input arg "options" will display minimization output
        
        num_fit_points = 2*length(xdata);

        gaus_range = 10*gaus_stdev;

        step = gaus_range/(num_fit_points);

        gaus_fit = zeros(num_fit_points,2);

        for l = 1:num_fit_points
            fit_x = l*step + gaus_avg - 0.5*gaus_range;
            gaus_fit(l,1) = fit_x;
            gaus_fit(l,2) = gaus_a*exp(-((fit_x - gaus_avg)/(sqrt(2)*gaus_stdev))^2);

        end

        %text box assignments
        inside_plot = [0.1375 0.73 0.5 0.19]; % edges: x y width height
        
        if save_figs == 1
            
            if not_optimized == 1
                
                title_string_full =sprintf('%s %d / %d (not optimized)',title_string,i,num_chunk);
                
              %  save_file_path = fullfile(savepath,sprintf('%s_%d_no_opt.png',plotname,i));
                
               % save_file_path_discharge = fullfile(savepath,sprintf('%s_%d_discharges_no_opt.png',plotname,i));
                
                save_file_path = fullfile(savepath,sprintf('%s_%d_no_opt',plotname,i));
                
                save_file_path_discharge = fullfile(savepath,sprintf('%s_%d_discharges_no_opt',plotname,i));
                
            else

                title_string_full =sprintf('%s %d / %d',title_string,i,num_chunk);
                
                %save_file_path = fullfile(savepath,sprintf('%s_%d.png',plotname,i));
                save_file_path = fullfile(savepath,sprintf('%s_%d',plotname,i));
                %save_file_path = fullfile(savepath,sprintf('%s_%d.emf',plotname,i));
            
                %save_file_path_discharge = fullfile(savepath,sprintf('%s_%d_discharges.png',plotname,i));
                save_file_path_discharge = fullfile(savepath,sprintf('%s_%d_discharges',plotname,i));
            
                
            end
            
            smallest_font_size = 13;
            
            figure1 = figure('visible','off');    
            
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 5.5 3.];
            plotstyle = 'linear';
%             plotstyle = 'linlog';
            
            if plotstyle == 'linear'
            
            h2 = histogram(xdata,'BinWidth',current_bval,'FaceColor',[0.3010 0.7450 0.9330]); hold on;
            
            plot(gaus_fit(:,1),gaus_fit(:,2),'r-','LineWidth',2.0);
            
            max_h2 = max(h2.Values);
            
            % find the biggest value on the histogram so we can scale correctly
            graph_max = max([max_count max_h2 gaus_a]);
            
            % round_to rounds the x axis bounds
            round_to = 20;
            if ((gaus_avg - 6*gaus_stdev)/5.0 > 1.0)
                
%                 axis ([ min([subset_xdata(1) 5.0*floor((gaus_avg - 6*gaus_stdev)/5.0)]) max([5.0*ceil((gaus_avg + 6*gaus_stdev)/5.0) subset_xdata(end)]) 0 5.*ceil(1.1*(graph_max)/5.)]);
                axis ([ round_to*min(floor([subset_xdata(1)/round_to (gaus_avg - 6*gaus_stdev)/round_to]))  round_to*ceil(ordered_xdata(end)/round_to) 0 5.*ceil(1.1*(graph_max)/5.)]);

            else               
                
                axis ([ min([subset_xdata(1) floor(gaus_avg - 6*gaus_stdev)]) max([ceil(gaus_avg + 6*gaus_stdev) subset_xdata(end)]) 0 5.*ceil(1.1*(graph_max)/5.)]);

            end
            
            
            ax = gca;
            ax.TickDir = 'out'; % make ticks point out
            
            ax.FontSize = smallest_font_size;
            outerpos = ax.OuterPosition;
            
%            ax.XScale='log';
            
            ti = ax.TightInset;
            
            left = outerpos(1);
            bottom = outerpos(2) + 0.25*ti(2);
            ax_width = outerpos(3);
            ax_height = outerpos(4) - 0.25*ti(2);
            
            ax.OuterPosition = [left bottom ax_width ax_height];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% try to split into linear and log
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% section

            elseif plotstyle == 'linlog'

                ax1 = subplot(121); %the linear part of the axis
                h2 = histogram(xdata(xdata<=gaus_fit(end,1)),'BinWidth',...
                    current_bval,'FaceColor',[0.3010 0.7450 0.9330]); hold on;
                max_h2 = max(h2.Values); hold on;

                % find the biggest value on the histogram so we can scale correctly
                graph_max = max([max_count max_h2 gaus_a]);

                plot(gaus_fit(:,1),gaus_fit(:,2),'r-','LineWidth',2.0);

                ax2 = subplot(122); %log part of the axis
                histogram(xdata(xdata>gaus_fit(end,1)),'BinWidth',...
                    current_bval,'FaceColor',[0.9330 0.3010 0.7450]); hold on;


                set(ax1,'units','normalized','position',[0.1 0.1 0.4 0.8]);
                set(ax2,'units','normalized','position',[0.5 0.1 0.4 0.8]);

                set(ax2,'xscale','log','xlim',[gaus_fit(end,1) max(xdata)],'yticklabel','');
                %set([ax1 ax2],'ylim',[-1.5 1.5],'ytick',-1.5:0.5:1.5,'box','off');
                set([ax1 ax2],'ylim',[0 5.*ceil(1.1*(graph_max)/5.)]);

                set(ax2,'yticklabel','');
                uistack(ax1,'top');
    %             grid(ax2,'on');

                ax1.TickDir = 'out'; % make ticks point out

                ax1.FontSize = smallest_font_size;

                ax2.TickDir = 'out'; % make ticks point out

                ax2.FontSize = smallest_font_size;
                
            end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             annotation('textbox',inside_plot,'String',...
%                 sprintf('%.1f \\pm %.1f pA',gaus_avg,gaus_stdev),...
%                 'FontSize',smallest_font_size,'BackgroundColor',[1 1 1],'FitBoxToText','on');
            annotation('textbox',inside_plot,'String',...
                sprintf('%.1f \\pm %.1f pA\nN = %i + %i',gaus_avg,gaus_stdev,length(xdata)-count_discharges,count_discharges),...
                'FontSize',smallest_font_size,'BackgroundColor','none',...
                'FitBoxToText','on','EdgeColor','none');

            
            title(title_string_full,'FontSize',smallest_font_size)
            xlabel(xlabel_string,'FontSize',smallest_font_size)
            ylabel(ylabel_string,'FontSize',smallest_font_size)            
            
            
            %print (save_file_path,'-dmeta','-r600');
            print(save_file_path,'-dpng');    
            %print(sprintf('%s.svg',save_file_path),'-dsvg')
            
        end

        % check if we had zero discharges this chunk
        if count_discharges == 0        

            discharge_times = [discharge_times [time_this_chunk(end) ; i]];
            
            discharge_vals = [discharge_vals [0.0 ; i]];
            
            discharge_stdevs = [discharge_stdevs [0.0 ; i]];

            discharges_per_chunk = [ discharges_per_chunk [0 ; i]];
        
        else   
            
            discharge_times = [discharge_times time_of_discharge];

            discharge_vals = [discharge_vals discharge_values];

            discharge_stdevs = [discharge_stdevs stdev_discharge_values];

            discharges_per_chunk = [ discharges_per_chunk [count_discharges ; i]];
        
        end
        
        if count_discharges_cutoff == 0        

            discharge_times_cutoff = [discharge_times_cutoff [time_this_chunk(end) ; i]];
            
            discharge_vals_cutoff = [discharge_vals_cutoff [0.0 ; i]];
            
            discharge_stdevs_cutoff = [discharge_stdevs_cutoff [0.0 ; i]];

            discharges_per_chunk_cutoff = [ discharges_per_chunk_cutoff [0 ; i]];            
        
        else   
            
            discharge_times_cutoff = [discharge_times_cutoff time_of_discharge_cutoff];

            discharge_vals_cutoff = [discharge_vals_cutoff discharge_values_cutoff];

            discharge_stdevs_cutoff = [discharge_stdevs_cutoff stdev_discharge_values_cutoff];

            discharges_per_chunk_cutoff = ...
                [ discharges_per_chunk_cutoff [count_discharges_cutoff ; i ]];
        
        end
        
        time_length_of_chunk(i) = time_this_chunk(end) - time_this_chunk(1);
        
        start_stop_chunk_times = [start_stop_chunk_times; [time_this_chunk(1),time_this_chunk(end)]];
        
        if not_optimized == 0
            
            xavg_history = [xavg_history gaus_avg];

            xstd_history = [xstd_history gaus_stdev];

            xamp_history = [xamp_history gaus_a];
        
        end
        
        if length(xavg_history) > 2
            
            xavg_history = xavg_history(end-2:end);
            
            xstd_history = xstd_history(end-2:end);
            
        end       

    gaus_avg_list = [gaus_avg_list gaus_avg];
    
    gaus_stdev_list = [gaus_stdev_list gaus_stdev/sqrt(in_this_chunk)];
        
    end 
%     discharges_per_chunk
%     fprintf('%i discharges this chunk',discharges_per_chunk(end));
    
    if num_not_optimized > 0
        
        fprintf('---- find_discarges.m: %i chunk distributions not optimized in %s \n',...
            num_not_optimized,plotname);
        
    end