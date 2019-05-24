function [discharge_times,discharge_times_cutoff, ...
    discharge_vals, discharge_vals_cutoff,...
    discharge_stdevs,discharge_stdevs_cutoff,...
    discharges_per_chunk,discharges_per_chunk_cutoff,...
    time_length_of_chunk] = ...
    find_discharges_final(data_raw,...
    data_weight_raw,time,x_avg_per_chunk,...
    x_stdev_per_chunk,num_chunk,chunk_array,bval,minimize,...
    scale,tail,stdev_discharge_cutoff,title_string,plotname,save_figs,savepath)
    
% function [discharge_times, discharge_vals, ...
%     discharge_stdevs,discharges_per_chunk,dph,dph_std,median_discharge_size,...
%     overall_dph,overall_dph_std,overall_median_discharge_size] = ...
%     find_discharges_final(data_raw,...
%     data_weight_raw,time,x_avg_per_chunk,...
%     x_stdev_per_chunk,num_chunk,chunk_array,bval,minimize,...
%     scale,tail,stdev_discharge_cutoff,title_string,plotname,save_figs,savepath)
    
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
    
    time_length_of_chunk = zeros(1,num_chunk);

    % the ramping sorting is imperfect and there is usually a point or two that
    % is an intermediate point that gets grouped in. We don't want to count
    % this as a discharge so we trim the ends of the segments by a few points.
    trim = 5;
    fprintf('---- find_discarges.m: scanning discharges in: %s \n',plotname);
    
%     if (minimize == 1)
%         disp('numerical optimization will vary Gaussian amplitude, keeping average and stdev fixed.');
% 
%     elseif (minimize == 2)
%         disp('numerical optimization will vary Gaussian amplitude and average, keeping stdev fixed.')
% 
%     elseif (minimize == 3)
%         disp('numerical optimization will vary Gaussian amplitude, average, and stdev.')
% 
%     end
    
%     fprintf('number of input arguments = %d \n',nargin);
%     fprintf('number of output arguments = %d \n',nargout);

    for i = 1:num_chunk

        close(gcf);

        % we're interested in the number of discharges, the values of the
        % discharges, and the time at which the discharges happened
        count_discharges = 0;
        count_discharges_cutoff = 0;

%{        
        discharge_values_pass = zeros(1,1);
        discharge_values_cutoff_pass = zeros(1,1);

        time_of_discharge_pass = zeros(1,1);
        time_of_discharge_cutoff_pass = zeros(1,1);

        stdev_discharge_values_pass = zeros(1,1);
        stdev_discharge_values_cutoff_pass = zeros(1,1);
%}        
        
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
        
        most_xdata = floor(0.99*length(xdata));
        ordered_xdata = sort(xdata);
        plot_subset_xdata = ordered_xdata(1:most_xdata);
        
      %  xavg = 15.;
        
        
%        xavg = x_avg_per_chunk(i);
        xavg = mode(floor(ordered_xdata(1:most_xdata)));

        xmode_index = find(floor(ordered_xdata(1:most_xdata)) == xavg);
        xmode_middle = xmode_index(floor(length(xmode_index)/2.));
        
        
        one_std = floor(0.34*most_xdata);
        
        if most_xdata - xmode_middle + 1 >= one_std
            
            xstd = ordered_xdata(xmode_middle + one_std) - xavg;
            
        else
            
            
            xstd = xavg - ordered_xdata(xmode_middle - one_std);
            
        end
        
        %xstd = (ordered_xdata(floor(0.95*most_xdata)) - ordered_xdata(1))/4.0;
         
      %  xstd = x_stdev_per_chunk(i);
        
        
        
        %figure0 = figure;
        figure0 = figure('visible','off');
        
        h1 = histogram(xdata,'BinWidth',bval);
%       h1 = histcounts(xdata,'BinWidth',bval);
        nbins = h1.NumBins;
        [counts,edges] = histcounts(xdata,nbins);
        
%         fprintf('maximum counts = %d \n',max(counts));        
%         fprintf('x avg for chunk %d = %f \n',i,xavg);
%         fprintf('x std for chunk %d = %f \n',i,xstd);
            


        x_hist_data = zeros(1,nbins);

%         for k = 1:nbins
% 
%             x_hist_data(k) = min_x+bval*(k-0.5);
% 
%         end
        
        for k = 1:length(edges)-1

            x_hist_data(k) = 0.5*edges(k) + 0.5*edges(k+1);

        end


        % indices of bins that we want to fit our Gaussian to:

      %  [in_range] = find( x_hist_data < xavg + 5*xstd & x_hist_data > xavg - 5*xstd);
        [in_range] = find( x_hist_data < xavg + 10*xstd & x_hist_data > xavg - 10*xstd);


        % the counts in the bins of interest:

        opt_counts = [];
        opt_x_hist_data = [];

        for j = 1:length(in_range)

         %   opt_counts = counts(min(in_range):max(in_range));

            opt_counts = [opt_counts counts(in_range(j))];

            % the x-values of the bins of interest:
            %opt_x_hist_data = x_hist_data(min(in_range):max(in_range));

            opt_x_hist_data = [opt_x_hist_data x_hist_data(in_range(j))];

        end

%         disp('opt_x_hist_data:')
%         disp(opt_x_hist_data)
%         disp('opt_counts: ')
%         disp(opt_counts)


        % view minimization output
      %  options = optimset('Display','iter','FunValCheck','on','TolX',1e-5,'TolFun',1e-5);
        options = optimset('FunValCheck','on','TolX',1e-4,'TolFun',1e-4);

        max_count = find(max(counts));


        if (minimize == 1)

            % optimize the amplitude of the Gaussian for fixed average and stdev:
        %    fun = @(x)gaus_min_amp(x,opt_x_hist_data,opt_counts,xavg,xstd);
            fun = @(x)gaus_min_amp(x,opt_x_hist_data,opt_counts,xavg,xstd);

            x0 = 0.75*max(counts);
          %  x0 = max(counts);

            bestx = fminsearch(fun,x0,options);

            gaus_a = bestx(1);
           % gaus_a = 500.;

            gaus_avg = xavg;

            gaus_stdev = xstd;
         %   gaus_stdev = 1.5;

        elseif (minimize == 2)
            

            % optimize the amplitude and average of the Gaussian for fixed stdev:    
            fun = @(x)gaus_min_amp_avg(x,opt_x_hist_data,opt_counts,xavg,xstd);

%            x0 = [0.75*max(counts) mean(opt_x_hist_data)];
            x0 = [max(opt_counts) xavg];

            bestx = fminsearch(fun,x0,options);

            gaus_a = bestx(1);

            gaus_avg = bestx(2);

            gaus_stdev = xstd;
            
        elseif (minimize == 3)

            % optimize the amplitude and average of the Gaussian for fixed stdev:    
            fun = @(x)gaus_min_amp_avg_stdev(x,opt_x_hist_data,opt_counts);

         %   x0 = [sqrt(0.75*max(counts)) xavg xstd];
         
            x0 = [sqrt(0.75*max(counts)) xavg xstd];
            
          %  x0 = [max(counts) xavg 1.5];

            bestx = fminsearch(fun,x0,options);

            gaus_a = bestx(1)^2;

            gaus_avg = bestx(2);

            gaus_stdev = abs(bestx(3));
            
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

            else

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

%         fprintf('chunk no: %d \n',i);
%         fprintf('# discharges: %d \n',count_discharges);
%         fprintf('%.2f pA at %.2f m \n',discharge_values,time_of_discharge);


        %%%%% Find discharges with stdev cutoff test
        
        for j = (chunk_array(:,1,i)+1)+trim:chunk_array(:,2,i) - trim

            %one tailed test. we only count values greater than xavg + 5 stdev
            if (tail == 1) 

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

            end

        end               

        discharge_values_cutoff = discharge_values_cutoff_pass(:,2:end);
        time_of_discharge_cutoff = time_of_discharge_cutoff_pass(:,2:end);
        stdev_discharge_values_cutoff = stdev_discharge_values_cutoff_pass(:,2:end);
        
        %%%%% cutoff test end

        %input arg "options" will display minimization output
    %    bestx = fminsearch(fun,x0,options);

        num_fit_points = 2*length(xdata);

        gaus_range = 10*gaus_stdev;

        step = gaus_range/(num_fit_points);

        gaus_fit = zeros(num_fit_points,2);

        for l = 1:num_fit_points
            fit_x = l*step + gaus_avg - 0.5*gaus_range;
            gaus_fit(l,1) = fit_x;
            gaus_fit(l,2) = gaus_a*exp(-((fit_x - gaus_avg)/(sqrt(2)*gaus_stdev))^2);

        end


%         fprintf('loop: %d \n',i);
%         fprintf('gaussian average: %.1f \n',gaus_avg);
%         fprintf('xstd: %.1f \n',gaus_stdev);
%         fprintf('number of discharges: %d \n',count_discharges);


        %%%%%%%%%%%%%%%%%%%%%%% gaussian fit of ramp data %%%%%%%%%%%%%%%%%%%%%

        %edges = [-70000 -700:.8:-550 0];
        %histogram(lcm1_avg_ramp_up,edges); hold on;


        % Do log10 but keep sign
        %xlog = sign(x).*log10(abs(x));
        % Just to get axis limits
        %plot(xlog,y,'o')
        histogram(xdata);
        %ax1=gca;
        %ax1.XScale='log';
        % Get limits
        lims = xlim;
        wdth = diff(lims);
        % Wrap negative data around to positive side
        %xlog(xlog<0) = xlog(xlog<0) + wdth;

        %x(x<0) = x(x<0) + wdth;
        %x = x + wdth;
        
        if save_figs == 1
            % Plot
            %plot(xlog,y,'o')

            title_string_full =sprintf('%s %d / %d',title_string,i,num_chunk);                      
            
            figure1 = figure('visible','off');                       
%            figure1 = figure;                        
            
            h2 = histogram(xdata,'BinWidth',bval); hold on;
            max_h2 = max(h2.Values);
%             fprintf('maximum count from plotted histogram = %d \n',max_h2);
                                            
            
            
          %  histogram(xdata,nbins); hold on;
           % histogram(xdata,'BinWidth',bval); hold on;
            %histogram(opt_x_hist_data,'BinWidth',bval); hold on;
            
            plot(gaus_fit(:,1),gaus_fit(:,2),'r-','LineWidth',2.0);
            
%             fprintf('chunk # %d \n', i);
%             fprintf('gaus amp = %f \n',gaus_a);
%             fprintf('gaus avg = %f \n',gaus_avg);
%             fprintf('gaus stdev = %f \n',gaus_stdev);
%             fprintf('x range min = %f \n',floor(gaus_avg - 6*gaus_stdev));
%             fprintf('x range max = %f \n',ceil(gaus_avg + 6*gaus_stdev));
%             fprintf(' y range max = %f \n',5.*ceil(1.1*max(counts)/5.));

            if ((gaus_avg - 6*gaus_stdev)/5.0 > 1.0)

                axis ([ 5.0*floor((gaus_avg - 6*gaus_stdev)/5.0) 5.0*ceil((gaus_avg + 6*gaus_stdev)/5.0) 0 5.*ceil(1.1*(max_h2)/5.)]);

            else
               
                axis ([ floor(gaus_avg - 6*gaus_stdev) ceil(gaus_avg + 6*gaus_stdev) 0 5.*ceil(1.1*(max_h2)/5.)]);

            end

            ax = gca;
            ax.TickDir = 'out'; % make ticks point out
            ax.FontSize = 14;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            
            left = outerpos(1);
            bottom = outerpos(2) + 0.25*ti(2);
            ax_width = outerpos(3);
            ax_height = outerpos(4) - 0.25*ti(2);
            
            ax.OuterPosition = [left bottom ax_width ax_height];
            
            title(title_string_full,'FontSize',16)
            xlabel(xlabel_string,'FontSize',16)
            ylabel(ylabel_string,'FontSize',16)
            %ax.XScale='log';

            % Mess with ticks
            %disp('original tick labels:')
            %tck = get(gca,'XTick')'

            % Shift those that were wrapped from negative to positive (above) back 
            % to their original values
            %tck(tck>lims(2)) = tck(tck>lims(2)) - wdth;
            %tck = tck - wdth;

            % Convert to string, then remove any midpoint
            %tcklbl = num2str(tck);
            %tcklbl(tck==lims(2),:) = ' ';
            %disp('new tick labels:')
            %tcklbl
            % Update tick labels
            %set(gca,'XTickLabel',tcklbl)
            
            save_file_path = fullfile(savepath,sprintf('%s_%d.png',plotname,i));
            %print (sprintf('%s_%d',plotname,i),'-dpng');
            print (save_file_path,'-dpng');
            
            
        end

        % check if we had zero discharges this chunk
        if isempty(time_of_discharge)
%{            
            discharge_times = [discharge_times time_this_chunk(end)];
            discharge_times_cutoff = [discharge_times_cutoff time_this_chunk(end)];
            
            discharge_vals = [discharge_vals 0.0];
            discharge_vals_cutoff = [discharge_vals_cutoff 0.0];
            
            discharge_stdevs = [discharge_stdevs 0.0];
            discharge_stdevs_cutoff = [discharge_stdevs_cutoff 0.0];

            discharges_per_chunk = [ discharges_per_chunk 0];
            discharges_per_chunk_cutoff = [ discharges_per_chunk_cutoff 0];
%}            

            discharge_times = [discharge_times [time_this_chunk(end) ; i]];
            discharge_times_cutoff = [discharge_times_cutoff [time_this_chunk(end) ; i]];
            
            discharge_vals = [discharge_vals [0.0 ; i]];
            discharge_vals_cutoff = [discharge_vals_cutoff [0.0 ; i]];
            
            discharge_stdevs = [discharge_stdevs [0.0 ; i]];
            discharge_stdevs_cutoff = [discharge_stdevs_cutoff [0.0 ; i]];

            discharges_per_chunk = [ discharges_per_chunk [0 ; i]];
            discharges_per_chunk_cutoff = [ discharges_per_chunk_cutoff [0 ; i]];            
        
        else   
            
            discharge_times = [discharge_times time_of_discharge];
            discharge_times_cutoff = [discharge_times_cutoff time_of_discharge_cutoff];

            discharge_vals = [discharge_vals discharge_values];
            discharge_vals_cutoff = [discharge_vals_cutoff discharge_values_cutoff];

            discharge_stdevs = [discharge_stdevs stdev_discharge_values];
            discharge_stdevs_cutoff = [discharge_stdevs_cutoff stdev_discharge_values_cutoff];

            discharges_per_chunk = [ discharges_per_chunk [count_discharges ; i]];
            discharges_per_chunk_cutoff = ...
                [ discharges_per_chunk_cutoff [count_discharges_cutoff ; i ]];
        
        end
        
        time_length_of_chunk(i) = time_this_chunk(end) - time_this_chunk(1);

    end 