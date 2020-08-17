function f = hv_plot_xy_errors(plot_title,x_axis_title,y_axis_title,...
    input_style,annote,legend_show,legend_names,colorbar_values,...
    save_fig,savepath,plotname,bounds,...
    xdata_1,xdata_stdev_1,ydata_1,ydata_stdev_1,varargin)

%{

Properties of plot_x_y_errors:

hv_plot_xy_errors(annote,legend_show,...
   xdata_1,xdata_stdev_1,ydata_1,ydata_stdev_1,varargin)

annote: set to "1" to generate an annotation box

legend_show: set to "1" to generate a legend

x_data_i, xdata_stdev_i,ydata_i,ydata_stdev_i: arrays of x and y data and
their standard deviations. Enter up to three sets of data

This plots up to three full error bar (x and y errors) datasets, assuming 
the error is symmetric for x and y. Error bars for x (or y, for that 
matter)can be ignored by setting x_data-stdev_i to an array of zeros.

%}

    voltages = abs(colorbar_values);
    
    legend_string = '';
% Determine how many datasets to plot. 
    num_total_args = nargin;
    num_var_args = length(varargin);
    
%     generic_names = 0;
%     
%     if isempty(legend_names)
%         
%         generic_names = 1;
%         legend_names = {'plot 1'};
%         
%     end

%    disp('Starting hv_plot_xy_errors');

%     fprintf('User provided %d total arguments. \n',num_total_args);
%     fprintf('User provided %d optional arguments. \n',num_var_args);


    overplot = 0;
    
    num_series = 1;
    
    if input_style == 1
    
        if num_var_args == 0
            disp('User provided one dataset to plot.');
        elseif num_var_args == 3
    %        disp('User provided two datasets to plot.');
            num_series = 2;
            xdata_2 = varargin{1};
            xdata_stdev_2 = varargin{2};
            ydata_2 = varargin{3};
            ydata_stdev_2 = zeros(length(ydata_2));
            overplot = 1;

        elseif num_var_args == 4
    %        disp('User provided two datasets to plot.');
            num_series = 2;
            xdata_2 = varargin{1};
            xdata_stdev_2 = varargin{2};
            ydata_2 = varargin{3};
            ydata_stdev_2 = varargin{4};
            overplot = 1;

        elseif num_var_args == 7
    %        disp('User provided three datasets to plot.');
            num_series = 3;
            xdata_2 = varargin{1};
            xdata_stdev_2 = varargin{2};
            ydata_2 = varargin{3};
            ydata_stdev_2 = varargin{4};
            xdata_3 = varargin{5};
            xdata_stdev_3 = varargin{6};
            ydata_3 = varargin{7};
            ydata_stdev_3 = zeros(length(ydata_3));
            overplot = 2;

        elseif num_var_args == 8
    %        disp('User provided three datasets to plot.');
            num_series = 3;
            xdata_2 = varargin{1};
            xdata_stdev_2 = varargin{2};
            ydata_2 = varargin{3};
            ydata_stdev_2 = varargin{4};
            xdata_3 = varargin{5};
            xdata_stdev_3 = varargin{6};
            ydata_3 = varargin{7};
            ydata_stdev_3 = varargin{8};
            overplot = 2;

        end

        if num_var_args >= 12
            % 4 datasets
            num_series = 4;
            xdata_2 = varargin{1};
            xdata_stdev_2 = varargin{2};
            ydata_2 = varargin{3};
            ydata_stdev_2 = varargin{4};
            xdata_3 = varargin{5};
            xdata_stdev_3 = varargin{6};
            ydata_3 = varargin{7};
            ydata_stdev_3 = varargin{8};
            xdata_4 = varargin{9};
            xdata_stdev_4 = varargin{10};
            ydata_4 = varargin{11};
            ydata_stdev_4 = varargin{12};
            overplot = 3;

        end

        if num_var_args >= 16
            % 5 datasets
            num_series = 5;
            xdata_5 = varargin{13};
            xdata_stdev_5 = varargin{14};
            ydata_5 = varargin{15};
            ydata_stdev_5 = varargin{16};
            overplot = 4;

        end

        if num_var_args >= 20
            % 6 datasets
            num_series = 6;
            xdata_6 = varargin{17};
            xdata_stdev_6 = varargin{18};
            ydata_6 = varargin{19};
            ydata_stdev_6 = varargin{20};
            overplot = 5;

        end

        if num_var_args >= 24
            % 7 datasets
            num_series = 7;
            xdata_7 = varargin{21};
            xdata_stdev_7 = varargin{22};
            ydata_7 = varargin{23};
            ydata_stdev_7 = varargin{24};
            overplot = 6;

        end

        if num_var_args >= 28
            % 8 datasets
            num_series = 8;
            xdata_8 = varargin{25};
            xdata_stdev_8 = varargin{26};
            ydata_8 = varargin{27};
            ydata_stdev_8 = varargin{28};
            overplot = 7;

        end

        if num_var_args >= 32
            % 9 datasets
            num_series = 9;
            xdata_9 = varargin{29};
            xdata_stdev_9 = varargin{30};
            ydata_9 = varargin{31};
            ydata_stdev_9 = varargin{32};
            overplot = 8;

        end
    
    elseif input_style == 2
        
        this_set_indices = varargin{1};
        
        num_series = length(this_set_indices);
        

       

        if length(colorbar_values) > 0
            
            ordered_voltages = sort(voltages);
            ordered_color_code_by_voltage = [1];
            color_code_by_voltage = [];

            color_ht = 1;

            for j = 2:length(ordered_voltages)

                if ordered_voltages(j) > ordered_voltages(j-1)

                    color_ht = color_ht + 1;

                end

                ordered_color_code_by_voltage = [ordered_color_code_by_voltage color_ht];

            end


            for j = 1:length(voltages)

                this_color_ht_index = find(ordered_voltages == voltages(j),1);

                color_code_by_voltage = ...
                    [color_code_by_voltage ordered_color_code_by_voltage(this_color_ht_index)];

            end
        end
        

        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        xrange = max(xdata_1) - min(xdata_1);
        yrange = max(ydata_1) - min(ydata_1);
        
        xmin = 0.0;
        xrange_padded = 5*ceil((xmin+xrange)*1.15/5);
        
        xmax = xmin + xrange_padded;
        
        ymin = -2;
        yrange_padded = 20*ceil((yrange)*1.2/20);
        
        ymax = yrange_padded;
                
    end % end input_style == 1


                % marker pattern
    marker_pattern_one_sequence = ['o','x','^','s','>'];
    marker_pattern = marker_pattern_one_sequence;

    for j = 1:num_series

        marker_pattern = [marker_pattern marker_pattern_one_sequence];

    end
    
    if isempty(legend_names)
        
        legend_names = {'data 1'};
            
        for i = 2:num_series

            legend_names(end+1) = {sprintf('data %d',i)};

        end

    end
    
    
    annote_string = annote;
    
    xlabel_string = x_axis_title;
    
    ylabel_string = y_axis_title;
    
    title_string = plot_title;
    
    num_sets = 1;
    
    %text box assignments
    inside_plot = [0.4 0.73 0.31 0.19]; % edges: x y width height
    outside_plot = [0.73 0.45 0.21 0.21];

    grid_box = zeros(num_sets,4);

    for i = 1:num_sets
        grid_text(i,:) = [ 0.4 0.73 0.31 0.19 ];
    end

% generate colormap of red, blue, and green with intermediate intensities
% determined by the amount of datasets being plotted (i.e. more datasets
% will result in more intermediate intensities).
    
    redhsv = [0 1 1];     % red HSV
    red_list_hsv = zeros(num_sets + 1,3);
    red_list_rgb = zeros(num_sets + 1,3);
    red_step = round(0.75/(num_sets + 1),4);
    
    for i = 1:num_sets+1
        sat = 1 + red_step*(i - num_sets - 1);
        val = 1 + red_step*(i - num_sets - 1)/3.;
        red_list_hsv(i,1) = redhsv(1);
        red_list_hsv(i,2) = sat;
        red_list_hsv(i,3) = val;
        red_list_rgb(i,:,:) = hsv2rgb(red_list_hsv(i,:,:));
    end

    bluehsv = [0.667 1 1];     % blue HSV
    blue_list_hsv = zeros(num_sets + 1,3);
    blue_list_rgb = zeros(num_sets + 1,3);
    blue_step = round(0.75/(num_sets + 1),4);
    
    for i = 1:num_sets+1
        sat = 1 + blue_step*(i - num_sets - 1);
        val = 1 + blue_step*(i - num_sets - 1)/3.;
        blue_list_hsv(i,1) = bluehsv(1);
        blue_list_hsv(i,2) = sat;
        blue_list_hsv(i,3) = val;
        blue_list_rgb(i,:,:) = hsv2rgb(blue_list_hsv(i,:,:));
    end

    grnhsv = [0.333 1 1];     % green HSV
    grn_list_hsv = zeros(num_sets + 1,3);
    grn_list_rgb = zeros(num_sets + 1,3);
    grn_step = round(0.75/(num_sets + 1),4);
    
    for i = 1:num_sets+1
        sat = 1 + grn_step*(i - num_sets - 1);
        val = 1 + grn_step*(i - num_sets - 1)/3.;
        grn_list_hsv(i,1) = grnhsv(1);
        grn_list_hsv(i,2) = sat;
        grn_list_hsv(i,3) = val;
        grn_list_rgb(i,:,:) = hsv2rgb(grn_list_hsv(i,:,:));
    end

    num_colors_rgb = 3*length(redhsv(:,1));

    join_list_rgb = zeros(num_colors_rgb,3);
    
    

    for i =1:num_sets + 1
        join_list_rgb(i,:,:) = red_list_rgb(i,:,:);
        join_list_rgb(i+num_sets + 1,:,:) = grn_list_rgb(i,:,:);
        join_list_rgb(i+2*num_sets+2,:,:) = blue_list_rgb(i,:,:);

    end

    cmap = colormap(join_list_rgb);
    
    % remove middle third of jet colors to get rid of the lighter blues and
    % yellows, which don't show up well.
    
    if input_style == 1
        
        num_colors = num_series;
 
    elseif input_style == 2
        
        num_colors = max(color_code_by_voltage);
    end
    
    a_third = ceil(num_colors /3);
    
%     a_fourth = ceil(num_colors/4);
%     extra = 4 - mod(num_colors + num_colors,4)
%     a_third = num_colors + 1 - mod(num_colors,4);
    
    
    middle_ish = floor((num_colors+0.5*a_third)/2);
    
    jet_set = jet(num_colors+a_third);
    
    jet_subset = zeros(num_colors,3);
    
    jet_subset(1:middle_ish,:) = jet_set(1:middle_ish,:);
    
    jet_subset(middle_ish+1:end,:) = jet_set(middle_ish+a_third+1:end,:);
    
    red_blue_green_hsv = zeros(3,3);
    
    red_blue_green_hsv(1,:) = [0 1 1];
    red_blue_green_hsv(2,:) = [0.667 1 1];
    red_blue_green_hsv(3,:) = [0.333 1 1];
    red_blue_green_rgb = hsv2rgb(red_blue_green_hsv);
    
    if input_style == 1 && num_series <= 3
        
        color_palette = red_blue_green_rgb;
        
    elseif num_series == 4
        
        jet_sub = jet(5);
        
        color_palette = [[0. 0. 1.];jet_sub(1,:);jet_sub(5,:);[1. 0 0]];
        
        
        
    else
        
        color_palette = jet_subset;
%         color_palette = [jet(num_colors);
        
    end
    


    figure1 = figure('Units','normalized');
    
    fig = gcf;
    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 8 6];


%     fig.PaperPosition = [0 0 5.33 2.5]; % 2 x 2 grid
%     fig.PaperPosition = [0 0 5.33 3.50];
    fig.PaperPosition = [0 0 4.67 2.67];
    
    if length(bounds) == 4
        
        axis(bounds)
        
    end
    
    % this sequence is only for generating the legend icons without the
    % error bars.
    
    if input_style == 1
        
        markersize = 5.0;
        linewidth = 0.875;

        for i = 1:overplot+1
            
            hold on;
            
            var_xdata = eval(sprintf('xdata_%i',i));
            var_ydata = eval(sprintf('ydata_%i',i));
            
            plot(var_xdata(1),var_ydata(1),...
            marker_pattern(i),'Color', color_palette(i,:),...
            'MarkerSize', markersize, 'LineWidth', linewidth);
        
        end
     
    elseif input_style == 2
        
        markersize = 2.0;
%         linewidth = 1.125;
        linewidth = 0.5;
        
        for i = 1:length(this_set_indices)
        
            hold on;
%             plot(xdata_1(this_set_indices(i,1)),...
%                 ydata_1(this_set_indices(i,1)),...
%                 marker_pattern(i),'Color', color_palette(color_code_by_voltage(i),:),...
%                 'MarkerSize',markersize,...
%                 'LineWidth', linewidth);
            
            plot(xdata_1(this_set_indices(i,1)),...
                ydata_1(this_set_indices(i,1)),...
                '-','Color', color_palette(color_code_by_voltage(i),:),...
                'LineWidth', linewidth);
            
        end
            
    end
     
%%%

    

%      pbaspect([1.33 1 1])
    ax = gca; % current axes
    ax.FontSize = 9;
    ax.TickDir = 'out'; % make ticks point out
    %title(title_string,'FontSize',20)
    title(title_string,'FontSize',12)
%    xlabel(xlabel_string,'FontSize',16)
%    ylabel(ylabel_string,'FontSize',16)
    xlabel(xlabel_string,'FontSize',10)
    ylabel(ylabel_string,'FontSize',10)
    legend_font = 8;
    
    % get axis position, move up the y co-ordinate by just little bit and 
    % offset the height by the same amount
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)+0.05 pos(3) pos(4)-0.05]);
    set(gca,'fontWeight','bold');
    
    
    if legend_show == 1
        
        if input_style == 1
                      
            l = legend('show'); 

            l.String = legend_names; 

            l.FontSize = legend_font;
%             l.Location = 'northeast outside';  
            l.Location = 'northwest'; 
            
        elseif input_style == 2
        
            colormap(cmap);
            colormap(color_palette);
            c = colorbar;
            c.FontSize = legend_font;
            
            tick_locations = linspace(0.5/max(color_code_by_voltage),...
                (max(color_code_by_voltage)-0.5)/(max(color_code_by_voltage)),...
                max(color_code_by_voltage));
            
            c.Ticks = tick_locations;
            

            colorbar_ticklabel_string = cell(1,color_code_by_voltage(end));
            
            unique_voltages = [];
            
            for i =1:length(voltages)
                
                duplicates = (unique_voltages == voltages(i));
                
                if sum(duplicates) == 0
                
                    unique_voltages(end+1) = voltages(i);
                    
                end
                    
            end
            
            
            ordered_voltages = sort(unique_voltages);
            
            

            for i = 1:max(color_code_by_voltage)
                
                
                colorbar_ticklabel_string{i} = sprintf('%.1f kV',ordered_voltages(i));

            end
            
            c.TickLabels = colorbar_ticklabel_string;
            
        end

    end
    
    % this section actually plots the data
    
    if input_style == 1

        for i = 1:overplot+1
            
            hold on;
            
            var_xdata = eval(sprintf('xdata_%i',i));
            var_xdata_stdev = eval(sprintf('xdata_stdev_%i',i));
            var_ydata = eval(sprintf('ydata_%i',i));
            var_ydata_stdev = eval(sprintf('ydata_stdev_%i',i));
        
            [zero_error] = find(var_ydata_stdev == 0.0);

            [nonzero_error] = find(var_ydata_stdev ~= 0.0);

            if length(nonzero_error) > 0

                errorbar(var_xdata(nonzero_error),var_ydata(nonzero_error),...
                    var_ydata_stdev(nonzero_error),var_ydata_stdev(nonzero_error),...
                    var_xdata_stdev(nonzero_error),var_xdata_stdev(nonzero_error),...
                    marker_pattern(i),'Color', color_palette(i,:),...
                    'MarkerSize', markersize, 'LineWidth', linewidth);

            end

            if length(zero_error) > 0

                plot(var_xdata(zero_error),var_ydata(zero_error),...
                    marker_pattern(i),'Color', color_palette(i,:),...
                    'MarkerSize', markersize, 'LineWidth', linewidth);

            end
         
        end
     
    elseif input_style == 2
        
        % option for adding a solid black line to connect the data points.
%         hold on;
%         
%         plot(xdata_1(this_set_indices(1,1):this_set_indices(end,2)),...
%                 ydata_1(this_set_indices(1,1):this_set_indices(end,2)),...
%                 '-k','LineWidth', 0.5);
        
        last_set_xdata_temp = [];
        
        last_set_xdata_stdev_temp = [];
        
        last_set_ydata_temp = [];
        
        last_set_ydata_stdev_temp = [];

        for i = 1:length(this_set_indices)
            
%         for i = 1:1
        
            hold on;
            
            xdata_temp = transpose([last_set_xdata_temp xdata_1(this_set_indices(i,1):this_set_indices(i,2))]);
            
            xdata_stdev_temp = transpose([last_set_xdata_stdev_temp xdata_stdev_1(this_set_indices(i,1):this_set_indices(i,2))]);
            
            ydata_temp = transpose([last_set_ydata_temp transpose(ydata_1(this_set_indices(i,1):this_set_indices(i,2)))]);
            
            ydata_stdev_temp = transpose([last_set_ydata_stdev_temp transpose(ydata_stdev_1(this_set_indices(i,1):this_set_indices(i,2)))]);
                        
            [zero_error] = find(ydata_stdev_temp == 0.0);
           
            [nonzero_error] = find(ydata_stdev_temp ~= 0.0);
            
            % for the purpose of drawing the fill plot, we assign a small
            % width to uncertainties of 0.0. This avoids breaking the
            % program.
            if length(nonzero_error) > 0 && length(zero_error) > 0
                
                ydata_stdev_temp(zero_error) = 0.5*min(ydata_stdev_temp(nonzero_error));
                
            end
            
            
            % need this to connect the lines of multiple series on the
            % plot!
            last_set_xdata_temp = xdata_temp(end);
            
            last_set_xdata_stdev_temp = xdata_stdev_temp(end); 
            
            last_set_ydata_temp = ydata_temp(end);
            
            last_set_ydata_stdev_temp = ydata_stdev_temp(end);
            
            
            if length(nonzero_error) > 0 && length(zero_error) < length(ydata_stdev_temp)
            
%                 fill([xdata_temp(nonzero_error);flipud(xdata_temp(nonzero_error))],...
%                     [ydata_temp(nonzero_error)-ydata_stdev_temp(nonzero_error);flipud(ydata_temp(nonzero_error)+ydata_stdev_temp(nonzero_error))],...
%                     color_palette(color_code_by_voltage(i),:),'linestyle','-','EdgeColor', color_palette(color_code_by_voltage(i),:),...
%                     'LineWidth', linewidth);
                
                fill([xdata_temp;flipud(xdata_temp)],...
                    [ydata_temp-ydata_stdev_temp;flipud(ydata_temp+ydata_stdev_temp)],...
                    color_palette(color_code_by_voltage(i),:),'linestyle','-','EdgeColor', color_palette(color_code_by_voltage(i),:),...
                    'LineWidth', linewidth);
 
            elseif length(zero_error) == length(ydata_stdev_temp)
                
%                 fill([xdata_temp(nonzero_error);flipud(xdata_temp(nonzero_error))],...
%                     [ydata_temp(nonzero_error)-ydata_stdev_temp(nonzero_error);flipud(ydata_temp(nonzero_error)+ydata_stdev_temp(nonzero_error))],...
%                     color_palette(color_code_by_voltage(i),:),'linestyle','-','EdgeColor', color_palette(color_code_by_voltage(i),:),...
%                     'LineWidth', linewidth);
                
%                 marker_pattern(i),'Color', color_palette(color_code_by_voltage(i),:),...
                plot(xdata_temp(zero_error),...
                ydata_temp(zero_error),...
                'o-','Color', color_palette(color_code_by_voltage(i),:),...
                'MarkerSize', markersize,...
                'MarkerFaceColor',color_palette(color_code_by_voltage(i),:),...
                'LineWidth', 2*linewidth);
            
            end
            
            
            
        end
        
    end
    

    if (length(annote) > 0)
        annotation(figure1,'textbox',outside_plot,'String',annote_string,...
            'FontSize',10,'BackgroundColor',[1 1 1]);
    end
    
    if save_fig >= 1
        
%         fig = gcf;
%         fig.PaperUnits = 'inches';
%         fig.PaperPosition = [0 0 12 9]; 

%         title(title_string,'FontSize',16);
%         
%         xlabel(xlabel_string,'FontSize',14);
%         
%         ylabel(ylabel_string,'FontSize',14);
    
%         if legend_show == 1
%             
%             l.FontSize = 10;   
% 
%         end

%          axis([0 120 0 1000]);
        ax.Box = 'on';
        
        save_file_path = fullfile(savepath,sprintf('%s.png',plotname));
        print (save_file_path,'-dpng');
        
        if save_fig == 2
        
            save_file_path = fullfile(savepath,sprintf('%s.emf',plotname));
            print (save_file_path,'-dmeta','-r600');
            
            % PDF and EPS looks like shit in MATLAB.
%             save_file_path = fullfile(savepath,sprintf('%s.pdf',plotname));
%             print (fig,save_file_path,'-dpdf');
%             
%             save_file_path = fullfile(savepath,sprintf('%s.eps',plotname));
%             print (save_file_path);
            
        end
        
    end
    
end