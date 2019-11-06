function f = hv_plot_xy_errors(plot_title,x_axis_title,y_axis_title,...
    input_style,annote,legend_show,...
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


    legend_string = '';
% Determine how many datasets to plot. 
    num_total_args = nargin;
    num_var_args = length(varargin);
    legend_names = {'plot 1'};

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
            legend_names(end+1) = {'plot 2'};

        elseif num_var_args == 4
    %        disp('User provided two datasets to plot.');
            num_series = 2;
            xdata_2 = varargin{1};
            xdata_stdev_2 = varargin{2};
            ydata_2 = varargin{3};
            ydata_stdev_2 = varargin{4};
            overplot = 1;
            legend_names(end+1) = {'plot 2'};

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
            legend_names(end+1) = {'plot 2'};
            legend_names(end+1) = {'plot 3'};

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
            legend_names(end+1) = {'plot 2'};
            legend_names(end+1) = {'plot 3'};

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
            legend_names(end+1) = {'plot 2'};
            legend_names(end+1) = {'plot 3'};
            legend_names(end+1) = {'plot 4'};

        end

        if num_var_args >= 16
            % 5 datasets
            num_series = 5;
            xdata_5 = varargin{13};
            xdata_stdev_5 = varargin{14};
            ydata_5 = varargin{15};
            ydata_stdev_5 = varargin{16};
            overplot = 4;
            legend_names(end+1) = {'plot 5'};

        end

        if num_var_args >= 20
            % 6 datasets
            num_series = 6;
            xdata_6 = varargin{17};
            xdata_stdev_6 = varargin{18};
            ydata_6 = varargin{19};
            ydata_stdev_6 = varargin{20};
            overplot = 5;
            legend_names(end+1) = {'plot 6'};

        end

        if num_var_args >= 24
            % 7 datasets
            num_series = 7;
            xdata_7 = varargin{21};
            xdata_stdev_7 = varargin{22};
            ydata_7 = varargin{23};
            ydata_stdev_7 = varargin{24};
            overplot = 6;
            legend_names(end+1) = {'plot 7'};

        end

        if num_var_args >= 28
            % 8 datasets
            num_series = 8;
            xdata_8 = varargin{25};
            xdata_stdev_8 = varargin{26};
            ydata_8 = varargin{27};
            ydata_stdev_8 = varargin{28};
            overplot = 7;
            legend_names(end+1) = {'plot 8'};

        end

        if num_var_args >= 32
            % 9 datasets
            num_series = 9;
            xdata_9 = varargin{29};
            xdata_stdev_9 = varargin{30};
            ydata_9 = varargin{31};
            ydata_stdev_9 = varargin{32};
            overplot = 8;
            legend_names(end+1) = {'plot 9'};

        end
    
    elseif input_style == 2
        
        this_set_indices = varargin{1};
        
        num_series = length(this_set_indices);
         
        for i = 2:num_series
            
            legend_names(end+1) = {sprintf('plot %d',i)};
             
        end

% Ti13 legend entries and color code        
        legend_names = {'14.9 kV','17.8 kV','19.4 kV','17.8','0.7 kV','9.8 kV',...
            '17.8 kV','17.8 kV','16.1 kV','17.9 kV','19.5 kV',...
            '20.5 kV','20.5 kV','21.9 kV','21.9 kV','21.9 kV',...
            '23.3 kV','24.8 kV','26.2 kV','26.2 kV','26.2 kV',...
            '+6.2 kV','27.6 kV','27.6 kV','29.1 kV','14.7 kV'};
        
        color_code_by_voltage = [4 6 8 6 1 2 6 6 5 7 9 10 10 11 11 11 12 13 14 14 14 14 15 15 16 3];

% Nb56 legend entries and color code

%         legend_names = {'12 kV','13 kV','14 kV','15 kV','16 kV',...
%             '17 kV','18 kV','19 kV','20 kV'};
        
%         color_code_by_voltage = [1 2 3 4 5 6 7 8 9];
        
        xrange = max(xdata_1) - min(xdata_1);
        yrange = max(ydata_1) - min(ydata_1);
        
        xmin = 0.0;
        xrange_padded = 5*ceil((xmin+xrange)*1.15/5);
        
        xmax = xmin + xrange_padded;
        
        ymin = -2;
        yrange_padded = 20*ceil((yrange)*1.2/20);
        
        ymax = yrange_padded;
        
            % marker pattern
        marker_pattern_one_sequence = ['o','x','^','s','>'];
        marker_pattern = marker_pattern_one_sequence;

        for j = 1:ceil(length(this_set_indices)/length(marker_pattern))

            marker_pattern = [marker_pattern marker_pattern_one_sequence];

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
    % yellows
    
%    num_colors = num_series;
 
    %for the Ti13 12/2018 -- 5/2019 datasets
      num_colors = 16;
    
    a_third = ceil(num_colors /3);
    
    middle_ish = floor((num_colors+0.5*a_third)/2);
    
    jet_set = jet(num_colors+a_third);
    
    jet_subset = zeros(num_colors,3);
    
    jet_subset(1:middle_ish,:) = jet_set(1:middle_ish,:);
    
    jet_subset(middle_ish+1:end,:) = jet_set(middle_ish+a_third+1:end,:);
        
    


    %titles and shit

    approximate_time_points = 4800;
    x_label = 'x label';
    y_label = 'y label';
    xtick_numbers = [ 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75];
    ytick_numbers = [-125 -100 -75 -50 -25 0 25 50 75 100 125];
    


%     ymin_right = -30;
%     ymax_right = 30;
%     plot_bounds_right = [xmin xmax ymin_right ymax_right];
% 
%     % y = 0 line
%     zero_line = zeros(2,2);
%     zero_line(1,1) = xmin;
%     zero_line(2,1) = xmax;


%%%%%%%%%%%%%%%%%%%% psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This plots up to three full error bar (x and y errors) datasets, assuming 
% the error is symmetric for x and y. Error bars for x (or y, for that 
% matter)can be ignored by setting x_data-stdev_i to an array of zeros.

    figure1 = figure('Units','normalized');
    
    fig = gcf;
    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 8 6];


%     fig.PaperPosition = [0 0 5.33 2.5]; % 2 x 2 grid
    fig.PaperPosition = [0 0 5.33 3.50];
    
    if length(bounds) == 4
        
        axis(bounds)
        
    end
    
%     markersize = 7.0;
    markersize = 5.0;
    linewidth = 1.125;

%%%    
    
    % this sequence is only for generating the legend icons without the
    % error bars.
    
    if input_style == 1

        if num_series <= 3

        plot(xdata_1(1),ydata_1(1),...
            'o','Color', 'red','MarkerSize', markersize, 'LineWidth', linewidth);

        else

            plot(xdata_1(1),ydata_1(1),...
            'o','Color', jet_subset(1,:),'MarkerSize', markersize, 'LineWidth', linewidth);

        end


         if overplot > 0
             hold on;

             if num_series <= 3

                 plot(xdata_2(1),ydata_2(1),...
                     'x','Color', 'blue','MarkerSize', markersize, 'LineWidth', linewidth);

             else             

                 plot(xdata_2(1),ydata_2(1),...
                     'x','Color', jet_subset(2,:),'MarkerSize', markersize, 'LineWidth', linewidth);

             end

         end

         if overplot > 1

             hold on;

             if num_series <= 3

                 plot(xdata_3(1),ydata_3(1),...
                     '^','Color', 'green','MarkerSize', markersize, 'LineWidth', linewidth);

             else

                 plot(xdata_3(1),ydata_3(1),...
                '^','Color', jet_subset(3,:),'MarkerSize', markersize, 'LineWidth', linewidth);

             end

         end

         if overplot > 2

             hold on;
            plot(xdata_4(1),ydata_4(1),...
                's','Color', jet_subset(4,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 3

             hold on;
            plot(xdata_5(1),ydata_5(1),...
                '>','Color', jet_subset(5,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 4

             hold on;
            plot(xdata_6(1),ydata_6(1),...
                'x','Color', jet_subset(6,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 5

             hold on;
            plot(xdata_7(1),ydata_7(1),...
                '^','Color', jet_subset(7,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 6

             hold on;
            plot(xdata_8(1),ydata_8(1),...
                's','Color', jet_subset(8,:),'MarkerSize', markersize, 'LineWidth', linewidth);

         end

         if overplot > 7

             hold on;
            plot(xdata_9(1),ydata_9(1),...
                '>','Color', jet_subset(9,:),'MarkerSize', markersize, 'LineWidth', linewidth);

         end
     
    elseif input_style == 2
        
        for i = 1:length(this_set_indices)
        
            hold on;
            plot(xdata_1(this_set_indices(i,1)),...
                ydata_1(this_set_indices(i,1)),...
                marker_pattern(i),'Color', jet_subset(color_code_by_voltage(i),:),...
                'MarkerSize',markersize,...
                'LineWidth', 1.0);
        end
            
    end
     
%%%

    

%      pbaspect([1.33 1 1])
    ax = gca; % current axes
    ax.FontSize = 8;
    ax.TickDir = 'out'; % make ticks point out
    %title(title_string,'FontSize',20)
    title(title_string,'FontSize',12)
%    xlabel(xlabel_string,'FontSize',16)
%    ylabel(ylabel_string,'FontSize',16)
    xlabel(xlabel_string,'FontSize',10)
    ylabel(ylabel_string,'FontSize',10)
    
    
    if legend_show == 1
                      
        l = legend('show'); 
%         l.String = legend_string; 
        l.String = legend_names; 
%        l.FontSize = 16;
%        l.FontSize = 12;
        l.FontSize = 7;
        l.Location = 'northeast outside';    

    end
    
    % this section actually plots the data
    
    if input_style == 1
        
        if num_series <= 3

            errorbar(xdata_1,ydata_1,...
                ydata_stdev_1,ydata_stdev_1,xdata_stdev_1,xdata_stdev_1,...
                'o','Color', 'red','MarkerSize', markersize, 'LineWidth', linewidth);

        else

            errorbar(xdata_1,ydata_1,...
                ydata_stdev_1,ydata_stdev_1,xdata_stdev_1,xdata_stdev_1,...
                'o','Color', jet_subset(1,:),'MarkerSize', markersize, 'LineWidth', linewidth);

        end


         if overplot > 0

             hold on;

             if num_series <= 3

                 errorbar(xdata_2,ydata_2,...
                     ydata_stdev_2,ydata_stdev_2,xdata_stdev_2,xdata_stdev_2,...
                     'x','Color','blue' ,'MarkerSize', markersize, 'LineWidth', linewidth);

             else

                 errorbar(xdata_2,ydata_2,...
                     ydata_stdev_2,ydata_stdev_2,xdata_stdev_2,xdata_stdev_2,...
                     'x','Color',jet_subset(2,:),'MarkerSize', markersize, 'LineWidth', linewidth);

             end

         end

         if overplot > 1

             hold on;

             if num_series <= 3

                 errorbar(xdata_3,ydata_3,...
                     ydata_stdev_3,ydata_stdev_3,xdata_stdev_3,xdata_stdev_3,...
                     '^','Color', 'green','MarkerSize', markersize, 'LineWidth', linewidth);

             else

                 errorbar(xdata_3,ydata_3,...
                     ydata_stdev_3,ydata_stdev_3,xdata_stdev_3,xdata_stdev_3,...
                     '^','Color', jet_subset(3,:),'MarkerSize', markersize, 'LineWidth', linewidth);

             end

         end

         if overplot > 2
             hold on;
            errorbar(xdata_4,ydata_4,...
                ydata_stdev_4,ydata_stdev_4,xdata_stdev_4,xdata_stdev_4,...
                 's','Color', jet_subset(4,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 3
             hold on;
            errorbar(xdata_5,ydata_5,...
                ydata_stdev_5,ydata_stdev_5,xdata_stdev_5,xdata_stdev_5,...
                 '>','Color', jet_subset(5,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 4
             hold on;
            errorbar(xdata_6,ydata_6,...
                ydata_stdev_6,ydata_stdev_6,xdata_stdev_6,xdata_stdev_6,...
                 'x','Color', jet_subset(6,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 5
             hold on;
            errorbar(xdata_7,ydata_7,...
                ydata_stdev_7,ydata_stdev_7,xdata_stdev_7,xdata_stdev_7,...
                 '^','Color', jet_subset(7,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 6
             hold on;
            errorbar(xdata_8,ydata_8,...
                ydata_stdev_8,ydata_stdev_8,xdata_stdev_8,xdata_stdev_8,...
                 's','Color', jet_subset(8,:),'MarkerSize', markersize, 'LineWidth', linewidth);
         end

         if overplot > 7
             hold on;
            errorbar(xdata_9,ydata_9,...
                ydata_stdev_9,ydata_stdev_9,xdata_stdev_9,xdata_stdev_9,...
                '>','Color', jet_subset(9,:),'MarkerSize', markersize, 'LineWidth', linewidth);

         end
     
    elseif input_style == 2
        
        % option for adding a solid black line to connect the data points.
%         hold on;
%         
%         plot(xdata_1(this_set_indices(1,1):this_set_indices(end,2)),...
%                 ydata_1(this_set_indices(1,1):this_set_indices(end,2)),...
%                 '-k','LineWidth', 0.5);
        
        for i = 1:length(this_set_indices)
        
            hold on;
            
            errorbar(xdata_1(this_set_indices(i,1):this_set_indices(i,2)),...
                ydata_1(this_set_indices(i,1):this_set_indices(i,2)),...
                ydata_stdev_1(this_set_indices(i,1):this_set_indices(i,2)),...
                ydata_stdev_1(this_set_indices(i,1):this_set_indices(i,2)),...
                xdata_stdev_1(this_set_indices(i,1):this_set_indices(i,2)),...
                xdata_stdev_1(this_set_indices(i,1):this_set_indices(i,2)),...
                marker_pattern(i),'Color', jet_subset(color_code_by_voltage(i),:),...
                'MarkerSize', markersize,...
                'LineWidth', linewidth);
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
            print (save_file_path,'-dmeta');
            
        end
        
    end
    
end