function f = hv_plot_xy_errors(plot_title,annote,legend_show,...
    save_fig,savepath,plotname,...
    xdata_1,xdata_stdev_1,ydata_1,ydata_stdev_1,varargin)

% Properties of plot_x_y_errors:
%
% hv_plot_xy_errors(annote,legend_show,...
%    xdata_1,xdata_stdev_1,ydata_1,ydata_stdev_1,varargin)
%
% annote: set to "1" to generate an annotation box
%
% legend_show: set to "1" to generate a legend
%
% x_data_i, xdata_stdev_i,ydata_i,ydata_stdev_i: arrays of x and y data and
% their standard deviations. Enter up to three sets of data
%
% This plots up to three full error bar (x and y errors) datasets, assuming 
% the error is symmetric for x and y. Error bars for x (or y, for that 
% matter)can be ignored by setting x_data-stdev_i to an array of zeros.
%
%

    legend_string = '';
% Determine how many datasets to plot. 
    num_total_args = nargin;
    num_var_args = length(varargin);
    num_total_args_str = sprintf('%d',num_total_args);
    num_var_args_str = sprintf('%d',num_var_args);
    mes1 = ['User provided ',num_total_args_str,' total arguments.'];
    mes2 = ['User provided ',num_var_args_str,' optional arguments.'];
    disp('Starting hv_plot_xy_errors');
    disp(mes1);
    disp(mes2);

    overplot = 0;
    
    if num_var_args == 0
        disp('User provided one dataset to plot.');
    elseif num_var_args == 3
        disp('User provided two datasets to plot.');
        xdata_2 = varargin{1};
        xdata_stdev_2 = varargin{2};
        ydata_2 = varargin{3};
        ydata_stdev_2 = zeros(length(ydata_2));
        overplot = 1;
        
    elseif num_var_args == 4
        disp('User provided two datasets to plot.');
        xdata_2 = varargin{1};
        xdata_stdev_2 = varargin{2};
        ydata_2 = varargin{3};
        ydata_stdev_2 = varargin{4};
        overplot = 1;
        
    elseif num_var_args == 7
        disp('User provided three datasets to plot.');
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
        disp('User provided three datasets to plot.');
        xdata_2 = varargin{1};
        xdata_stdev_2 = varargin{2};
        ydata_2 = varargin{3};
        ydata_stdev_2 = varargin{4};
        xdata_3 = varargin{5};
        xdata_stdev_3 = varargin{6};
        ydata_3 = varargin{7};
        ydata_stdev_3 = varargin{8};
        overplot = 2;
    
    else
        disp('User provided incorrect number of optional arguments. Ignoring all optional arguments...');
    end

    annote_string = annote;
    xlabel_string = 'xlabel string';
    ylabel_string = 'ylabel string';
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

    num_colors = 3*length(redhsv(:,1));
    join_list_rgb = zeros(num_colors,3);

    for i =1:num_sets + 1
        join_list_rgb(i,:,:) = red_list_rgb(i,:,:);
        join_list_rgb(i+num_sets + 1,:,:) = grn_list_rgb(i,:,:);
        join_list_rgb(i+2*num_sets+2,:,:) = blue_list_rgb(i,:,:);

    end

    cmap = colormap(join_list_rgb);

    %titles and shit

    approximate_time_points = 4800;
    x_label = 'x label';
    y_label = 'y label';
    xtick_numbers = [ 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75];
    ytick_numbers = [-125 -100 -75 -50 -25 0 25 50 75 100 125];
    xmin = 0.0;
    xmax = +40.0;
    ymin = 0.5;
    ymax = +5.5;
    plot_bounds = [xmin xmax ymin ymax];

    ymin_right = -30;
    ymax_right = 30;
    plot_bounds_right = [xmin xmax ymin_right ymax_right];

    % y = 0 line
    zero_line = zeros(2,2);
    zero_line(1,1) = xmin;
    zero_line(2,1) = xmax;


%%%%%%%%%%%%%%%%%%%% psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This plots up to three full error bar (x and y errors) datasets, assuming 
% the error is symmetric for x and y. Error bars for x (or y, for that 
% matter)can be ignored by setting x_data-stdev_i to an array of zeros.

    figure1 = figure('Units','normalized');
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 6]; 
    
    errorbar(xdata_1,ydata_1,...
         ydata_stdev_1,ydata_stdev_1,xdata_stdev_1,xdata_stdev_1,...
         'o','Color', 'red','MarkerSize', 10, 'LineWidth', 2.0);
     if overplot > 0
         hold on;
        errorbar(xdata_2,ydata_2,...
            ydata_stdev_2,ydata_stdev_2,xdata_stdev_2,xdata_stdev_2,...
             'x','Color', 'blue','MarkerSize', 10, 'LineWidth', 2.0);
     end

     if overplot > 1
         hold on;
        errorbar(xdata_3,ydata_3,...
            ydata_stdev_3,ydata_stdev_3,xdata_stdev_3,xdata_stdev_3,...
             '^','Color', 'green','MarkerSize', 10, 'LineWidth', 2.0);
     end

    pbaspect([1.33 1 1])
    ax = gca; % current axes
    ax.TickDir = 'out'; % make ticks point out
    title(title_string,'FontSize',20)
    xlabel(xlabel_string,'FontSize',16)
    ylabel(ylabel_string,'FontSize',16)
    if legend_show == 1
              
        l = legend('show'); l.String = legend_string; 
        l.FontSize = 16; l.Location = 'northeast outside';
    
    end

    if (length(annote) > 0)
        annotation(figure1,'textbox',outside_plot,'String',annote_string,...
            'FontSize',16,'BackgroundColor',[1 1 1]);
    end
    
    if save_fig == 1
        
%         fig = gcf;
%         fig.PaperUnits = 'inches';
%         fig.PaperPosition = [0 0 12 9]; 
        
        save_file_path = fullfile(savepath,sprintf('%s.png',plotname));
        print (save_file_path,'-dpng');
    end
        
    
%     left_color = cmap(num_sets+1+0*num_sets+2,:);
%     right_color = cmap(num_sets+1+2*num_sets+2,:);
    % set(fig,'defaultAxesColorOrder',[left_color; right_color]);

%     figure2= figure('Units','normalized')
%     plot(1:1:35,1:1:35,...
%         '^','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0);
    

end