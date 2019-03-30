function f = hv_grid_one_axis(annote,legend_show,...
    xdata_1,ydata_1,varargin)

%Properties of plot:
%
% plot = 1: 
% This plots up to three full error bar (x and y errors) datasets, assuming 
% the error is symmetric for x and y. Error bars for x (or y, for that 
% matter)can be ignored by setting x_data-stdev_i to an array of zeros.
%
% plot = 2: 

% This simply plots x vs. y. With 3 overplots selected, data_1 and data_2
% use the same symbol to represent one feature their sets share, while 
% data3 and data4 share a different symbol to represent their shared
% feature. data1 and data3 are the same color, as are data2 and data4, to
% indicate that they are from the same data set.
%
% plot = 3:
% This makes a grid of plots of x vs. y with up to two datasets
%
% plot = 4:
%
% note: this won't run in versions previous to 2016a, use plotyy instead
%
% this plots a grid of plots with two vertical axes. Up to two datasets can
% be plotted on the left vertical axis. One dataset can be plotted on the right
% vertical axis. This was originally intended for plotting the log of
% positive (one dataset) and negative (the other dataset)
%  on the left vertical axis and a linear plot of the whole dataset on the
%  right axis.

    num_total_args = nargin;
    num_var_args = length(varargin);
    num_total_args_str = sprintf('%d',num_total_args);
    num_var_args_str = sprintf('%d',num_var_args);
    mes1 = ['User provided ',num_total_args_str,' total arguments.'];
    mes2 = ['User provided ',num_var_args_str,' optional arguments.'];
    disp('Starting hv_grid_one_axis');
    disp(mes1);
    disp(mes2);

    overplot = 0;
    
    if num_var_args == 0
        
      overplot = 0;
      disp('User provided one dataset to plot.');
        
    elseif num_var_args == 2
        
        overplot = 1;
        disp('User provided two datasets to plot.');
        xdata_2 = varargin{1};
        ydata_2 = varargin{2};
        
    else
        
        disp('Seriously you fucked up the input arguments. Check the documentation dumby.');
        
    end
    
    
    legend_string = '';
    annote_string = 'annote string';
    xlabel_string = 'xlabel string';
    ylabel_string = 'ylabel string';
    title_string = 'title string';
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
    title_string = 'plot title';
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
    
    for i = 1:num_sets
        
%%%%%%%%% grid of field v. time with corresponding leakage current %%%%%%%%

% This makes a grid of plots of x vs. y with up to two datasets

        figure3= figure('Units','normalized')

        if i > 1
            hold on;
        end

       subplot(2*num_sets,1,2*num_sets + 1 - i)
       plot(xdata_1(i,:),ydata_1(i,:),...
           '-','Color', cmap(1+i+ 0*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0); 
       pbaspect([4 1 1]); ax = gca; ax.FontSize = 16; ax.TickDir = 'out';

       if overplot > 0
           hold on;
           subplot(2*num_sets,1,2*num_sets - i)
           plot(xdata_2(i,:),ydata_2(i,:) ,...
               '-','Color', cmap(1+i+ 1*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0);
       end
       pbaspect([4 1 1])
       ax = gca; % current axes
        text('Position',[xmax ymax],'String',sprintf('%d1.0',i),...
            'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
       ax.FontSize = 16;
       ax.TickDir = 'out'; % make ticks point out

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%% Grid Super Axis Labels %%%%%%%%%%%%%%%%%%%%%%%%%%


    subplot(num_sets,1,num_sets)
    xlabel(xlabel_string,'FontSize',32)
    subplot(num_sets,1,num_sets-1)
    ylabel(ylabel_string,'FontSize',32)
    subplot(num_sets,1,1)
    title(title_string,'FontSize',40)

    if legend_show == 1
        l = legend('show'); 
        l.FontSize = 32; 
        l.Location = 'northeast outside';
    end

    ax = gca; % current axes
    ax.FontSize = 32;
    ax.TickDir = 'out'; % make ticks point out

end