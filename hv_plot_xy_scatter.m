function f = hv_plot_xy_scatter(title_string,annote,legend_show,...
    xdata_1,ydata_1,varargin)

%Properties of hv_plot_xy_scatter:
%
%function f = hv_plot_xy_scatter(annote,legend_show,...
%    xdata_1,ydata_1,varargin)
% This simply plots x vs. y. With 3 overplots selected, data_1 and data_2
% use the same symbol to represent one feature their sets share, while 
% data3 and data4 share a different symbol to represent their shared
% feature. data1 and data3 are the same color, as are data2 and data4, to
% indicate that they are from the same data set.
    
    legend_string = '';
    
    % Determine how many datasets to plot. 
    num_total_args = nargin;
    num_var_args = length(varargin);
    num_total_args_str = sprintf('%d',num_total_args);
    num_var_args_str = sprintf('%d',num_var_args);
    mes1 = ['User provided ',num_total_args_str,' total arguments.'];
    mes2 = ['User provided ',num_var_args_str,' optional arguments.'];
    disp('Starting hv_plot_xy_scatter');
    disp(mes1);
    disp(mes2);
    
    overplot = 0;
    
    if num_var_args ==0
    
        disp('User provided one dataset to plot.');
        
    elseif num_var_args == 2
        
        overplot = 1;
        disp('User provided two datasets to plot.');
        xdata_2 = varargin{1};
        ydata_2 = varargin{2};
        
    elseif num_var_args == 4
        
        overplot = 2;
        disp('User provided three datasets to plot.');
        xdata_2 = varargin{1};
        ydata_2 = varargin{2};
        xdata_3 = varargin{3};
        ydata_3 = varargin{4};
        
    elseif num_var_args == 6
    
        overplot = 3;
        disp('User provided four datasets to plot.');
        xdata_2 = varargin{1};
        ydata_2 = varargin{2};
        xdata_3 = varargin{3};
        ydata_3 = varargin{4};
        xdata_4 = varargin{5};
        ydata_4 = varargin{6};
        
    else 
        
        disp('You hath fucked up royally. Check your arguments and the documentation.')
        
    end
    
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
    
    annote_string = 'annote string';
    
    xlabel_string = 'xlabel string';
    
    ylabel_string = 'ylabel string';

    
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
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHUNK TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% ramp data plot of mean up chunk leakage current v chunk # %%%%%%%

% This simply plots x vs. y. With 3 overplots selected, data_1 and data_2
% use the same symbol to represent one feature their sets share, while 
% data3 and data4 share a different symbol to represent their shared
% feature. data1 and data3 are the same color, as are data2 and data4, to
% indicate that they are from the same data set.



        figure2= figure('Units','normalized')
        plot(xdata_1(i,:),ydata_1(i,:),...
            '^','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0);
    %cmap(1+i+ 2*(num_sets+1),:)
        if overplot > 0
            hold on;
            plot(xdata_2(i,:),ydata_2(i,:),...
                '^','Color', cmap(1+i+ 0*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0);
        end

        if overplot > 1
            hold on;
            plot(xdata_3(i,:),ydata_3(i,:),...
                'v','Color', cmap(1+i+ 2*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0);
        end

        if overplot > 2
            hold on;
            plot(xdata_4(i,:),ydata_4(i,:),...
                'v','Color', cmap(1+i+ 0*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0);
        end

       if legend_show == 1
           l = legend('show'); l.String = legend_string; 
           l.FontSize = 32; l.Location = 'northeast outside';
       end

       if annote == 1
           annotation(figure2,'textbox',outside_plot,'String',...
               annote_string,'FontSize',32,'BackgroundColor',[1 1 1]);
       end

       pbaspect([1.33 1 1])           
       ax = gca; % current axes
       ax.TickDir = 'out'; % make ticks point out
       ax.FontSize = 32;
       title(title_string,'FontSize',40)
       xlabel(xlabel_string,'FontSize',32)
       ylabel(ylabel_string,'FontSize',32)


        

    end
    
end