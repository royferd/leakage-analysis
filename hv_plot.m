function f = hv_plot(plot, overplot,annote,legend_show,...
    xdata_1,xdata_stdev_1,ydata_1,ydata_stdev_1,xdata_2,...
    xdata_stdev_2,ydata_2,ydata_stdev_2,xdata_3,xdata_stdev_3,ydata_3,...
    ydata_stdev_3,xdata_4,xdata_stdev_4,ydata_4,ydata_stdev_4)

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

    annote_string = 'annote_string';
    xlabel_string = 'xlabel_string';
    ylabel_string = 'ylabel_string';
    title_string = 'title_string';
    num_sets = length(xdata_1(:,1));
    
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


%%%%%%%%%%%%%%%%%%%% psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This plots up to three full error bar (x and y errors) datasets, assuming 
% the error is symmetric for x and y. Error bars for x (or y, for that 
% matter)can be ignored by setting x_data-stdev_i to an array of zeros.
    if plot == 1

        figure1 =figure('Units','normalized')
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
        title(title_string,'FontSize',40)
        xlabel(xlabel_string,'FontSize',32)
        ylabel(ylabel_string,'FontSize',32)
        if legend_show == 1
            l = legend('show'); l.String = legend_string; 
            l.FontSize = 32; l.Location = 'northeast outside';
        end
        
        if annote == 1
            annotation(figure1,'textbox',outside_plot,'String',annote_string,...
                'FontSize',32,'BackgroundColor',[1 1 1]);
        end
    end
    
    left_color = cmap(num_sets+1+0*num_sets+2,:);
    right_color = cmap(num_sets+1+2*num_sets+2,:);
    % set(fig,'defaultAxesColorOrder',[left_color; right_color]);

    for i = 1:num_sets
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHUNK TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% ramp data plot of mean up chunk leakage current v chunk # %%%%%%%

% This simply plots x vs. y. With 3 overplots selected, data_1 and data_2
% use the same symbol to represent one feature their sets share, while 
% data3 and data4 share a different symbol to represent their shared
% feature. data1 and data3 are the same color, as are data2 and data4, to
% indicate that they are from the same data set.

        if plot == 2

            figure2= figure('Units','normalized')
            plot(xdata_1(i,:),ydata_1(i,:),...
                '^','Color', cmap(1+i+ 2*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0);

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
        
%%%%%%%%% grid of field v. time with corresponding leakage current %%%%%%%%

% This makes a grid of plots of x vs. y with up to two datasets

        if plot == 3
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
    end

    %%%%%%% %grid of log leakage current & viktage v. time  %%%%%%%%%%%%%%%%%%%

% note: this won't run in versions previous to 2016a, use plotyy instead

% this plots a grid of plots with two vertical axes. Up to two datasets can
% be plotted on the left vertical axis. One dataset can be plotted on the right
% vertical axis. This was originally intended for plotting the log of
% positive (one dataset) and negative (the other dataset)
%  on the left vertical axis and a linear plot of the whole dataset on the
%  right axis.

        if plot == 4
             figure4= figure('Units','normalized')
            subplot(num_sets,1,num_sets + 1 - i)
            ax = gca; % current axes

            text(...
               'Position',[xmax ymax_right],...
               'String',sprintf('%d1.0',i),...
               'HorizontalAlignment','right','VerticalAlignment','top',...
               'FontSize',14);
            ax.TickDir = 'out'; % make ticks point out

            yyaxis left
            plot(xdata_1(i,:),ydata_1(i,:),...
                'x','Color', cmap(1+i+ 1*(num_sets+1),:),'MarkerSize', 4,...
                'LineWidth', 1.0); hold on;

            if overplot > 0
                plot(xdata_2(i,:),ydata_2(i,:),'o','Color',...
                    cmap(1+i+ 2*(num_sets+1),:),'MarkerSize', 4,'LineWidth', 1.0);
            end
            axis(plot_bounds)
            
            yyaxis right
            plot(xdata_3(i,:),ydata_3(i,:),...
               '-','Color', cmap(1+i+ 3*(num_sets+1),:),'MarkerSize', 3,...
               'LineWidth', 2.0);
            ylim([ymin_right ymax_right])
            pbaspect([7 1 1])
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%% Grid Super Axis Labels %%%%%%%%%%%%%%%%%%%%%%%%%%

    if plot == 3 || plot == 4
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
end