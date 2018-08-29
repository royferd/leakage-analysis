function f = hv_plot(plot, overplot,annote,annote_string,legend_show,legend_string,xlabel_string,...
    ylabel_string,title_string,xdata_1,xdata_stdev_1,ydata_1,ydata_stdev_1,xdata_2,...
    xdata_stdev_2,ydata_2,ydata_stdev_2,xdata_3,xdata_stdev_3,ydata_3,...
    ydata_stdev_3,xdata_4,xdata_stdev_4,ydata_4,ydata_stdev_4,num_sets)

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
    legend_titles = [filenames];
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


%%%%%%%%%%%%%%%%%%% sample set plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot == 2

        figure2 =figure('Units','normalized')
        ramp_up_point = 58;
        errorbar(xdata_1,ydata_1,...
             0.5*ydata_stdev_1,... 
             0.5*ydata_stdev_1,...
             0.5*xdata_stdev_1,...
             0.5*xdata_stdev_1,...
             'o','Color', 'red','MarkerSize', 10, 'LineWidth', 2.0);
         if overplot > 0
             hold on;
            errorbar(xdata_2,ydata_2,...
                 0.5*ydata_stdev_2,...
                 0.5*ydata_stdev_2,...
                 0.5*xdata_stdev_2,...
                 0.5*xdata_stdev_2,...
                 'o','Color', 'blue','MarkerSize', 10, 'LineWidth', 2.0);
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
    end
    
%%%%%%%%%%%%%% plot leakage current v. meas. test # %%%%%%%%%%%%%%%%%%%%%%    
    if plot == 3

        figure3 = figure('Units','normalized')
        errorbar(xdata_1,ydata_1,ydata_stdev_1,...
            's','Color', 'blue','MarkerSize', 12, 'LineWidth', 1.0); hold on; %nA
        
        if overplot > 0
            errorbar(xdata_2,ydata_2,ydata_stdev_2,...
                'g^','MarkerSize', 8, 'LineWidth', 1.0); hold on; %nA
        end
        
        if overplot > 1
             errorbar(xdata_3,ydata_3,ydata_stdev_3,...
                'rx','MarkerSize', 8, 'LineWidth', 1.0); hold on; %nA
        end

        pbaspect([1.33 1 1])
        ax = gca; % current axes
        ax.TickDir = 'out'; % make ticks point out
        title(title_string,'FontSize',40)
        xlabel(xlabel_string,'FontSize',32)
        ylabel(ylabel_string,'FontSize',32)
        
        if annote == 1
            annotation(figure3,'textbox',outside_plot,'String',annote_string,...
                'FontSize',32,'BackgroundColor',[1 1 1]);
        end
    end
    
    left_color = cmap(num_sets+1+0*num_sets+2,:);
    right_color = cmap(num_sets+1+2*num_sets+2,:);
    % set(fig,'defaultAxesColorOrder',[left_color; right_color]);

    for i = 1:num_sets
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHUNK TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% ramp data plot of mean up chunk leakage current v chunk # %%%%%%%
        if plot == 8

            figure5= figure('Units','normalized')
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
               annotation(figure5,'textbox',outside_plot,'String',...
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

%%%%%%% %grid of log leakage current & viktage v. time  %%%%%%%%%%%%%%%%%%%

% note: this won't run in versions previous to 2016a, use plotyy instead
        if plot == 5
            subplot(num_sets,1,num_sets + 1 - i)
        %    figure
            ax = gca; % current axes

            text(...
               'Position',[xmax ymax_right],...
               'String',legend_titles(i),...
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
        
%%%%%%%%% grid of field v. time with corresponding leakage current %%%%%%%%
        if plot == 6
            if i > 1
                hold on;
            end
           subplot(2*num_sets,1,2*num_sets + 1 - i)
           plot(xdata_1(i,:),ydata_1(i,:),...
               '-','Color', cmap(1+i+ 0*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0); 
           pbaspect([4 1 1]);    ax = gca; ax.FontSize = 16; ax.TickDir = 'out';
           
           if overplot > 0
               hold on;
               subplot(2*num_sets,1,2*num_sets - i)
               plot(xdata_2(i,:),ydata_2(i,:) ,...
                   '-','Color', cmap(1+i+ 1*(num_sets+1),:),'MarkerSize', 8, 'LineWidth', 2.0);
           end
           pbaspect([4 1 1])
           ax = gca; % current axes
            text('Position',[xmax ymax],'String',legend_titles(i),...
                'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
           ax.FontSize = 16;
           ax.TickDir = 'out'; % make ticks point out
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%% Grid Super Axis Labels %%%%%%%%%%%%%%%%%%%%%%%%%%

    if plot == 5 || plot == 6
        subplot(num_sets,1,num_sets)
        xlabel(xlabel_string,'FontSize',32)
        subplot(num_sets,1,num_sets-1)
        ylabel(ylabel_string,'FontSize',32)
        subplot(num_sets,1,1)
        title(title_string,'FontSize',40)

        if legend_show == 1
            l = legend('show'); 
            l.String = legend_titles; 
            l.FontSize = 32; 
            l.Location = 'northeast outside';
        end

        ax = gca; % current axes
        ax.FontSize = 32;
        ax.TickDir = 'out'; % make ticks point out
    end
end