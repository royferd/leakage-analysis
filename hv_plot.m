function f = hv_plot(xdata_1,xdata_stdev_1,ydata_1,ydata_stdev_1,xdata_2,...
    xdata_stdev_2,ydata_2,ydata_stdev_2,xdata_3,xdata_stdev_3,ydata_3,...
    ydata_stdev_3,xdata_4,xdata_stdev_4,ydata_4,ydata_stdev_4,num_sets)

    %start_point(5) = 50; %lost 30 min o fdata for 2-10-2017
    %time(5,:) = time(5,:) + 33.0; % lost 30 min of data for 2-10-2017

    for i = 1:num_sets
        str1 = 'average LO leakage current (pA): ';
        %str2 = ' I_{leak} = a + bV_{PS}';
        str2 = sprintf('%.0f',lcm1_avg_ramp_up_avg(i) - lcm1_avg_offset(i));
        str3 = ' \pm ';
        str4 = sprintf('%.0f',lcm1_avg_ramp_up_stdev(i));
    %    str5 = 'average stdev HI leakage current (pA): ';
    %    str6 = sprintf('%.0f',lcm1_avg_ramp_up_stdev_avg_chunk(i));
    %    str7 = ' \pm ';
    %    str8 = sprintf('%.0f',lcm1_avg_ramp_up_stdev_stdev_chunk);
        %str3 = 'a: ';
        %str4 = sprintf('%.0f', a);
        %str5 = ' \pm ';
        %str6 = sprintf('%.0f', a_stdev);
        %str7 = ' pA';
        %str8 = 'b: ';
        %str9 = sprintf('%.1f', b);
        %str10 = ' \pm ';
        %str11 = sprintf('%.1f', b_stdev);
        %str12 = ' pA/kV';
        %str7 = 'resistance: ';
        %str8 = num2str(round(resistance,2));
        %str9 = ' +/- ';
        %str10 = num2str(round(resistance_stdev,1));
        %str11 = ' P\Omega';
        %str = {[ str1], [str2 str3 str4 str5 str6], [str7 str8 str9 str10 str11]};
        str = {[ str1], [str2 str3 str4]};
    end
    
        % str1 = 'y_{fit}(pA) = a + bI_{src}';
    % str2 = 'leakage offset (pA): ';
    % str3 = sprintf('%.0f',a(1,1)*1e3);
    % str4 = ' \pm ';
    % str5 = sprintf('%.0f',a(1,2)*1e3);
    % str6 = 'slope: ';
    % str7 = sprintf('%.3f',b(1,1));
    % str8 = ' \pm ';
    % str9 = sprintf('%.3f',b(1,2));
    % str = {[ str1], [str2 str3 str4 str5],[str6 str7 str8 str9]};

    %text box assignments
    inside_plot = [0.4 0.73 0.31 0.19]; % edges: x y width height
    outside_plot = [0.73 0.45 0.21 0.21];

    grid_box = zeros(num_files,4);

    for i = 1:num_files
        grid_text(i,:) = [ 0.4 0.73 0.31 0.19 ];
    end

    
    redhsv = [0 1 1];     % red HSV
    red_list_hsv = zeros(num_files + 1,3);
    red_list_rgb = zeros(num_files + 1,3);
    red_step = round(0.75/(num_files + 1),4);
    for i = 1:num_files+1
        sat = 1 + red_step*(i - num_files - 1);
        val = 1 + red_step*(i - num_files - 1)/3.;
        red_list_hsv(i,1) = redhsv(1);
        red_list_hsv(i,2) = sat;
        red_list_hsv(i,3) = val;
        red_list_rgb(i,:,:) = hsv2rgb(red_list_hsv(i,:,:));
    end


    bluehsv = [0.667 1 1];     % blue HSV
    blue_list_hsv = zeros(num_files + 1,3);
    blue_list_rgb = zeros(num_files + 1,3);
    blue_step = round(0.75/(num_files + 1),4);
    for i = 1:num_files+1
        sat = 1 + blue_step*(i - num_files - 1);
        val = 1 + blue_step*(i - num_files - 1)/3.;
        blue_list_hsv(i,1) = bluehsv(1);
        blue_list_hsv(i,2) = sat;
        blue_list_hsv(i,3) = val;
        blue_list_rgb(i,:,:) = hsv2rgb(blue_list_hsv(i,:,:));
    end

    grnhsv = [0.333 1 1];     % green HSV
    grn_list_hsv = zeros(num_files + 1,3);
    grn_list_rgb = zeros(num_files + 1,3);
    grn_step = round(0.75/(num_files + 1),4);
    for i = 1:num_files+1
        sat = 1 + grn_step*(i - num_files - 1);
        val = 1 + grn_step*(i - num_files - 1)/3.;
        grn_list_hsv(i,1) = grnhsv(1);
        grn_list_hsv(i,2) = sat;
        grn_list_hsv(i,3) = val;
        grn_list_rgb(i,:,:) = hsv2rgb(grn_list_hsv(i,:,:));
    end

    num_colors = 3*length(redhsv(:,1));
    join_list_rgb = zeros(num_colors,3);

    for i =1:num_files + 1
        join_list_rgb(i,:,:) = red_list_rgb(i,:,:);
        join_list_rgb(i+num_files + 1,:,:) = grn_list_rgb(i,:,:);
        join_list_rgb(i+2*num_files+2,:,:) = blue_list_rgb(i,:,:);

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

    if plot == 1
    %%%%%%%%%%%%%%%%%%%%%%% of max field v. gap size %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure1 = figure('Units','normalized')
    plot(xdata_1, ydata_1,...
       'o','Color', cmap(num_sets+1,:),'MarkerSize', 6, 'LineWidth', 2.0); %nA
    pbaspect([1.33 1 1])
    ax = gca; % current axes
    ax.TickDir = 'out'; % make ticks point out
    title(title_string,'FontSize',40)
    xlabel(x_label,'FontSize',32)
    ylabel(y_label,'FontSize',32)
    end

    %%%%%%%%%%%%%%%%%%% sample set plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plot == 2
    %%%%%%%%%%%%%%%%%%%% psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure2 =figure('Units','normalized')
        ramp_up_point = 58;
        errorbar(xdata_1,ydata_1,...
             0.5*ydata_stdev_1,... 
             0.5*ydata_stdev_1,...
             0.5*xdata_stdev_1,...
             0.5*xdata_stdev_1,...
             'o','Color', 'red','MarkerSize', 10, 'LineWidth', 2.0); hold on;
        errorbar(xdata_2,ydata_2,...
             0.5*ydata_stdev_2,...
             0.5*ydata_stdev_2,...
             0.5*xdata_stdev_2,...
             0.5*xdata_stdev_2,...
             'o','Color', 'blue','MarkerSize', 10, 'LineWidth', 2.0); hold on; 
        %axis(plot_bounds)
        pbaspect([1.33 1 1])
        ax = gca; % current axes
        ax.TickDir = 'out'; % make ticks point out
        title(filenames,'FontSize',40)
        xlabel('imon','FontSize',32)
        ylabel('vmon','FontSize',32)
        l = legend('show'); l.String = [{'ramp up'}, {'ramp down'}]; l.FontSize = 32; l.Location = 'northeast outside';
    end
    
    if plot == 3
    %%%%%%%% plot leakage current v. meas. test # %%%%%%%%%%%%%%%%%%%%%%
        figure1 = figure('Units','normalized')
        errorbar(xdata_1,ydata_1,ydata_stdev_1,...
            's','Color', 'blue','MarkerSize', 12, 'LineWidth', 1.0); hold on; %nA
        if overplot == 1
            errorbar(xdata_2,ydata_2,ydata_stdev_2,...
                'g^','MarkerSize', 8, 'LineWidth', 1.0); hold on; %nA
        end
        if overplot > 1
             errorbar(xdata_3,ydata_3,ydata_stdev_3,...
                'rx','MarkerSize', 8, 'LineWidth', 1.0); hold on; %nA
        end

        %     %axis(plot_bounds)
         pbaspect([1.33 1 1])
         ax = gca; % current axes
         ax.TickDir = 'out'; % make ticks point out
         title(title_string,'FontSize',40)
        xlabel(x_label,'FontSize',32)
        ylabel(y_label,'FontSize',32)
         annotation(figure1,'textbox',...
           outside_plot,'String',str,'FontSize',32,'BackgroundColor',[1 1 1]);
    end
    
    left_color = cmap(num_sets+1,:);
    right_color = cmap(num_sets+1+2*num_files+2,:);
    % set(fig,'defaultAxesColorOrder',[left_color; right_color]);

    for i = 1:num_sets


        if plot == 4
         %%%%%%% grid of psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%
            figure1 =figure('Units','normalized')
            subplot(num_sets,1,num_sets + 1 - i)
            plot(xdata_1(i,:),ydata_1(i,:),...
                'o','Color', cmap(i+1,:),'MarkerSize', 2, 'LineWidth', 2.0); %nA
            %axis(plot_bounds)
            pbaspect([1.33 1 1])
            ax = gca; % current axes
            ax.TickDir = 'out'; % make ticks point out
        end
        
%%%%%% %grid of log leakage current & viktage v. time  %%%%%%%%%%%%%%%%%%%
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
                'x','Color', cmap(i+1+num_sets+1,:),'MarkerSize', 4,...
                'LineWidth', 1.0); hold on;
            if overplot > 0
                plot(xdata_2(i,:) - xdata_offset_2(i),...
                    ydata_2(i,:) - ydata_offset_2(i),...
                    'o','Color', cmap(i+1,:),'MarkerSize', 4,...
                    'LineWidth', 1.0);
            end
            axis(plot_bounds)
            
            yyaxis right
            plot(xdata_3(i,:),ydata_3(i,:),...
               '-','Color', cmap(i+1+2*num_sets+2,:),'MarkerSize', 3,...
               'LineWidth', 2.0);
            ylim([ymin_right ymax_right])
            pbaspect([7 1 1])
        end
    
    %%%%%%%%%%%%%%%%% grid of power ps voltage v. time %%%%%%%%%%%%%%%%%%%%%%%%
        if plot == 6
        %    subplot(num_sets,1,num_sets + 1 - i)
            plot(xdata_1(i,:),ydata_1(i,:),...
                '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
            pbaspect([4 1 1])
            ax = gca; % current axes
            text('Position',[xmax ymax],'String',legend_titles(i),...
                'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
            ax.FontSize = 16;
            ax.TickDir = 'out'; % make ticks point out
        end

    % %%%%%%%%% grid of field v. time with corresponding leakage current %%%%%%%%
    %    subplot(2*num_sets,1,2*num_sets + 1 - i)
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0); 
    %    pbaspect([4 1 1]);    ax = gca; ax.FontSize = 16; ax.TickDir = 'out'; hold on;
    %    subplot(2*num_sets,1,2*num_sets - i)
    %    plot(xdata_2(i,:) - xdata_offset_2(i),...
    %        ydata_2(i,:) - ydata_offset_2(i),...
    %        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    % %    axis(plot_bounds)
    %    pbaspect([4 1 1])
    %    ax = gca; % current axes
    %     text('Position',[xmax ymax],'String',legend_titles(i),...
    %         'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
    %    ax.FontSize = 16;
    %    ax.TickDir = 'out'; % make ticks point out

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHUNK TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % %%%%%%%%%%%%%% ramp data plot of ramp_up voltage v time %%%%%%%%%%%%%%%%%%%
    % %   subplot(num_sets,1,num_sets + 1 - i)
    %    %note: for 2017-10-23-175202-hv-1.txt, time_ramp_up and vmon_avg_ramp_up
    %    %have different lengths... possible bugs
    % 
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %    'o','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
    %    %axis(plot_bounds)
    %    pbaspect([1.33 1 1])
    %    ax = gca; % current axes
    % %    ax.FontSize = 16;
    %    ax.TickDir = 'out'; % make ticks point out

       %%%%%%%%%%%%%% ramp data plot of ramp voltage v time %%%%%%%%%%%%%%%%%%
    %    figure
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %        'x','Color', 'black','MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %    plot(xdata_2(i,:),ydata_2(i,:),...
    %        'o','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %    plot(xdata_3(i,:),ydata_3(i,:),...
    %        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %    plot(xdata_4(i,:),ydata_4(i,:),...
    %        'o','Color', cmap(1+i+2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0);
    %    %axis(plot_bounds)
    %    pbaspect([1.33 1 1])
    %    ax = gca; % current axes
    % %    ax.FontSize = 16;
    %    ax.TickDir = 'out'; % make ticks point out
    %    
    %       %%%%%%%%%%%%%% ramp data plot of ramp leakage v time %%%%%%%%%%%%%%%%%%
    %    figure
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %        'x','Color', 'black','MarkerSize', 8, 'LineWidth', 2.0); hold on;
    % %    plot(xdata_1(i,:),ydata_1(i,:),...
    % %        's','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    % 
    % % plot(xdata_2(i,:),ydata_2(i,:),...
    % %        '^','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %    plot(xdata_2(i,:),ydata_2(i,:),...
    %        'o','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    % 
    % %   plot(xdata_3(i,:),ydata_3(i,:),...
    % %        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    % %    plot(xdata_3(i,:), ydata_3(i,:),...
    % %        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %    plot(xdata_3(i,:),ydata_3(i,:),...
    %        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    % %    plot(xdata_4(i,:),ydata_4(i,:),...
    % %        '^','Color', cmap(1+i+2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    % 
    % %    plot(xdata_4(i,:),ydata_4(i,:),...
    % %        'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;   
    % %    plot(xdata_4(i,:),ydata_4(i,:),...
    % %        'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;    
    %    %axis(plot_bounds)
    %    pbaspect([1.33 1 1])
    %    ax = gca; % current axes
    % %    ax.FontSize = 16;
    %    ax.TickDir = 'out'; % make ticks point out

    %%%%%%%%%%%%%%%%% discharging summed leakage data v time %%%%%%%%%%%%%%%%%%
       figure
       plot(xdata_1(i,:),ydata_1(i,:),...
           'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
       plot(xdata_2(i,:),ydata_2(i,:),...
           'x','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
       l = legend('show'); l.String = [{'avg + discharge'},{'avg - discharge'}]; 
       l.FontSize = 32; l.Location = 'northeast outside';
       pbaspect([1.33 1 1])
       ax = gca; % current axes
       ax.FontSize = 32;
       ax.TickDir = 'out'; % make ticks point out
       title('Ramping current','FontSize',40)
       xlabel('time (s)','FontSize',32)
       ylabel('leakage current (pA)','FontSize',32)
       %axis(plot_bounds)
       pbaspect([1.33 1 1])
       ax = gca; % current axes
    %    ax.FontSize = 16;
       ax.TickDir = 'out'; % make ticks point out

    %%%%%%%%%%%%%%%% charging summed leakage data v time %%%%%%%%%%%%%%%%%%%%%%
       figure
       plot(xdata_1(i,:),ydata_1(i,:),...
           'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
       plot(xdata_2(i,:),ydata_2(i,:),...
           'x','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
       l = legend('show'); l.String = [{'avg + charge'},{'avg - charge'}]; 
       l.FontSize = 32; l.Location = 'northeast outside';
       pbaspect([1.33 1 1])
       ax = gca; % current axes
       ax.FontSize = 32;
       ax.TickDir = 'out'; % make ticks point out
       title('Ramping current','FontSize',40)
       xlabel('time (s)','FontSize',32)
       ylabel('leakage current (pA)','FontSize',32)
       %axis(plot_bounds)
       pbaspect([1.33 1 1])
       ax = gca; % current axes
    %    ax.FontSize = 16;
       ax.TickDir = 'out'; % make ticks point out  

    % % %%%%%%%%%% ramp data plot of mean up chunk voltage v chunk # %%%%%%%
    %     figure1 = figure('Units','normalized')
    %     plot(xdata_1(i,:), ydata_1(i,:),...
    %     'x','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0);
    %     %errorbar(xdata_1(i,:), ydata_1(i,:),ydata_stdev_1(i,:),'o');
    % %    axis(plot_bounds)
    %     pbaspect([1.33 1 1])
    %     ax = gca; % current axes
    %     ax.FontSize = 32;
    %     ax.TickDir = 'out'; % make ticks point out
    %     title('hvps HI ramp means','FontSize',40)
    %     xlabel('chunk #','FontSize',32)
    %     ylabel('high voltage (-kV)','FontSize',32)
    %     annotation(figure1,'textbox',...
    %     outside_plot,'String',{['avg HI (-kV):'],...
    %     [sprintf('%.3f',vmon_avg_ramp_up_avg(i) - vmon_avg_offset(i)) ' \pm ' sprintf('%.3f',vmon_avg_ramp_up_stdev(i))]},...
    %     'FontSize',32,'BackgroundColor',[1 1 1]);
    % 
    % %%%%%%%%%%% ramp data plot of mean down chunk voltage v chunk # %%%%%%%%%%
    %     figure2 = figure('Units','normalized')
    %     plot(xdata_1(i,:),ydata_1(i,:),...
    %     'o','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0);
    %     %errorbar(xdata_1(i,:),ydata_1(i,:),ydata_stdev_1(i,:),'o');
    % %    axis(plot_bounds)
    %     pbaspect([1.33 1 1])
    %     ax = gca; % current axes
    %     ax.FontSize = 32;
    %     ax.TickDir = 'out'; % make ticks point out
    %     title('hvps LO ramp means','FontSize',40)
    %     xlabel('chunk #','FontSize',32)
    %     ylabel('high voltage (-kV)','FontSize',32)
    %     annotation(figure2,'textbox',...
    %     outside_plot,'String',{['avg LO (-kV):'],...
    %     [sprintf('%.3f',vmon_avg_ramp_down_avg(i) - vmon_avg_offset(i)) ' \pm ' sprintf('%.3f',vmon_avg_ramp_down_stdev(i))]},...
    %     'FontSize',32,'BackgroundColor',[1 1 1]);
    % 
    %%%%%%%%% ramp data plot of stdev hi/lo chunk ps voltage v chunk # %%%%%%%
    %     figure3 = figure('Units','normalized')
    %     plot(xdata_1,(i,:),ydata_1(i,:),...
    %     'x','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %     plot(xdata_2(i,:),ydata_2(i,:),...
    %     'o','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0);
    %    %errorbar(xdata_1(i,:),ydata_1(i,:),ydata_stdev_1(i,:),'o');
    %     l = legend('show'); l.String = [{'HI'},{'LO'}]; l.FontSize = 32; l.Location = 'northeast outside';
    % %    axis(plot_bounds)
    %     pbaspect([1.33 1 1])
    %     ax = gca; % current axes
    %     ax.FontSize = 32;
    %     ax.TickDir = 'out'; % make ticks point out
    %     title('hvps ramp stdevs','FontSize',40)
    %     xlabel('chunk #','FontSize',32)
    %     ylabel('high voltage (-kV)','FontSize',32)
    %     annotation(figure3,'textbox',...
    %     outside_plot,'String',{['avg HI (kV):'],...
    %     [sprintf('%.3f',vmon_avg_ramp_up_stdev_avg_chunk(i)) ' \pm ' sprintf('%.3f',vmon_avg_ramp_up_stdev_stdev_chunk(i))],...
    %     ['avg LO (kV):'],...
    %     [sprintf('%.3f',vmon_avg_ramp_down_stdev_avg_chunk(i)) ' \pm ' sprintf('%.3f',vmon_avg_ramp_up_stdev_stdev_chunk(i))]},...
    %     'FontSize',32,'BackgroundColor',[1 1 1]);
    % 
    %%%%%%%%%% ramp data plot of mean down chunk leakage current v chunk # %%%%%%%
    %     figure4= figure('Units','normalized')
    %     plot(xdata_1(i,:),ydata_1(i,:),...
    %     'x','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0); 
    %     if inclusive_data == 1
    %         hold on;
    %         plot(xdata_2(i,:),ydata_2(i,:),...
    %             'o','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0);
    %     end
    %     pbaspect([1.33 1 1])
    %     ax = gca; % current axes
    %     ax.FontSize = 32;
    %     ax.TickDir = 'out'; % make ticks point out
    %     title('lcm1 LO ramp means','FontSize',40)
    %     xlabel('chunk #','FontSize',32)
    %     ylabel('leakage current (pA)','FontSize',32)
    %     annotation(figure4,'textbox',...
    %     outside_plot,'String',{['avg LO inc (pA):'],...
    %     [sprintf('%.1f',lcm1_avg_ramp_down_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_down_stdev(i))],...
    %     ['avg LO exc (pA):'],...
    %     [sprintf('%.1f',lcm1_avg_ramp_down_inc_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_down_inc_stdev(i))]},...
    %     'FontSize',32,'BackgroundColor',[1 1 1]);
    % % 
    %%%%%%%%%% ramp data plot of mean up chunk leakage current v chunk # %%%%%%%
        figure5= figure('Units','normalized')
        plot(xdata_1(i,:),ydata_1(i,:),...
            '^','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0); hold on;
        plot(xdata_2(i,:),ydata_2(i,:),...
            '^','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0); hold on;
        plot(xdata_3(i,:),ydata_3(i,:),...
            'v','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0); hold on;
        plot(xdata_4(i,:),ydata_4(i,:),...
            'v','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0);
            l = legend('show'); l.String = [{'exclusive HI'},{'inclusive HI'},{'exclusive LO'},{'inclusive LO'}]; l.FontSize = 32; l.Location = 'northeast outside';
        pbaspect([1.33 1 1])
        ax = gca; % current axes
        ax.FontSize = 32;
        ax.TickDir = 'out'; % make ticks point out
        title('lcm1 HI ramp means','FontSize',40)
        xlabel('chunk #','FontSize',32)
        ylabel('leakage current (pA)','FontSize',32)
        annotation(figure5,'textbox',...
        outside_plot,'String',{...
        ['avg HI ex. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_up_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_up_stdev(i))],...
        ['avg HI in. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_up_inc_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_up_inc_stdev(i))],...
        ['avg LO ex. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_down_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_down_stdev(i))],...
        ['avg LO in. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_down_inc_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_down_inc_stdev(i))]},...
        'FontSize',32,'BackgroundColor',[1 1 1]);

    %%%%%%%%% ramp data plot of HI stdev inclusive/exclusive chunk leakage current v chunk # %%%%%%%
        figure6= figure('Units','normalized')
        plot(xdata_1(i,:),ydata_1(i,:),...
        '^','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
        plot(xdata_2(i,:),ydata_2(i,:),...
        '^','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
        plot(xdata_3(i,:),ydata_3(i,:),...
        'v','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
        plot(xdata_4(i,:),ydata_4(i,:),...
        'v','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
        l = legend('show'); l.String = [{'exclusive HI'},{'inclusive HI'},{'exclusive LO'},{'inclusive LO'}]; l.FontSize = 32; l.Location = 'northeast outside';
    %     l = legend('show'); l.String = [{'exclusive'},{'inclusive'}]; l.FontSize = 32; l.Location = 'northeast outside';
    %    axis(plot_bounds)
        pbaspect([1.33 1 1])
        ax = gca; % current axes
        ax.FontSize = 32;
        ax.TickDir = 'out'; % make ticks point out
        title('lcm1 ramp stdevs','FontSize',40)
        xlabel('chunk #','FontSize',32)
        ylabel('leakage current (pA)','FontSize',32)
        annotation(figure6,'textbox',...
        outside_plot,'String',{...
        ['avg HI ex. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_up_stdev_avg_chunk(i)) ' \pm ' ...
        sprintf('%.1f',lcm1_avg_ramp_up_stdev_stdev_chunk(i))],...
        ['avg HI in. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_up_inc_stdev_avg_chunk(i)) ' \pm ' ...
        sprintf('%.1f',lcm1_avg_ramp_up_inc_stdev_stdev_chunk(i))],...
        ['avg LO ex. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_down_stdev_avg_chunk(i)) ' \pm ' ...
        sprintf('%.1f',lcm1_avg_ramp_down_stdev_stdev_chunk(i))],...
        ['avg LO in. (pA):'],...
        [sprintf('%.1f',lcm1_avg_ramp_down_inc_stdev_avg_chunk(i)) ' \pm ' ...
        sprintf('%.1f',lcm1_avg_ramp_down_inc_stdev_stdev_chunk(i))]},...
        'FontSize',32,'BackgroundColor',[1 1 1]);

    % %%%%%%%%% ramp data plot of LO stdev inclusive/exclusive chunk leakage current v chunk # %%%%%%%
    %     figure6= figure('Units','normalized')
    %     plot(xdata_1(i,:),ydata_1(i,:),...
    %     'x','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
    %     plot(xdata_2(i,:),ydata_2(i,:),...
    %     'o','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
    %     l = legend('show'); l.String = [{'exclusive'},{'inclusive'}]; l.FontSize = 32; l.Location = 'northeast outside';
    % %    axis(plot_bounds)
    %     pbaspect([1.33 1 1])
    %     ax = gca; % current axes
    %     ax.FontSize = 32;
    %     ax.TickDir = 'out'; % make ticks point out
    %     title('lcm1 ramp stdevs','FontSize',40)
    %     xlabel('chunk #','FontSize',32)
    %     ylabel('leakage current (pA)','FontSize',32)
    %     annotation(figure6,'textbox',...
    %     outside_plot,'String',{['avg LO ex. (pA):'],...
    %     [sprintf('%.1f',lcm1_avg_ramp_down_stdev_avg_chunk(i)) ' \pm ' ...
    %     sprintf('%.1f',lcm1_avg_ramp_down_stdev_stdev_chunk(i))],...
    %     ['avg LO in. (pA):'],...
    %     [sprintf('%.1f',lcm1_avg_ramp_down_inc_stdev_avg_chunk(i)) ' \pm ' ...
    %     sprintf('%.1f',lcm1_avg_ramp_down_inc_stdev_stdev_chunk(i))]},...
    %     'FontSize',32,'BackgroundColor',[1 1 1]);
    % 
    % %%%%%%%%% ramp data plot of stdev hi/lo chunk leakage current v chunk # %%%%%%%
    %     figure6= figure('Units','normalized')
    %     plot(xdata_1(i,:),ydata_1(i,:),...
    %     'x','Color', cmap(2,:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
    %     plot(xdata_2(i,:),ydata_2(i,:),...
    %     'o','Color', 'blue','MarkerSize', 8, 'LineWidth', 1.5);
    %     %errorbar(xdata_1(i,:),ydata_1(i,:),ydata_stdev_1(i,:),'o');
    %     l = legend('show'); l.String = [{'HI'},{'LO'}]; l.FontSize = 32; l.Location = 'northeast outside';
    % %    axis(plot_bounds)
    %     pbaspect([1.33 1 1])
    %     ax = gca; % current axes
    %     ax.FontSize = 32;
    %     ax.TickDir = 'out'; % make ticks point out
    %     title('lcm1 ramp stdevs','FontSize',40)
    %     xlabel('chunk #','FontSize',32)
    %     ylabel('leakage current (pA)','FontSize',32)
    %     annotation(figure6,'textbox',...
    %     outside_plot,'String',{['avg HI (pA):'],...
    %     [sprintf('%.1f',lcm1_avg_ramp_up_stdev_avg_chunk(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_up_stdev_stdev_chunk(i))],...
    %     ['avg LO (pA):'],...
    %     [sprintf('%.1f',lcm1_avg_ramp_down_stdev_avg_chunk(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_down_stdev_stdev_chunk(i))]},...
    %     'FontSize',32,'BackgroundColor',[1 1 1]);

    %%%%%%%%%% ramp data plot of stdev hi leakage current v chunk # %%%%%%%
    %     figure7= figure('Units','normalized')
    %     plot(xdata_1(i,:),ydata_1(i,:),...
    %     'x','Color', cmap(2,:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
    %     %errorbar(xdata_1(i,:),ydata_1(i,:),ydata_stdev_1(i,:)'o');
    % %    l = legend('show'); l.String = [{'HI'},{'LO'}]; l.FontSize = 32; l.Location = 'northeast outside';
    % %    axis(plot_bounds)
    %     pbaspect([1.33 1 1])
    %     ax = gca; % current axes
    %     ax.FontSize = 32;
    %     ax.TickDir = 'out'; % make ticks point out
    %     title('lcm1 ramp stdevs','FontSize',40)
    %     xlabel('chunk #','FontSize',32)
    %     ylabel('leakage current (pA)','FontSize',32)
    %     annotation(figure7,'textbox',...
    %     outside_plot,'String',{['avg HI (pA):'],...
    %     [sprintf('%.1f',lcm1_avg_ramp_up_stdev_avg_chunk(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_up_stdev_stdev_chunk(i))]
    %     },...
    %     'FontSize',32,'BackgroundColor',[1 1 1]);
    % 
    % %%%%%%% test chunk counting algorithm psvoltage v time %%%%%%%%%%%%%%%%%%%
    %    figure33= figure('Units','normalized')
    % %    subplot(num_files,1,num_files + 1 - i)
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %    'o','Color', cmap(1+i,:),'MarkerSize', 6, 'LineWidth', 1.0); hold on;
    % 
%        plot(xdata_2(i,:),ydata_2(i,:),...
%        'x','Color', 'green','MarkerSize', 6, 'LineWidth', 2.0); hold on;

%        plot(xdata_3(i,:),ydata_3(i,:),...
%       'x','Color', 'blue','MarkerSize', 6, 'LineWidth', 2.0); hold on;
%    
%        plot(xdata_4(i,:),ydata_4(i,:),...
%       'x','Color', 'cyan','MarkerSize', 6, 'LineWidth', 2.0); hold on;
%    
    %    %axis(plot_bounds)
    %    pbaspect([1.33 1 1])
    %    ax = gca; % current axes
    % %    ax.XDir = 'reverse'; % voltage decreases left to right
    %     ax.FontSize = 32;
    %    ax.TickDir = 'out'; % make ticks point out

    %%%%%%%%%%%%%% ramp data plot of lcm1 v time %%%%%%%%%%%%%%%%%%%%%%%%%
    %   subplot(num_sets,1,num_sets + 1 - i)
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %    '.','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
    %    axis(plot_bounds)
    %    pbaspect([4 1 1])
    %    ax = gca; % current axes
    %%    ax.XDir = 'reverse'; % voltage decreases left to right
    %%    ax.FontSize = 32;
    %    ax.TickDir = 'out'; % make ticks point out

    %%%%%%%%%%%%%%%%% grid of imon v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    subplot(2, 1 ,1)
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %        '-','Color', 'green','MarkerSize', 8, 'LineWidth', 2.0);
    %    subplot(2,1,2)
    %    plot(xdata_2(i,:),ydata_2(i,:),...
    %        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
        %axis(plot_bounds)
        %pbaspect([4 1 1])
    %    ax = gca; % current axes
    %    ax.FontSize = 16;
    %    ax.TickDir = 'out'; % make ticks point out

    %%%%%%%%%%%%%%%%%%%% grid of pressure v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %convert to mV
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %        '-','Color', 'green','MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %    plot(xdata_2(i,:),ydata_2(i,:),...
    %        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
    %    axis(plot_bounds)
    %    pbaspect([1.33 1 1])
    %    ax = gca; % current axes
    %%    ax.FontSize = 16;
    %    ax.TickDir = 'out'; % make ticks point out
    %    legend_titles = {'pressure' ,'leakage current'};

    %%%%%%%%%%%%%%%%%%% grid of power ps current v. time %%%%%%%%%%%%%%%%%%%%%%
    %    subplot(num_sets,1,num_sets + 1 - i)
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %        '-','Color',cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
    %%    plot(xdata_2,ydata_2(i,:),...
    %%        'k--','LineWidth',2.0);
    %    %axis(plot_bounds)
    %    pbaspect([1.33 1 1])
    %    ax = gca; % current axes
    %%    ax.FontSize = 32;
    %    ax.TickDir = 'out'; % make ticks point out

    %%%%%%%%%%%%%%%%% grid of resistance v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    subplot(num_sets,1,num_sets + 1 - i)
    %    plot(xdata_1(i,:),ydata_1(i,:),...
    %        '-','Color',cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
    %    axis(plot_bounds)
    %    pbaspect([3 1 1])
    %    ax = gca; % current axes
    %    ax.FontSize = 16;
    %    ax.TickDir = 'out'; % make ticks point out


    end

    %%%%%%%%%%%%%%%%%%%%%% Grid Super Axis Labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %subplot(num_files,1,num_files)
    %xlabel(x_label,'FontSize',32)
    %subplot(num_files,1,num_files-1)
    %ylabel(y_label,'FontSize',32)
    %subplot(num_files,1,1)
    %title(title_string,'FontSize',40)

    % title(title_string,'FontSize',40)
    % xlabel(x_label,'FontSize',32)
    % ylabel(y_label,'FontSize',32)
    % l = legend('show'); 
    %l.String = legend_titles; 
    l.FontSize = 32; 
    l.Location = 'northeast outside';
    %ax = gca; % current axes
    %ax.XDir = 'reverse' % voltage decreases left to right
    %ax.FontSize = 32;
    %ax.TickDir = 'out'; % make ticks point out

    %%%% bounds and tick labels. comment out to let matlab autoscale %%%%%%%%%%

    %ax.XTick = xtick_numbers;
    %ax.YTick = ytick_numbers;
    %axis(plot_bounds)

    %%%%                                                             %%%%%%%%%%

    %pbaspect([1.33 1 1])
    % annotation(figure1,'textbox',...
    %    outside_plot,'String',str,'FontSize',32,'BackgroundColor',[1 1 1]);
end