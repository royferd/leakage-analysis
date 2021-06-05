function plot_discharge_hist(discharge_vals,...
    title_string,hour_number,save_file_path,xlabel_string,ylabel_string)

% take a vector of values (positive and/or negative OK) and plot both on a
% log-log histogram.

    discharge_vals_this_hour = discharge_vals;

    [counts,edges] = histcounts(discharge_vals_this_hour);

    current_bval = ceil(abs(edges(2) - edges(1)));

    pos_indices = find(discharge_vals_this_hour > 0);

    pos_nan = nan(1,length(pos_indices));

    neg_indices = find(discharge_vals_this_hour < 0);

    neg_nan = nan(1,length(neg_indices));


    pos_values_this_hour = [discharge_vals_this_hour(pos_indices)];

    neg_values_this_hour = [abs(discharge_vals_this_hour(neg_indices))];

%     title_string_full =sprintf('%s hr %d (%d discharges)',title_string,hour_number,...
%         length([pos_values_this_hour neg_values_this_hour]));
    title_string_full =sprintf('%s (hr %d)',title_string,hour_number);


    figure2 = figure('visible','off'); 

    [~,edges] = histcounts(log10([pos_values_this_hour neg_values_this_hour]));

    min_discharge_this_hour = min([pos_values_this_hour neg_values_this_hour]);

    edge_max = ceil(10.^edges(end));

    num_edges = floor(log2(10.^edges(end))) -1;

%                         edge_max/(2.^(num_edges-1))

%                         min_discharge_this_hour
                        
    if min_discharge_this_hour > edge_max/(2.^(num_edges-2))
        
        while min_discharge_this_hour > edge_max/(2.^(num_edges-2))
            
            num_edges = num_edges - 1;
            
        end                        

    end
        
    if edge_max/(2.^(num_edges-1)) > min_discharge_this_hour

         while edge_max/(2.^(num_edges-1)) > min_discharge_this_hour

            num_edges = num_edges+1;

         end

    end

    edges_base_two = zeros(1,num_edges);

    for i = 1:num_edges

        edges_base_two(i) = edge_max/(2.^(num_edges-i));

    end


%     edges_base_two
    
    bin_width_two = edges_base_two(2) - edges_base_two(1);


    log_base =edges_base_two(2)/edges_base_two(1);

    legend_string = {};

    if isempty(pos_values_this_hour) == 0

        h3 = histogram(pos_values_this_hour,...
            edges_base_two,'LineWidth',1.5,'FaceAlpha',0.3,'FaceColor','blue',...
            'EdgeColor','blue'); hold on;

        %legend_string{end+1} = '+discharge';
        legend_string{end+1} = sprintf('+ve (%d)',length(pos_values_this_hour));

    end

    if isempty(neg_values_this_hour) == 0

        h4 = histogram(neg_values_this_hour,...
            edges_base_two,'LineWidth',1.5,'FaceAlpha',0.3,'FaceColor','red',...
            'EdgeColor','red'); hold on;

        %legend_string{end+1} = '-discharge';
        legend_string{end+1} = sprintf('-ve (%d)',length(neg_values_this_hour));

    end

    
    fontsize = 12;

    if isempty(legend_string) == 0
        
        l = legend('show'); 

        l.String = legend_string;

        l.FontSize = fontsize;
  %  l.Location = 'northeast outside';  

    end

    fig = gcf;
    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0 0 8 6];


%     fig.PaperPosition = [0 0 5.33 2.5]; % 2 x 2 grid
%     fig.PaperPosition = [0 0 5.33 3.50];
    %fig.PaperPosition = [0 0 4.67 2.67];
    fig.PaperPosition = [0 0 3.5 3.5];
    

    ax = gca;
    ax.TickDir = 'out'; % make ticks point out
    ax.FontSize = fontsize;
    ax.XScale='log';
    ax.YScale='log';
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;

%                         annotation('textbox',inside_plot,'String',...
%                             sprintf('bin width %.1f \\pm %.1f pA',current_bval),...
%                             'FontSize',16,'BackgroundColor',[1 1 1],'FitBoxToText','on');  



    left = outerpos(1);
    bottom = outerpos(2) + 0.25*ti(2);
    ax_width = outerpos(3);
    ax_height = outerpos(4) - 0.25*ti(2);

    ax.OuterPosition = [left bottom ax_width ax_height];

    title(title_string_full,'FontSize',fontsize+1)
    
    xlabel(xlabel_string,'FontSize',fontsize+1)

    ylabel(sprintf('discharges / %.0f \\times 2^{n} pA',bin_width_two),'FontSize',fontsize+1)

    print(save_file_path,'-dpng');
    print(sprintf('%s.svg',save_file_path),'-dsvg');