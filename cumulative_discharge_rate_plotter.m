clearvars 
close all

program = which('cumulative_discharge_rate_plotter','-all');
% % create a folder to place all saved images in

% electrodes = 'Nb56';
% electrodes = 'Ti13';
% electrodes = 'Nb78';
% electrodes = 'Nb23';
electrodes = 'NET';

% sim_dates = '12/13/2018--5/1/2019';

analysis_file_patterns = {'*discharge-rate-pos*','*discharge-rate-neg*',...
    '*discharge-rate-zero*'};

% polarity of HV data identified in title and filename
plot_polarity_name = {'+HV', '-HV', '0HV'};

file_polarity_name = {'pos','neg','zero'};

% analysis_folder_name_parent = 'nb56-all-simulation-analysis';
% analysis_folder_name_parent = '2019-07-11-all-ti13-sim-analysis-files';
% analysis_folder_name_parent = '2019-07-11-nb78-all-simulation-analysis';
% analysis_folder_name_parent = '2020-03-12-nb23-analysis-files';
% analysis_folder_name_parent = '2020-06-17-NET-analysis-files';
% analysis_folder_name_parent = '2020-06-24-nb78-all-analysis-again';

% analysis_folder_name_parent = '2020-07-02-Nb78-analysis';
% analysis_folder_name_parent = '2020-07-02-Nb56-analysis';
% analysis_folder_name_parent = '2020-07-02-Ti13-analysis';
analysis_folder_name_parent = '2020-07-02-NET-analysis';
% analysis_folder_name_parent = '2020-07-02-Nb23-analysis';

% analysis_folder_name = '2019-07-11'; % latest nb56 analysis
analysis_folder_name = ''; % latest ti13 analysis

current_time = datetime('now','Format','yyyy-MM-dd-HHmmss');

newfolder = sprintf('%s-cumulative-discharge-rate-plotter',current_time(1));

current_directory = sprintf('%s',pwd);

path_to_analysis_files = fullfile(current_directory,analysis_folder_name_parent,analysis_folder_name);

fullpath = fullfile(current_directory,newfolder);

sizeset = [6 Inf];

formatSpec = '%f %f %f %f %f %f';

mkdir(fullpath);

% Nb56 bounds    
%     discharge_rate_bounds = [0 35 -100 350];
%     
%     discharge_size_bounds = [0 35 -30 70];
    
%     discharge_rate_bounds = [];
%     
%     discharge_size_bounds = [];


% Ti13 bounds
%     discharge_rate_bounds = [-5 120 -1000 8000];
%     
%     discharge_size_bounds = [-5 120 -500 200];
    
%     discharge_rate_bounds = [];
%     
%     discharge_size_bounds = [];
    

    
% Nb78 bounds
%     discharge_rate_bounds = [0 31 -500 3000];
%     
%     discharge_size_bounds = [0 31 -50 450];
    
%     discharge_rate_bounds = [];
%     
%     discharge_size_bounds = [];
    
% Nb23 bounds
%     discharge_rate_bounds = [-5 130 -200 3500];
%     
%     discharge_size_bounds = [-5 130 -100 2000];
    
%     discharge_rate_bounds = [];
%     
%     discharge_size_bounds = [];

    
% NET bounds

    discharge_rate_bounds = [0 30 -200 7000];
    
    discharge_size_bounds = [0 30 -100 2500];
    
%     discharge_rate_bounds = [];
%     
%     discharge_size_bounds = [];

    
    
    

    
    
% for i = 1:length(analysis_file_patterns)
for i = 1:2
    
    pname_sig_rate = sprintf('%s 5\\sigma discharge rates (%s)',electrodes,...
        string(plot_polarity_name(i)));
    
    pname_sig_size = sprintf('%s 5\\sigma discharge sizes (%s)',electrodes,...
        string(plot_polarity_name(i)));
    
    pname_cut_rate = sprintf('%s 100 pA cutoff discharge rates (%s)',electrodes,...
        string(plot_polarity_name(i)));
    
    pname_cut_size = sprintf('%s 100 pA cutoff discharge sizes (%s)',electrodes,...
        string(plot_polarity_name(i)));
    
    cd(path_to_analysis_files);

%     sim_filenames = dir('*discharge-rate-pos*');
    sim_filenames = dir(cell2mat(analysis_file_patterns(i)));

    num_files = length(sim_filenames);

    one_big_data_set = [];

    this_set_indices = [];

    number_rows = 0;
    
    number_rows_first_set = 0;

    for j=1:num_files

        fileID = fopen(fullfile(current_directory,analysis_folder_name_parent,...
            analysis_folder_name,sprintf('%s',sim_filenames(j).name)),'r');

        data_set = fscanf(fileID, formatSpec, sizeset);

        data_set = data_set';

        if j == 1
            
            number_rows_first_set = length(data_set(:,1));
            
        end
        
        number_rows_this_set = length(data_set(:,1));

        for k = 1:number_rows_this_set

            number_rows = number_rows + 1;

            one_big_data_set(end+1,:) = [number_rows j data_set(k,:)];

        end

       this_set_indices(end+1,:) = [(number_rows - number_rows_this_set + 1) number_rows];

       fclose(fileID);

    end
    
    dph_sigma = zeros(number_rows,2);
    dph_cutoff = zeros(number_rows,2);
    median_sigma = zeros(number_rows,1);
    median_cutoff = zeros(number_rows,1);

    dph_sigma_baseline = mean(one_big_data_set(1:number_rows_first_set,3));
    
    dph_cutoff_baseline = mean(one_big_data_set(1:number_rows_first_set,6));

    size_sigma_baseline = mean(one_big_data_set(1:number_rows_first_set,5));
    
    size_cutoff_baseline = mean(one_big_data_set(1:number_rows_first_set,8));
    
%     dph_sigma_baseline = 0.0;
%     
%     dph_cutoff_baseline = 0.0;
% 
%     size_sigma_baseline = 0.0;
%     
%     size_cutoff_baseline = 0.0;

    %     yname_rate = 'discharge rate - baseline (hr^{-1})';
%     
%     yname_size = 'median discharge size - baseline (pA)';
    
    yname_rate_sigma = sprintf('discharge rate - %.0f  (hr^{-1})',dph_sigma_baseline);
    
    yname_rate_cutoff = sprintf('discharge rate - %.0f  (hr^{-1})',dph_cutoff_baseline);
    
    yname_size_sigma = sprintf('median discharge size - %.0f  (pA)',size_sigma_baseline);
    
    yname_size_cutoff = sprintf('median discharge size - %.0f  (pA)',size_cutoff_baseline);

    for j = 1:number_rows

        dph_sigma(j,1) = one_big_data_set(j,3) - dph_sigma_baseline;

        dph_sigma(j,2) = one_big_data_set(j,4);
        
        median_sigma(j) = one_big_data_set(j,5) - size_sigma_baseline;
        
        dph_cutoff(j,1) = one_big_data_set(j,6) - dph_cutoff_baseline;
        
        dph_cutoff(j,2) = one_big_data_set(j,7);
        
        median_cutoff(j) = one_big_data_set(j,8) - size_cutoff_baseline;

    end
    
    cd(current_directory);


    
%     hv_plot_xy_errors(sprintf('%s 5\\sigma discharge rates (%s)',electrodes,string(plot_polarity_name(i))),...
%         'simulation time (hr)','discharges  per  hour',...
%         2,'',1,2,2,fullpath,sprintf('%s-discharge-rates-sigma-%s',electrodes,string(file_polarity_name(i))),...
%         discharge_rate_bounds,1:number_rows,zeros(length(1:number_rows),1),...
%         dph_sigma(:,1),dph_sigma(:,2),this_set_indices);
%     
%     hv_plot_xy_errors(sprintf('%s 5\\sigma discharge sizes (%s)',electrodes,string(plot_polarity_name(i))),...
%         'simulation time (hr)','median discharge size (pA)',...
%         2,'',1,2,2,fullpath,sprintf('%s-median-discharges-sigma-%s',electrodes,string(file_polarity_name(i))),...
%         discharge_size_bounds,1:number_rows,zeros(length(1:number_rows),1),...
%         median_sigma,zeros(length(median_sigma),1),this_set_indices);
%     
%     hv_plot_xy_errors(sprintf('%s 100 pA cutoff discharge rates (%s)',electrodes,string(plot_polarity_name(i))),...
%         'simulation time (hr)','discharges  per  hour',...
%         2,'',1,2,2,fullpath,sprintf('%s-discharge-rates-cutoff-%s',electrodes,string(file_polarity_name(i))),...
%         discharge_rate_bounds,1:number_rows,zeros(length(1:number_rows),1),...
%         dph_cutoff(:,1),dph_cutoff(:,2),this_set_indices);
%     
%     hv_plot_xy_errors(sprintf('%s 100 pA cutoff discharge sizes (%s)',electrodes,string(plot_polarity_name(i))),...
%         'simulation time (hr)','median discharge size (pA)',...
%         2,'',1,2,2,fullpath,sprintf('%s-median-discharges-cutoff-%s',electrodes,string(file_polarity_name(i))),...
%         discharge_size_bounds,1:number_rows,zeros(length(1:number_rows),1),...
%         median_cutoff,zeros(length(median_cutoff),1),this_set_indices);

    hv_plot_xy_errors(pname_sig_rate,'conditioning time (hr)',yname_rate_sigma,...
        2,'',1,2,2,fullpath,...
        sprintf('%s-discharge-rates-sigma-%s',electrodes,string(file_polarity_name(i))),...
        discharge_rate_bounds,linspace(1,number_rows,number_rows),...
        linspace(0,0,number_rows),dph_sigma(:,1),dph_sigma(:,2),this_set_indices);
    
    hv_plot_xy_errors(pname_sig_size,'conditioning time (hr)',yname_size_sigma,...
        2,'',1,2,2,fullpath,sprintf('%s-median-discharges-sigma-%s',electrodes,...
        string(file_polarity_name(i))),...
        discharge_size_bounds,linspace(1,number_rows,number_rows),linspace(0,0,number_rows),...
        median_sigma,zeros(length(median_sigma),1),this_set_indices);
    
    hv_plot_xy_errors(pname_cut_rate,'conditioning time (hr)',yname_rate_cutoff,...
        2,'',1,2,2,fullpath,sprintf('%s-discharge-rates-cutoff-%s',electrodes,...
        string(file_polarity_name(i))),...
        discharge_rate_bounds,linspace(1,number_rows,number_rows),linspace(0,0,number_rows),...
        dph_cutoff(:,1),dph_cutoff(:,2),this_set_indices);
    
    hv_plot_xy_errors(pname_cut_size,'conditioning time (hr)',yname_size_cutoff,...
        2,'',1,2,2,fullpath,sprintf('%s-median-discharges-cutoff-%s',electrodes,...
        string(file_polarity_name(i))),...
        discharge_size_bounds,linspace(1,number_rows,number_rows),linspace(0,0,number_rows),...
        median_cutoff,zeros(length(median_cutoff),1),this_set_indices);
    
    % create a file with the discharge sizes, rates, and voltages
    text_parameters_path = fullfile(fullpath,sprintf('%s-conditioning-timeline-%s-polarity.txt',...
        current_time,string(file_polarity_name(i))));

    fileID = fopen(text_parameters_path,'w');
    fprintf(fileID,'%s cumulative conditioning timeline for %s \n',string(plot_polarity_name(i)),electrodes);
    fprintf(fileID,'analysis ran at: %s \n',current_time);
    fprintf(fileID,'analysis ran by %s \n\n', program{:});

    fprintf(fileID,'hour # \t discharges per hour \t median discharge size (pA)\n\n');

    fprintf(fileID,'******************** 5 sigma data ******************** \n\n');
    for j = 1:number_rows

        fprintf(fileID,'%d: \t %.1f +/- %.1f \t %.1f \n',j,one_big_data_set(j,3),...
            one_big_data_set(j,4),one_big_data_set(j,5));

    end

    fprintf(fileID,'\n');

    fprintf(fileID,'******************** cutoff data ******************** \n\n');

    for j = 1:number_rows

        fprintf(fileID,'%d: \t %.1f +/- %.1f \t %.1f \n',j,one_big_data_set(j,6),...
            one_big_data_set(j,7),one_big_data_set(j,8));

    end

    fprintf(fileID,'\n');
    fprintf(fileID,'dph 5sigma baseline = %.1f discharges / hr\n',dph_sigma_baseline);
    fprintf(fileID,'median size 5sigma baseline = %.1f pA\n',size_sigma_baseline);
    fprintf(fileID,'dph cutoff baseline = %.1f discharges / hr\n',dph_cutoff_baseline);
    fprintf(fileID,'median size cutoff baseline = %.1f pA\n',size_cutoff_baseline);

    fclose(fileID);
    
end