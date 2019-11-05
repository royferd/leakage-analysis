% % create a folder to place all saved images in
electrodes = 'Nb56';
sim_dates = '12/13/2018--5/1/2019';

analysis_file_patterns = {'*discharge-rate-pos*','*discharge-rate-neg*',...
    '*discharge-rate-zero*'};

% polarity of HV data identified in title and filename
plot_polarity_name = {'+HV', '-HV', '0HV'};

file_polarity_name = {'pos','neg','zero'};

analysis_folder_name_parent = 'nb56-all-simulation-analysis';

analysis_folder_name = '2019-07-11'; % latest nb56 analysis

current_time = datetime('now','Format','yyyy-MM-dd-HHmmss');

newfolder = sprintf('%s-cumulative-discharge-rate-plotter',current_time(1));

current_directory = sprintf('%s',pwd);

path_to_analysis_files = fullfile(current_directory,analysis_folder_name_parent,analysis_folder_name);

fullpath = fullfile(current_directory,newfolder);

sizeset = [6 Inf];

formatSpec = '%f %f %f %f %f %f';

mkdir(fullpath);

% for i = 1:length(analysis_file_patterns)
for i = 2:2
    cd(path_to_analysis_files);

%     sim_filenames = dir('*discharge-rate-pos*');
    sim_filenames = dir(cell2mat(analysis_file_patterns(i)));

    num_files = length(sim_filenames);

    one_big_data_set = [];

    this_set_indices = [];

    number_rows = 0;

    for j=1:num_files

        fileID = fopen(fullfile(current_directory,analysis_folder_name_parent,...
            analysis_folder_name,sprintf('%s',sim_filenames(j).name)),'r');

        data_set = fscanf(fileID, formatSpec, sizeset);

        data_set = data_set';

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

    for j = 1:number_rows

        dph_sigma(j,1) = one_big_data_set(j,3);

        dph_sigma(j,2) = one_big_data_set(j,4);
        
        median_sigma(j) = one_big_data_set(j,5);
        
        dph_cutoff(j,1) = one_big_data_set(j,6);
        
        dph_cutoff(j,2) = one_big_data_set(j,7);
        
        median_cutoff(j) = one_big_data_set(j,8);

    end
    
    cd(current_directory);

    hv_plot_xy_errors(sprintf('%s 5\\sigma discharge rates (%s)',electrodes,string(plot_polarity_name(i))),...
        'simulation time (hr)','discharges per hour',...
        2,'',1,2,fullpath,sprintf('%s-discharge-rates-sigma-%s',electrodes,string(file_polarity_name(i))),...
        [0 35 0 500],1:number_rows,zeros(length(1:number_rows),1),...
        dph_sigma(:,1),dph_sigma(:,2),this_set_indices);
    
    hv_plot_xy_errors(sprintf('%s 5\\sigma discharge sizes (%s)',electrodes,string(plot_polarity_name(i))),...
        'simulation time (hr)','median discharge size (pA)',...
        2,'',1,2,fullpath,sprintf('%s-median-discharges-sigma-%s',electrodes,string(file_polarity_name(i))),...
        [0 35 0 100],1:number_rows,zeros(length(1:number_rows),1),...
        median_sigma,zeros(length(median_sigma),1),this_set_indices);
    
    hv_plot_xy_errors(sprintf('%s 100 pA cutoff discharge rates (%s)',electrodes,string(plot_polarity_name(i))),...
        'simulation time (hr)','discharges per hour',...
        2,'',1,2,fullpath,sprintf('%s-discharge-rates-cutoff-%s',electrodes,string(file_polarity_name(i))),...
        [0 35 -3 70],1:number_rows,zeros(length(1:number_rows),1),...
        dph_cutoff(:,1),dph_cutoff(:,2),this_set_indices);
    
    hv_plot_xy_errors(sprintf('%s 100 pA cutoff discharge sizes (%s)',electrodes,string(plot_polarity_name(i))),...
        'simulation time (hr)','median discharge size (pA)',...
        2,'',1,2,fullpath,sprintf('%s-median-discharges-cutoff-%s',electrodes,string(file_polarity_name(i))),...
        [0 35 -40 800],1:number_rows,zeros(length(1:number_rows),1),...
        median_cutoff,zeros(length(median_cutoff),1),this_set_indices);

end

% row 3 = dph (5 sigma)
% row 4 = dph_std (5 sigma)
% row 5 = median discharge size (5 sigma)
% row 6 = dph (cutoff)
% row 7 = dph_std (cutoff)
% row 8 = median discharge size (cutoff)

% value, stdev
