clearvars 
close all

program = which('steady_state_plotter','-all');

analysis_file_pattern = '*steady-state-vs-voltage*.txt';

analysis_folder_name_parent ='2020-09-11-steady-state-analysis';

analysis_folder_name = ''; % latest ti13 analysis

current_time = datetime('now','Format','yyyy-MM-dd-HHmmss');

newfolder = sprintf('%s-steady-state-rate-plotter',current_time(1));

current_directory = sprintf('%s',pwd);

path_to_analysis_files = fullfile(current_directory,analysis_folder_name_parent,analysis_folder_name);

fullpath = fullfile(current_directory,newfolder);

%number of columns in the files with discharge rates. 

sizeset = [3 Inf];

formatSpec = '%f %f %f';

mkdir(fullpath);

cd(path_to_analysis_files);

sim_filenames = dir(analysis_file_pattern);
% sim_filenames = dir(cell2mat(analysis_file_patterns(i)));

sim_filenames(1).name
length(sim_filenames)
sprintf('%s',sim_filenames(1).name)

num_files = length(sim_filenames);

one_big_data_set = [];


number_rows = 0;

number_rows_first_set = 0;

nb56_voltage = [];
nb56_steady_state_leakage = [];
nb56_steady_state_stdev = [];

nb23_voltage = [];
nb23_steady_state_leakage = [];
nb23_steady_state_stdev = [];

nb78_voltage = [];
nb78_steady_state_leakage = [];
nb78_steady_state_stdev = [];

ti13_voltage = [];
ti13_steady_state_leakage = [];
ti13_steady_state_stdev = [];



for j=1:num_files

    fileID = fopen(fullfile(current_directory,analysis_folder_name_parent,...
        analysis_folder_name,sprintf('%s',sim_filenames(j).name)),'r');
    el_label = fgetl(fileID);

    data_set = fscanf(fileID, formatSpec, sizeset);

    data_set = data_set';

    if j == 1

        number_rows_first_set = length(data_set(:,1));

    end

    number_rows_this_set = length(data_set(:,1));
    
    voltage = data_set(:,1);
    
    steady_state_leakage = data_set(:,2);
    
    steady_state_stdev = data_set(:,3);

    fclose(fileID);
    
    if el_label == 'Nb56'
        
        nb56_voltage = voltage;
        
        nb56_steady_state_leakage = steady_state_leakage;
        
        nb56_steady_state_stdev = steady_state_stdev;
        
    elseif el_label == 'Nb23'
        
        nb23_voltage = voltage;
        
        nb23_steady_state_leakage = steady_state_leakage;
        
        nb23_steady_state_stdev = steady_state_stdev;
        
    elseif el_label == 'Nb78'
        
        nb78_voltage = voltage;
        
        nb78_steady_state_leakage = steady_state_leakage;
        
        nb78_steady_state_stdev = steady_state_stdev;
        
    elseif el_label == 'Ti13'
        
        ti13_voltage = voltage;
        
        ti13_steady_state_leakage = steady_state_leakage;
        
        ti13_steady_state_stdev = steady_state_stdev;
        
    end


end  
    
cd(current_directory)

labels = {'Nb56','Nb78','Ti13','Nb23'};
% labels = {'Nb56','Nb78'};

hv_plot_xy_errors(sprintf('steady-state leakage current'),'voltage (kV)','leakage current (pA)',...
1,'',1,labels,[],...
2,fullpath,'steady-state-vs-voltage-global',[-30,30,-20,10],...
nb56_voltage,zeros(1,length(nb56_voltage)),nb56_steady_state_leakage,zeros(1,length(nb56_voltage)));


    
% hv_plot_xy_errors(sprintf('steady-state leakage current'),'voltage (kV)','leakage current (pA)',...
% 1,'',1,labels,[],...
% 2,fullpath,'steady-state-vs-voltage-global',[],...
% nb56_voltage,zeros(1,length(nb56_voltage)),nb56_steady_state_leakage,nb56_steady_state_stdev,...
% nb78_voltage,zeros(1,length(nb78_voltage)),nb78_steady_state_leakage,nb78_steady_state_stdev);
