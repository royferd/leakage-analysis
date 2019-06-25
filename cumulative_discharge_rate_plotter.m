% % create a folder to place all saved images in

current_time = datetime('now','Format','yyyy-MM-dd-HHmmss');

newfolder = sprintf('%s-cumulative-discharge-rate-plotter',current_time(1));

current_directory = sprintf('%s',pwd);

fullpath = fullfile(current_directory,newfolder);

sim_filenames = [string('2019-06-24-153425-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-153706-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-153948-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-154151-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-154417-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-154657-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-154920-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-155201-discharge-rate-neg-5-sig.txt');
    string('2019-06-24-155415-discharge-rate-neg-5-sig.txt')];

num_files = length(sim_filenames);

sizeset = [6 Inf];

formatSpec = '%f %f %f %f %f %f';

one_big_data_set = [];

this_set_indices = [];

number_rows = 0;

for i=1:num_files
    
    fileID = fopen(fullfile(current_directory,'nb56-discharge-rates-all-sims',...
        sprintf('%s',sim_filenames(i))),'r');

    data_set = fscanf(fileID, formatSpec, sizeset);

    data_set = data_set';
    
    number_rows_this_set = length(data_set(:,1));
    
    for j = 1:number_rows_this_set
        
        number_rows = number_rows + 1;
        
        one_big_data_set(end+1,:) = [number_rows i data_set(j,:)];
        
    end
    
   this_set_indices(end+1,:) = [(number_rows - number_rows_this_set + 1) number_rows];
    
end


hv_plot_xy_errors('ramp simulation (dis)charging voltage plotted over all','',...
    1,0,fullpath,sprintf('nb56-discharge-rates-all-sims'),...
    one_big_data_set(this_set_indices(1,1):this_set_indices(1,2),1),...
    zeros(length(one_big_data_set(this_set_indices(1,1):this_set_indices(1,2),1)),1),...
    one_big_data_set(this_set_indices(1,1):this_set_indices(1,2),3),...
    one_big_data_set(this_set_indices(1,1):this_set_indices(1,2),4),...
    one_big_data_set(this_set_indices(2,1):this_set_indices(2,2),1),...
    zeros(length(one_big_data_set(this_set_indices(2,1):this_set_indices(2,2),1)),1),...
    one_big_data_set(this_set_indices(2,1):this_set_indices(2,2),3),...
    one_big_data_set(this_set_indices(2,1):this_set_indices(2,2),4),...
    one_big_data_set(this_set_indices(3,1):this_set_indices(3,2),1),...
    zeros(length(one_big_data_set(this_set_indices(3,1):this_set_indices(3,2),1)),1),...
    one_big_data_set(this_set_indices(3,1):this_set_indices(3,2),3),...
    one_big_data_set(this_set_indices(3,1):this_set_indices(3,2),4),...
    one_big_data_set(this_set_indices(4,1):this_set_indices(4,2),1),...
    zeros(length(one_big_data_set(this_set_indices(4,1):this_set_indices(4,2),1)),1),...
    one_big_data_set(this_set_indices(4,1):this_set_indices(4,2),3),...
    one_big_data_set(this_set_indices(4,1):this_set_indices(4,2),4),...
    one_big_data_set(this_set_indices(5,1):this_set_indices(5,2),1),...
    zeros(length(one_big_data_set(this_set_indices(5,1):this_set_indices(5,2),1)),1),...
    one_big_data_set(this_set_indices(5,1):this_set_indices(5,2),3),...
    one_big_data_set(this_set_indices(5,1):this_set_indices(5,2),4),...
    one_big_data_set(this_set_indices(6,1):this_set_indices(6,2),1),...
    zeros(length(one_big_data_set(this_set_indices(6,1):this_set_indices(6,2),1)),1),...
    one_big_data_set(this_set_indices(6,1):this_set_indices(6,2),3),...
    one_big_data_set(this_set_indices(6,1):this_set_indices(6,2),4),...
    one_big_data_set(this_set_indices(7,1):this_set_indices(7,2),1),...
    zeros(length(one_big_data_set(this_set_indices(7,1):this_set_indices(7,2),1)),1),...
    one_big_data_set(this_set_indices(7,1):this_set_indices(7,2),3),...
    one_big_data_set(this_set_indices(7,1):this_set_indices(7,2),4),...
    one_big_data_set(this_set_indices(8,1):this_set_indices(8,2),1),...
    zeros(length(one_big_data_set(this_set_indices(8,1):this_set_indices(8,2),1)),1),...
    one_big_data_set(this_set_indices(8,1):this_set_indices(8,2),3),...
    one_big_data_set(this_set_indices(8,1):this_set_indices(8,2),4),...
    one_big_data_set(this_set_indices(9,1):this_set_indices(9,2),1),...
    zeros(length(one_big_data_set(this_set_indices(9,1):this_set_indices(9,2),1)),1),...
    one_big_data_set(this_set_indices(9,1):this_set_indices(9,2),3),...
    one_big_data_set(this_set_indices(9,1):this_set_indices(9,2),4));