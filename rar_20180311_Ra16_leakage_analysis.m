clearvars -except masterlist binsize
close all
%clc

% initialize 3D array. Needs to be as larges as the biggest data set.
% if num_cols different, then sample set files have different formats
% and should not be compared with this program. program is for runs taken
% from February 1st 2017 onward, for which the picoammeter
% is manually set to 2 uA range and the LabVIEW VI has ion gauge output in the 
% data files.
% 08-24-2017 changing this to read data from jhv8-development onward, which
%  is compatible with the new bipolar power supply. There is a 7th column
%  that defines the polarity of Vmon. Also updating the pressure
%  conversion for the new all-range model.
% 11-09-2017 updated to be compatible with bipolar and unipolar EDM ramp
%  simulations. Tested that code is consistent with code used for 7/18 and
%  7/05 simulations. CHUNK can now correctly parse HI (negative) and LO
%  (positive) voltages. Voltage polarity is properly propagated for both
%  polarities. Different sample rates are now configurable by user. 
% 11-19-2017 changing adding module to plot current from 11/13 leakage
% current precision test. 
% 11-26-2017 CHUNK code now collects the data between the ramp up and ramp
% down intervals, affectionally referred to as the "trash" data. This
% includes data where the voltage is ramping up/down or is at 0. The 0 data
% will be useful for comparing the leakage when the voltage is HI/LO vs.
% off. The ramping data will be useful for comparing the symmetry of the
% charging up/down currents.

set(0, 'DefaultFigureRenderer', 'Painters');

EDM_sim = 1;        % 0 = not an EDM sim
                    % 1 = EDM sim. CHUNK stuff applies

inclusive_data = 1; % 0 = do not incorporate ramp data points into ramp CHUNKS (exclusive)
%                      1 = incorporate ramp data (inclusive)
power_supply = 1; % 0 = unipolar Acopian 
%                   1 = bipolar AK (installed and tested 8-18-2017)
%                   2 = Ra EDM Spellman power supply for 2016 run

pressure_gauge = 1; % 0 = ion gauge 
%                     1 = all-range gauge (installed 7-18-2017)

sample_size_setting = 1; % # samples taken per measurement
                            % 0: Spinlab sample set, 10 samples by default
                            % 1: fast data, 1 
                            %sample = # samples / sampling frequency.

sample_rate = 1 ; % 0 = data saved every 0.02 min (8192 samples / 8 kHz)
                  % 1 = data saved every 0.01 min (8192 samples / 16 kHz)
                  % 2 = data saved every 0.056 min ( guessing 1024 samples / 18 kHz)
leakage_sensitivity_test = 0; % 0 = not a sensitivity test
                              % 1 = leakage sensitivity test

file_struct = dir('*hv-1.txt');
num_files = length(file_struct); 

num_cols = 0;
num_rows = 0;
filenames = cell(num_files,1);

if sample_size_setting == 0
    sample_size = 10;
else 
    sample_size = 1;
end
                  
if power_supply == 0 || power_supply == 1

    for i = 1:num_files
        set = dlmread(file_struct(i).name);
        if length(set(:,1)) > num_rows
            num_rows = length(set(:,1));
        end
        if length(set(1,:)) > num_cols
            num_cols = length(set(1,:));
        end
        filenames{i} = file_struct(i).name(1:10);
    end
    
    data = zeros(num_files,num_rows,num_cols);    

    for i = 1:num_files
        set = dlmread(file_struct(i).name);
        set_rows = length(set(:,1));
        set_cols = length(set(1,:));
        data(i,1:set_rows,1:set_cols) = dlmread(file_struct(i).name);
    end
    
elseif power_supply == 2
    
    sizeset = [4 Inf];
    formatSpec = '{%f, %f, %f, %f}\n';

    for i = 1:num_files
        fileID = fopen(file_struct(i).name,'r');
        set = fscanf(fileID, formatSpec, sizeset);
        set = set';
        
        if length(set(:,1)) > num_rows
            num_rows = length(set(:,1));
        end
        if length(set(1,:)) > num_cols
            num_cols = length(set(1,:));
        end

        filenames{i} = file_struct(i).name(1:10);
        fclose(fileID);
    end
    
    
    data = zeros(num_files,num_rows,num_cols); 
    
    
    for i = 1:num_files
        fileID = fopen(file_struct(i).name,'r');
        set = fscanf(fileID, formatSpec, sizeset);
        set = set';
        
        for j = 1:length(set(:,1))
            data(i,j,1) = set(j,1);
            data(i,j,2) = set(j,2);
            data(i,j,3) = set(j,3);
            data(i,j,4) = set(j,4);
            
        end
    fclose(fileID);
    end    
end
        



%data = dlmread('2017-02-01-163817-hv-sample-set-data.txt');
%numpoints = length(data);




gap_size = zeros(num_files,1);
%gap_size(1) = .19; %gap size in cm
%gap_size(2) = .15; 
%gap_size(3) = .12;
%gap_size(4) = .1;
%gap_size(5) = .08;
%gap_size(6) = .06;
%gap_size(7) = .05;
%gap_size(8) = .04;

%preallocate arrays

time = zeros(num_files,num_rows,1);
time_raw = zeros(num_files,num_rows,1);

vmon_avg = zeros(num_files,num_rows,1);
vmon_avg_mag = zeros(num_files,num_rows,1);
vmon_avg_raw = zeros(num_files,num_rows,1);
vmon_avg_raw_mag = zeros(num_files,num_rows,1);
vmon_avg_stdev = zeros(num_files,num_rows,1);
vmon_avg_stdev_raw = zeros(num_files,num_rows,1);
vmon_weight_raw = zeros(num_files,num_rows,1);
vmon_weight = zeros(num_files,num_rows,1);
vmon_avg_wt = zeros(num_files,num_rows,1);
vmon_avg_wt_raw = zeros(num_files,num_rows,1);
vmon_avg_wt_raw_mag = zeros(num_files,num_rows,1);

imon_avg = zeros(num_files,num_rows,1);
imon_avg_raw = zeros(num_files,num_rows,1);
imon_avg_stdev = zeros(num_files,num_rows,1);
imon_avg_stdev_raw = zeros(num_files,num_rows,1);
imon_weight_raw = zeros(num_files,num_rows,1);
imon_weight = zeros(num_files,num_rows,1);
imon_avg_wt = zeros(num_files,num_rows,1);
imon_avg_wt_raw = zeros(num_files,num_rows,1);

lcm1_avg = zeros(num_files,num_rows,1);
lcm1_avg_log_neg = zeros(num_files,num_rows,2);
lcm1_avg_log_pos = zeros(num_files,num_rows,2);
lcm1_avg_raw = zeros(num_files,num_rows,1);
lcm1_avg_stdev = zeros(num_files,num_rows,1);
lcm1_avg_stdev_raw = zeros(num_files,num_rows,1);
lcm1_avg_wt = zeros(num_files,num_rows,1);
lcm1_avg_wt_raw = zeros(num_files,num_rows,1);

pressure_avg_raw = zeros(num_files,num_rows,1);
pressure_avg_stdev_raw = zeros(num_files,num_rows,1);
pressure_avg = zeros(num_files,num_rows,1);
pressure_avg_stdev = zeros(num_files,num_rows,1);
polarity_avg_raw = zeros(num_files,num_rows,1);
lcm1_weight = zeros(num_files,num_rows,1);
lcm1_weight_raw = zeros(num_files,num_rows,1);
field_avg_raw = zeros(num_files,num_rows,1);
field_avg = zeros(num_files,num_rows,1);
ohm_avg = zeros(num_files,num_rows,1);
ohm_avg_raw = zeros(num_files,num_rows,1);

polarity_sign = ones(num_files,num_rows);

%linear fit arrays
wxsq = zeros(num_files, num_rows, 1);
wxy = zeros(num_files, num_rows, 1); 
wy = zeros(num_files, num_rows, 1);
wx = zeros(num_files, num_rows, 1);
wx_sq = zeros(num_files, num_rows, 1);

%linearly fitted parameter arrays
del = zeros(num_files, 1);
a = zeros(num_files, 1);
b = zeros(num_files, 1);
a_stdev = zeros(num_files, 1);
b_stdev = zeros(num_files, 1);
resistance = zeros(num_files, 1);
resistance_stdev = zeros(num_files, 1);
num_fitpoints = zeros(num_files, 1);
fitline = zeros(num_files, 2, num_rows*20);

%residual array
residual = zeros(num_files,num_rows,1);

numpoints = zeros(num_files,1);

%I add an extra 0.01 minutes

if sample_rate == 0
    sampling_time = 0.03;
elseif sample_rate == 1
    sampling_time = 0.02;
elseif sample_rate == 2
    sampling_time = 0.001;
end

if power_supply == 0
    disp('Using -30 kV Acopian power supply')
    for i = 1:num_files
        numpoints(i) = num_rows;

        for j =1:numpoints(i)
           if data(i,j,2) == 0.0
               numpoints(i) = j-1;
               break;
           end
        end

        for j =1:numpoints(i)
            time_raw(i,j) = data(i,j,1);
            vmon_avg_raw_mag(i,j) = data(i,j,3);
            imon_avg_raw(i,j) = data(i,j,4);
            lcm1_avg_raw(i,j) = data(i,j,5);
            pressure_avg_raw(i,j) = data(i,j,7);
            polarity_avg_raw(i,j) = 3.4;
            vmon_avg_stdev_raw(i,j) = data(i,j,9)/sqrt(sample_size_setting);
            vmon_weight_raw(i,j) = abs((vmon_avg_stdev_raw(i,j))^(-2));
            imon_avg_stdev_raw(i,j) = data(i,j,10)/sqrt(sample_size_setting);
            imon_weight_raw(i,j) = abs((imon_avg_stdev_raw(i,j))^(-2));
            lcm1_avg_stdev_raw(i,j) = data(i,j,11)/sqrt(sample_size_setting);
            pressure_avg_stdev_raw(i,j) = data(i,j,12);
            lcm1_weight_raw(i,j) = abs(1/(lcm1_avg_stdev_raw(i,j)^2));
            ohm_avg_raw(i,j) = vmon_avg(i,j)/imon_avg(i,j)*1e3; % Mohm
            vmon_avg_wt_raw_mag(i,j) = vmon_avg_raw_mag(i,j)*vmon_weight_raw(i,j);
            imon_avg_wt_raw(i,j) = imon_avg_raw(i,j)*imon_weight_raw(i,j);
            field_avg_raw(i,j) = vmon_avg_raw_mag(i,j)/gap_size(i); % (raw vmon / cm)
        end
    time(i,:) = (time_raw(i,:) - time_raw(i,1))/60.0; %min
    end
elseif power_supply == 1
    disp('Using +/- 30 kV Applied Kilovolts power supply')
    for i = 1:num_files
        numpoints(i) = num_rows;

        for j =1:numpoints(i)
           if data(i,j,2) == 0.0
               numpoints(i) = j-1;
               break;
           end
        end

        for j =1:numpoints(i)
            time_raw(i,j) = data(i,j,1);
            vmon_avg_raw_mag(i,j) = data(i,j,3);
            imon_avg_raw(i,j) = data(i,j,4);
            lcm1_avg_raw(i,j) = data(i,j,5);
            pressure_avg_raw(i,j) = data(i,j,7);
            polarity_avg_raw(i,j) = data(i,j,8); % ~20 mV Low (positive), ~ 3.4 V high (negative)
            vmon_avg_stdev_raw(i,j) = data(i,j,10)/sqrt(sample_size_setting); %don't think reduced sqrt is true for fast data
            vmon_weight_raw(i,j) = abs((vmon_avg_stdev_raw(i,j))^(-2));
            imon_avg_stdev_raw(i,j) = data(i,j,11)/sqrt(sample_size_setting);
            imon_weight_raw(i,j) = abs((imon_avg_stdev_raw(i,j))^(-2));
            lcm1_avg_stdev_raw(i,j) = data(i,j,12)/sqrt(sample_size_setting);
            pressure_avg_stdev_raw(i,j) = data(i,j,14);
            lcm1_weight_raw(i,j) = abs(1/(lcm1_avg_stdev_raw(i,j)^2));
            lcm1_avg_wt_raw(i,j) = lcm1_avg_raw(i,j)*lcm1_weight_raw(i,j);
            
            ohm_avg_raw(i,j) = vmon_avg(i,j)/imon_avg(i,j)*1e3; % Mohm
            vmon_avg_wt_raw_mag(i,j) = vmon_avg_raw_mag(i,j)*vmon_weight_raw(i,j);
            imon_avg_wt_raw(i,j) = imon_avg_raw(i,j)*imon_weight_raw(i,j);
            field_avg_raw(i,j) = vmon_avg_raw_mag(i,j)/gap_size(i); % (raw vmon / cm)
        end
    time(i,:) = (time_raw(i,:) - time_raw(i,1))/60.0; %min
    end
elseif power_supply == 2
    disp('Using Ra EDM Spellman power supply')
    for i = 1:num_files
        numpoints(i) = num_rows;

        for j =2:numpoints(i)
           if data(i,j,2) == 0.0
               numpoints(i) = j-1;
               break;
           end
        end

        for j =1:numpoints(i)
            time_raw(i,j) = data(i,j,2); % ms, s, min?
            lcm1_avg_raw(i,j) = data(i,j,3);
            lcm1_avg_stdev_raw(i,j) = data(i,j,4); %sample size?
            lcm1_weight_raw(i,j) = lcm1_avg_stdev_raw(i,j)^(-2);
            lcm1_avg_wt_raw(i,j) = lcm1_avg_raw(i,j)*lcm1_weight_raw(i,j);
        end
    time(i,:) = (time_raw(i,:) - time_raw(i,1))/60.; %min
    end
end

for i = 1:num_files
    for j =1:numpoints(i)
        if polarity_avg_raw(i,j) > 1.5
            polarity_sign(i,j) = -1;
        else polarity_sign(i,j) = +1;
        end
            
    end
end
%scale factors
if power_supply == 0
    disp('Using -30 kV Acopian analog output scaling factors')
    vmon_avg_scale = 1e1; %(kV)
    imon_avg_scale = 1e3; % (uA)
    lcm1_avg_scale = (-2)*1e4; %(pA)
    ohm_avg_scale = 1e3; % Mohm
elseif power_supply == 1
    disp('Using +/-30 kV Applied Kilovolts analog output scaling factors')
    vmon_avg_scale = 30/10; % (kV)
    imon_avg_scale = 250/10; %(uA)
    lcm1_avg_scale = (-2)*1e4; %(pA)
    ohm_avg_scale = 1e3; % Mohm
elseif power_supply == 2
    disp('Using Ra EDM transimpedance amplifier')
    lcm1_avg_scale = (1/0.4762)*1e3; %(pA)
    vmon_avg_scale = 0;
    imon_avg_scale = 0;
    ohm_avg_scale = 0;
end

% for i = 1:num_files
%     for j= 1:numpoints
%         vmon_avg_raw(i,j) = vmon_avg_raw_mag(i,j)*polarity_sign(i,j);
%         vmon_avg_wt_raw(i,j) = vmon_avg_wt_raw_mag(i,j)*polarity_sign(i,j);
%     end
% end

if pressure_gauge == 0
    disp('Using ion gauge scaling factors')
    for i = 1:num_files
        for j = 1:numpoints(i)
            pressure_avg(i,j) = 10^(pressure_avg_raw(i,j) - 10); %(Torr)
            pressure_avg_stdev(i,j) = 10^(pressure_avg_stdev_raw(i,j) - 10); %(Torr)
        end
    end
    
elseif pressure_gauge == 1
    disp('Using all-range pressure gauge scaling factors')
    for i = 1:num_files
        for j = 1:numpoints(i)
            pressure_avg(i,j) = 10^(2*pressure_avg_raw(i,j) - 11); %(Torr)
            pressure_avg_stdev(i,j) = 10^(2*pressure_avg_stdev_raw(i,j) - 11); %(Torr)
        end
    end
end

%average time step in seconds
time_step = zeros(num_files,1);
time_difference_raw = zeros(num_files,num_rows);

for i = 1:num_files
    for j = 1:numpoints(i)
        vmon_avg_raw(i,j) = vmon_avg_raw_mag(i,j)*polarity_sign(i,j);
        vmon_avg_wt_raw(i,j) = vmon_avg_wt_raw_mag(i,j)*polarity_sign(i,j);
        vmon_avg(i,j) = vmon_avg_raw(i,j)*vmon_avg_scale;
        vmon_avg_mag(i,j) = vmon_avg_raw_mag(i,j)*vmon_avg_scale;
        vmon_avg_stdev(i,j) = vmon_avg_stdev_raw(i,j)*vmon_avg_scale;
        vmon_weight(i,j) = (vmon_avg_stdev(i,j))^-2;
        
        imon_avg(i,j) = imon_avg_raw(i,j)*imon_avg_scale;
        imon_avg_stdev(i,j) = imon_avg_stdev_raw(i,j)*imon_avg_scale;
        imon_weight(i,j) = (vmon_avg_stdev(i,j))^-2;
        
        lcm1_avg(i,j) = lcm1_avg_raw(i,j)*lcm1_avg_scale;
%        lcm1_avg_log(i,j) = lcm1_avg_log(i,j)*lcm1_avg_scale;
        lcm1_avg_stdev(i,j) = lcm1_avg_stdev_raw(i,j)*abs(lcm1_avg_scale);
        lcm1_weight(i,j) = 1/(lcm1_avg_stdev(i,j)^2);
        lcm1_avg_wt(i,j) = lcm1_avg_wt_raw(i,j)*lcm1_avg_scale;
        
        ohm_avg(i,j) = vmon_avg(i,j)/imon_avg(i,j)*ohm_avg_scale;
        vmon_avg_wt(i,j) = vmon_avg_wt_raw(i,j)*vmon_avg_scale;
        imon_avg_wt(i,j) = imon_avg_wt_raw(i,j)*imon_avg_scale;
        field_avg(i,j) = field_avg_raw(i,j)*vmon_avg_scale; %(kV/cm)
    end
    
    for j = round(0.4 * numpoints(i),0):round(0.6 * numpoints(i),0)
        time_difference_raw(i,j) = time_raw(i,j) - time_raw(i,j-1);
    end
    time_step(i) = mean(time_difference_raw(i,round(0.4 * numpoints(i),0):round(0.6 * numpoints(i),0)));
end

current_source_avg = zeros(num_files,35); %(in pA)
current_source_avg_stdev = zeros(num_files,35); %(in pA)
for i =1:num_files
    current_source_avg(i,1) = 0.0000;
    current_source_avg(i,2) = -2.0000;
    current_source_avg(i,3) = -1.6384;
    current_source_avg(i,4) = -0.8192;
    current_source_avg(i,5) = -0.4096;
    current_source_avg(i,6) = -0.2048;
    current_source_avg(i,7) = -0.1024;
    current_source_avg(i,8) = -0.0512;
    current_source_avg(i,9) = -0.0256;
    current_source_avg(i,10) = -0.0128;
    current_source_avg(i,11) = -0.0064;
    current_source_avg(i,12) = -0.0032;
    current_source_avg(i,13) = -0.0016;
    current_source_avg(i,14) = -0.0008;
    current_source_avg(i,15) = -0.0004;
    current_source_avg(i,16) = -0.0002;
    current_source_avg(i,17) = -0.0001;
    current_source_avg(i,18) = 0.0000;
    current_source_avg(i,19) = +0.0001;
    current_source_avg(i,20) = +0.0002;
    current_source_avg(i,21) = +0.0004;
    current_source_avg(i,22) = +0.0008;
    current_source_avg(i,23) = +0.0016;
    current_source_avg(i,24) = +0.0032;
    current_source_avg(i,25) = +0.0064;
    current_source_avg(i,26) = +0.0128;
    current_source_avg(i,27) = +0.0256;
    current_source_avg(i,28) = +0.0512;
    current_source_avg(i,29) = +0.1024;
    current_source_avg(i,30) = +0.2048;
    current_source_avg(i,31) = +0.4096;
    current_source_avg(i,32) = +0.8192;
    current_source_avg(i,33) = +1.6384;
    current_source_avg(i,34) = +2.0000;
    current_source_avg(i,35) = 0.0000;
    
    for j =1:35
        current_source_avg_stdev(i,j) = abs(current_source_avg(i,j))*0.004 + 0.002;
    end
end
start_point = ones(num_files,1);
end_point = zeros(num_files,1);
start_offset = zeros(num_files,2,1);
end_offset = zeros(num_files,2,1);
offset_length = 10;


if (power_supply == 0) || (power_supply == 1)
    %define offset chunks arbitrarily in case we don't have a lot of points
    for i =1:num_files
        start_offset(i,1) = 1;
        start_offset(i,2) = 2;
    end

    %define offset chunks to be 10 points at beginning and end of run
    for i = 1:num_files
        for j = (offset_length+1):numpoints(i)
    % look for 100 V bump indicating the HV switch is on. Note: PS resistance 
    % fluctuates wildly under 1.1 kV
           if  abs(vmon_avg(i,j)) > 0.1 && numpoints(i)>1000
               start_point(i) = j;  % supply is turned on
               start_offset(i,1,1) = j - offset_length;
               start_offset(i,2,1) = j - 1;
               break;
           end
        end
    end

    %find_error = zeros(2,num_files);
    for i = 1:num_files
        for j = 1:numpoints(i)
           if  (numpoints(i) - j > 0 && vmon_avg_mag(i,numpoints(i)-j) ... 
                   < 0.1 && abs(lcm1_avg(i,numpoints(i)-j)) < 2*1e3 && ... 
                   time(i,numpoints(i))  - time(i,numpoints(i)-j) < 10 && ... 
                   numpoints(i) > 1000)
               end_point(i) = numpoints(i) - j;
               end_offset(i,1,1) = end_point(i) - offset_length;
               end_offset(i,2,1) = end_point(i) - 1;
               break;
           else
               end_point(i) = numpoints(i);
               end_offset(i,1,1) = start_offset(i,1,1);
               end_offset(i,2,1) = start_offset(i,2,1);
           end
        end
    end
elseif power_supply == 2
    %define offset chunk to start 12 s after file starts
    for i =1:num_files
        a = find(time(i,:) < 0.2);
        start_offset(i,1) = a(end);
        start_offset(i,2) = a(end) + offset_length -1;
    end


    for i = 1:num_files
        for j = 1:numpoints(i)
           if  (numpoints(i) - j > 0 && vmon_avg_mag(i,numpoints(i)-j) ... 
                   < 0.1 && abs(lcm1_avg(i,numpoints(i)-j)) < 2*1e3 && ... 
                   time(i,numpoints(i))  - time(i,numpoints(i)-j) < 10 && ... 
                   numpoints(i) > 1000)
               end_point(i) = numpoints(i) - j;
               end_offset(i,1,1) = end_point(i) - offset_length;
               end_offset(i,2,1) = end_point(i) - 1;
               break;
           else
               end_point(i) = numpoints(i);
               end_offset(i,1,1) = start_offset(i,1,1);
               end_offset(i,2,1) = start_offset(i,2,1);
           end
        end
    end
end


%subtract off the offsets, i.e. when power supply voltage = 0
vmon_avg_offset_raw = zeros(num_files,1);
vmon_avg_offset = zeros(num_files,1);
vmon_avg_offset_raw_wt_sum = zeros(num_files,1);
vmon_weight_offset_sum = zeros(num_files,1);
vmon_avg_offset_raw_end_wt = zeros(num_files,offset_length);
vmon_avg_offset_raw_start_wt = zeros(num_files,offset_length);
vmon_avg_offset_raw_start = zeros(num_files,1);
vmon_avg_offset_raw_end = zeros(num_files,1);
vmon_weight_offset_raw_start = zeros(num_files,1);
vmon_weight_offset_raw_end = zeros(num_files,1);

imon_avg_offset_raw = zeros(num_files,1);
imon_avg_offset = zeros(num_files,1);
imon_avg_offset_raw_wt_sum = zeros(num_files,1);
imon_weight_offset_sum = zeros(num_files,1);
imon_avg_offset_raw_end_wt = zeros(num_files,offset_length);
imon_avg_offset_raw_start_wt = zeros(num_files,offset_length);
imon_avg_offset_raw_start = zeros(num_files,1);
imon_avg_offset_raw_end = zeros(num_files,1);
imon_weight_offset_raw_start = zeros(num_files,1);
imon_weight_offset_raw_end = zeros(num_files,1);

lcm1_avg_offset_raw = zeros(num_files,1);
lcm1_avg_offset_raw_log = zeros(num_files,1);
lcm1_acc_test_offset_raw = zeros(num_files,1);
lcm1_avg_offset = zeros(num_files,1);
lcm1_avg_offset_log = zeros(num_files,1);
lcm1_acc_test_offset = zeros(num_files,1);
lcm1_avg_offset_raw_start_wt = zeros(num_files,offset_length);
lcm1_avg_offset_raw_start = zeros(num_files,1);
lcm1_avg_offset_raw_end_wt = zeros(num_files,offset_length);
lcm1_avg_offset_raw_end = zeros(num_files,1);
lcm1_weight_offset_raw_start = zeros(num_files,1);
lcm1_weight_offset_raw_end = zeros(num_files,1);
lcm1_avg_offset_raw_wt_sum = zeros(num_files,1);
lcm1_weight_offset_sum = zeros(num_files,1);
lcm1_avg_offset_raw = zeros(num_files,1);
lcm1_avg_offset = zeros(num_files,1);
field_avg_offset = zeros(num_files,1);

for i = 1:num_files
    for j = start_offset(i,1):start_offset(i,2);
        vmon_avg_offset_raw_start_wt(i,j-start_offset(i,1)+1) = vmon_avg_raw(i,j)*vmon_weight_raw(i,j);
        imon_avg_offset_raw_start_wt(i,j-start_offset(i,1)+1) = imon_avg_raw(i,j)*imon_weight_raw(i,j);
        lcm1_avg_offset_raw_start_wt(i,j-start_offset(i,1)+1) = lcm1_avg_raw(i,j)*lcm1_weight_raw(i,j);
    end
end

for i = 1:num_files
    vmon_avg_offset_raw_start(i) = sum(vmon_avg_offset_raw_start_wt(i,:))/sum(vmon_weight_raw(i,start_offset(i,1):start_offset(i,2)));
    vmon_weight_offset_raw_start(i) = (std(vmon_avg_raw(i,start_offset(i,1):start_offset(i,2)),vmon_weight_raw(i,start_offset(i,1):start_offset(i,2))))^(-2);
    
    imon_avg_offset_raw_start(i) = sum(imon_avg_offset_raw_start_wt(i,:))/sum(imon_weight_raw(i,start_offset(i,1):start_offset(i,2)));
    imon_weight_offset_raw_start(i) = (std(imon_avg_raw(i,start_offset(i,1):start_offset(i,2)),imon_weight_raw(i,start_offset(i,1):start_offset(i,2))))^(-2);    
    
    lcm1_avg_offset_raw_start(i) = sum(lcm1_avg_offset_raw_start_wt(i,:))/sum(lcm1_weight_raw(i,start_offset(i,1):start_offset(i,2)));
    lcm1_weight_offset_raw_start(i) = (std(lcm1_avg_raw(i,start_offset(i,1):start_offset(i,2)),lcm1_weight_raw(i,start_offset(i,1):start_offset(i,2))))^(-2);
end

for i = 1:num_files
    for j = end_offset(i,1):end_offset(i,2);
        vmon_avg_offset_raw_end_wt(i,j-end_offset(i,1)+1) = vmon_avg_raw(i,j)*vmon_weight_raw(i,j);
        imon_avg_offset_raw_end_wt(i,j-end_offset(i,1)+1) = imon_avg_raw(i,j)*imon_weight_raw(i,j);
        lcm1_avg_offset_raw_end_wt(i,j-end_offset(i,1)+1) = lcm1_avg_raw(i,j)*lcm1_weight_raw(i,j);
    end
end

for i = 1:num_files
    vmon_avg_offset_raw_end(i) = sum(vmon_avg_offset_raw_end_wt(i,:))/sum(vmon_weight_raw(i,end_offset(i,1):end_offset(i,2)));
    vmon_weight_offset_raw_end(i) = (std(vmon_avg_raw(i,end_offset(i,1):end_offset(i,2)),vmon_weight_raw(i,end_offset(i,1):end_offset(i,2))))^(-2);

    imon_avg_offset_raw_end(i) = sum(imon_avg_offset_raw_end_wt(i,:))/sum(imon_weight_raw(i,end_offset(i,1):end_offset(i,2)));
    imon_weight_offset_raw_end(i) = (std(imon_avg_raw(i,end_offset(i,1):end_offset(i,2)),imon_weight_raw(i,end_offset(i,1):end_offset(i,2))))^(-2);
    
    lcm1_avg_offset_raw_end(i) = sum(lcm1_avg_offset_raw_end_wt(i,:))/sum(lcm1_weight_raw(i,end_offset(i,1):end_offset(i,2)));
    lcm1_weight_offset_raw_end(i) = (std(lcm1_avg_raw(i,end_offset(i,1):end_offset(i,2)),lcm1_weight_raw(i,end_offset(i,1):end_offset(i,2))))^(-2);
end

for i = 1:num_files
    vmon_avg_offset_raw_wt_sum(i) = vmon_avg_offset_raw_start(i)*vmon_weight_offset_raw_start(i) + vmon_avg_offset_raw_end(i)*vmon_weight_offset_raw_end(i);
    vmon_weight_offset_sum(i) = vmon_weight_offset_raw_start(i)+vmon_weight_offset_raw_end(i);

    imon_avg_offset_raw_wt_sum(i) = imon_avg_offset_raw_start(i)*imon_weight_offset_raw_start(i) + imon_avg_offset_raw_end(i)*imon_weight_offset_raw_end(i);
    imon_weight_offset_sum(i) = imon_weight_offset_raw_start(i)+imon_weight_offset_raw_end(i);
    
    lcm1_avg_offset_raw_wt_sum(i) = lcm1_avg_offset_raw_start(i)*lcm1_weight_offset_raw_start(i) + lcm1_avg_offset_raw_end(i)*lcm1_weight_offset_raw_end(i);
    lcm1_weight_offset_sum(i) = lcm1_weight_offset_raw_start(i)+lcm1_weight_offset_raw_end(i);
end

for i =1:num_files
    vmon_avg_offset_raw(i) = vmon_avg_offset_raw_wt_sum(i)/vmon_weight_offset_sum(i);
    vmon_avg_offset(i) = vmon_avg_offset_raw(i)*vmon_avg_scale;
    
    imon_avg_offset_raw(i) = imon_avg_offset_raw_wt_sum(i)/imon_weight_offset_sum(i);
    imon_avg_offset(i) = imon_avg_offset_raw(i)*imon_avg_scale;
    
    field_avg_offset(i) = vmon_avg_offset(i) / gap_size(i);
    
    lcm1_avg_offset_raw(i) = lcm1_avg_offset_raw_wt_sum(i)/lcm1_weight_offset_sum(i);
    lcm1_acc_test_offset_raw(i) = ((lcm1_avg_raw(i,1)*lcm1_weight_raw(i,1)+...
        lcm1_avg_raw(i,18)*lcm1_weight_raw(i,18)+...
        lcm1_avg_raw(i,35)*lcm1_weight_raw(i,35))/(lcm1_weight_raw(i,1)+...
        lcm1_weight_raw(i,18)+lcm1_weight_raw(i,35)));
    
    lcm1_avg_offset(i) = lcm1_avg_offset_raw(i)*lcm1_avg_scale;
    lcm1_acc_test_offset(i) = lcm1_acc_test_offset_raw(i)*lcm1_avg_scale;
end

% calculate log10s of leakage current. For negative leakage, put sign
% in front of the log of the magnitude of the leakage.

for i = 1:num_files
    
    for j = 1:numpoints(i)
        
        if lcm1_avg(i,j) - lcm1_avg_offset(i) < 0
            lcm1_avg_log_neg(i,j,1) = time(i,j);            
            lcm1_avg_log_neg(i,j,2) = log10(abs(lcm1_avg(i,j) - lcm1_avg_offset(i)));
        else
            lcm1_avg_log_pos(i,j,1) = time(i,j);
            lcm1_avg_log_pos(i,j,2) = log10(lcm1_avg(i,j) - lcm1_avg_offset(i));
        end
        
        
    end
    
    for j = 1:numpoints(i)
         
        if lcm1_avg_log_neg(i,j,2) < 1
            lcm1_avg_log_neg(i,j,1) = NaN;
            lcm1_avg_log_neg(i,j,2) = NaN;
        end
        
        if lcm1_avg_log_pos(i,j,2) < 1
            lcm1_avg_log_pos(i,j,1) = NaN;
            lcm1_avg_log_pos(i,j,2) = NaN;
        end
        
    end
end


field_max = zeros(num_files,1);
leakage_max = zeros(num_files,1);
voltage_max = zeros(num_files,1);
for i = 1:num_files
    field_max(i) = max(abs(field_avg(i,start_point(i):end_point(i)))) - field_avg_offset(i);
    voltage_max(i) = max(abs(vmon_avg(i,start_point(i):end_point(i)))) - vmon_avg_offset(i);
    [M,I] = max(abs(lcm1_avg(i,start_point(i):end_point(i))));
    leakage_max(i) = lcm1_avg(i,I+ start_point(i) - 1) - lcm1_avg_offset(i);
end
    

%%%%%%%%%%%%%%%%%%%%%%%% ramp test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%08-28-2017 put this into github. implement identical code, both
%%polarities.


vmon_avg_avg = zeros(num_files,1);
vmon_avg_total_stdev = zeros(num_files,1);
vmon_ramp_high = zeros(num_files,1);
vmon_ramp_low = zeros(num_files,1);
vmon_ramp_high_deviation = zeros(num_files,1);
vmon_ramp_low_deviation = zeros(num_files,1);

lcm1_avg_avg = zeros(num_files,1);
lcm1_avg_total_stdev = zeros(num_files,1);
lcm1_ramp_high = zeros(num_files,1);
lcm1_ramp_low = zeros(num_files,1);
lcm1_ramp_high_deviation = zeros(num_files,1);
lcm1_ramp_low_deviation = zeros(num_files,1);

% avg around 0. To separate HI/LO values from 0 values, ignore voltages
% within pm 1 kV of vmon_avg_avg

for i = 1:num_files    
    vmon_avg_avg(i) = mean(vmon_avg(i,start_point(i):end_point(i)));
    vmon_avg_total_stdev(i) = std(vmon_avg(i,start_point(i):end_point(i)));
    vmon_hi_vals = find(vmon_avg(i,:) < vmon_avg_avg(i) - 5);
    vmon_lo_vals = find(vmon_avg(i,:) > vmon_avg_avg(i) + 5);
    vmon_ramp_high_deviation(i) = std(vmon_avg(i,vmon_hi_vals),vmon_weight(i,vmon_hi_vals));
    vmon_ramp_low_deviation(i) = std(vmon_avg(i,vmon_lo_vals),vmon_weight(i,vmon_lo_vals));
    vmon_ramp_high(i) = sum(vmon_avg_wt(i,vmon_hi_vals))/sum(vmon_weight_raw(i,vmon_hi_vals));
    vmon_ramp_low(i) = sum(vmon_avg_wt(i,vmon_lo_vals))/sum(vmon_weight_raw(i,vmon_lo_vals));    
    
    lcm1_avg_avg(i) = mean(lcm1_avg(i,start_point(i):end_point(i)));
    lcm1_avg_total_stdev(i) = std(lcm1_avg(i,start_point(i):end_point(i)),lcm1_weight(i,start_point(i):end_point(i)));
    lcm1_hi_vals = find(lcm1_avg(i,:) < lcm1_avg_avg(i) - 5);
    lcm1_lo_vals = find(lcm1_avg(i,:) > lcm1_avg_avg(i) + 5);
    lcm1_ramp_high_deviation(i) = std(lcm1_avg(i,lcm1_hi_vals),lcm1_weight(i,lcm1_hi_vals));
    lcm1_ramp_low_deviation(i) = std(lcm1_avg(i,lcm1_lo_vals),lcm1_weight(i,lcm1_lo_vals));
    lcm1_ramp_high(i) = sum(lcm1_avg_wt(i,lcm1_hi_vals))/sum(lcm1_weight_raw(i,lcm1_hi_vals));
    lcm1_ramp_low(i) = sum(lcm1_avg_wt(i,lcm1_lo_vals))/sum(lcm1_weight_raw(i,lcm1_lo_vals));
end

%ramp_deviation = 0.05;

vmon_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
vmon_weight_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);

vmon_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);
vmon_weight_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);

vmon_avg_trash_raw_pass = zeros(num_files,num_rows,1);
vmon_weight_avg_trash_raw_pass = zeros(num_files,num_rows,1);

time_ramp_up_pass = zeros(num_files,num_rows,1);
time_ramp_down_pass = zeros(num_files,num_rows,1);
time_trash_pass = zeros(num_files,num_rows,1);

lcm1_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
lcm1_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);
lcm1_avg_trash_raw_pass = zeros(num_files,num_rows,1);

lcm1_weight_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
lcm1_weight_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);
lcm1_weight_avg_trash_raw_pass = zeros(num_files,num_rows,1);



%Chunk array(file number, chunk number, chunk boundary element): For each 
%'chunk,' I define a chunk number (first, second, etc) and the elements 
%that mark the boundaries of each chunk. For example, the second chunk for
%the first file is given as up_chunk_array(1,1,2). The boundaries of the
%second chunk would be up_chunk_array(1,2,2) + 1 and up_chunk_array(1,2,3).

up_chunk_array_pass = zeros(num_files,2,num_rows); 
down_chunk_array_pass = zeros(num_files,2,num_rows);
trash_chunk_array_pass = zeros(num_files,2,num_rows);

up_array_count = zeros(num_files,1);
down_array_count = zeros(num_files,1);
trash_array_count = zeros(num_files,1);

%I start them at one out of necessity and subtract one after the sorting
num_ramp_up_points = ones(num_files,1);
num_ramp_down_points = ones(num_files,1);
num_trash_points = ones(num_files,1);

count_up_chunks = zeros(num_files,1);
count_down_chunks = zeros(num_files,1);
count_trash_chunks = zeros(num_files,1);

index_up_chunk_data_pass = zeros(num_files,1);
index_down_chunk_data_pass = zeros(num_files,1);
index_trash_data_pass = zeros(num_files,1);

if (power_supply == 0) || (power_supply == 1)
    for i = 1:num_files
        count_up_chunks(i) = 0;
        count_down_chunks(i) = 0;
        count_trash_chunks(i) = 0;
        
        for j = start_point(i):numpoints(i)           
            if abs(vmon_avg(i,j) - vmon_ramp_high(i)) < 3*vmon_ramp_high_deviation(i) && vmon_avg(i,j) < vmon_ramp_high(i) + vmon_ramp_high_deviation(i)
                index_up_chunk_data_pass(i,end+1) = j;
                num_ramp_up_points(i) = num_ramp_up_points(i) + 1;
                vmon_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = vmon_avg_raw(i,j);
                vmon_weight_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = vmon_weight_raw(i,j);
                lcm1_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_avg_raw(i,j);
                lcm1_weight_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_weight_raw(i,j);
                time_ramp_up_pass(i,num_ramp_up_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                up_array_count(i) = up_array_count(i)+1;
                if time_ramp_up_pass(i,num_ramp_up_points(i)) - time_ramp_up_pass(i,num_ramp_up_points(i)-1) > sampling_time                
                    count_up_chunks(i) = count_up_chunks(i) + 1;
                    up_chunk_array_pass(i,1, count_up_chunks(i)) = count_up_chunks(i);
                    up_chunk_array_pass(i,2, count_up_chunks(i)) = up_array_count(i) - 1;
                end
            elseif abs(vmon_avg(i,j) - vmon_ramp_low(i)) < 3*vmon_ramp_low_deviation(i) && vmon_avg(i,j) > vmon_ramp_low(i) - vmon_ramp_low_deviation(i)
                index_down_chunk_data_pass(i,end+1) = j;
                num_ramp_down_points(i) = num_ramp_down_points(i) + 1;
                vmon_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = vmon_avg_raw(i,j);
                vmon_weight_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = vmon_weight_raw(i,j);
                lcm1_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = lcm1_avg_raw(i,j);
                lcm1_weight_avg_ramp_down_raw_pass(i,num_ramp_down_points(i)) = lcm1_weight_raw(i,j);
                time_ramp_down_pass(i,num_ramp_down_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                down_array_count(i) = down_array_count(i) + 1;
                if time_ramp_down_pass(i,num_ramp_down_points(i)) - time_ramp_down_pass(i,num_ramp_down_points(i)-1) > sampling_time                
                    count_down_chunks(i) = count_down_chunks(i) + 1;
                    down_chunk_array_pass(i,1, count_down_chunks(i)) = count_down_chunks(i);
                    down_chunk_array_pass(i,2, count_down_chunks(i)) = down_array_count(i)-1;
                end
            else
                index_trash_data_pass(i,end+1) = j;
                num_trash_points(i) = num_trash_points(i) + 1;
                vmon_avg_trash_raw_pass(i,num_trash_points(i)) = vmon_avg_raw(i,j);
                vmon_weight_avg_trash_raw_pass(i,num_trash_points(i)) = vmon_weight_raw(i,j);
                lcm1_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_avg_raw(i,j);
                lcm1_weight_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_weight_raw(i,j);
                time_trash_pass(i,num_trash_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                trash_array_count(i) = trash_array_count(i) + 1;
                if time_trash_pass(i,num_trash_points(i)) - time_trash_pass(i,num_trash_points(i)-1) > sampling_time                
                    count_trash_chunks(i) = count_trash_chunks(i) + 1;
                    trash_chunk_array_pass(i,1, count_trash_chunks(i)) = count_trash_chunks(i);
                    trash_chunk_array_pass(i,2, count_trash_chunks(i)) = trash_array_count(i)-1;
                end
            end
        end
    end
elseif power_supply == 2
    for i = 1:num_files
        count_up_chunks(i) = 0;
        count_down_chunks(i) = 0;
        for j = start_point(i):numpoints(i)
            if (lcm1_avg(i,j) - lcm1_avg_offset(i) < +6*1e3) &&  (lcm1_avg(i,j) - lcm1_avg_offset(i) > -6*1e3)
                num_ramp_up_points(i) = num_ramp_up_points(i) + 1;
                vmon_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = vmon_avg_raw(i,j);
                vmon_weight_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = vmon_weight_raw(i,j);
                lcm1_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_avg_raw(i,j);
                lcm1_weight_avg_ramp_up_raw_pass(i,num_ramp_up_points(i)) = lcm1_weight_raw(i,j);
                time_ramp_up_pass(i,num_ramp_up_points(i)) = time(i,j); 
                up_array_count(i) = up_array_count(i)+1;
                if time_ramp_up_pass(i,num_ramp_up_points(i)) - time_ramp_up_pass(i,num_ramp_up_points(i)-1) > 5*sampling_time          
                    count_up_chunks(i) = count_up_chunks(i) + 1;
                    up_chunk_array_pass(i,1, count_up_chunks(i)) = count_up_chunks(i);
                    up_chunk_array_pass(i,2, count_up_chunks(i)) = up_array_count(i) - 1;
                end
            else
                index_trash_data_pass(i,end+1) = j;
                num_trash_points(i) = num_trash_points(i) + 1;
                lcm1_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_avg_raw(i,j);
                lcm1_weight_avg_trash_raw_pass(i,num_trash_points(i)) = lcm1_weight_raw(i,j);
                time_trash_pass(i,num_trash_points(i)) = time(i,j); % removed time(i,start_point(i)) subtraction
                trash_array_count(i) = trash_array_count(i) + 1;
                if time_trash_pass(i,num_trash_points(i)) - time_trash_pass(i,num_trash_points(i)-1) > 5*sampling_time                
                    count_trash_chunks(i) = count_trash_chunks(i) + 1;
                    trash_chunk_array_pass(i,1, count_trash_chunks(i)) = count_trash_chunks(i);
                    trash_chunk_array_pass(i,2, count_trash_chunks(i)) = trash_array_count(i)-1;
                end
            end
        end
    end
end

for i = 1:num_files
    num_ramp_up_points(i) = num_ramp_up_points(i) - 1;
    num_ramp_down_points(i) = num_ramp_down_points(i) -1;
    num_trash_points(i) = num_trash_points(i) -1;
end    

num_up_chunks = zeros(num_files,1);
num_down_chunks = zeros(num_files,1);
num_trash_chunks = zeros(num_files,1);

for i = 1:num_files
    %count to third to last trash chunk
    num_trash_chunks(i) = uint8(count_trash_chunks(i) -4);
    
    %count to second to last down chunk
    num_down_chunks(i) = uint8(count_down_chunks(i) - 2);

    %count up to 2nd to last up chunk
    num_up_chunks(i) = uint8(count_up_chunks(i) - 2);

    
    %this isn't counting correctly and I don't know why.
    if EDM_sim == 1
        
        truth = logical(num_trash_chunks(i) == (num_down_chunks(i) + num_up_chunks(i)));
        
        if truth == 1
            disp('num trash chunks = num_down_chunks + num_up_chunks... yay!')
        elseif truth == 0
            chunk_warn = ['num_trash_chunks != num_down_chunks + num_up_chunks...'...
                ' check your chunks, something may be wrong with chunk'...
                ' identifier algorithm. Inclusive chunk analysis may not work correctly.'];
            disp(chunk_warn)
        end                
    end
  
end

    
up_chunk_array = zeros(num_files,2,num_up_chunks);
down_chunk_array = zeros(num_files,2,num_down_chunks);
trash_chunk_array = zeros(num_files,2,num_trash_chunks);

time_ramp_up = zeros(num_files,num_ramp_up_points,1);
time_ramp_down = zeros(num_files,num_ramp_down_points,1);
time_trash = zeros(num_files,num_trash_points,1);

lcm1_avg_ramp_up_raw = zeros(num_files,num_ramp_up_points,1);
lcm1_avg_ramp_up = zeros(num_files,num_ramp_up_points,1);
lcm1_avg_ramp_up_raw_wt = zeros(num_files,num_ramp_up_points,1);

lcm1_avg_ramp_down_raw = zeros(num_files,num_ramp_down_points,1);
lcm1_avg_ramp_down = zeros(num_files,num_ramp_down_points,1);
lcm1_avg_ramp_down_raw_wt = zeros(num_files,num_ramp_down_points,1);

lcm1_avg_trash_raw = zeros(num_files,num_trash_points,1);
lcm1_avg_trash = zeros(num_files,num_trash_points,1);
lcm1_avg_trash_raw_wt = zeros(num_files,num_trash_points,1);

vmon_weight_avg_ramp_up_raw = zeros(num_files,num_ramp_up_points,1);
vmon_weight_avg_ramp_up = zeros(num_files,num_ramp_up_points,1);

vmon_weight_avg_ramp_down_raw = zeros(num_files,num_ramp_down_points,1);
vmon_weight_avg_ramp_down = zeros(num_files,num_ramp_down_points,1);

vmon_weight_avg_trash_raw = zeros(num_files,num_trash_points,1);
vmon_weight_avg_trash = zeros(num_files,num_trash_points,1);

vmon_avg_ramp_up_raw = zeros(num_files,num_ramp_up_points,1);
vmon_avg_ramp_up = zeros(num_files,num_ramp_up_points,1);

vmon_avg_ramp_down_raw = zeros(num_files,num_ramp_down_points,1);
vmon_avg_ramp_down = zeros(num_files,num_ramp_down_points,1);

vmon_avg_trash_raw = zeros(num_files,num_trash_points,1);
vmon_avg_trash = zeros(num_files,num_trash_points,1);

lcm1_avg_charge_neg_pass = zeros(num_files,1);
lcm1_stdev_charge_neg_pass = zeros (num_files,1);
lcm1_weight_charge_neg_raw_pass = zeros (num_files,1);
lcm1_avg_charge_neg_raw_pass = zeros(num_files,1);
lcm1_stdev_charge_neg_raw_pass = zeros (num_files,1);
lcm1_charge_neg_pass_time = zeros(num_files,1);

lcm1_avg_charge_pos_pass = zeros(num_files,1);
lcm1_stdev_charge_pos_pass = zeros (num_files,1);
lcm1_weight_charge_pos_raw_pass = zeros (num_files,1);
lcm1_avg_charge_pos_raw_pass = zeros(num_files,1);
lcm1_stdev_charge_pos_raw_pass = zeros (num_files,1);
lcm1_charge_pos_pass_time = zeros(num_files,1);

lcm1_avg_discharge_pos_pass = zeros(num_files,1);
lcm1_stdev_discharge_pos_pass = zeros (num_files,1);
lcm1_avg_discharge_pos_raw_pass = zeros(num_files,1);
lcm1_stdev_discharge_pos_raw_pass = zeros (num_files,1);
lcm1_weight_discharge_pos_raw_pass = zeros (num_files,1);
lcm1_discharge_pos_pass_time = zeros(num_files,1);

lcm1_avg_discharge_neg_pass = zeros(num_files,1);
lcm1_stdev_discharge_neg_pass = zeros (num_files,1);
lcm1_weight_discharge_neg_raw_pass = zeros (num_files,1);
lcm1_avg_discharge_neg_raw_pass = zeros(num_files,1);
lcm1_stdev_discharge_neg_raw_pass = zeros (num_files,1);
lcm1_discharge_neg_pass_time = zeros(num_files,1);

% 3/11/2018 inv_weight is just stdev. I fucked up variable names so I have to use
%this nomenclature to avoid a lot of work I don't want to do right now.

lcm1_weight_avg_ramp_up_raw = zeros(num_files,num_ramp_up_points,1);
lcm1_inv_weight_avg_ramp_up = zeros(num_files,num_ramp_up_points,1);

lcm1_weight_avg_ramp_down_raw = zeros(num_files,num_ramp_down_points,1);
lcm1_inv_weight_avg_ramp_down = zeros(num_files,num_ramp_down_points,1);

lcm1_weight_avg_trash_raw = zeros(num_files,num_trash_points,1);
lcm1_inv_weight_avg_trash = zeros(num_files,num_trash_points,1);
lcm1_inv_weight_avg_trash_raw = zeros(num_files,num_trash_points,1);
    
% radius of points to sample for steady-state average and stdev
lcm1_avg_trash_steady_avg = zeros(num_files,num_trash_chunks);
lcm1_avg_trash_steady_stdev = zeros(num_files,num_trash_chunks);
steady_state_range = 10;
%track index of discharge currents and charging currents in trash data.
%column 1: positive/negative (dis)charge
%column 2: end point of previous chunk
%column 3: end point of this chunk
discharge_index_pass = zeros(num_files,num_trash_chunks,3);
charge_index_pass = zeros(num_files,num_trash_chunks,3);



for i = 1:num_files
    for j = 1:num_ramp_up_points(i)
        lcm1_avg_ramp_up_raw(i,j) = lcm1_avg_ramp_up_raw_pass(i,j+1);
        lcm1_weight_avg_ramp_up_raw(i,j) = lcm1_weight_avg_ramp_up_raw_pass(i,j+1);
        vmon_avg_ramp_up_raw(i,j) = vmon_avg_ramp_up_raw_pass(i,j+1);
        vmon_weight_avg_ramp_up_raw(i,j) = vmon_weight_avg_ramp_up_raw_pass(i,j+1);
        time_ramp_up(i,j) = time_ramp_up_pass(i,j+1);
        
        lcm1_avg_ramp_up(i,j) = lcm1_avg_ramp_up_raw(i,j)*lcm1_avg_scale;
        %lcm1_weight_avg_ramp_up(i,j) = lcm1_weight_avg_ramp_up_raw(i,j)*lcm1_avg_scale;
        lcm1_inv_weight_avg_ramp_up(i,j) = (lcm1_weight_avg_ramp_up_raw(i,j))^(-1/2)*lcm1_avg_scale;
        vmon_avg_ramp_up(i,j) = vmon_avg_ramp_up_raw(i,j)*vmon_avg_scale;
        vmon_weight_avg_ramp_up(i,j) = vmon_weight_avg_ramp_up_raw(i,j)*vmon_avg_scale;
    end

    for j = 1:num_ramp_down_points(i)
        lcm1_avg_ramp_down_raw(i,j) = lcm1_avg_ramp_down_raw_pass(i,j+1);
        lcm1_weight_avg_ramp_down_raw(i,j) = lcm1_weight_avg_ramp_down_raw_pass(i,j+1);
        vmon_avg_ramp_down_raw(i,j) = vmon_avg_ramp_down_raw_pass(i,j+1);
        vmon_weight_avg_ramp_down_raw(i,j) = vmon_weight_avg_ramp_down_raw_pass(i,j+1);
        time_ramp_down(i,j) = time_ramp_down_pass(i,j+1);
        
        lcm1_avg_ramp_down(i,j) = lcm1_avg_ramp_down_raw(i,j)*lcm1_avg_scale;
        %lcm1_weight_avg_ramp_down(i,j) = lcm1_weight_avg_ramp_down_raw(i,j)*lcm1_avg_scale;
        lcm1_inv_weight_avg_ramp_down(i,j) = (lcm1_weight_avg_ramp_down_raw(i,j))^(-1/2)*lcm1_avg_scale;
        vmon_avg_ramp_down(i,j) = vmon_avg_ramp_down_raw(i,j)*vmon_avg_scale;
        vmon_weight_avg_ramp_down(i,j) = vmon_weight_avg_ramp_down_raw(i,j)*vmon_avg_scale;
    end

    for j = 1:num_trash_points(i)
        lcm1_avg_trash_raw(i,j) = lcm1_avg_trash_raw_pass(i,j+1);
        lcm1_weight_avg_trash_raw(i,j) = lcm1_weight_avg_trash_raw_pass(i,j+1);
        vmon_avg_trash_raw(i,j) = vmon_avg_trash_raw_pass(i,j+1);
        vmon_weight_avg_trash_raw(i,j) = vmon_weight_avg_trash_raw_pass(i,j+1);
        time_trash(i,j) = time_trash_pass(i,j+1);
        
        lcm1_avg_trash(i,j) = lcm1_avg_trash_raw(i,j)*lcm1_avg_scale;
        lcm1_inv_weight_avg_trash_raw(i,j) = (lcm1_weight_avg_trash_raw(i,j))^(-1/2);
        lcm1_inv_weight_avg_trash(i,j) = lcm1_inv_weight_avg_trash_raw(i,j)*abs(lcm1_avg_scale);
        vmon_avg_trash(i,j) = vmon_avg_trash_raw(i,j)*vmon_avg_scale;
        vmon_weight_avg_trash(i,j) = vmon_weight_avg_trash_raw(i,j)*vmon_avg_scale;
        
    end
    
end

for i =1:num_files
    
    %start at third trash chunk
    for j =1:num_trash_chunks(i) + 1
        trash_chunk_array(i,1,j) = trash_chunk_array_pass(i,2,j+1);
        trash_chunk_array(i,2,j) = trash_chunk_array_pass(i,2,j+2);
    end
    
    %start at 2nd down chunk
    for j =1:num_down_chunks(i) + 1
%        if time_ramp_down(i,down_chunk_array_pass(i,2,j+1)) > time_trash(i,trash_chunk_array(1,2,1))
        down_chunk_array(i,1,j) = down_chunk_array_pass(i,2,j);
        down_chunk_array(i,2,j) = down_chunk_array_pass(i,2,j+1);
%        end
    end
    
    %start at 2nd up chunk
    for j =1:num_up_chunks(i) + 1
        up_chunk_array(i,1,j) = up_chunk_array_pass(i,2,j);
        up_chunk_array(i,2,j) = up_chunk_array_pass(i,2,j+1);
    end
    
    for j = 1:num_trash_chunks(i)
        chunk_begin = trash_chunk_array(i,2,j);
        chunk_end = trash_chunk_array(i,2,j+1);
        middle = round((chunk_end - chunk_begin)/2,0);
        middle_lo = chunk_begin + middle - steady_state_range;
        middle_hi = chunk_begin + middle + steady_state_range;
        lcm1_avg_trash_steady_avg(i,j) = mean(lcm1_avg_trash(i,middle_lo:middle_hi));
        lcm1_avg_trash_steady_stdev(i,j) = mean(lcm1_inv_weight_avg_trash(i,middle_lo:middle_hi));
    end

    %pick out the ramping currents. +V -> 0, 0 -> -V, -V -> 0, 0 -> +V.
    max_length_discharge = 0;
    for j = 1:num_trash_chunks(i)
        %first point might be the start of the ramp, which could be a small
        %stdev. Increment by 1 so we know we're starting on a ramp point
        chunk_begin = trash_chunk_array(i,2,j)+1;
        chunk_end = trash_chunk_array(i,2,j+1);
        middle = round((chunk_end - chunk_begin)/2,0);
        for k = chunk_begin:chunk_begin+ middle
            if lcm1_inv_weight_avg_trash(i,k+1) < 1.5 * lcm1_avg_trash_steady_stdev(i,j) && lcm1_inv_weight_avg_trash(i,k+1) > lcm1_inv_weight_avg_trash(i,k+2)                
                if max(lcm1_avg_trash(i,chunk_begin:k)) < lcm1_avg_trash_steady_avg(i,j) && abs(max(lcm1_avg_trash(i,chunk_begin:k))) > 3 * abs(lcm1_avg_trash_steady_avg(i,j))
                    %discharge from +V
                    discharge_index_pass(i,j,1) = 1;
                elseif max(lcm1_avg_trash(i,chunk_begin:k)) > lcm1_avg_trash_steady_avg(i,j) && abs(max(lcm1_avg_trash(i,chunk_begin:k))) > 3 * abs(lcm1_avg_trash_steady_avg(i,j))
                    %discharge from -V
                    discharge_index_pass(i,j,1) = -1;
                end

                discharge_index_pass(i,j,2) = chunk_begin;
                discharge_index_pass(i,j,3) = k+2;
                
                break            
            end                    
        end
    % if the condition for the discharging current for a segment is not 
    %met, I simply define the segment to start from chunk_begin and have
    %thesame # of points as the previous segment
    if discharge_index_pass(i,j,1) == 0
%            disp(j)
            discharge_index_pass(i,j,1) = - discharge_index_pass(i,j-1,1);
            discharge_index_pass(i,j,2) = chunk_begin;
            discharge_index_pass(i,j,3) = chunk_begin + discharge_index_pass(i,j-1,3) - discharge_index_pass(i,j-1,2);                     
    end
    
    if max_length_discharge < discharge_index_pass(i,j,3) - discharge_index_pass(i,j,2)
        max_length_discharge = discharge_index_pass(i,j,3) - discharge_index_pass(i,j,2);
    end
    
    end
    
    max_length_charge = 0;
    for j = 1:num_trash_chunks
        chunk_begin = trash_chunk_array(i,2,j)+1;
        chunk_end = trash_chunk_array(i,2,j+1);
        middle = round((chunk_end - chunk_begin)/2,0);
        for k = chunk_begin+ middle:chunk_end                       
            if lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k) < 1.5 * lcm1_avg_trash_steady_stdev(i,j) && lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k) > lcm1_inv_weight_avg_trash(i,chunk_end + chunk_begin + middle - k - 1)
%                 if lcm1_avg_trash(i,chunk_end + middle - k) > lcm1_avg_trash(i,chunk_end + middle - k - 1)
                if max(lcm1_avg_trash(i,k:chunk_end)) >  lcm1_avg_trash_steady_avg(i,j) && abs(max(lcm1_avg_trash(i,k:chunk_end))) > 3*abs(lcm1_avg_trash_steady_avg(i,j))
                    % charging to +V
                    charge_index_pass(i,j,1) = 1;
%                 elseif lcm1_avg_trash(i,chunk_end + middle - k) < lcm1_avg_trash(i,chunk_end + middle - k - 1)
                elseif max(lcm1_avg_trash(i,k:chunk_end)) <  lcm1_avg_trash_steady_avg(i,j)   && abs(max(lcm1_avg_trash(i,k:chunk_end))) > 3*abs(lcm1_avg_trash_steady_avg(i,j))                
                    % charging to -V
                    charge_index_pass(i,j,1) = -1;
                end

                charge_index_pass(i,j,2) = chunk_end + chunk_begin + middle - k -1;
                charge_index_pass(i,j,3) = chunk_end;                
                
                break                                        
            end
        end
    % if the condition for detecting the charging current is not satisfied
    % for a chunk, I simplly define the charging segment to be the same #
    % of points as the previous segment
        if charge_index_pass(i,j,1) == 0
    %        disp(j)
            charge_index_pass(i,j,1) = -charge_index_pass(i,j-1,1);
            charge_index_pass(i,j,2) = chunk_end - (charge_index_pass(i,j-1,3) - charge_index_pass(i,j-1,2));
            charge_index_pass(i,j,3) = chunk_end;
        end

        if max_length_charge < charge_index_pass(i,j,3) - charge_index_pass(i,j,2)
            max_length_charge = charge_index_pass(i,j,3) - charge_index_pass(i,j,2);
%             disp(j) 
%             disp(max_length_charge)
        end
        
    end
    
    %want to make the charge/discharge length long enough so that the last
    %points are almost steady-state
    max_length_charge = max_length_charge+round(0.09*max_length_charge+1,0);
    max_length_discharge = max_length_discharge+round(0.09*max_length_discharge+1,0);
    
    %make all segments the same length which is determined by the maximum
    %segment size
    for j = 1:num_trash_chunks
        
        charge_index_pass(i,j,2) = charge_index_pass(i,j,3) - max_length_charge;
        discharge_index_pass(i,j,3) = discharge_index_pass(i,j,2) + max_length_discharge;
        
    end
end

%We want all but the first discharge_index and all but the last
%charge_index
charge_index = charge_index_pass(:,1:end-1,:);
discharge_index = discharge_index_pass(:,2:end,:);

charge_neg_index_count = 0;
charge_neg_index = zeros(num_files,num_trash_chunks);

discharge_neg_index_count = 0;
discharge_neg_index = zeros(num_files,num_trash_chunks);

charge_pos_index_count = 0;
charge_pos_index = zeros(num_files,num_trash_chunks);

discharge_pos_index_count = 0;
discharge_pos_index = zeros(num_files,num_trash_chunks);


%create vectors containing time, leakage avg, and leakage stdev data 
%for + V -> 0, -V -> 0, 0 -> +V, and 0 -> -V
for i = 1:num_files
    for j = 1:num_trash_chunks -1
        if discharge_index(i,j,1) == 1
%             lcm1_discharge_pos_pass(i,2,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = lcm1_avg_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));
%             lcm1_discharge_pos_pass(i,3,end - (discharge_index(i,j,3)-discharge_index(i,j,2)):end) = lcm1_inv_weight_avg_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));
%             lcm1_discharge_pos_pass(i,1,end - (discharge_index(i,j,3)-discharge_index(i,j,2)):end) = time_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));
            discharge_pos_index_count = discharge_pos_index_count+1;
            discharge_pos_index(i,discharge_pos_index_count+1) = discharge_pos_index(i,discharge_pos_index_count) + discharge_index(i,j,3) - discharge_index(i,j,2)+1;            
            
            lcm1_avg_discharge_pos_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = lcm1_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));
            lcm1_stdev_discharge_pos_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) =lcm1_inv_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));
            lcm1_weight_discharge_pos_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) =lcm1_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));
            lcm1_discharge_pos_pass_time(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = time_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));
        
        elseif discharge_index(i,j,1) == -1
%             lcm1_discharge_neg_pass(i,2,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = lcm1_avg_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));
%             lcm1_discharge_neg_pass(i,3,end - (discharge_index(i,j,3)-discharge_index(i,j,2)):end) = lcm1_inv_weight_avg_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));
%             lcm1_discharge_neg_pass(i,1,end - (discharge_index(i,j,3)-discharge_index(i,j,2)):end) = time_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));
            discharge_neg_index_count = discharge_neg_index_count+1;
            discharge_neg_index(i,discharge_neg_index_count+1) = discharge_neg_index(i,discharge_neg_index_count) + discharge_index(i,j,3) - discharge_index(i,j,2)+1;
            lcm1_avg_discharge_neg_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = lcm1_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));
            lcm1_stdev_discharge_neg_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) =lcm1_inv_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));
            lcm1_weight_discharge_neg_raw_pass(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) =lcm1_weight_avg_trash_raw(i,discharge_index(i,j,2):discharge_index(i,j,3));
            lcm1_discharge_neg_pass_time(i,end+1:end+1+discharge_index(i,j,3)-discharge_index(i,j,2)) = time_trash(i,discharge_index(i,j,2):discharge_index(i,j,3));

        end
        
        if charge_index(i,j,1) == 1
%             lcm1_charge_pos_pass(i,2,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = lcm1_avg_trash(i,charge_index(i,j,2):charge_index(i,j,3));
%             lcm1_charge_pos_pass(i,3,end - (charge_index(i,j,3)-charge_index(i,j,2)):end) = lcm1_inv_weight_avg_trash(i,charge_index(i,j,2):charge_index(i,j,3));
%             lcm1_charge_pos_pass(i,1,end - (charge_index(i,j,3)-charge_index(i,j,2)):end) = time_trash(i,charge_index(i,j,2):charge_index(i,j,3));
            charge_pos_index_count = charge_pos_index_count+1; 
            lcm1_avg_charge_pos_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = lcm1_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));
            lcm1_weight_charge_pos_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) =lcm1_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));
            lcm1_stdev_charge_pos_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) =lcm1_inv_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));
            lcm1_charge_pos_pass_time(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = time_trash(i,charge_index(i,j,2):charge_index(i,j,3));
            charge_pos_index(i,charge_pos_index_count+1) = charge_pos_index(i,charge_pos_index_count) + charge_index(i,j,3) - charge_index(i,j,2)+1;
            
        elseif charge_index(i,j,1) == -1
            charge_neg_index_count = charge_neg_index_count+1;
            lcm1_avg_charge_neg_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = lcm1_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));
            lcm1_stdev_charge_neg_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) =lcm1_inv_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));
            lcm1_weight_charge_neg_raw_pass(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) =lcm1_weight_avg_trash_raw(i,charge_index(i,j,2):charge_index(i,j,3));
            lcm1_charge_neg_pass_time(i,end+1:end+1+charge_index(i,j,3)-charge_index(i,j,2)) = time_trash(i,charge_index(i,j,2):charge_index(i,j,3));
            charge_neg_index(i,charge_neg_index_count+1) = charge_neg_index(i,charge_neg_index_count) + charge_index(i,j,3) - charge_index(i,j,2)+1;
        end
    end
end

%get rid of those pesky zeros in front
% lcm1_discharge_pos = lcm1_discharge_pos_pass(:,:,2:end);
% lcm1_discharge_neg = lcm1_discharge_neg_pass(:,:,2:end);
% lcm1_charge_pos = lcm1_charge_pos_pass(:,:,2:end);

lcm1_avg_charge_neg_raw = lcm1_avg_charge_neg_raw_pass(:,2:end);
lcm1_stdev_charge_neg_raw = lcm1_stdev_charge_neg_raw_pass(:,2:end);
lcm1_weight_charge_neg_raw = lcm1_weight_charge_neg_raw_pass(:,2:end);
lcm1_charge_neg_time = lcm1_charge_neg_pass_time(:,2:end);

lcm1_avg_charge_pos_raw = lcm1_avg_charge_pos_raw_pass(:,2:end);
lcm1_stdev_charge_pos_raw = lcm1_stdev_charge_pos_raw_pass(:,2:end);
lcm1_weight_charge_pos_raw = lcm1_weight_charge_pos_raw_pass(:,2:end);
lcm1_charge_pos_time = lcm1_charge_pos_pass_time(:,2:end);

lcm1_avg_discharge_pos_raw = lcm1_avg_discharge_pos_raw_pass(:,2:end);
lcm1_stdev_discharge_pos_raw = lcm1_stdev_discharge_pos_raw_pass(:,2:end);
lcm1_weight_discharge_pos_raw = lcm1_weight_discharge_pos_raw_pass(:,2:end);
lcm1_discharge_pos_time = lcm1_discharge_pos_pass_time(:,2:end);

lcm1_avg_discharge_neg_raw = lcm1_avg_discharge_neg_raw_pass(:,2:end);
lcm1_stdev_discharge_neg_raw = lcm1_stdev_discharge_neg_raw_pass(:,2:end);
lcm1_weight_discharge_neg_raw = lcm1_weight_discharge_neg_raw_pass(:,2:end);
lcm1_discharge_neg_time = lcm1_discharge_neg_pass_time(:,2:end);


lcm1_avg_charge_neg = lcm1_avg_charge_neg_raw*lcm1_avg_scale;
lcm1_stdev_charge_neg = lcm1_stdev_charge_neg_raw * abs(lcm1_avg_scale);


lcm1_avg_charge_pos = lcm1_avg_charge_pos_raw*lcm1_avg_scale;
lcm1_stdev_charge_pos = lcm1_stdev_charge_pos_raw*abs(lcm1_avg_scale);


lcm1_avg_discharge_pos = lcm1_avg_discharge_pos_raw* lcm1_avg_scale;
lcm1_stdev_discharge_pos = lcm1_stdev_discharge_pos_raw * abs(lcm1_avg_scale);


lcm1_avg_discharge_neg = lcm1_avg_discharge_neg_raw * lcm1_avg_scale;
lcm1_stdev_discharge_neg = lcm1_stdev_discharge_neg_raw * abs(lcm1_avg_scale);



lcm1_avg_ramp_down_inc_raw_pass = zeros(num_files,1);
lcm1_weight_ramp_down_inc_raw_pass = zeros(num_files,1);
time_ramp_down_inc_pass = zeros(num_files,1);
down_chunk_inc_array_pass = zeros(num_files,2,num_rows);

lcm1_avg_discharge_pos_sum_current_raw = zeros(num_files,max_length_discharge);
lcm1_avg_discharge_pos_sum_weight_raw =  zeros(num_files,max_length_discharge);
lcm1_avg_charge_pos_sum_current_raw = zeros(num_files,max_length_charge);

lcm1_avg_ramp_up_inc_raw_pass = zeros(num_files,1);
lcm1_weight_ramp_up_inc_raw_pass = zeros(num_files,1);
time_ramp_up_inc_pass = zeros(num_files,1);
up_chunk_inc_array_pass = zeros(num_files,2,num_rows);

lcm1_avg_discharge_neg_sum_current_raw = zeros(num_files,max_length_discharge);
lcm1_avg_discharge_neg_sum_weight_raw = zeros(num_files,max_length_discharge);
lcm1_avg_charge_neg_sum_current_raw = zeros(num_files,max_length_charge);

%This builds the inclusive data sets. charging and discharging segments are
%attached to either end of the appropriate ramp segments.
if inclusive_data == 1
     for i = 1:num_files


       for j = 1:num_down_chunks(i)
           down_chunk_inc_array_pass(i,1,j) = length(lcm1_avg_ramp_down_inc_raw_pass(i,:));
           lcm1_avg_ramp_down_inc_raw_pass(i,end+1:end+1+charge_pos_index(i,j+1) - (charge_pos_index(i,j)+1) ) = lcm1_avg_charge_pos_raw(i,charge_pos_index(i,j)+1:charge_pos_index(i,j+1));
           lcm1_weight_ramp_down_inc_raw_pass(i,end+1:end+1+charge_pos_index(i,j+1) - (charge_pos_index(i,j)+1) ) = lcm1_weight_charge_pos_raw(i,charge_pos_index(i,j)+1:charge_pos_index(i,j+1));
           time_ramp_down_inc_pass(i,end+1:end+1+charge_pos_index(i,j+1) - (charge_pos_index(i,j)+1) ) = lcm1_charge_pos_time(i,charge_pos_index(i,j)+1:charge_pos_index(i,j+1));

           lcm1_avg_ramp_down_inc_raw_pass(i,end+1:end+1+down_chunk_array(i,2,j+1) - (down_chunk_array(i,2,j)+1) ) = lcm1_avg_ramp_down_raw(i,down_chunk_array(i,2,j)+1:down_chunk_array(i,2,j+1));
           lcm1_weight_ramp_down_inc_raw_pass(i,end+1:end+1+down_chunk_array(i,2,j+1) - (down_chunk_array(i,2,j)+1) ) = lcm1_weight_avg_ramp_down_raw(i,down_chunk_array(i,2,j)+1:down_chunk_array(i,2,j+1));           
           time_ramp_down_inc_pass(i,end+1:end+1+down_chunk_array(i,2,j+1) - (down_chunk_array(i,2,j)+1) ) = time_ramp_down(i,down_chunk_array(i,2,j)+1:down_chunk_array(i,2,j+1));

           lcm1_avg_ramp_down_inc_raw_pass(i,end+1:end+1+discharge_pos_index(i,j+1) - (discharge_pos_index(i,j)+1) ) =    lcm1_avg_discharge_pos_raw(i,discharge_pos_index(i,j)+1:discharge_pos_index(i,j+1));
           lcm1_weight_ramp_down_inc_raw_pass(i,end+1:end+1+discharge_pos_index(i,j+1) - (discharge_pos_index(i,j)+1) ) =    lcm1_weight_discharge_pos_raw(i,discharge_pos_index(i,j)+1:discharge_pos_index(i,j+1));
           time_ramp_down_inc_pass(i,end+1:end+1+discharge_pos_index(i,j+1) - (discharge_pos_index(i,j)+1) ) =    lcm1_discharge_pos_time(i,discharge_pos_index(i,j)+1:discharge_pos_index(i,j+1));
           down_chunk_inc_array_pass(i,2,j) = length(lcm1_avg_ramp_down_inc_raw_pass(i,:));
           
           for k = 1:max_length_discharge
               lcm1_avg_discharge_pos_sum_current_raw(i,k) = lcm1_avg_discharge_pos_sum_current_raw(i,k)+  lcm1_avg_discharge_pos_raw(i,discharge_pos_index(i,j)+1+k);
               lcm1_avg_discharge_pos_sum_weight_raw(i,k) = lcm1_avg_discharge_pos_sum_weight_raw(i,k)+ lcm1_weight_discharge_pos_raw(i,discharge_pos_index(i,j)+1+k);
           end
           
           for k = 1:max_length_charge
               lcm1_avg_charge_pos_sum_current_raw(i,k) = lcm1_avg_charge_pos_sum_current_raw(i,k)+ lcm1_avg_charge_pos_raw(i,charge_pos_index(i,j)+1+k);               
           end
       end

       for j = 1:num_up_chunks(i)
           up_chunk_inc_array_pass(i,1,j) = length(lcm1_avg_ramp_up_inc_raw_pass(i,:));
           lcm1_avg_ramp_up_inc_raw_pass(i,end+1:end+1+charge_neg_index(i,j+1) - (charge_neg_index(i,j)+1) ) = lcm1_avg_charge_neg_raw(i,charge_neg_index(i,j)+1:charge_neg_index(i,j+1));
           lcm1_weight_ramp_up_inc_raw_pass(i,end+1:end+1+charge_neg_index(i,j+1) - (charge_neg_index(i,j)+1) ) = lcm1_weight_charge_neg_raw(i,charge_neg_index(i,j)+1:charge_neg_index(i,j+1));
           time_ramp_up_inc_pass(i,end+1:end+1+charge_neg_index(i,j+1) - (charge_neg_index(i,j)+1) ) = lcm1_charge_neg_time(i,charge_neg_index(i,j)+1:charge_neg_index(i,j+1));

           lcm1_avg_ramp_up_inc_raw_pass(i,end+1:end+1+up_chunk_array(i,2,j+1) - (up_chunk_array(i,2,j)+1) ) = lcm1_avg_ramp_up_raw(i,up_chunk_array(i,2,j)+1:up_chunk_array(i,2,j+1));
           lcm1_weight_ramp_up_inc_raw_pass(i,end+1:end+1+up_chunk_array(i,2,j+1) - (up_chunk_array(i,2,j)+1) ) = lcm1_weight_avg_ramp_up_raw(i,up_chunk_array(i,2,j)+1:up_chunk_array(i,2,j+1));           
           time_ramp_up_inc_pass(i,end+1:end+1+up_chunk_array(i,2,j+1) - (up_chunk_array(i,2,j)+1) ) = time_ramp_up(i,up_chunk_array(i,2,j)+1:up_chunk_array(i,2,j+1));

           lcm1_avg_ramp_up_inc_raw_pass(i,end+1:end+1+discharge_neg_index(i,j+1) - (discharge_neg_index(i,j)+1) ) =    lcm1_avg_discharge_neg_raw(i,discharge_neg_index(i,j)+1:discharge_neg_index(i,j+1));
           lcm1_weight_ramp_up_inc_raw_pass(i,end+1:end+1+discharge_neg_index(i,j+1) - (discharge_neg_index(i,j)+1) ) =    lcm1_weight_discharge_neg_raw(i,discharge_neg_index(i,j)+1:discharge_neg_index(i,j+1));
           time_ramp_up_inc_pass(i,end+1:end+1+discharge_neg_index(i,j+1) - (discharge_neg_index(i,j)+1) ) =    lcm1_discharge_neg_time(i,discharge_neg_index(i,j)+1:discharge_neg_index(i,j+1));
           up_chunk_inc_array_pass(i,2,j) = length(lcm1_avg_ramp_up_inc_raw_pass(i,:));
           
           for k = 1:max_length_discharge
               lcm1_avg_discharge_neg_sum_current_raw(i,k) = lcm1_avg_discharge_neg_sum_current_raw(i,k)+ lcm1_avg_discharge_neg_raw(i,discharge_neg_index(i,j)+1+k);
               lcm1_avg_discharge_neg_sum_weight_raw(i,k) = lcm1_avg_discharge_neg_sum_weight_raw(i,k)+ lcm1_weight_discharge_neg_raw(i,discharge_neg_index(i,j)+1+k);
           end
           
           for k = 1:max_length_charge
               lcm1_avg_charge_neg_sum_current_raw(i,k) = lcm1_avg_charge_neg_sum_current_raw(i,k)+ lcm1_avg_charge_neg_raw(i,charge_neg_index(i,j)+1+k);               
           end
           
       end
    end    
end    


lcm1_avg_ramp_down_inc_raw = lcm1_avg_ramp_down_inc_raw_pass(:,2:end);
lcm1_weight_ramp_down_inc_raw = lcm1_weight_ramp_down_inc_raw_pass(:,2:end);

time_ramp_down_inc = time_ramp_down_inc_pass(:,2:end);
lcm1_avg_ramp_down_inc = lcm1_avg_ramp_down_inc_raw * lcm1_avg_scale;

down_chunk_inc_array = zeros(num_files,2,max(num_down_chunks(:)));
num_ramp_down_inc_points = zeros(num_files,1);

lcm1_avg_ramp_up_inc_raw = lcm1_avg_ramp_up_inc_raw_pass(:,2:end);
lcm1_weight_ramp_up_inc_raw = lcm1_weight_ramp_up_inc_raw_pass(:,2:end);

time_ramp_up_inc = time_ramp_up_inc_pass(:,2:end);
lcm1_avg_ramp_up_inc = lcm1_avg_ramp_up_inc_raw * lcm1_avg_scale;

up_chunk_inc_array = zeros(num_files,2,max(num_up_chunks(:)));
num_ramp_up_inc_points = zeros(num_files,1);

lcm1_avg_charge_neg_norm = zeros(num_files,max_length_charge);
lcm1_avg_charge_pos_norm = zeros(num_files,max_length_charge);
lcm1_avg_discharge_neg_norm = zeros(num_files,max_length_discharge);
lcm1_avg_discharge_pos_norm = zeros(num_files,max_length_discharge);

for i = 1:num_files
    for j = 1:max_length_charge
        lcm1_avg_charge_neg_norm(i,j) = (lcm1_avg_charge_neg_sum_current_raw(i,j)-lcm1_avg_offset_raw(i))/max_length_charge*lcm1_avg_scale;
        lcm1_avg_charge_pos_norm(i,j) = (lcm1_avg_charge_pos_sum_current_raw(i,j)-lcm1_avg_offset_raw(i))/max_length_charge*lcm1_avg_scale;
    end
    
    for j = 1:max_length_discharge
        lcm1_avg_discharge_neg_norm(i,j) = (lcm1_avg_discharge_neg_sum_current_raw(i,j)-lcm1_avg_offset_raw(i))/max_length_discharge *  lcm1_avg_scale;
        lcm1_avg_discharge_pos_norm(i,j) = (lcm1_avg_discharge_pos_sum_current_raw(i,j)-lcm1_avg_offset_raw(i))/max_length_discharge * lcm1_avg_scale;
    end
end
    
for i = 1: num_files
    num_ramp_down_inc_points = length(lcm1_avg_ramp_down_inc(i,:));
    for j = 1:num_down_chunks(i)
        
        down_chunk_inc_array(i,1,j) = down_chunk_inc_array_pass(i,1,j) - 1;
        down_chunk_inc_array(i,2,j) = down_chunk_inc_array_pass(i,2,j) - 1;
        
    end
    
    num_ramp_up_inc_points = length(lcm1_avg_ramp_up_inc(i,:));
    for j = 1:num_up_chunks(i)
        
        up_chunk_inc_array(i,1,j) = up_chunk_inc_array_pass(i,1,j) - 1;
        up_chunk_inc_array(i,2,j) = up_chunk_inc_array_pass(i,2,j) - 1;
        
    end
    
end


%Use function 'chunkify' to average each chunk dataset into a point.
%Calculate average of each chunk, weighted stdev of each chunk, average of
%all the chunks together, stdev of all the chunks together, stdev of the
%average of all the chunks together, and stdev of all the stdevs together
if EDM_sim ==1
    [vmon_avg_ramp_up_avg_chunk, vmon_avg_ramp_up_stdev_chunk,...
        vmon_avg_ramp_up_avg,vmon_avg_ramp_up_stdev,...
        vmon_avg_ramp_up_stdev_avg_chunk,vmon_avg_ramp_up_stdev_stdev_chunk] ...
        = chunkify(num_files,num_up_chunks,...
        num_ramp_up_points,up_chunk_array,vmon_avg_ramp_up_raw,...
        vmon_weight_avg_ramp_up_raw,vmon_avg_scale);

    [lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,...
        lcm1_avg_ramp_up_avg,lcm1_avg_ramp_up_stdev,...
        lcm1_avg_ramp_up_stdev_avg_chunk,lcm1_avg_ramp_up_stdev_stdev_chunk] ...
        = chunkify(num_files,num_up_chunks,...
        num_ramp_up_points,up_chunk_array,lcm1_avg_ramp_up_raw,...
        lcm1_weight_avg_ramp_up_raw,lcm1_avg_scale);

    [vmon_avg_ramp_down_avg_chunk, vmon_avg_ramp_down_stdev_chunk,...
        vmon_avg_ramp_down_avg,vmon_avg_ramp_down_stdev,...
        vmon_avg_ramp_down_stdev_avg_chunk,vmon_avg_ramp_down_stdev_stdev_chunk] ...
        = chunkify(num_files,num_down_chunks,...
        num_ramp_down_points,down_chunk_array,vmon_avg_ramp_down_raw,...
        vmon_weight_avg_ramp_down_raw,vmon_avg_scale);

    [lcm1_avg_ramp_down_avg_chunk, lcm1_avg_ramp_down_stdev_chunk,...
        lcm1_avg_ramp_down_avg,lcm1_avg_ramp_down_stdev,...
        lcm1_avg_ramp_down_stdev_avg_chunk,lcm1_avg_ramp_down_stdev_stdev_chunk] ...
        = chunkify(num_files,num_down_chunks,...
        num_ramp_down_points,down_chunk_array,lcm1_avg_ramp_down_raw,...
        lcm1_weight_avg_ramp_down_raw,lcm1_avg_scale);

    [vmon_avg_trash_avg_chunk, vmon_avg_trash_stdev_chunk,...
        vmon_avg_trash_avg,vmon_avg_trash_stdev,...
        vmon_avg_trash_stdev_avg_chunk,vmon_avg_trash_stdev_stdev_chunk] ...
        = chunkify(num_files,num_trash_chunks,...
        num_trash_points,trash_chunk_array,vmon_avg_trash_raw,...
        vmon_weight_avg_trash_raw,vmon_avg_scale);

    [lcm1_avg_trash_avg_chunk, lcm1_avg_trash_stdev_chunk,...
        lcm1_avg_trash_avg,lcm1_avg_trash_stdev,...
        lcm1_avg_trash_stdev_avg_chunk,lcm1_avg_trash_stdev_stdev_chunk] ...
        = chunkify(num_files,num_trash_chunks,...
        num_trash_points,trash_chunk_array,lcm1_avg_trash_raw,...
        lcm1_weight_avg_trash_raw,lcm1_avg_scale);
    
    if inclusive_data == 1

        [lcm1_avg_ramp_down_inc_chunk, lcm1_avg_ramp_down_inc_stdev_chunk,...
            lcm1_avg_ramp_down_inc_avg,lcm1_avg_ramp_down_inc_stdev,...
            lcm1_avg_ramp_down_inc_stdev_avg_chunk,lcm1_avg_ramp_down_inc_stdev_stdev_chunk] ...
            = chunkify(num_files,num_down_chunks,...
            num_ramp_down_inc_points,down_chunk_inc_array,lcm1_avg_ramp_down_inc_raw,...
            lcm1_weight_ramp_down_inc_raw,lcm1_avg_scale);

        [lcm1_avg_ramp_up_inc_chunk, lcm1_avg_ramp_up_inc_stdev_chunk,...
            lcm1_avg_ramp_up_inc_avg,lcm1_avg_ramp_up_inc_stdev,...
            lcm1_avg_ramp_up_inc_stdev_avg_chunk,lcm1_avg_ramp_up_inc_stdev_stdev_chunk] ...
            = chunkify(num_files,num_up_chunks,...
            num_ramp_up_inc_points,up_chunk_inc_array,lcm1_avg_ramp_up_inc_raw,...
            lcm1_weight_ramp_up_inc_raw,lcm1_avg_scale);
    end
        
end


gaus_a = zeros(num_files,1);
gaus_range = zeros(num_files,1);
num_fit_points = zeros(num_files,1);

for i = 1:num_files
    num_fit_points(i) = 2*num_ramp_up_points(i);
    gaus_a(i) = 1/(abs(lcm1_avg_ramp_up_stdev(i))*sqrt(2*pi));
    gaus_range(i) = 12*abs(lcm1_avg_ramp_up_stdev(i));
    step = gaus_range(i)/(num_fit_points(i));
    gaus_fit = zeros(num_ramp_up_points(i),2);
    for j = 1:num_fit_points(i)
        fit_x = j*step + lcm1_avg_ramp_up_avg(i) - 0.5*gaus_range(i);
        gaus_fit(j,1) = fit_x;
        gaus_fit(j,2) = gaus_a(i)*exp(-(fit_x - lcm1_avg_ramp_up_avg(i))^2/(2*lcm1_avg_ramp_up_stdev(i)^2));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%% ramp test code end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fit LSR slope to picoammeter vs. current source input

if leakage_sensitivity_test == 1
    fit_params = zeros(num_files,4);
    %a: offset, offset uncertainty
    a = zeros(num_files,2);
    %b: slope, slope uncertainty
    b = zeros(num_files,2);
    %calculate residuals and the line
    leakage_residual = zeros(num_files,num_rows,2);
    fit_line_x = zeros(num_files,2);
    fit_line_y = zeros(num_files,2);
    trim = 7; %artificially zoom in on the points in the ~ 100 pA range

    for i =1:num_files
        fit_params(i,:) = array_linfit(current_source_avg(i,:),lcm1_avg_raw(i,:),lcm1_avg_stdev_raw(i,:),numpoints(i));
        %convert y values to nA, slope will be overall unitless
        a(i,1) = fit_params(i,1)*lcm1_avg_scale*1e-3;
        a(i,2) = fit_params(i,2)*abs(lcm1_avg_scale)*1e-3;
        b(i,1) = fit_params(i,3)*lcm1_avg_scale*1e-3;
        b(i,2) = fit_params(i,4)*abs(lcm1_avg_scale)*1e-3;
        fit_line_x(i,1) = 1.06*current_source_avg(i,1+trim);
        fit_line_x(i,2) = 1.06*current_source_avg(i,end-trim);
        fit_line_y(i,1) = a(i,1)+b(i,1)*fit_line_x(i,1);
        fit_line_y(i,2) = a(i,1)+b(i,1)*fit_line_x(i,2);
        for j = 1:numpoints(i)
            leakage_residual(i,j,1) = lcm1_avg(i,j)*1e-3 - (a(i) + b(i)*current_source_avg(i,j));
        end
    end
end

% compare charging in +V dirxn vs. -V dirxn
%sum the leakage current when charging top electrode in the +V dirxn
lcm1_avg_trash_charge_pos = zeros(num_files,1);
lcm1_avg_trash_stdev_charge_pos = zeros(num_files,1);
%sum the leakage curent when charging top electrode in the -V dirxn
lcm1_avg_trash_charge_neg = zeros(num_files,1);
lcm1_avg_trash_stdev_charge_neg = zeros(num_files,1);


for i =1:num_files
    %only count an even # of up ramps and down ramps
    if mod(num_trash_chunks(i),2) == 1
        num_trash_chunks(i) = num_trash_chunks(i) - 1;
    end
    for j = 1:num_trash_chunks(i)
       chunk_begin = trash_chunk_array(i,2,j) +1;
       chunk_end = trash_chunk_array(i,2,j+1);
       if mod(j,2) == 1
           lcm1_avg_trash_charge_pos(i,end+1) = sum(lcm1_avg_trash(i,chunk_begin:chunk_end) - lcm1_avg_offset(i));
           lcm1_avg_trash_stdev_charge_pos(i,end+1) = sum(lcm1_inv_weight_avg_trash(i,chunk_begin:chunk_end));
       elseif mod(j,2) == 0
           lcm1_avg_trash_charge_neg(i,end+1) = sum(lcm1_avg_trash(i,chunk_begin:chunk_end) - lcm1_avg_offset(i));
           lcm1_avg_trash_stdev_charge_neg(i,end+1) = sum(lcm1_inv_weight_avg_trash(i,chunk_begin:chunk_end));
       end
    end
end






%start_point(5) = 50; %lost 30 min o fdata for 2-10-2017
%time(5,:) = time(5,:) + 33.0; % lost 30 min of data for 2-10-2017

for i = 1:num_files
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
title_string = 'leakage current comparison plot';
x_label = 'PS voltage (-kV)';
y_label = 'residual (pA)';
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





   
%%%%%%%%%%%%%%%%%%%%%%% of max field v. gap size %%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(gap_size, field_max,...
%    'o','Color', cmap(num_files+1,:),'MarkerSize', 6, 'LineWidth', 2.0); %nA
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out
   
%%%%%%%%%%%%%%%%%%%%%% of max leakage v. gap size %%%%%%%%%%%%%%%%%%%%%%%%%

% plot(gap_size, leakage_max*1e-3,... %convert to nA
%    'o','Color', cmap(i+1,:),'MarkerSize', 6, 'LineWidth', 2.0); %nA
% %axis(plot_bounds)
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out
   
   
%%%%%%%%%%%%%%%%%%%%%  of max field v. max voltage %%%%%%%%%%%%%%%%%%%%%%%%

% plot(field_max, voltage_max,... %convert to nA
%    'o','Color', cmap(i+1,:),'MarkerSize', 6, 'LineWidth', 2.0); %nA
% axis(plot_bounds)
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out

%%%%%%%%%%%%%%%%%%%%%  of max field v. max leakage %%%%%%%%%%%%%%%%%%%%%%%%
% figure4 = figure('Units','normalized')
% plot(field_max, leakage_max*1e-3,... %convert to nA
%    'o','Color', cmap(num_files+1,:),'MarkerSize', 6, 'LineWidth', 2.0); %nA
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out

%%%%%%%%%%%%%%%%%%%%%%%%%  of max voltage v. gap %%%%%%%%%%%%%%%%%%%%%%%%%%
% figure5 = figure('Units','normalized')
% plot(gap_size, voltage_max,... %convert to nA
%    'o','Color', cmap(num_files+1,:),'MarkerSize', 6, 'LineWidth', 2.0); %nA
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out


title(title_string,'FontSize',40)
xlabel(x_label,'FontSize',32)
ylabel(y_label,'FontSize',32)



%%%%%%%%%%%%%%%%%%% sample set plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure1 =figure('Units','normalized')
% ramp_up_point = 58;
% errorbar(imon_avg(1,1:ramp_up_point) - imon_avg_offset(1),...
%      vmon_avg(1,1:ramp_up_point) - vmon_avg_offset(1),...
%      0.5*vmon_avg_stdev(1,1:ramp_up_point),... 
%      0.5*vmon_avg_stdev(1,1:ramp_up_point),...
%      0.5*imon_avg_stdev(1,1:ramp_up_point),...
%      0.5*imon_avg_stdev(1,1:ramp_up_point),...
%      'o','Color', 'red','MarkerSize', 10, 'LineWidth', 2.0); hold on;
% errorbar(imon_avg(1,ramp_up_point+1:end) - imon_avg_offset(1),...
%      vmon_avg(1,ramp_up_point+1:end) - vmon_avg_offset(1),...
%      0.5*vmon_avg_stdev(1,ramp_up_point+1:end),...
%      0.5*vmon_avg_stdev(1,ramp_up_point+1:end),...
%      0.5*imon_avg_stdev(1,ramp_up_point+1:end),...
%      0.5*imon_avg_stdev(1,ramp_up_point+1:end),...
%      'o','Color', 'blue','MarkerSize', 10, 'LineWidth', 2.0); hold on; 
% %axis(plot_bounds)
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out
% title(filenames,'FontSize',40)
% xlabel('imon','FontSize',32)
% ylabel('vmon','FontSize',32)
% l = legend('show'); l.String = [{'ramp up'}, {'ramp down'}]; l.FontSize = 32; l.Location = 'northeast outside';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%% psvoltage v. HV Set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure2 =figure('Units','normalized')
% errorbar(time_raw(1,1:ramp_up_point),...
%     vmon_avg(1,1:ramp_up_point) - vmon_avg_offset(1),...
%     0.5*vmon_avg_stdev(1,1:ramp_up_point),... 
%     0.5*vmon_avg_stdev(1,1:ramp_up_point),...
%     'o','Color', 'red','MarkerSize', 10, 'LineWidth', 2.0); hold on;
% errorbar(time_raw(1,ramp_up_point+1:end),...
%     vmon_avg(1,ramp_up_point+1:end) - vmon_avg_offset(1),...
%     0.5*vmon_avg_stdev(1,ramp_up_point+1:end),...
%     0.5*vmon_avg_stdev(1,ramp_up_point+1:end),...
%     'x','Color', 'blue','MarkerSize', 10, 'LineWidth', 2.0); hold on; 
% %axis(plot_bounds)
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out
% title(filenames,'FontSize',40)
% xlabel('setv','FontSize',32)
% ylabel('vmon','FontSize',32)
% l = legend('show'); l.String = [{'ramp up'}, {'ramp down'}]; l.FontSize = 32; l.Location = 'northeast outside';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%% ps current v. HV Set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure3 =figure('Units','normalized')
% errorbar(time_raw(1,1:ramp_up_point),...
%     imon_avg(1,1:ramp_up_point) - imon_avg_offset(1),...
%     0.5*imon_avg_stdev(1,1:ramp_up_point),... 
%     0.5*imon_avg_stdev(1,1:ramp_up_point),...
%     'o','Color', 'red','MarkerSize', 10, 'LineWidth', 2.0); hold on;
% errorbar(time_raw(1,ramp_up_point+1:end),...
%     imon_avg(1,ramp_up_point+1:end) - imon_avg_offset(1),...
%     0.5*imon_avg_stdev(1,ramp_up_point+1:end),...
%     0.5*imon_avg_stdev(1,ramp_up_point+1:end),...
%     'x','Color', 'blue','MarkerSize', 10, 'LineWidth', 2.0); hold on; 
% %axis(plot_bounds)
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out
% title(filenames,'FontSize',40)
% xlabel('setv','FontSize',32)
% ylabel('imon','FontSize',32)
% l = legend('show'); l.String = [{'ramp up'}, {'ramp down'}]; l.FontSize = 32; l.Location = 'northeast outside';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%% ps current v. HV Set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure4 =figure('Units','normalized')
% 
% errorbar(time_raw(1,1:ramp_up_point),...
%     pressure_avg(1,1:ramp_up_point),...
%     0.5*pressure_avg_stdev(1,1:ramp_up_point),... 
%     0.5*pressure_avg_stdev(1,1:ramp_up_point),...
%     'o','Color', 'red','MarkerSize', 10, 'LineWidth', 2.0); hold on;
% errorbar(time_raw(1,ramp_up_point+1:end),...
%     pressure_avg(1,ramp_up_point+1:end),...
%     0.5*pressure_avg_stdev(1,ramp_up_point+1:end),...
%     0.5*pressure_avg_stdev(1,ramp_up_point+1:end),...
%     'x','Color', 'blue','MarkerSize', 10, 'LineWidth', 2.0); hold on; 
% %axis(plot_bounds)
% pbaspect([1.33 1 1])
% ax = gca; % current axes
% ax.TickDir = 'out'; % make ticks point out
% title(filenames,'FontSize',40)
% xlabel('setv','FontSize',32)
% ylabel('ig','FontSize',32)
% l = legend('show'); l.String = [{'ramp up'}, {'ramp down'}]; l.FontSize = 32; l.Location = 'northeast outside';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


left_color = cmap(num_files+1,:);
right_color = cmap(num_files+1+2*num_files+2,:);
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);

for i = 1:num_files
    
%    figure1 =figure('Units','normalized')

%%%%%%%% grid of psvoltage v. pscurrent %%%%%%%%%%%%%%%%%%%%%%%%%%
%   subplot(num_files,1,num_files + 1 - i)
%     plot(imon_avg(i,start_point(i):end_point(i)) - imon_avg_offset(i),...
%         vmon_avg(i,start_point(i):end_point(i)) - vmon_avg_offset(i),...
%         'o','Color', cmap(i+1,:),'MarkerSize', 2, 'LineWidth', 2.0); %nA
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
%    ax.TickDir = 'out'; % make ticks point out    

%%%%%%%% grid of field v. ps voltage %%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(num_files,1,num_files + 1 - i)
%     plot(vmon_avg(i,start_point(i):end_point(i)) - vmon_avg_offset(i),...
%         field_avg(i,start_point(i):end_point(i)) - field_avg_offset(i),...
%        'o','Color', cmap(i+1,:),'MarkerSize', 2, 'LineWidth', 2.0); %nA
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
%    ax.TickDir = 'out'; % make ticks point out

%%%%%%%% grid of leakage current v. ps voltage %%%%%%%%%%%%%%%%%%%%%%%%%%
%      plot(vmon_avg(i,start_point(i):end_point(i)) - vmon_avg_offset(i),...
%          (lcm1_avg(i,start_point(i):end_point(i)) - lcm1_avg_offset(i))*1e-3,...
%         'o','Color', cmap(i+1,:),'MarkerSize', 2, 'LineWidth', 2.0); %nA
% %     %axis(plot_bounds)
%      pbaspect([1.33 1 1])
%      ax = gca; % current axes
%      ax.TickDir = 'out'; % make ticks point out


%%%%%%%%%% grid of leakage current stdev v. ps voltage %%%%%%%%%%%%%%%%%%%%

% %    subplot(num_files,1,num_files + 1 - i)
%    plot(vmon_avg(i,start_point(i):end_point(i)),...
%        lcm1_avg_stdev(i,start_point(i):end_point(i))*1e0,...
%        'o','Color', cmap(i+1,:),'MarkerSize', 2, 'LineWidth', 2.0); %pA
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
%    ax.TickDir = 'out'; % make ticks point out


% %%%%%%%%%%%% %grid of leakage current v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    subplot(num_files,1,num_files + 1 - i)
%%note: this won't run in versions previous to 2016a, use plotyy instead
%    yyaxis left
%       plot(time(i,start_point(i):end_point(i)) - time(i,start_point(i)),...
%        vmon_avg(i,start_point(i):end_point(i))-vmon_avg_offset(i),...
%        '-','Color', cmap(i+1+2*num_files+2,:),'MarkerSize', 3,...
%        'LineWidth', 2.0);
%    axis(plot_bounds)
% 
%%note: this won't run in versions previous to 2016a, use plotyy instead
%    yyaxis right
%    plot(time(i,start_point(i):end_point(i)) - time(i,start_point(i)),...
%        lcm1_avg(i,start_point(i):end_point(i))-lcm1_avg_offset(i),...
%        '-','Color', cmap(i+1,:),'MarkerSize', 3,...
%        'LineWidth', 2.0);
% %   ylim([ymin_right ymax_right])
%    pbaspect([7 1 1])
%    ax = gca; % current axes
%    text(...
%        'Position',[xmax ymax_right],...
%        'String',legend_titles(i),...
%        'HorizontalAlignment','right','VerticalAlignment','top',...
%        'FontSize',14);
%    ax.TickDir = 'out'; % make ticks point out

%%%%%%% %grid of log leakage current & viktage v. time  %%%%%%%%%%%%%%%%%%%
%     subplot(num_files,1,num_files + 1 - i)
% %    figure
%     ax = gca; % current axes
% 
%     text(...
%        'Position',[xmax ymax_right],...
%        'String',legend_titles(i),...
%        'HorizontalAlignment','right','VerticalAlignment','top',...
%        'FontSize',14);
%     ax.TickDir = 'out'; % make ticks point out
% %%note: this won't run in versions previous to 2016a, use plotyy instead
%     yyaxis left
% 
%        plot(lcm1_avg_log_pos(i,start_point(i):approximate_time_points,1) - time(i,start_point(i)),...
%        lcm1_avg_log_pos(i,start_point(i):approximate_time_points,2),...
%        'x','Color', cmap(i+1+num_files+1,:),'MarkerSize', 4,...
%        'LineWidth', 1.0); hold on;
%     plot(lcm1_avg_log_neg(i,start_point(i):approximate_time_points,1) - time(i,start_point(i)),...
%        lcm1_avg_log_neg(i,start_point(i):approximate_time_points,2),...
%        'o','Color', cmap(i+1,:),'MarkerSize', 4,...
%        'LineWidth', 1.0);   
%    
%    axis(plot_bounds)
% %%note: this won't run in versions previous to 2016a, use plotyy instead
%     yyaxis right
%     plot(time(i,start_point(i):approximate_time_points) - time(i,start_point(i)),...
%        vmon_avg(i,start_point(i):approximate_time_points)-vmon_avg_offset(i),...
%        '-','Color', cmap(i+1+2*num_files+2,:),'MarkerSize', 3,...
%        'LineWidth', 2.0);
%     ylim([ymin_right ymax_right])
% 
% 
% 
%     pbaspect([7 1 1])



%%%%%%%%%%%%%%%%% grid of power ps voltage v. time %%%%%%%%%%%%%%%%%%%%%%%%
% %    subplot(num_files,1,num_files + 1 - i)
%    plot(time(i,start_offset(i):end_offset(i)) - time(i,start_offset(i)),...
%        (vmon_avg(i,start_offset(i):end_offset(i)) - vmon_avg_offset(i)),...
%        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
% %    axis(plot_bounds)
%    pbaspect([4 1 1])
%    ax = gca; % current axes
%     text('Position',[xmax ymax],'String',legend_titles(i),...
%         'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
%    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out



% %%%%%%%%% grid of field v. time with corresponding leakage current %%%%%%%%
%    subplot(2*num_files,1,2*num_files + 1 - i)
%    plot(time(i,start_offset(i):end_offset(i)) - time(i,start_offset(i)),...
%        (vmon_avg(i,start_offset(i):end_offset(i)) - vmon_avg_offset(i))/gap_size,...
%        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0); 
%    pbaspect([4 1 1]);    ax = gca; ax.FontSize = 16; ax.TickDir = 'out'; hold on;
%    subplot(2*num_files,1,2*num_files - i)
%    plot(time(i,start_offset(i):end_offset(i)) - time(i,start_offset(i)),...
%        lcm1_avg(i,start_offset(i):end_offset(i)) - lcm1_avg_offset(i),...
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
% %   subplot(num_files,1,num_files + 1 - i)
%    %note: for 2017-10-23-175202-hv-1.txt, time_ramp_up and vmon_avg_ramp_up
%    %have different lengths... possible bugs
% 
%    plot(time_ramp_up(i,start_point(i):num_ramp_up_points), ...
%        vmon_avg_ramp_up(i,start_point(i):num_ramp_up_points) - vmon_avg_offset(i),...
%    'o','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
% %    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out
   
   %%%%%%%%%%%%%% ramp data plot of ramp voltage v time %%%%%%%%%%%%%%%%%%
%    figure
%    plot(time(i,1:numpoints(i)),...
%        vmon_avg(i,1:numpoints(i)) - vmon_avg_offset(i),...
%        'x','Color', 'black','MarkerSize', 8, 'LineWidth', 2.0); hold on;
%    plot(time_trash(i,trash_chunk_array(i,2,1)+1:trash_chunk_array(i,2,end)),...
%        vmon_avg_trash(i,trash_chunk_array(i,2,1)+1:trash_chunk_array(i,2,end)) - vmon_avg_offset(i),...
%        'o','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
%    plot(time_ramp_up(i,up_chunk_array(i,2,1)+1:up_chunk_array(i,2,end)), ...
%        vmon_avg_ramp_up(i,up_chunk_array(i,2,1)+1:up_chunk_array(i,2,end)) - vmon_avg_offset(i),...
%        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
%    plot(time_ramp_down(i,down_chunk_array(i,2,1)+1:down_chunk_array(i,2,end)), ...
%        vmon_avg_ramp_down(i,down_chunk_array(i,2,1)+1:down_chunk_array(i,2,end)) - vmon_avg_offset(i),...
%        'o','Color', cmap(1+i+2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0);
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
% %    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out
%    
%       %%%%%%%%%%%%%% ramp data plot of ramp leakage v time %%%%%%%%%%%%%%%%%%
%    figure
%    plot(time(i,1:numpoints(i)),...
%        lcm1_avg(i,1:numpoints(i)) - lcm1_avg_offset(i),...
%        'x','Color', 'black','MarkerSize', 8, 'LineWidth', 2.0); hold on;
% %    plot(time_trash(i,trash_chunk_array(i,2,1)+1:trash_chunk_array(i,2,end)),...
% %        lcm1_avg_trash(i,trash_chunk_array(i,2,1)+1:trash_chunk_array(i,2,end)) - lcm1_avg_offset(i),...
% %        's','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
% 
% % plot(time_ramp_up(i,up_chunk_array(i,2,1)+1:up_chunk_array(i,2,end)), ...
% %        lcm1_avg_ramp_up(i,up_chunk_array(i,2,1)+1:up_chunk_array(i,2,end)) - lcm1_avg_offset(i),...
% %        '^','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
%    plot(time_ramp_up_inc(i,:), ...
%        lcm1_avg_ramp_up_inc(i,:) - lcm1_avg_offset(i),...
%        'o','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
% 
% %   plot(lcm1_charge_neg_time(i,:),lcm1_avg_charge_neg(i,:) - lcm1_avg_offset(i),...
% %        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
% %    plot(lcm1_discharge_neg_time(i,:),lcm1_avg_discharge_neg(i,:) - lcm1_avg_offset(i),...
% %        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
%    plot(time_ramp_down_inc(i,:), ...
%        lcm1_avg_ramp_down_inc(i,:) - lcm1_avg_offset(i),...
%        'o','Color', cmap(1+i+ 1*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
% %    plot(time_ramp_down(i,down_chunk_array(i,2,1)+1:down_chunk_array(i,2,end)), ...
% %        lcm1_avg_ramp_down(i,down_chunk_array(i,2,1)+1:down_chunk_array(i,2,end)) - lcm1_avg_offset(i),...
% %        '^','Color', cmap(1+i+2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
% 
% %    plot(lcm1_charge_pos_time(i,:),lcm1_avg_charge_pos(i,:) - lcm1_avg_offset(i),...
% %        'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;   
% %    plot(lcm1_discharge_pos_time(i,:),lcm1_avg_discharge_pos(i,:) - lcm1_avg_offset(i),...
% %        'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;    
% 
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
% %    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out
   
   
   
%%%%%%%%%%%%%%%%% discharging summed leakage data v time %%%%%%%%%%%%%%%%%%
   figure
   plot(1:time_step(i):(max_length_discharge+1)*time_step(i),...
       -(lcm1_avg_discharge_pos_norm(i,:)),...
       'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
   plot(1:time_step(i):(max_length_discharge+1)*time_step(i),...
       lcm1_avg_discharge_neg_norm(i,:),...
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
   plot(1:time_step(i):(max_length_charge+1)*time_step(i),...
       lcm1_avg_charge_pos_norm(i,:),...
       'o','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
   plot(1:time_step(i):(max_length_charge+1)*time_step(i),...
       -(lcm1_avg_charge_neg_norm(i,:)),...
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
%     plot([1:1:num_up_chunks(i)], vmon_avg_ramp_up_avg_chunk(i,1:num_up_chunks(i)) - vmon_avg_offset(i),...
%     'x','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0);
%     %errorbar([1:1:num_up_chunk_rows(i), vmon_avg_ramp_up_avg_chunk, vmon_avg_ramp_up_stdev_chunk,'o');
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
% 
% %%%%%%%%%%% ramp data plot of mean down chunk voltage v chunk # %%%%%%%%%%
%     figure2 = figure('Units','normalized')
%     plot([ 1:1:num_down_chunks(i) ], vmon_avg_ramp_down_avg_chunk(i,1:num_down_chunks(i)) - vmon_avg_offset(i),...
%     'o','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0);
%     %errorbar([1:1:num_up_chunk_rows(i), lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,'o');
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
% 
%%%%%%%%% ramp data plot of stdev hi/lo chunk ps voltage v chunk # %%%%%%%
%     figure3 = figure('Units','normalized')
%     plot(([1:1:num_up_chunks(i)] - 0.5)*2, vmon_avg_ramp_up_stdev_chunk(i,1:num_up_chunks(i)),...
%     'x','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0); hold on;
%     plot([ 1:1:num_down_chunks(i) ]*2, vmon_avg_ramp_down_stdev_chunk(i,1:num_down_chunks(i)),...
%     'o','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0);
%    %errorbar([1:1:num_up_chunk_rows(i), lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,'o');
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
%     plot([ 1:1:num_down_chunks(i) ], lcm1_avg_ramp_down_avg_chunk(i,1:num_down_chunks(i)) - lcm1_avg_offset(i),...
%     'x','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0); 
%     if inclusive_data == 1
%         hold on;
%         plot([ 1:1:num_down_chunks(i) ], lcm1_avg_ramp_down_inc_chunk(i,1:num_down_chunks(i)) - lcm1_avg_offset(i),...
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
    plot([1:1:num_up_chunks(i)], lcm1_avg_ramp_up_avg_chunk(i,1:num_up_chunks(i)) - lcm1_avg_offset(i),...
        '^','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0); hold on;
    plot([ 1:1:num_up_chunks(i) ], lcm1_avg_ramp_up_inc_chunk(i,1:num_up_chunks(i)) - lcm1_avg_offset(i),...
        '^','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0); hold on;
    plot([2:2:2*num_down_chunks(i)], lcm1_avg_ramp_down_avg_chunk(i,1:num_down_chunks(i)) - lcm1_avg_offset(i),...
        'v','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0); hold on;
    plot([ 2:2:2*num_down_chunks(i) ], lcm1_avg_ramp_down_inc_chunk(i,1:num_down_chunks(i)) - lcm1_avg_offset(i),...
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
    plot(([1:1:num_up_chunks(i)] - 0.5)*2, lcm1_avg_ramp_up_stdev_chunk(i,1:num_up_chunks(i)),...
    '^','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
    plot(([1:1:num_up_chunks(i)] - 0.5)*2, lcm1_avg_ramp_up_inc_stdev_chunk(i,1:num_up_chunks(i)),...
    '^','Color', cmap(1+i+ 0*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
    plot([2:2:2*num_up_chunks(i)], lcm1_avg_ramp_down_stdev_chunk(i,1:num_down_chunks(i)),...
    'v','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
    plot([2:2:2*num_up_chunks(i)], lcm1_avg_ramp_down_inc_stdev_chunk(i,1:num_down_chunks(i)),...
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
%     plot(([1:1:num_down_chunks(i)] - 0.5)*2, lcm1_avg_ramp_down_stdev_chunk(i,1:num_down_chunks(i)),...
%     'x','Color', cmap(1+i+ 2*(num_files+1),:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
%     plot(([1:1:num_down_chunks(i)] - 0.5)*2, lcm1_avg_ramp_down_inc_stdev_chunk(i,1:num_down_chunks(i)),...
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
% 
% %%%%%%%%% ramp data plot of stdev hi/lo chunk leakage current v chunk # %%%%%%%
%     figure6= figure('Units','normalized')
%     plot(([1:1:num_up_chunks(i)] - 0.5)*2, lcm1_avg_ramp_up_stdev_chunk(i,1:num_up_chunks(i)),...
%     'x','Color', cmap(2,:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
%     plot([ 1:1:num_down_chunks(i) ]*2, lcm1_avg_ramp_down_stdev_chunk(i,1:num_down_chunks(i)),...
%     'o','Color', 'blue','MarkerSize', 8, 'LineWidth', 1.5);
%     %errorbar([[1:1:num_up_chunk_rows(i)], lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,'o');
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
%     plot([1:1:num_up_chunks(i)], lcm1_avg_ramp_up_stdev_chunk(i,1:num_up_chunks(i)),...
%     'x','Color', cmap(2,:),'MarkerSize', 8, 'LineWidth', 1.5); hold on;
%     %errorbar([[1:1:num_up_chunk_rows(i)], lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,'o');
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
%    plot(time(i,start_offset(i):end_offset(i)), vmon_avg(i,start_offset(i):end_offset(i)) - vmon_avg_offset(i),...
%    'o','Color', cmap(1+i,:),'MarkerSize', 6, 'LineWidth', 1.0); hold on;
% 
%    for j = 1:num_up_chunks(i)
%        chunk_begin = up_chunk_array(i,2,j) +1;
%        chunk_end = up_chunk_array(i,2,j+1);
%        plot(time_ramp_up(i,chunk_begin:chunk_end), vmon_avg_ramp_up(i,chunk_begin:chunk_end) - vmon_avg_offset(i),...
%        'x','Color', 'green','MarkerSize', 6, 'LineWidth', 2.0); hold on;
%    end
%    for j = 1:num_down_chunks(i)
%        chunk_begin = down_chunk_array(i,2,j) +1;
%        chunk_end = down_chunk_array(i,2,j+1);
%        plot(time_ramp_down(i,chunk_begin:chunk_end), vmon_avg_ramp_down(i,chunk_begin:chunk_end) - vmon_avg_offset(i),...
%       'x','Color', 'blue','MarkerSize', 6, 'LineWidth', 2.0); hold on;
%    end
%    
%    for j = 1:num_trash_chunks(i)
%        chunk_begin = trash_chunk_array(i,2,j) +1;
%        chunk_end = trash_chunk_array(i,2,j+1);
%        plot(time_trash(i,chunk_begin:chunk_end), vmon_avg_trash(i,chunk_begin:chunk_end) - vmon_avg_offset(i),...
%       'x','Color', 'cyan','MarkerSize', 6, 'LineWidth', 2.0); hold on;
%    end
%    
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
% %    ax.XDir = 'reverse'; % voltage decreases left to right
%     ax.FontSize = 32;
%    ax.TickDir = 'out'; % make ticks point out

%%%%%%%%%%%%%% ramp data plot of lcm1 v time %%%%%%%%%%%%%%%%%%%%%%%%%
%   subplot(num_files,1,num_files + 1 - i)
%    plot(time(i,start_point(i):end_point(i))  - time(i,start_point(i)),...
%        lcm1_avg(i,start_point(i):end_point(i)) - lcm1_avg_offset(i),...
%    '.','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
%
%    axis(plot_bounds)
%    pbaspect([4 1 1])
%    ax = gca; % current axes
%%    ax.XDir = 'reverse'; % voltage decreases left to right
%%    ax.FontSize = 32;
%    ax.TickDir = 'out'; % make ticks point out

%%%%%%%%%%%%%%%%% grid of imon v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    subplot(2, 1 ,1)
%    plot(time(i,start_point(i):end_point(i)) - time(i,start_point(i)),...
%        pressure_avg_raw(i,start_point(i):end_point(i)),...
%        '-','Color', 'green','MarkerSize', 8, 'LineWidth', 2.0);
%    subplot(2,1,2)
%    plot(time(i,start_point(i):end_point(i)) - time(i,start_point(i)),...
%        lcm1_avg_stdev_raw(i,start_point(i):end_point(i)),...
%        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
    %axis(plot_bounds)
    %pbaspect([4 1 1])
%    ax = gca; % current axes
%    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out

%%%%%%%%%%%%%%%%%%%% grid of pressure v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert to mV
%    plot(time(i,start_point(i):end_point(i)) - time(i,start_point(i)),...
%        (pressure_avg_raw(i,start_point(i):end_point(i))-3.5)*1e3,...
%        '-','Color', 'green','MarkerSize', 8, 'LineWidth', 2.0); hold on;
%    plot(time(i,start_point(i):end_point(i)) - time(i,start_point(i)),...
%        lcm1_avg_raw(i,start_point(i):end_point(i))*1e3,...
%        '-','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
%    axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
%%    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out
%    legend_titles = {'pressure' ,'leakage current'};

%%%%%%%%%%%%%%%%%%% grid of power ps current v. time %%%%%%%%%%%%%%%%%%%%%%
%    subplot(num_files,1,num_files + 1 - i)
%    plot(time(i,start_point(i):numpoints(i)) - time(i,start_point(i)), imon_avg(i,start_point(i):numpoints(i)),...
%        '-','Color',cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0); hold on;
%%    plot([time(i,3184) - time(i,start_point(i)) time(i,3184) - time(i,start_point(i))],[min(imon_avg(i,:)) max(imon_avg(i,:))],...
%%        'k--','LineWidth',2.0);
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
%%    ax.FontSize = 32;
%    ax.TickDir = 'out'; % make ticks point out

%%%%%%%%%%%%%%%%% grid of resistance v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    subplot(num_files,1,num_files + 1 - i)
%    plot(time(i,start_point(i):numpoints(i)) - time(i,start_point(i)), ohm_avg(i,start_point(i):numpoints(i)),...
%        '-','Color',cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
%    axis(plot_bounds)
%    pbaspect([3 1 1])
%    ax = gca; % current axes
%    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out


end

% %%%%%%%% plot leakage current v. meas. test # %%%%%%%%%%%%%%%%%%%%%%
% figure1 = figure('Units','normalized')
% errorbar((1:1:35),...         
%      current_source_avg(1,:),current_source_avg_stdev(1,:),...
%     's','Color', 'blue','MarkerSize', 12, 'LineWidth', 1.0); hold on; %nA
% errorbar((1:1:35),...
%      (lcm1_avg(1,:) - lcm1_acc_test_offset(1))*1e-3,lcm1_avg_stdev(1,:)*1e-3,...
%     'g^','MarkerSize', 8, 'LineWidth', 1.0); hold on; %nA
%  errorbar((1:1:35),...
%      (lcm1_avg(2,:) - lcm1_acc_test_offset(2))*1e-3,lcm1_avg_stdev(2,:)*1e-3,...
%     'rx','MarkerSize', 8, 'LineWidth', 1.0); hold on; %nA
% 
% %     %axis(plot_bounds)
%  pbaspect([1.33 1 1])
%  ax = gca; % current axes
%  ax.TickDir = 'out'; % make ticks point out
 
%  %%%%%%% plot leakage sensitivity fit lines %%%%%%%%%%%%%%%%%%%%%%

% figure1 = figure('Units','normalized')
% 
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
% 
% plot(fit_line_x(1,:)*1e3,...
%      (fit_line_y(1,:)-a(1,1))*1e3,...
%     'k-', 'LineWidth', 1.0); hold on; %pA
% errorbar(current_source_avg(1,1+trim:end-trim)*1e3,...         
%      lcm1_avg(1,1+trim:end-trim)-a(1,1)*1e3,lcm1_avg_stdev(1,1+trim:end-trim),...
%     'o','Color', 'blue','MarkerSize', 12, 'LineWidth', 1.5); hold on; %pA
% %axis(plot_bounds)
%  pbaspect([1.33 1 1])
%  ax = gca; % current axes
%  ax.TickDir = 'out'; % make ticks point out
%  
%  title(title_string,'FontSize',40)
% xlabel(x_label,'FontSize',32)
% ylabel(y_label,'FontSize',32)
% % l = legend('show'); 
%  
%  annotation(figure1,'textbox',...
%    outside_plot,'String',str,'FontSize',32,'BackgroundColor',[1 1 1]);
%  
%   %%%%%%%% plot leakage sensitivity fit line residual %%%%%%%%%%%%%%%%%%%%%
% %convert to pA
% figure2 = figure('Units','normalized')
% 
% plot(fit_line_x(1,:)*1e3,...
%      [0 0]*1e3,...
%     'k--', 'LineWidth', 2.0); hold on; %nA
% errorbar(current_source_avg(1,1+trim:end-trim)*1e3,...         
%      leakage_residual(1,1+trim:end-trim,1)*1e3,lcm1_avg_stdev(1,1+trim:end-trim),...
%     'o','Color', 'blue','MarkerSize', 12, 'LineWidth', 1.5); hold on; %nA
% %axis(plot_bounds)
%  pbaspect([1.33 1 1])
%  ax = gca; % current axes
%  ax.TickDir = 'out'; % make ticks point out
%  
%  title(title_string,'FontSize',40)
% xlabel(x_label,'FontSize',32)
% ylabel(y_label,'FontSize',32)
% % l = legend('show'); 

%%%%%%%%%%%%%%%%%%%%%%% gaussian fit of ramp data %%%%%%%%%%%%%%%%%%%%%%%%%

%edges = [-70000 -700:.8:-550 0];
%histogram(lcm1_avg_ramp_up,edges); hold on;
%plot(gaus_fit(:,1),gaus_fit(:,2)*num_ramp_up_points,'r-','LineWidth',3.0);




%%%%%%%%%%%%%%%%%% gaussian fit of ramp data end %%%%%%%%%%%%%%%%%%%%%%%%%%


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