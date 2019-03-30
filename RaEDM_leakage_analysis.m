clearvars -except masterlist binsize
close all
%clc

%{

initialize 3D array. Needs to be as larges as the biggest data set.
if num_cols different, then sample set files have different formats
and should not be compared with this program. program is for runs taken
from February 1st 2017 onward, for which the picoammeter
is manually set to 2 uA range and the LabVIEW VI has ion gauge output in the 
data files.

08-24-2017 changing this to read data from jhv8-development onward, which
 is compatible with the new bipolar power supply. There is a 7th column
 that defines the polarity of Vmon. Also updating the pressure
 conversion for the new all-range model.

11-09-2017 updated to be compatible with bipolar and unipolar EDM ramp
 simulations. Tested that code is consistent with code used for 7/18 and
 7/05 simulations. CHUNK can now correctly parse HI (negative) and LO
 (positive) voltages. Voltage polarity is properly propagated for both
 polarities. Different sample rates are now configurable by user. 

11-19-2017 changing adding module to plot current from 11/13 leakage
current precision test. 

11-26-2017 CHUNK code now collects the data between the ramp up and ramp
down intervals, affectionally referred to as the "trash" data. This
includes data where the voltage is ramping up/down or is at 0. The 0 data
will be useful for comparing the leakage when the voltage is HI/LO vs.
off. The ramping data will be useful for comparing the symmetry of the
charging up/down currents.

08-21-2018 I'm separating all the plotting into a separate program. This
program will call it for plots

%}

set(0, 'DefaultFigureRenderer', 'Painters');

EDM_sim = 1;        % 0 = not an EDM sim
                    % 1 = EDM sim. CHUNK stuff applies

inclusive_data = 1; % 0 = do not incorporate ramp data points into ramp CHUNKS (exclusive)
%                      1 = incorporate ramp data (inclusive)

power_supply = 3; % 0 = unipolar Acopian 
%                   1 = bipolar AK (installed and tested 8-18-2017)
%                   2 = Ra EDM Spellman power supply for 2016 run
%                   3 = Ra EDM Spellman power supply for 2018/2019 data
%                   files

pressure_gauge = 1; % 0 = ion gauge 
%                     1 = all-range gauge (installed 7-18-2017)

sample_size_setting = 1; % # samples taken per measurement
                            % 0: Spinlab sample set, 10 samples by default
                            % 1: fast data, 1 
                            %sample = # samples / sampling frequency.

sample_rate = 2 ; % 0 = data saved every 0.02 min (8192 samples / 8 kHz)
                  % 1 = data saved every 0.01 min (8192 samples / 16 kHz)
                  % 2 = data saved every 0.056 min ( guessing 1024 samples
                  % / 18 kHz) Ra EDM system

leakage_sensitivity_test = 0; % 0 = not a sensitivity test
                              % 1 = leakage sensitivity test

%file_struct = dir('*hv-1.txt');
file_struct = dir('*21kv.txt');

num_files = length(file_struct); 

num_cols = 0;

num_rows = 0;

filenames = cell(num_files,1);

if sample_size_setting == 0
    
    sample_size = 10;
    
else 
    
    sample_size = 1;
    
end
                  
%if power_supply == 0 || power_supply == 1
if power_supply == 0 || power_supply == 1 || power_supply == 3

    for i = 1:num_files
        
        set = dlmread(file_struct(i).name);
        
        if length(set(:,1)) > num_rows
            
            num_rows = length(set(:,1));
            
        end
        
        if length(set(1,:)) > num_cols
            
            num_cols = length(set(1,:));
            
        end
        
        filenames{i} = file_struct(i).name(1:17);
        
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

        filenames{i} = file_struct(i).name(1:17);
        
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

vhvd_avg_raw = zeros(num_files,num_rows,1);
vhvd_avg_stdev_raw = zeros(num_files,num_rows,1);
lc_pkpk_raw = zeros(num_files,num_rows,1);
lc_pkpk = zeros(num_files,num_rows,1);

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
          
            % ~20 mV Low (positive), ~ 3.4 V high (negative)
            polarity_avg_raw(i,j) = data(i,j,8); 
         
            %don't think reduced sqrt is true for fast data
            vmon_avg_stdev_raw(i,j) = data(i,j,10)/sqrt(sample_size_setting); 
        
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

elseif power_supply == 3
    
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
            
            time_raw(i,j) = data(i,j,2); % s
            
            lcm1_avg_raw(i,j) = data(i,j,3);
            
            imon_avg_raw(i,j) = data(i,j,4);
            
            vmon_avg_raw_mag(i,j) = data(i,j,5);
                        
            lcm1_avg_stdev_raw(i,j) = data(i,j,6); %sample size?
            
            lcm1_weight_raw(i,j) = lcm1_avg_stdev_raw(i,j)^(-2);
            
            lcm1_avg_wt_raw(i,j) = lcm1_avg_raw(i,j)*lcm1_weight_raw(i,j);
            
            imon_avg_stdev_raw(i,j) = data(i,j,7);
            
            vmon_avg_stdev_raw(i,j) = data(i,j,8);
            
            vmon_weight_raw(i,j) = abs((vmon_avg_stdev_raw(i,j))^(-2));
            
            vmon_avg_wt_raw_mag(i,j) = vmon_avg_raw_mag(i,j)*vmon_weight_raw(i,j);
            
            polarity_avg_raw(i,j) = sign(data(i,j,9));
            
            vhvd_avg_raw(i,j) = data(i,j,9);
            
            vhvd_avg_stdev_raw(i,j) = data(i,j,10);
            
            lc_pkpk_raw(i,j) = data(i,j,11);
            
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

elseif power_supply == 3
    
    disp('Using Ra EDM transimpedance amplifier')
    
    lcm1_avg_scale = 60.*1e3; % 60 nA/V = 6e4 pA / V
    
    vmon_avg_scale = 30./10.; % 30 kV / 10 V
    
    imon_avg_scale = 300./10.; % 300 uA / 10 V
    
    ohm_avg_scale = 0;    
    
end



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
        
        lcm1_avg_stdev(i,j) = lcm1_avg_stdev_raw(i,j)*abs(lcm1_avg_scale);
        
        lcm1_weight(i,j) = 1/(lcm1_avg_stdev(i,j)^2);
        
        lcm1_avg_wt(i,j) = lcm1_avg_wt_raw(i,j)*lcm1_avg_scale;
        
        lc_pkpk(i,j) = lc_pkpk(i,j)*lcm1_avg_scale;
        
        
        ohm_avg(i,j) = vmon_avg(i,j)/imon_avg(i,j)*ohm_avg_scale;
        
        vmon_avg_wt(i,j) = vmon_avg_wt_raw(i,j)*vmon_avg_scale;
        
        imon_avg_wt(i,j) = imon_avg_wt_raw(i,j)*imon_avg_scale;
        
        field_avg(i,j) = field_avg_raw(i,j)*vmon_avg_scale; %(kV/cm)
 
    end
    
    for j = round(0.4 * numpoints(i),0):round(0.6 * numpoints(i),0)
    
        time_difference_raw(i,j) = time_raw(i,j) - time_raw(i,j-1);
  
    end
    
    time_step(i) = ...
        mean(time_difference_raw(i,...
        round(0.4 * numpoints(i),0):round(0.6 * numpoints(i),0)));

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
    
%elseif power_supply == 2
elseif power_supply == 2 || power_supply == 3
    
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

% calculate offset values, i.e. when high voltage is off.

lcm1_avg_offset_raw = calc_offset(lcm1_avg_raw,...
    lcm1_weight_raw,start_offset,end_offset,offset_length);

imon_avg_offset_raw = calc_offset(imon_avg_raw,imon_weight_raw,...
    start_offset,end_offset,offset_length);

vmon_avg_offset_raw = calc_offset(vmon_avg_raw,vmon_weight_raw,...
    start_offset,end_offset,offset_length);

for i = 1:num_files
   
    lcm1_avg_offset(i) = lcm1_avg_offset_raw(i)*lcm1_avg_scale;
    
    imon_avg_offset(i) = imon_avg_offset_raw(i)*imon_avg_scale;
    
    vmon_avg_offset(i) = vmon_avg_offset_raw(i)*vmon_avg_scale;
    
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

% separate voltage and current into high negative voltage, high positive
% voltage, and everything else (trash). E.g. for "ramp up" (negative high
% voltage):
%
% lcm1_avg_ramp_up_raw = all the leakage current values (raw voltage) at
%   negative high voltage in one array.
% up_chunk_array = the indices of the start/stop indices of each "chunk,"
%   or interval of time at negative high voltage
% num_ramp_up_points = the # of chunks. It's a few less than the length of
%   up_chunk_array because we don't count the first/last chunks

[lcm1_avg_trash_raw,lcm1_weight_avg_trash_raw,...
    lcm1_avg_trash_stdev_raw,time_trash,num_trash_points,...
    num_trash_chunks,lcm1_avg_ramp_down_raw,...
    lcm1_weight_avg_ramp_down_raw,lcm1_avg_ramp_down_stdev_raw,...
    time_ramp_down,num_ramp_down_points,num_down_chunks,...
    lcm1_avg_ramp_up_raw,lcm1_weight_avg_ramp_up_raw,...
    lcm1_avg_ramp_up_stdev_raw,time_ramp_up,num_ramp_up_points,...
    num_up_chunks,up_chunk_array,...
    down_chunk_array,trash_chunk_array,index_up_chunk_data,index_trash_chunk_data]...
    = sort_state(vmon_avg,vmon_weight,vmon_avg_wt,vmon_weight_raw,...
    start_point,end_point,lcm1_avg_raw,lcm1_weight_raw,time,num_rows,...
    numpoints,sampling_time,EDM_sim,power_supply);

[vmon_avg_trash_raw,vmon_weight_avg_trash_raw,...
    vmon_avg_trash_stdev_raw,time_trash,num_trash_points,...
    num_trash_chunks,vmon_avg_ramp_down_raw,...
    vmon_weight_avg_ramp_down_raw,vmon_avg_ramp_down_stdev,...
    time_ramp_down,num_ramp_down_points,num_ramp_down_chunks,...
    vmon_avg_ramp_up_raw,vmon_weight_avg_ramp_up_raw,...
    vmon_avg_ramp_up_stdev_raw,time_ramp_up,num_ramp_up_points,...
    num_ramp_up_chunks,up_chunk_array,...
    down_chunk_array,trash_chunk_array,index_up_chunk_data,index_trash_chunk_data]...
    = sort_state(vmon_avg,vmon_weight,vmon_avg_wt,vmon_weight_raw,...
    start_point,end_point,vmon_avg_raw,vmon_weight_raw,time,num_rows,...
    numpoints,sampling_time,EDM_sim,power_supply);

lcm1_avg_ramp_down = lcm1_avg_ramp_down_raw*lcm1_avg_scale;

lcm1_avg_ramp_up_stdev = lcm1_avg_ramp_up_stdev_raw*abs(lcm1_avg_scale);

lcm1_avg_ramp_up = lcm1_avg_ramp_up_raw*lcm1_avg_scale;

lcm1_avg_ramp_down_stdev = lcm1_avg_ramp_down_stdev_raw*abs(lcm1_avg_scale);

lcm1_avg_trash_stdev = lcm1_avg_trash_stdev_raw*abs(lcm1_avg_scale);

lcm1_avg_trash = lcm1_avg_trash_raw * lcm1_avg_scale;

% Separate the ramping data from the 'trash data', i.e. all the data that
% is not either high negative or high positive voltage
% [lcm1_avg_charge_neg_raw, lcm1_stdev_charge_neg_raw,...
%     lcm1_weight_charge_neg_raw, lcm1_charge_neg_time,lcm1_avg_charge_pos_raw,...
%     lcm1_stdev_charge_pos_raw,lcm1_weight_charge_pos_raw,...
%     lcm1_charge_pos_time,lcm1_avg_discharge_pos_raw,...
%     lcm1_stdev_discharge_pos_raw,lcm1_weight_discharge_pos_raw,...
%     lcm1_discharge_pos_time,lcm1_avg_discharge_neg_raw,...
%     lcm1_stdev_discharge_neg_raw,lcm1_weight_discharge_neg_raw,...
%     lcm1_discharge_neg_time,max_length_discharge,max_length_charge,...
%     charge_pos_index,charge_neg_index,discharge_pos_index,...
%     discharge_neg_index,discharge_index,charge_index,lcm1_avg_zero_raw,...
%     lcm1_stdev_zero_raw,lcm1_weight_avg_zero_raw,zero_time] = ...
%     find_charging_data(num_trash_chunks,trash_chunk_array,lcm1_inv_weight_avg_trash,...
%     lcm1_avg_trash,lcm1_avg_trash_raw,lcm1_inv_weight_avg_trash_raw,...
%     lcm1_weight_avg_trash_raw,time_trash);
% 
% [vmon_avg_charge_neg_raw, vmon_stdev_charge_neg_raw,...
%     vmon_weight_charge_neg_raw, vmon_charge_neg_time,vmon_avg_charge_pos_raw,...
%     vmon_stdev_charge_pos_raw,vmon_weight_charge_pos_raw,...
%     vmon_charge_pos_time,vmon_avg_discharge_pos_raw,...
%     vmon_stdev_discharge_pos_raw,vmon_weight_discharge_pos_raw,...
%     vmon_discharge_pos_time,vmon_avg_discharge_neg_raw,...
%     vmon_stdev_discharge_neg_raw,vmon_weight_discharge_neg_raw,...
%     vmon_discharge_neg_time,max_length_discharge,max_length_charge,...
%     charge_pos_index,charge_neg_index,discharge_pos_index,...
%     discharge_neg_index,discharge_index,charge_index,vmon_avg_zero_raw,...
%     vmon_stdev_zero_raw,vmon_weight_avg_zero_raw,zero_time] = ...
%     find_charging_data(num_trash_chunks,trash_chunk_array,lcm1_inv_weight_avg_trash,...
%     lcm1_avg_trash,lcm1_avg_trash_raw,lcm1_inv_weight_avg_trash_raw,...
%     lcm1_weight_avg_trash_raw,time_trash);

[lcm1_avg_charge_neg_raw, lcm1_stdev_charge_neg_raw,...
    lcm1_weight_charge_neg_raw, lcm1_charge_neg_time,lcm1_avg_charge_pos_raw,...
    lcm1_stdev_charge_pos_raw,lcm1_weight_charge_pos_raw,...
    lcm1_charge_pos_time,lcm1_avg_discharge_pos_raw,...
    lcm1_stdev_discharge_pos_raw,lcm1_weight_discharge_pos_raw,...
    lcm1_discharge_pos_time,lcm1_avg_discharge_neg_raw,...
    lcm1_stdev_discharge_neg_raw,lcm1_weight_discharge_neg_raw,...
    lcm1_discharge_neg_time,max_length_discharge,max_length_charge,...
    charge_pos_index,charge_neg_index,discharge_pos_index,...
    discharge_neg_index,discharge_index,charge_index,lcm1_avg_zero_raw,...
    lcm1_stdev_zero_raw,lcm1_weight_avg_zero_raw,...
    zero_time,num_zero_chunks,num_zero_points,zero_chunk_array] = ...
    find_charging_data(num_trash_chunks,trash_chunk_array,lcm1_avg_trash_stdev,...
    lcm1_avg_trash,lcm1_avg_trash_raw,...
    lcm1_avg_trash_stdev_raw,...
    lcm1_weight_avg_trash_raw,time_trash);

[vmon_avg_charge_neg_raw, vmon_stdev_charge_neg_raw,...
    vmon_weight_charge_neg_raw, vmon_charge_neg_time,vmon_avg_charge_pos_raw,...
    vmon_stdev_charge_pos_raw,vmon_weight_charge_pos_raw,...
    vmon_charge_pos_time,vmon_avg_discharge_pos_raw,...
    vmon_stdev_discharge_pos_raw,vmon_weight_discharge_pos_raw,...
    vmon_discharge_pos_time,vmon_avg_discharge_neg_raw,...
    vmon_stdev_discharge_neg_raw,vmon_weight_discharge_neg_raw,...
    vmon_discharge_neg_time,max_length_discharge,max_length_charge,...
    charge_pos_index,charge_neg_index,discharge_pos_index,...
    discharge_neg_index,discharge_index,charge_index,vmon_avg_zero_raw,...
    vmon_stdev_zero_raw,vmon_weight_avg_zero_raw,...
    zero_time,num_zero_chunks,num_zero_points,zero_chunk_array] = ...
    find_charging_data(num_trash_chunks,trash_chunk_array,lcm1_avg_trash_stdev,...
    lcm1_avg_trash,lcm1_avg_trash_raw,...
    lcm1_avg_trash_stdev_raw,...
    lcm1_weight_avg_trash_raw,time_trash);


lcm1_avg_zero = lcm1_avg_zero_raw*lcm1_avg_scale;
lcm1_stdev_zero = lcm1_stdev_zero_raw*abs(lcm1_avg_scale);

[num_ramp_up_inc_points,num_ramp_down_inc_points,...
    up_chunk_inc_array,down_chunk_inc_array,lcm1_avg_ramp_up_inc_raw,...
    lcm1_weight_ramp_down_inc_raw,lcm1_weight_ramp_up_inc_raw,...
    lcm1_avg_ramp_down_inc_raw,time_ramp_up_inc,time_ramp_down_inc,...
    lcm1_avg_charge_neg_inc_sum_current_raw,lcm1_avg_charge_pos_inc_sum_current_raw,...
    lcm1_avg_discharge_neg_inc_sum_current_raw,lcm1_avg_discharge_pos_inc_sum_current_raw,...
    lcm1_avg_discharge_neg_inc_sum_weight_raw,lcm1_avg_discharge_pos_inc_sum_weight_raw] = ...
    include_charging_data(num_rows,num_up_chunks,num_down_chunks,up_chunk_array,...
    down_chunk_array,lcm1_discharge_pos_time,lcm1_discharge_neg_time,...
    lcm1_weight_discharge_pos_raw,lcm1_weight_discharge_neg_raw,...
    lcm1_avg_discharge_pos_raw,lcm1_avg_discharge_neg_raw,charge_pos_index,...
    discharge_pos_index,charge_neg_index,...
    discharge_neg_index,time_ramp_up,time_ramp_down,...
    lcm1_weight_avg_ramp_up_raw,lcm1_weight_avg_ramp_down_raw,lcm1_avg_ramp_down_raw,...
    lcm1_avg_ramp_up_raw,lcm1_charge_pos_time,...
    lcm1_charge_neg_time,lcm1_weight_charge_pos_raw,lcm1_weight_charge_neg_raw,...
    lcm1_avg_charge_pos_raw,lcm1_avg_charge_neg_raw,...
    max_length_charge,max_length_discharge);


lcm1_avg_ramp_down_inc = lcm1_avg_ramp_down_inc_raw * lcm1_avg_scale;

lcm1_avg_ramp_up_inc = lcm1_avg_ramp_up_inc_raw * lcm1_avg_scale;



% Use function 'chunkify' to average each chunk dataset into a point.
% Calculate average of each chunk, weighted stdev of each chunk, average of
% all the chunks together, stdev of all the chunks together, stdev of the
% average of all the chunks together, and stdev of all the stdevs together
if EDM_sim ==1
    [vmon_avg_ramp_up_avg_chunk, vmon_avg_ramp_up_stdev_chunk,...
        vmon_avg_ramp_up_avg,vmon_avg_ramp_up_total_stdev,...
        vmon_avg_ramp_up_stdev_avg_chunk,vmon_avg_ramp_up_stdev_stdev_chunk,...
         vmon_ramp_up_stdev_of_stdev_of_chunk,vmon_ramp_up_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_up_chunks,...
        num_ramp_up_points,up_chunk_array,vmon_avg_ramp_up_raw,...
        vmon_weight_avg_ramp_up_raw,vmon_avg_scale);

    [lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,...
        lcm1_avg_ramp_up_avg,lcm1_avg_ramp_up_total_stdev,...
        lcm1_avg_ramp_up_stdev_avg_chunk,lcm1_avg_ramp_up_stdev_stdev_chunk,...
        lcm1_ramp_up_stdev_of_stdev_of_chunk,lcm1_ramp_up_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_up_chunks,...
        num_ramp_up_points,up_chunk_array,lcm1_avg_ramp_up_raw,...
        lcm1_weight_avg_ramp_up_raw,lcm1_avg_scale);

    [vmon_avg_ramp_down_avg_chunk, vmon_avg_ramp_down_stdev_chunk,...
        vmon_avg_ramp_down_avg,vmon_avg_ramp_down_total_stdev,...
        vmon_avg_ramp_down_stdev_avg_chunk,vmon_avg_ramp_down_stdev_stdev_chunk,...
        vmon_ramp_down_stdev_of_stdev_of_chunk,vmon_ramp_down_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_down_chunks,...
        num_ramp_down_points,down_chunk_array,vmon_avg_ramp_down_raw,...
        vmon_weight_avg_ramp_down_raw,vmon_avg_scale);

    [lcm1_avg_ramp_down_avg_chunk, lcm1_avg_ramp_down_stdev_chunk,...
        lcm1_avg_ramp_down_avg,lcm1_avg_ramp_down_total_stdev,...
        lcm1_avg_ramp_down_stdev_avg_chunk,lcm1_avg_ramp_down_stdev_stdev_chunk,...
        lcm1_ramp_down_stdev_of_stdev_of_chunk,lcm1_ramp_down_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_down_chunks,...
        num_ramp_down_points,down_chunk_array,lcm1_avg_ramp_down_raw,...
        lcm1_weight_avg_ramp_down_raw,lcm1_avg_scale);

    [vmon_avg_trash_avg_chunk, vmon_avg_trash_stdev_chunk,...
        vmon_avg_trash_avg,vmon_avg_trash_total_stdev,...
        vmon_avg_trash_stdev_avg_chunk,vmon_avg_trash_stdev_stdev_chunk,...
        vmon_trash_stdev_of_stdev_of_chunk,vmon_trash_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_trash_chunks,...
        num_trash_points,trash_chunk_array,vmon_avg_trash_raw,...
        vmon_weight_avg_trash_raw,vmon_avg_scale);

    [lcm1_avg_trash_avg_chunk, lcm1_avg_trash_stdev_chunk,...
        lcm1_avg_trash_avg,lcm1_avg_trash_total_stdev,...
        lcm1_avg_trash_stdev_avg_chunk,lcm1_avg_trash_stdev_stdev_chunk,...
        lcm1_trash_stdev_of_stdev_of_chunk,lcm1_trash_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_trash_chunks,...
        num_trash_points,trash_chunk_array,lcm1_avg_trash_raw,...
        lcm1_weight_avg_trash_raw,lcm1_avg_scale);
    
    [lcm1_avg_zero_avg_chunk, lcm1_avg_zero_stdev_chunk,...
        lcm1_avg_zero_avg,lcm1_avg_zero_total_stdev,...
        lcm1_avg_zero_stdev_avg_chunk,lcm1_avg_zero_stdev_stdev_chunk,...
        lcm1_zero_stdev_of_stdev_of_chunk,lcm1_zero_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_zero_chunks,...
        num_zero_points,zero_chunk_array,lcm1_avg_zero_raw,...
        lcm1_weight_avg_zero_raw,lcm1_avg_scale);
    
    [vmon_avg_zero_avg_chunk, vmon_avg_zero_stdev_chunk,...
        vmon_avg_zero_avg,vmon_avg_zero_total_stdev,...
        vmon_avg_zero_stdev_avg_chunk,vmon_avg_zero_stdev_stdev_chunk,...
        vmon_zero_stdev_of_stdev_of_chunk,vmon_zero_avg_stdev_of_chunk] ...
        = chunkify(num_files,num_zero_chunks,...
        num_zero_points,zero_chunk_array,vmon_avg_zero_raw,...
        vmon_weight_avg_zero_raw,vmon_avg_scale);
    
    if inclusive_data == 1

        [lcm1_avg_ramp_down_inc_chunk, lcm1_avg_ramp_down_inc_stdev_chunk,...
            lcm1_avg_ramp_down_inc_avg,lcm1_avg_ramp_down_inc_total_stdev,...
            lcm1_avg_ramp_down_inc_stdev_avg_chunk,lcm1_avg_ramp_down_inc_stdev_stdev_chunk,...
            lcm1_ramp_down_inc_stdev_of_stdev_of_chunk,lcm1_ramp_down_inc_avg_stdev_of_chunk] ...
            = chunkify(num_files,num_down_chunks,...
            num_ramp_down_inc_points,down_chunk_inc_array,lcm1_avg_ramp_down_inc_raw,...
            lcm1_weight_ramp_down_inc_raw,lcm1_avg_scale);

        [lcm1_avg_ramp_up_inc_chunk, lcm1_avg_ramp_up_inc_stdev_chunk,...
            lcm1_avg_ramp_up_inc_avg,lcm1_avg_ramp_up_inc_total_stdev,...
            lcm1_avg_ramp_up_inc_stdev_avg_chunk,lcm1_avg_ramp_up_inc_stdev_stdev_chunk,...
            lcm1_ramp_up_inc_stdev_of_stdev_of_chunk,lcm1_ramp_up_inc_avg_stdev_of_chunk] ...
            = chunkify(num_files,num_up_chunks,...
            num_ramp_up_inc_points,up_chunk_inc_array,lcm1_avg_ramp_up_inc_raw,...
            lcm1_weight_ramp_up_inc_raw,lcm1_avg_scale);
    end
        
end

% % % delete this shite
% % % gaus_a = zeros(num_files,1);
% % % gaus_range = zeros(num_files,1);
% % % num_fit_points = zeros(num_files,1);
% % % 
% % % for i = 1:num_files
% % %     
% % %     num_fit_points(i) = 2*num_ramp_up_points(i);
% % %     
% % %     gaus_a(i) = 1/(abs(lcm1_avg_ramp_up_total_stdev(i))*sqrt(2*pi));
% % %     
% % %     gaus_range(i) = 12*abs(lcm1_avg_ramp_up_total_stdev(i));
% % %     
% % %     step = gaus_range(i)/(num_fit_points(i));
% % %     
% % %     gaus_fit = zeros(num_ramp_up_points(i),2);
% % %     
% % %     for j = 1:num_fit_points(i)
% % %         
% % %         fit_x = j*step + lcm1_avg_ramp_up_avg(i) - 0.5*gaus_range(i);
% % %         
% % %         gaus_fit(j,1) = fit_x;
% % %         
% % %         gaus_fit(j,2) = gaus_a(i)*exp(-(fit_x - ...
% % %             lcm1_avg_ramp_up_avg(i))^2/(2*lcm1_avg_ramp_up_total_stdev(i)^2));
% % %         
% % %     end
% % %     
% % % end


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
        
        fit_params(i,:) = array_linfit(current_source_avg(i,:),...
            lcm1_avg_raw(i,:),lcm1_avg_stdev_raw(i,:),numpoints(i));
        
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
        
            leakage_residual(i,j,1) = (lcm1_avg(i,j)*1e-3 - ...
                (a(i) + b(i)*current_source_avg(i,j)));
        
        end
        
    end
    
end

% compare charging in +V dirxn vs. -V dirxn
% sum the leakage current when charging top electrode in the +V dirxn
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
           
           lcm1_avg_trash_charge_pos(i,end+1) = ...
               sum(lcm1_avg_trash(i,chunk_begin:chunk_end) - lcm1_avg_offset(i));
       
           lcm1_avg_trash_stdev_charge_pos(i,end+1) = ...
               sum(lcm1_avg_trash_stdev(i,chunk_begin:chunk_end));
     
       elseif mod(j,2) == 0
       
           lcm1_avg_trash_charge_neg(i,end+1) = ...
               sum(lcm1_avg_trash(i,chunk_begin:chunk_end) - lcm1_avg_offset(i));
         
           lcm1_avg_trash_stdev_charge_neg(i,end+1) = ...
               sum(lcm1_avg_trash_stdev(i,chunk_begin:chunk_end));
    
       end
       
    end
    
end


[ramp_up_avg_discharge_times, ramp_up_avg_discharge_vals, ramp_up_avg_discharge_stdevs,...
    ramp_up_avg_count_discharges,ramp_up_avg_dph,ramp_up_avg_dph_std] = ...
    find_discharges_final(lcm1_avg_ramp_up_raw,...
    lcm1_weight_avg_ramp_up_raw,time_ramp_up,lcm1_avg_ramp_up_avg_chunk,...
    lcm1_avg_ramp_up_stdev_chunk,num_up_chunks,...
    up_chunk_array,5.0,3,lcm1_avg_scale,2,sprintf('averages %.0f kV chunk ',...
    vmon_avg_ramp_up_avg),...
    'neg_voltage_avg_discharges',1);

one_big_up_chunk_array = zeros(1,2,1);
one_big_up_chunk_array(1,1,1) = up_chunk_array(1,1,1);
one_big_up_chunk_array(1,2,1) = up_chunk_array(1,2,end);

[ramp_up_avg_discharge_times_sum, ramp_up_avg_discharge_vals_sum, ramp_up_avg_discharge_stdevs_sum,...
    ramp_up_avg_count_discharges_sum,ramp_up_avg_dph_sum,ramp_up_avg_dph_std_sum] = ...
    find_discharges_final(lcm1_avg_ramp_up_raw,...
    lcm1_weight_avg_ramp_up_raw,time_ramp_up,lcm1_avg_ramp_up_avg,...
    lcm1_avg_ramp_up_total_stdev,1,...
    one_big_up_chunk_array,5.0,3,lcm1_avg_scale,2,sprintf('averages at %.0f kV chunk ',...
    vmon_avg_ramp_up_avg),...
    'neg_voltage_avg_discharges_summed',1);

one_big_down_chunk_array = zeros(1,2,1);
one_big_down_chunk_array(1,1,1) = down_chunk_array(1,1,1);
one_big_down_chunk_array(1,2,1) = down_chunk_array(1,2,end);

[ramp_down_avg_discharge_times, ramp_down_avg_discharge_vals, ramp_down_avg_discharge_stdevs,...
    ramp_down_avg_count_discharges,ramp_down_avg_dph,ramp_down_avg_dph_std] = ...
    find_discharges_final(lcm1_avg_ramp_down_raw,...
    lcm1_weight_avg_ramp_down_raw,time_ramp_down,lcm1_avg_ramp_down_avg_chunk,...
    lcm1_avg_ramp_down_stdev_chunk,num_down_chunks,...
    down_chunk_array,5.0,3,lcm1_avg_scale,2,sprintf('averages at %.0f kV chunk ',...
    vmon_avg_ramp_down_avg),...
    'pos_voltage_avg_discharges',1);

[ramp_down_avg_discharge_times_sum, ramp_down_avg_discharge_vals_sum, ...
    ramp_down_avg_discharge_stdevs_sum,...
    ramp_down_avg_count_discharges_sum,ramp_down_avg_dph_sum,ramp_down_avg_dph_std_sum] = ...
    find_discharges_final(lcm1_avg_ramp_down_raw,...
    lcm1_weight_avg_ramp_down_raw,time_ramp_down,lcm1_avg_ramp_down_avg,...
    lcm1_avg_ramp_down_total_stdev,1,...
    one_big_down_chunk_array,5.0,3,lcm1_avg_scale,2,sprintf('averages at %.0f kV chunk ',...
    vmon_avg_ramp_down_avg),...
    'pos_voltage_avg_discharges_summed',1);

one_big_zero_chunk_array = zeros(1,2,1);
one_big_zero_chunk_array(1,1,1) = zero_chunk_array(1,1,1);
one_big_zero_chunk_array(1,2,1) = zero_chunk_array(1,2,end);

[zero_avg_discharge_times, zero_avg_discharge_vals, zero_avg_discharge_stdevs,...
    zero_avg_count_discharges,zero_avg_dph,zero_avg_dph_std] = ...
    find_discharges_final(lcm1_avg_zero_raw,...
    lcm1_weight_avg_zero_raw,zero_time,lcm1_avg_zero_avg_chunk,...
    lcm1_avg_zero_stdev_chunk,num_zero_chunks,...
    zero_chunk_array,5.0,3,lcm1_avg_scale,2,sprintf('averages at %.0f kV chunk ',...
    vmon_avg_zero_avg),...
    'zero_voltage_avg_discharges',1);

[zero_avg_discharge_times_sum, zero_avg_discharge_vals_sum, zero_avg_discharge_stdevs_sum,...
    zero_avg_count_discharges_sum,zero_avg_dph_sum,zero_avg_dph_std_sum] = ...
    find_discharges_final(lcm1_avg_zero_raw,...
    lcm1_weight_avg_zero_raw,zero_time,lcm1_avg_zero_avg,...
    lcm1_avg_zero_total_stdev,1,...
    one_big_zero_chunk_array,5.0,3,lcm1_avg_scale,2,sprintf('averages at %.0f kV chunk ',...
    vmon_avg_zero_avg),...
    'zero_voltage_avg_discharges_summed',1);

[zero_stdev_discharge_times, zero_stdev_discharge_vals, zero_stdev_discharge_stdevs,...
    zero_stdev_count_discharges,zero_stdev_dph,zero_stdev_dph_std] = ...
    find_discharges_final(lcm1_stdev_zero_raw,...
    lcm1_weight_avg_zero_raw,zero_time,lcm1_zero_avg_stdev_of_chunk,...
    lcm1_avg_zero_stdev_stdev_chunk*ones(1,num_zero_chunks),num_zero_chunks,...
    zero_chunk_array,1.0,3,abs(lcm1_avg_scale),1,sprintf('stdevs at %.0f kV chunk ',...
    vmon_avg_zero_avg),...
    'zero_voltage_stdev_discharges',1);

[zero_stdev_discharge_times_sum, zero_stdev_discharge_vals_sum, zero_stdev_discharge_stdevs_sum,...
    zero_stdev_count_discharges_sum,zero_stdev_dph_sum,zero_stdev_dph_std_sum] = ...
    find_discharges_final(lcm1_stdev_zero_raw,...
    lcm1_weight_avg_zero_raw,zero_time,lcm1_avg_zero_total_stdev,...
    lcm1_avg_zero_stdev_stdev_chunk,1,...
    one_big_zero_chunk_array,0.25,3,abs(lcm1_avg_scale),1,sprintf('stdevs at %.0f kV chunk ',...
    vmon_avg_zero_avg),...
    'zero_voltage_stdev_discharges_summed',1);

[ramp_down_stdev_discharge_times, ramp_down_stdev_discharge_vals, ramp_down_stdev_discharge_stdevs,...
    ramp_down_stdev_count_discharges,ramp_down_stdev_dph,ramp_down_stdev_dph_std] = ...
    find_discharges_final(lcm1_avg_ramp_down_stdev_raw,...
    lcm1_weight_avg_ramp_down_raw,time_ramp_down,lcm1_ramp_down_avg_stdev_of_chunk,...
    lcm1_avg_ramp_down_stdev_stdev_chunk*ones(1,num_down_chunks),num_ramp_down_chunks,...
    down_chunk_array,1.0,3,abs(lcm1_avg_scale),1,sprintf('stdevs at %.0f kV chunk ',...
    vmon_avg_ramp_down_avg),...
    'pos_voltage_stdev_discharges',1);

[ramp_down_stdev_discharge_times_sum, ramp_down_stdev_discharge_vals_sum,...
    ramp_down_stdev_discharge_stdevs_sum,...
    ramp_down_stdev_count_discharges_sum,ramp_down_stdev_dph_sum,ramp_down_stdev_dph_std_sum] = ...
    find_discharges_final(lcm1_avg_ramp_down_stdev_raw,...
    lcm1_weight_avg_ramp_down_raw,time_ramp_down,lcm1_avg_ramp_down_total_stdev,...
    lcm1_avg_ramp_down_stdev_stdev_chunk,1,...
    one_big_down_chunk_array,0.25,3,abs(lcm1_avg_scale),1,sprintf('stdevs at %.0f kV chunk ',...
    vmon_avg_ramp_down_avg),...
    'pos_voltage_stdev_discharges_summed',1);

[ramp_up_stdev_discharge_times, ramp_up_stdev_discharge_vals, ramp_up_stdev_discharge_stdevs,...
    ramp_up_stdev_count_discharges,ramp_up_stdev_dph,ramp_up_stdev_dph_std] = ...
    find_discharges_final(lcm1_avg_ramp_up_stdev_raw,...
    lcm1_weight_avg_ramp_up_raw,time_ramp_up,lcm1_ramp_up_avg_stdev_of_chunk,...
    lcm1_avg_ramp_up_stdev_stdev_chunk*ones(1,num_up_chunks),num_ramp_up_chunks,...
    up_chunk_array,1.0,3,abs(lcm1_avg_scale),1,sprintf('stdevs at %.0f kV chunk ',...
    vmon_avg_ramp_up_avg),...
    'neg_voltage_stdev_discharges',1);

[ramp_up_stdev_discharge_times_sum, ramp_up_stdev_discharge_vals_sum,...
    ramp_up_stdev_discharge_stdevs_sum,...
    ramp_up_stdev_count_discharges_sum,ramp_up_stdev_dph_sum,ramp_up_stdev_dph_std_sum] = ...
    find_discharges_final(lcm1_avg_ramp_up_stdev_raw,...
    lcm1_weight_avg_ramp_up_raw,time_ramp_up,lcm1_avg_ramp_up_total_stdev,...
    lcm1_avg_ramp_up_stdev_stdev_chunk,1,...
    one_big_up_chunk_array,0.25,3,abs(lcm1_avg_scale),1,sprintf('stdevs at %.0f kV chunk ',...
    vmon_avg_ramp_up_avg),...
    'neg_voltage_stdev_discharges_summed',1);

% % % delete
% % % lcm1_avg_ramp_up = lcm1_avg_ramp_up_raw * lcm1_avg_scale;
% % % lcm1_avg_ramp_up_inv_weight = 1./sqrt(lcm1_weight_avg_ramp_up_raw) * abs(lcm1_avg_scale);

    
%A = [lcm1_avg_ramp_up_one; lcm1_avg_ramp_up_one_stdev];

% fileID = fopen('lcm1_ramp_data.txt','w');
% fprintf(fileID,'%.5e %.5e\n',A);
% fclose(fileID);


%{
lcm1_avg_charge_neg = lcm1_avg_charge_neg_raw*lcm1_avg_scale;
lcm1_stdev_charge_neg = lcm1_stdev_charge_neg_raw * abs(lcm1_avg_scale);


lcm1_avg_charge_pos = lcm1_avg_charge_pos_raw*lcm1_avg_scale;
lcm1_stdev_charge_pos = lcm1_stdev_charge_pos_raw*abs(lcm1_avg_scale);


lcm1_avg_discharge_pos = lcm1_avg_discharge_pos_raw* lcm1_avg_scale;
lcm1_stdev_discharge_pos = lcm1_stdev_discharge_pos_raw * abs(lcm1_avg_scale);


lcm1_avg_discharge_neg = lcm1_avg_discharge_neg_raw * lcm1_avg_scale;
lcm1_stdev_discharge_neg = lcm1_stdev_discharge_neg_raw * abs(lcm1_avg_scale);
%}

%{
lcm1_avg_charge_neg_norm = zeros(num_files,max_length_charge);

lcm1_avg_charge_pos_norm = zeros(num_files,max_length_charge);

lcm1_avg_discharge_neg_norm = zeros(num_files,max_length_discharge);

lcm1_avg_discharge_pos_norm = zeros(num_files,max_length_discharge);


for i = 1:num_files

    for j = 1:max_length_charge

        lcm1_avg_charge_neg_norm(i,j) = ...
            ((lcm1_avg_charge_neg_sum_current_raw(i,j) - ...
            lcm1_avg_offset_raw(i))/max_length_charge*lcm1_avg_scale);

        lcm1_avg_charge_pos_norm(i,j) = ...
            ((lcm1_avg_charge_pos_sum_current_raw(i,j) - ...
            lcm1_avg_offset_raw(i))/max_length_charge*lcm1_avg_scale);

    end

    for j = 1:max_length_discharge

        lcm1_avg_discharge_neg_norm(i,j) = ...
            ((lcm1_avg_discharge_neg_sum_current_raw(i,j) - ...
            lcm1_avg_offset_raw(i))/max_length_discharge *  lcm1_avg_scale);

        lcm1_avg_discharge_pos_norm(i,j) = ...
            ((lcm1_avg_discharge_pos_sum_current_raw(i,j) - ...
            lcm1_avg_offset_raw(i))/max_length_discharge * lcm1_avg_scale);

    end
    
end
%}

% all avg discharges
hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
    sum(ramp_up_avg_count_discharges),sum(ramp_down_avg_count_discharges), ...
    sum(zero_avg_count_discharges)),...
    1,ramp_up_avg_discharge_times,zeros(1,length(ramp_up_avg_discharge_times)),...
    ramp_up_avg_discharge_vals,ramp_up_avg_discharge_stdevs,ramp_down_avg_discharge_times,...
    zeros(1,length(ramp_down_avg_discharge_times)),...
    ramp_down_avg_discharge_vals,ramp_down_avg_discharge_stdevs,zero_avg_discharge_times, ...
    zeros(1,length(zero_avg_discharge_times)),zero_avg_discharge_vals,...
    zero_avg_discharge_stdevs);

% +/- V avg discharges
hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
    sum(ramp_up_avg_count_discharges),sum(ramp_down_avg_count_discharges),...
    sum(zero_avg_count_discharges)),...
    1,ramp_up_avg_discharge_times,zeros(1,length(ramp_up_avg_discharge_times)),...
    ramp_up_avg_discharge_vals,ramp_up_avg_discharge_stdevs,ramp_down_avg_discharge_times,...
    zeros(1,length(ramp_down_avg_discharge_times)),...
    ramp_down_avg_discharge_vals,ramp_down_avg_discharge_stdevs);

% all stdev discharges
hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
    sum(ramp_up_stdev_count_discharges),sum(ramp_down_stdev_count_discharges),...
    sum(zero_stdev_count_discharges)),...
    1,ramp_up_stdev_discharge_times,zeros(1,length(ramp_up_stdev_discharge_times)),...
    ramp_up_stdev_discharge_vals,zeros(1,length(ramp_up_stdev_discharge_vals)),...
    ramp_down_stdev_discharge_times,...
    zeros(1,length(ramp_down_stdev_discharge_times)),...
    ramp_down_stdev_discharge_vals,zeros(1,length(ramp_down_stdev_discharge_vals)),...
    zero_stdev_discharge_times, ...
    zeros(1,length(zero_stdev_discharge_times)),zero_stdev_discharge_vals,...
    zeros(1,length(zero_stdev_discharge_vals)));

% +/- V stdev discharges
hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
    sum(ramp_up_stdev_count_discharges),sum(ramp_down_stdev_count_discharges),...
    sum(zero_stdev_count_discharges)),...
    1,ramp_up_stdev_discharge_times,zeros(1,length(ramp_up_stdev_discharge_times)),...
    ramp_up_stdev_discharge_vals,zeros(1,length(ramp_up_stdev_discharge_vals)),...
    ramp_down_stdev_discharge_times,...
    zeros(1,length(ramp_down_stdev_discharge_times)),...
    ramp_down_stdev_discharge_vals,zeros(1,length(ramp_down_stdev_discharge_vals)));

hv_plot_xy_errors('',...
    1,1:length(ramp_up_stdev_dph),zeros(1,length(ramp_up_stdev_dph)),...
    ramp_up_stdev_dph,ramp_up_stdev_dph_std,...
    1:length(ramp_down_stdev_dph),zeros(1,length(ramp_down_stdev_dph)),...
    ramp_down_stdev_dph,ramp_down_stdev_dph_std);

% % +V stdev discharges plotted over all +V stdev points
% hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
%     sum(ramp_up_stdev_count_discharges),sum(ramp_down_stdev_count_discharges),...
%     sum(zero_stdev_count_discharges)),...
%     1,time_ramp_down,zeros(1,num_ramp_down_points),lcm1_avg_ramp_down_stdev,...
%     zeros(1,num_ramp_down_points),ramp_down_stdev_discharge_times, ...
%     zeros(1,length(ramp_down_stdev_discharge_times)),ramp_down_stdev_discharge_vals,...
%     zeros(1,length(ramp_down_stdev_discharge_vals)));
% 
% %%% - V stdev discharges plotted over all -V stdev points
% hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
%     sum(ramp_up_stdev_count_discharges),sum(ramp_down_stdev_count_discharges),...
%     sum(zero_stdev_count_discharges)),...
%     1,time_ramp_up,zeros(1,num_ramp_up_points),lcm1_avg_ramp_up_stdev,...
%     zeros(1,num_ramp_up_points),ramp_up_stdev_discharge_times, ...
%     zeros(1,length(ramp_up_stdev_discharge_times)),ramp_up_stdev_discharge_vals,...
%     zeros(1,length(ramp_up_stdev_discharge_vals)));
% 
% %%% 0 V stdev discharges plotted over all 0 V stdev points
% hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
%     sum(ramp_up_stdev_count_discharges),sum(ramp_down_stdev_count_discharges),...
%     sum(zero_stdev_count_discharges)),...
%     1,zero_time,zeros(1,num_zero_points),lcm1_stdev_zero,...
%     zeros(1,num_zero_points),zero_stdev_discharge_times, ...
%     zeros(1,length(zero_stdev_discharge_times)),zero_stdev_discharge_vals,...
%     zeros(1,length(zero_stdev_discharge_vals)));
%
%%%% +V avg discharges plotted over all -V avg points
% hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
%     ramp_up_avg_count_discharges,ramp_down_avg_count_discharges, zero_avg_count_discharges),...
%     1,time_ramp_down,zeros(1,num_ramp_down_points),lcm1_avg_ramp_down,...
%     zeros(1,num_ramp_down_points),ramp_down_avg_discharge_times, ...
%     zeros(1,length(ramp_down_avg_discharge_times)),ramp_down_avg_discharge_vals,...
%     zeros(1,length(ramp_down_avg_discharge_vals)));
% 
%%%% -V avg discharges plotted over all -V avg points
% hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
%     ramp_up_avg_count_discharges,ramp_down_avg_count_discharges, zero_avg_count_discharges),...
%     1,time_ramp_up,zeros(1,num_ramp_up_points),lcm1_avg_ramp_up,...
%     zeros(1,num_ramp_up_points),ramp_up_avg_discharge_times, ...
%     zeros(1,length(ramp_up_avg_discharge_times)),ramp_up_avg_discharge_vals,...
%     zeros(1,length(ramp_up_avg_discharge_vals)));
% 
%
%%% 0 V avg discharges plotted over all 0 V avg points
% hv_plot_xy_errors(sprintf('%d discharges at -20 kV \n %d discharges at +20 kV \n %d discharges at 0 kV',...
%     ramp_up_avg_count_discharges,ramp_down_avg_count_discharges, zero_avg_count_discharges),...
%     1,zero_time,zeros(1,num_zero_points),lcm1_avg_zero,...
%     zeros(1,num_zero_points),zero_avg_discharge_times, ...
%     zeros(1,length(zero_avg_discharge_times)),zero_avg_discharge_vals,...
%     zeros(1,length(zero_avg_discharge_vals)));

