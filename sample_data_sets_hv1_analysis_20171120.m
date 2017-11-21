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
% current precision test. README:
% date: 11-13-2017 (Monday afternoon)
% on shift: Roy Ready, Daniel Coulter
% README author: Roy Ready, Daniel Coulter
% Logbook #8 pg. 145-146
% 
% files:	2017-11-13-151308-hv-sample-set-data-1.txt
% 		2017-11-13-151330-hv-1.txt
% 		2017-11-13-151330-hv-avg-FFT-rootHz-1.txt
% 		2017-11-13-165144-hv-sample-set-data-1.txt
% 		2017-11-13-165214-hv-1.txt
% 		2017-11-13-165214-hv-avg-FFT-rootHz-1.txt
% 
% summary: Today we were testing the sensitivity of our leakage current 
%  measurements. We performed two tests, one with the current source directly 
%  connected to the picoammeter, and one with the current source connected to a 
%  mock setup of our leakage current circuit (current source -> Kapton -> feedthru 
%  -> HV cable -> protxn circuit -> picoammeter). We adjusted the current on the 
%  current source and measured the leakage current for both setups. We used the 
%  following values: 
% 
%  picoammeter: 200 nA fixed range, both channels
%  current source: 210 nA range for all measurements
% 
%  meas. no.	current source (nA)
%  1			 000.00
%  2			-200.00
%  3			-131.07
%  4			-065.54
%  5			-032.77
%  6			-016.38
%  7			-008.19
%  8			-004.10
%  9			-002.05
%  10			-001.02
%  11			-000.51
%  12			-000.26
%  13			-000.13
%  14			-000.06
%  15			-000.03
%  16			-000.02
%  17			-000.01
%  18			 000.00
%  19			+000.01
%  20			+000.02
%  21			+000.03
%  22			+000.06
%  23			+000.13
%  24			+000.26
%  25			+000.51
%  26			+001.02
%  27			+002.05
%  28			+004.10
%  29			+008.19
%  30			+016.38
%  31			+032.77
%  32			+065.54
%  33			+131.07
%  34			+200.00
%  35			 000.00
% 
% 
%  For the mock setup we could not get the background leakage current below 150 pA, 
%  likely due to vibrations from the roughing pump. We tried placing the mock setup 
%  on the cart and above the coffin with various types of padding (books, foam) to 
%  isolate the signal, but neither was significantly better than the other. In the 
%  end, we decided to put the mock setup on the cart, with foam and kimwipes to 
%  cushion the components. 
% 
% comments:
%  * significant vibration at the feedthru-cable interface, possibly due to the 2
%    devices not being bolted in as they normally would be on the vacuum chamber
%  * today's labors have convinced me that we need rubber vibration-damping mats
%    for the HV cart and server rack wheels. 
% 
% 11/21/2017 repeated the leakage sensitivity measurement with the current source
%  set to 2 nA range setting (lower bound accuracy of 2 pA). This program will 
%  plot picoammeter reading v. current input and a straight line will be fit to 
%  remove offset.
% README:
% date: 11/20/2017 (Monday afternoon)
% on shift: Daniel Coulter, Roy Ready (author)
% Logbook #9 pg. 41-45
% 
% files: 2017-11-20-153425-hv-sample-set-data-1.txt - current source calibration file
% 	        * picoammeter: 200 nA fixed, both channels
% 			* current source: 2 nA
%        2017-11-20-155710-hv-sample-set-data-1.txt - mock GND circuit measurement
% 	        * picoammeter: 200 nA fixed, both channels
% 			* current source: 2 nA
% 	   2017-11-20-163402-hv-sample-set-data-1.txt - trip protxn circuit test 1
% 	        * picoammeter: 20 mA fixed, both channels
% 			* current source: 2 mA
% 	   2017-11-20-165910-hv-sample-set-data-1.txt - trip protxn circuit test 2
% 	        * picoammeter: 20 mA fixed, both channels
% 			* current source: 20 mA
% 	   
% summary: performed leakage current sensitivity current similar to 11/13 test 
%  (ELOG entry RaEDM 572) but on a smaller range setting for the current source.
%  The smaller range will allow us to quantify the limit of our leakage current
%  measurements. 
%  
%  We tried to intentionally trip the protxn circuit to measure the current diverted
%  from the picoammeter. We were unsuccessful in our attempt to do so with the 
%  current source. The TVS diodes in the protection circuit have a voltage breakdown
%  that may be higher than the voltage from the current source (10 V < V_source < 20 V).
%  The current source can send up to 200 mA but the picoammeter is only rated for 
%  20 mA max input.
% 
% 
% 
%  For current source calibration and mock GND circuit measurement files, we used
%  the following current inputs:
% 
%  meas. no.	current source (nA)
%  1			 0.0000
%  2			-2.0000
%  3			-1.6384
%  4			-0.8192
%  5			-0.4096
%  6			-0.2048
%  7			-0.1024
%  8			-0.0512
%  9			-0.0256
%  10			-0.0128
%  11			-0.0064
%  12			-0.0032
%  13			-0.0016
%  14			-0.0008
%  15			-0.0004
%  16			-0.0002
%  17			-0.0001
%  18			 0.0000
%  19			+0.0001
%  20			+0.0002
%  21			+0.0004
%  22			+0.0008
%  23			+0.0016
%  24			+0.0032
%  25			+0.0064
%  26			+0.0128
%  27			+0.0256
%  28			+0.0512
%  29			+0.1024
%  30			+0.2048
%  31			+0.4096
%  32			+0.8192
%  33			+1.6384
%  34			+2.0000
%  35			 000.00
%  
%  For trip protxn circuit test 1, we used the following current inputs:
%  
%  meas. no.	current source (mA)
%  1			 0.0000
%  2			+0.0001
%  3			+0.0002
%  4			+0.0004
%  5			+0.0008
%  6			+0.0016
%  7			+0.0032
%  8			+0.0064
%  9			+0.0128
%  10			+0.0256
%  11			+0.0512
%  12			+0.1024
%  13			+0.2048
%  14			+0.4096
%  15			+0.8192
%  16			+1.6384
%  17			+2.0000
%  
%  For trip protxn circuit test 2 , we used the following current inputs:
%  
%  meas. no.	current source (mA)
%  1			 00.000
%  2			+00.512
%  3			+01.024
%  4			+02.048
%  5			+04.096
%  6			+08.192
%  7			+16.384
%  8			+20.000
%  
% comments:
%  * to trip the protxn circuit I believe we'd need to pulse a high (>20 V) voltage
%    across the diodes. The current source was not pulsed, but rather gradually
%    increased. The only evidence we have that the protxn circuit diverts current 
%    during a surge is that the picoammeter has not exploded since we installed 
%    the protxn circuit (it has exploded before that, however).


file_struct = dir('*hv*.txt');
num_files = length(file_struct); 

num_cols = 0;
num_rows = 0;
filenames = cell(num_files,1);
for i = 1:num_files
    set = dlmread(file_struct(i).name);
    if length(set(:,1)) > num_rows
        num_rows = length(set(:,1));
    end
    if length(set(1,:)) > num_cols
        num_cols = length(set(1,:));
    end
    %filenames = [filenames ' ' file_struct(i).name(1:10)];
    filenames{i} = file_struct(i).name(1:10);
end
        
data = zeros(num_files,num_rows,num_cols);    

for i = 1:num_files
    set = dlmread(file_struct(i).name);
    set_rows = length(set(:,1));
    set_cols = length(set(1,:));
    data(i,1:set_rows,1:set_cols) = dlmread(file_struct(i).name);
end


%data = dlmread('2017-02-01-163817-hv-sample-set-data.txt');
%numpoints = length(data);

power_supply = 1; % 0 = unipolar Acopian 
%                   1 = bipolar AK (installed and tested 8-18-2017)

pressure_gauge = 1; % 0 = ion gauge 
%                     1 = all-range gauge (installed 7-18-2017)

sample_size = 10; % # samples taken per measurement

sample_rate = 1 ; % 0 = data saved every 0.02 min (8192 samples / 8 kHz)
                  % 1 = data saved every 0.01 min (8192 samples / 16 kHz)


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
lcm1_avg_raw = zeros(num_files,num_rows,1);
lcm1_avg_stdev = zeros(num_files,num_rows,1);
lcm1_avg_stdev_raw = zeros(num_files,num_rows,1);
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
elseif sample_rate == 1;
    sampling_time = 0.02;
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
            vmon_avg_stdev_raw(i,j) = data(i,j,9)/sqrt(sample_size);
            vmon_weight_raw(i,j) = abs((vmon_avg_stdev_raw(i,j))^(-2));
            imon_avg_stdev_raw(i,j) = data(i,j,10)/sqrt(sample_size);
            imon_weight_raw(i,j) = abs((imon_avg_stdev_raw(i,j))^(-2));
            lcm1_avg_stdev_raw(i,j) = data(i,j,11)/sqrt(sample_size);
            pressure_avg_stdev_raw(i,j) = data(i,j,12);
            lcm1_weight_raw(i,j) = abs(1/(lcm1_avg_stdev_raw(i,j)^2));
            ohm_avg_raw(i,j) = vmon_avg(i,j)/imon_avg(i,j)*1e3; % Mohm
            vmon_avg_wt_raw_mag(i,j) = vmon_avg_raw_mag(i,j)*vmon_weight_raw(i,j);
            imon_avg_wt_raw(i,j) = imon_avg_raw(i,j)*imon_weight_raw(i,j);
            field_avg_raw(i,j) = vmon_avg_raw_mag(i,j)/gap_size(i); % (raw vmon / cm)
        end
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
            vmon_avg_stdev_raw(i,j) = data(i,j,10)/sqrt(sample_size);
            vmon_weight_raw(i,j) = abs((vmon_avg_stdev_raw(i,j))^(-2));
            imon_avg_stdev_raw(i,j) = data(i,j,11)/sqrt(sample_size);
            imon_weight_raw(i,j) = abs((imon_avg_stdev_raw(i,j))^(-2));
            lcm1_avg_stdev_raw(i,j) = data(i,j,12)/sqrt(sample_size);
            pressure_avg_stdev_raw(i,j) = data(i,j,14);
            lcm1_weight_raw(i,j) = abs(1/(lcm1_avg_stdev_raw(i,j)^2));
            ohm_avg_raw(i,j) = vmon_avg(i,j)/imon_avg(i,j)*1e3; % Mohm
            vmon_avg_wt_raw_mag(i,j) = vmon_avg_raw_mag(i,j)*vmon_weight_raw(i,j);
            imon_avg_wt_raw(i,j) = imon_avg_raw(i,j)*imon_weight_raw(i,j);
            field_avg_raw(i,j) = vmon_avg_raw_mag(i,j)/gap_size(i); % (raw vmon / cm)
        end
    end
end

for i = 1:num_files
    for j =1:numpoints(i)
        if polarity_avg_raw(i,j) > 1.5;
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
        
        ohm_avg(i,j) = vmon_avg(i,j)/imon_avg(i,j)*ohm_avg_scale;
        vmon_avg_wt(i,j) = vmon_avg_wt_raw(i,j)*vmon_avg_scale;
        imon_avg_wt(i,j) = imon_avg_wt_raw(i,j)*imon_avg_scale;
        field_avg(i,j) = field_avg_raw(i,j)*vmon_avg_scale; %(kV/cm)
    end
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

%define offset chunks arbitrarily in case we don't have a lot of points
for i =1:num_files
    start_offset(i,1) = 1;
    start_offset(i,2) = 2;
end

%define offset chunks to be 10 points at beginning and end of run
for i = 1:num_files
    time(i,:) = (time_raw(i,:) - time_raw(i,1))/60.0; %time starts at 0. in minutes
    for j = (offset_length+1):numpoints(i)
% look for 100 V bump indicating the HV switch is on. Note: PS resistance 
% fluctuates wildly under 1.1 kV
       if  abs(vmon_avg(i,j)) > 0.1 && numpoints(i)>1000;
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
               numpoints(i) > 1000);
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
lcm1_acc_test_offset_raw = zeros(num_files,1);
lcm1_avg_offset = zeros(num_files,1);
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
ramp_high_voltage = zeros(num_files,1);
ramp_low_voltage = zeros(num_files,1);
ramp_high_devition = zeros(num_files,1);
ramp_low_deviation = zeros(num_files,1);

% avg around 0. To separate HI/LO values from 0 values, ignore voltages
% within pm 1 kV of vmon_avg_avg

for i = 1:num_files    
    vmon_avg_avg(i) = mean(vmon_avg(i,start_point(i):end_point(i)));
    hi_vals = find(vmon_avg(i,:) < vmon_avg_avg(i) - 5);
    lo_vals = find(vmon_avg(i,:) > vmon_avg_avg(i) + 5);
    ramp_high_deviation(i) = std(vmon_avg(i,hi_vals),vmon_weight(i,hi_vals));
    ramp_low_deviation(i) = std(vmon_avg(i,lo_vals),vmon_weight(i,lo_vals));
    ramp_high_voltage(i) = sum(vmon_avg_wt(i,hi_vals))/sum(vmon_weight_raw(i,hi_vals));
    ramp_low_voltage(i) = sum(vmon_avg_wt(i,lo_vals))/sum(vmon_weight_raw(i,lo_vals));    
end

%ramp_deviation = 0.05;

vmon_avg_ramp_up_raw = zeros(num_files,num_rows,1);
vmon_avg_ramp_down_raw = zeros(num_files,num_rows,1);
time_ramp_up = zeros(num_files,num_rows,1);
time_ramp_down = zeros(num_files,num_rows,1);
vmon_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
vmon_weight_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
vmon_weight_avg_ramp_up_raw = zeros(num_files,num_rows,1);
vmon_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);
vmon_weight_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);
vmon_weight_avg_ramp_down_raw = zeros(num_files,num_rows,1);
time_ramp_up_pass = zeros(num_files,num_rows,1);
time_ramp_down_pass = zeros(num_files,num_rows,1);
lcm1_avg_ramp_up_raw = zeros(num_files,num_rows,1);
lcm1_avg_ramp_up_raw_wt = zeros(num_files,num_rows,1);
lcm1_avg_ramp_down_raw = zeros(num_files,num_rows,1);
lcm1_avg_ramp_down_raw_wt = zeros(num_files,num_rows,1);
lcm1_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
lcm1_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);
lcm1_weight_avg_ramp_up_raw_pass = zeros(num_files,num_rows,1);
lcm1_weight_avg_ramp_down_raw_pass = zeros(num_files,num_rows,1);
lcm1_weight_avg_ramp_up_raw = zeros(num_files,num_rows,1);
lcm1_weight_avg_ramp_down_raw = zeros(num_files,num_rows,1);

%store chunk # and index
up_chunk_array = zeros(num_files,2,num_rows);
up_chunk_array_pass = zeros(num_files,2,num_rows); 
down_chunk_array = zeros(num_files,2,num_rows);
down_chunk_array_pass = zeros(num_files,2,num_rows);



up_array_count = zeros(num_files,1);
down_array_count = zeros(num_files,1);
num_ramp_up_points = ones(num_files,1);
num_ramp_down_points = ones(num_files,1);
count_up_chunks = zeros(num_files,1);
count_down_chunks = zeros(num_files,1);

for i = 1:num_files
    count_up_chunks(i) = 0;
    count_down_chunks(i) = 0;
    for j = start_point(i):numpoints(i)           
        if abs(vmon_avg(i,j) - ramp_high_voltage(i)) < 3*ramp_high_deviation(i) && vmon_avg(i,j) < ramp_high_voltage(i) + ramp_high_deviation(i)
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
        elseif abs(vmon_avg(i,j) - ramp_low_voltage(i)) < 3*ramp_low_deviation(i) && vmon_avg(i,j) > ramp_low_voltage(i) - ramp_low_deviation(i)
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
        end
    end
end



num_up_chunks = zeros(num_files,1);
num_down_chunks = zeros(num_files,1);
index_up_chunks = zeros(num_files,1);
index_down_chunks = zeros(num_files,1);

%count the number of ramp up points and ramp down points for HI/LO sets
%for i = 1:num_files
%    for j = start_point(i):numpoints(i)
%        if vmon_avg_ramp_up_raw_pass(i,j) == 0.0
%            num_ramp_up_points(i) = j;
%            break;
%        end
%    end
%end

%for i = 1:num_files
%    for j = start_point(i):numpoints(i)
%        if vmon_avg_ramp_down_raw_pass(i,j) == 0.0
%            num_ramp_down_points(i) = j;
%            break;
%        end
%    end
%end


for i = 1:num_files
    %count up to 2nd to last up chunk
    index_up_chunks(i) = count_up_chunks(i) - 1;
    num_up_chunks(i) = index_up_chunks(i) - 1;
    %count to 3rd to last down chunk
    index_down_chunks(i) = count_down_chunks(i) - 3;
    num_down_chunks(i) = index_down_chunks(i) - 1;
end

num_down_chunk_rows = 0;
num_up_chunk_rows = 0;
for i =1:num_files
    if num_up_chunks(i) > num_up_chunk_rows;
        num_up_chunk_rows = num_up_chunks(i);
    end
    if num_down_chunks(i) > num_down_chunk_rows;
        num_down_chunk_rows = num_down_chunks(i);
    end
end
    

for i = 1:num_files
    for j = 1:num_ramp_up_points(i)
        lcm1_avg_ramp_up_raw(i,j) = lcm1_avg_ramp_up_raw_pass(i,j+1);
        lcm1_weight_avg_ramp_up_raw(i,j) = lcm1_weight_avg_ramp_up_raw_pass(i,j+1);
        vmon_avg_ramp_up_raw(i,j) = vmon_avg_ramp_up_raw_pass(i,j+1);
        vmon_weight_avg_ramp_up_raw(i,j) = vmon_weight_avg_ramp_up_raw_pass(i,j+1);
        time_ramp_up(i,j) = time_ramp_up_pass(i,j+1);
    end
end

for i = 1:num_files
    for j = 1:num_ramp_down_points(i)
        lcm1_avg_ramp_down_raw(i,j) = lcm1_avg_ramp_down_raw_pass(i,j+1);
        lcm1_weight_avg_ramp_down_raw(i,j) = lcm1_weight_avg_ramp_down_raw_pass(i,j+1);
        vmon_avg_ramp_down_raw(i,j) = vmon_avg_ramp_down_raw_pass(i,j+1);
        vmon_weight_avg_ramp_down_raw(i,j) = vmon_weight_avg_ramp_down_raw_pass(i,j+1);
        time_ramp_down(i,j) = time_ramp_down_pass(i,j+1);
    end
end

for i =1:num_files
    %start at 2nd up chunk
    for j =1:index_up_chunks(i)
        up_chunk_array(i,1,j) = up_chunk_array_pass(i,1,j+1);
        up_chunk_array(i,2,j) = up_chunk_array_pass(i,2,j+1);
    end
end

for i =1:num_files
    %start at 2nd down chunk
    for j =1:index_down_chunks(i)
        down_chunk_array(i,1,j) = down_chunk_array_pass(i,1,j+2);
        down_chunk_array(i,2,j) = down_chunk_array_pass(i,2,j+2);
    end
end





%7/11/2017 can't think of a way to store chunk averaging for more than 1 data set

%calculate averages
vmon_avg_ramp_up_raw_avg = zeros(num_files,1);
vmon_avg_ramp_up_raw_avg_chunk = zeros(num_files,num_up_chunk_rows,1);
vmon_avg_ramp_down_raw_avg = zeros(num_files,1);
vmon_avg_ramp_down_raw_avg_chunk = zeros(num_files,num_down_chunk_rows,1);
vmon_avg_ramp_up_raw_stdev = zeros(num_files,1);
vmon_avg_ramp_down_raw_stdev = zeros(num_files,1);
lcm1_avg_ramp_up_raw_avg = zeros(num_files,1);
lcm1_avg_ramp_up_raw_avg_chunk = zeros(num_files,num_up_chunk_rows,1);
lcm1_avg_ramp_down_raw_avg = zeros(num_files,1);
lcm1_avg_ramp_down_raw_avg_chunk = zeros(num_files,num_down_chunk_rows,1);
lcm1_avg_ramp_up_raw_stdev = zeros(num_files,1);
lcm1_avg_ramp_down_raw_stdev = zeros(num_files,1);


vmon_avg_ramp_up_avg = zeros(num_files,1);
vmon_avg_ramp_up_avg_chunk = zeros(num_files,num_up_chunk_rows,1);
vmon_avg_ramp_down_avg = zeros(num_files,1);
vmon_avg_ramp_down_avg_chunk = zeros(num_files,num_down_chunk_rows,1);
vmon_avg_ramp_up_stdev = zeros(num_files,1);
vmon_avg_ramp_up_stdev_chunk = zeros(num_files,num_up_chunk_rows,1);
vmon_avg_ramp_down_stdev = zeros(num_files,1);
vmon_avg_ramp_down_stdev_chunk = zeros(num_files,num_down_chunk_rows,1);
lcm1_avg_ramp_up_avg = zeros(num_files,1);
lcm1_avg_ramp_up_avg_chunk = zeros(num_files,num_up_chunk_rows,1);
lcm1_avg_ramp_down_avg = zeros(num_files,1);
lcm1_avg_ramp_down_avg_chunk = zeros(num_files,1);
lcm1_avg_ramp_up_stdev = zeros(num_files,1);
lcm1_avg_ramp_up_stdev_chunk = zeros(num_files,num_up_chunk_rows,1);
lcm1_avg_ramp_up_chunk = zeros(num_files,1);
lcm1_avg_ramp_down_stdev = zeros(num_files,1);
lcm1_avg_ramp_down_stdev_chunk = zeros(num_files,1);

vmon_avg_ramp_up_raw_avg_raw = zeros(num_files,1);
vmon_avg_ramp_up_raw_avg_chunk_raw = zeros(num_files,num_up_chunk_rows, 1);
vmon_avg_ramp_down_raw_avg_raw = zeros(num_files,1);
vmon_avg_ramp_down_raw_avg_chunk_raw = zeros(num_files,num_down_chunk_rows, 1);
vmon_avg_ramp_up_raw_stdev_raw = zeros(num_files,1);
vmon_avg_ramp_up_raw_stdev_raw_chunk = zeros(num_files,num_up_chunk_rows,1);
vmon_avg_ramp_down_raw_stdev_raw = zeros(num_files,1);
vmon_avg_ramp_down_raw_stdev_raw_chunk = zeros(num_files,num_down_chunk_rows,1);
lcm1_avg_ramp_up_raw_avg_raw = zeros(num_files,1);
lcm1_avg_ramp_up_raw_avg_chunk_raw = zeros(num_files,num_up_chunk_rows,1);
lcm1_avg_ramp_down_raw_avg_raw = zeros(num_files,1);
lcm1_avg_ramp_down_raw_avg_chunk_raw = zeros(num_files,num_down_chunk_rows,1);
lcm1_avg_ramp_up_raw_stdev_raw = zeros(num_files,1);
lcm1_avg_ramp_up_raw_stdev_raw_chunk = zeros(num_files,num_up_chunk_rows,1);
lcm1_avg_ramp_down_raw_stdev_raw = zeros(num_files,1);
lcm1_avg_ramp_down_raw_stdev_raw_chunk = zeros(num_files,num_down_chunk_rows,1);

lcm1_avg_ramp_down_raw_avg_wt = zeros(num_files,1);
vmon_avg_ramp_down_raw_avg_wt = zeros(num_files,1);
vmon_avg_ramp_up_raw_avg_wt = zeros(num_files,1);
lcm1_avg_ramp_up_raw_avg_wt = zeros(num_files,1);

for i = 1:num_files
    %vmon_avg_ramp_up_raw_avg(i) = mean(vmon_avg_ramp_up_raw(i,:));
    %vmon_avg_ramp_down_raw_avg(i) = mean(vmon_avg_ramp_down_raw(i,:));
    vmon_avg_ramp_up_raw_stdev(i) = std(vmon_avg_ramp_up_raw(i,:),vmon_weight_avg_ramp_up_raw(i,:));
    vmon_avg_ramp_down_raw_stdev(i) = std(vmon_avg_ramp_down_raw(i,:),vmon_weight_avg_ramp_down_raw(i,:));
    lcm1_avg_ramp_up_raw_stdev(i) = std(lcm1_avg_ramp_up_raw(i,:),lcm1_weight_avg_ramp_up_raw(i,:));
    lcm1_avg_ramp_down_raw_stdev(i) = std(lcm1_avg_ramp_down_raw(i,:),lcm1_weight_avg_ramp_down_raw(i,:));

end

for i = 1:num_files
    for j = 1:num_ramp_up_points(i)
        vmon_avg_ramp_up_raw_avg_wt(i,j) = vmon_avg_ramp_up_raw(i,j)*vmon_weight_avg_ramp_up_raw(i,j);
        lcm1_avg_ramp_up_raw_avg_wt(i,j) = lcm1_avg_ramp_up_raw(i,j)*lcm1_weight_avg_ramp_up_raw(i,j);
                
    end
end

for i = 1:num_files
    for j = 1:num_ramp_down_points(i)
        
        lcm1_avg_ramp_down_raw_avg_wt(i,j) = lcm1_avg_ramp_down_raw(i,j)*lcm1_weight_avg_ramp_down_raw(i,j);
        vmon_avg_ramp_down_raw_avg_wt(i,j) = vmon_avg_ramp_down_raw(i,j)*vmon_weight_avg_ramp_down_raw(i,j);
                
    end
end

for i = 1:num_files
    lcm1_avg_ramp_up_raw_avg(i) = sum(lcm1_avg_ramp_up_raw_avg_wt(i,:)) / sum(lcm1_weight_avg_ramp_up_raw(i,:));
    lcm1_avg_ramp_down_raw_avg(i) = sum(lcm1_avg_ramp_down_raw_avg_wt(i,:)) / sum(lcm1_weight_avg_ramp_down_raw(i,:));
    vmon_avg_ramp_down_raw_avg(i) = sum(vmon_avg_ramp_down_raw_avg_wt(i,:))/sum(vmon_weight_avg_ramp_down_raw(i,:));
    vmon_avg_ramp_up_raw_avg(i) = sum(vmon_avg_ramp_up_raw_avg_wt(i,:))/sum(vmon_weight_avg_ramp_up_raw(i,:));
end

vmon_weight_ramp_down_raw_stdev_raw_chunk = zeros(num_files,num_down_chunk_rows);
lcm1_weight_ramp_down_raw_stdev_raw_chunk = zeros(num_files,num_down_chunk_rows);
vmon_avg_ramp_down_raw_stdev_avg_raw_chunk_wt = zeros(num_files,num_down_chunk_rows,1);
lcm1_avg_ramp_down_raw_stdev_avg_raw_chunk_wt = zeros(num_files,num_down_chunk_rows,1);
vmon_weight_ramp_up_raw_stdev_raw_chunk = zeros(num_files,num_up_chunk_rows);
lcm1_weight_ramp_up_raw_stdev_raw_chunk = zeros(num_files,num_up_chunk_rows);
vmon_avg_ramp_up_raw_stdev_avg_raw_chunk_wt = zeros(num_files,num_up_chunk_rows,1);
lcm1_avg_ramp_up_raw_stdev_avg_raw_chunk_wt = zeros(num_files,num_up_chunk_rows,1);

vmon_avg_ramp_down_raw_stdev_stdev_raw_chunk = zeros(num_files,1);
vmon_avg_ramp_down_stdev_stdev_chunk = zeros(num_files,1);
lcm1_avg_ramp_down_raw_stdev_stdev_raw_chunk = zeros(num_files,1);
lcm1_avg_ramp_down_stdev_stdev_chunk = zeros(num_files,1);
vmon_avg_ramp_up_raw_stdev_stdev_raw_chunk = zeros(num_files,1);
vmon_avg_ramp_up_stdev_stdev_chunk = zeros(num_files,1);
lcm1_avg_ramp_up_raw_stdev_stdev_raw_chunk = zeros(num_files,1);
lcm1_avg_ramp_up_stdev_stdev_chunk = zeros(num_files,1);

vmon_avg_ramp_down_raw_stdev_avg_raw_chunk = zeros(num_files,1);
vmon_avg_ramp_down_stdev_avg_chunk = zeros(num_files,1);
lcm1_avg_ramp_down_raw_stdev_avg_raw_chunk = zeros(num_files,1);
lcm1_avg_ramp_down_stdev_avg_chunk = zeros(num_files,1);
vmon_avg_ramp_up_raw_stdev_avg_raw_chunk = zeros(num_files,1);
vmon_avg_ramp_up_stdev_avg_chunk = zeros(num_files,1);
lcm1_avg_ramp_up_raw_stdev_avg_raw_chunk = zeros(num_files,1);
lcm1_avg_ramp_up_stdev_avg_chunk = zeros(num_files,1);

%7/12 change chunk code when I figure out how to handle multiple data sets...
for i = 1:num_files
    for j = 1:num_down_chunks(i)
        chunk_begin = down_chunk_array(i,2,j) +1;
        chunk_end = down_chunk_array(i,2,j+1);
        for k = chunk_begin:chunk_end
            vmon_avg_ramp_down_raw_stdev_raw_chunk(i,j) =( std(vmon_avg_ramp_down_raw(i,chunk_begin:chunk_end),...
                vmon_weight_avg_ramp_down_raw(i,chunk_begin:chunk_end)));
            vmon_weight_ramp_down_raw_stdev_raw_chunk(i,j) = (vmon_avg_ramp_down_raw_stdev_raw_chunk(i,j))^(-2);
            lcm1_avg_ramp_down_raw_stdev_raw_chunk(i,j) =( std(lcm1_avg_ramp_down_raw(i,chunk_begin:chunk_end),...
                lcm1_weight_avg_ramp_down_raw(i,chunk_begin:chunk_end)));
            lcm1_weight_ramp_down_raw_stdev_raw_chunk(i,j) = (lcm1_avg_ramp_down_raw_stdev_raw_chunk(i,j))^(-2);
        end
        vmon_avg_ramp_down_raw_stdev_avg_raw_chunk_wt(i,j) = vmon_avg_ramp_down_raw_stdev_raw_chunk(i,j)*vmon_weight_ramp_down_raw_stdev_raw_chunk(i,j);
        lcm1_avg_ramp_down_raw_stdev_avg_raw_chunk_wt(i,j) = lcm1_avg_ramp_down_raw_stdev_raw_chunk(i,j)*lcm1_weight_ramp_down_raw_stdev_raw_chunk(i,j);
    end
end

    

for i = 1:num_files
        vmon_avg_ramp_down_raw_stdev_stdev_raw_chunk(i) =1/sqrt(sum(vmon_weight_ramp_down_raw_stdev_raw_chunk(:)));
        vmon_avg_ramp_down_raw_stdev_avg_raw_chunk(i) = sum(vmon_avg_ramp_down_raw_stdev_avg_raw_chunk_wt(:))/sum(vmon_weight_ramp_down_raw_stdev_raw_chunk(:));
        lcm1_avg_ramp_down_raw_stdev_stdev_raw_chunk(i) =1/sqrt(sum(lcm1_weight_ramp_down_raw_stdev_raw_chunk(:)));
        lcm1_avg_ramp_down_raw_stdev_avg_raw_chunk(i) = sum(lcm1_avg_ramp_down_raw_stdev_avg_raw_chunk_wt(:))/sum(lcm1_weight_ramp_down_raw_stdev_raw_chunk(:));
end

for i = 1:num_files
    for j = 1:num_up_chunks(i)
        chunk_begin = up_chunk_array(i,2,j) +1;
        chunk_end = up_chunk_array(i,2,j+1);
        for k = chunk_begin:chunk_end
            vmon_avg_ramp_up_raw_stdev_raw_chunk(i,j) =( std(vmon_avg_ramp_up_raw(i,chunk_begin:chunk_end),...
                vmon_weight_avg_ramp_up_raw(i,chunk_begin:chunk_end)));
            vmon_weight_ramp_up_raw_stdev_raw_chunk(i,j) = (vmon_avg_ramp_up_raw_stdev_raw_chunk(i,j))^(-2);
            lcm1_avg_ramp_up_raw_stdev_raw_chunk(i,j) =( std(lcm1_avg_ramp_up_raw(i,chunk_begin:chunk_end),...
                lcm1_weight_avg_ramp_up_raw(i,chunk_begin:chunk_end)));
            lcm1_weight_ramp_up_raw_stdev_raw_chunk(i,j) = (lcm1_avg_ramp_up_raw_stdev_raw_chunk(i,j))^(-2);
        end
        vmon_avg_ramp_up_raw_stdev_avg_raw_chunk_wt(i,j) = vmon_avg_ramp_up_raw_stdev_raw_chunk(i,j)*vmon_weight_ramp_up_raw_stdev_raw_chunk(i,j);
        lcm1_avg_ramp_up_raw_stdev_avg_raw_chunk_wt(i,j) = lcm1_avg_ramp_up_raw_stdev_raw_chunk(i,j)*lcm1_weight_ramp_up_raw_stdev_raw_chunk(i,j);
    end
end

    

for i = 1:num_files
        vmon_avg_ramp_up_raw_stdev_stdev_raw_chunk(i) =1/sqrt(sum(vmon_weight_ramp_up_raw_stdev_raw_chunk(i,:)));
        vmon_avg_ramp_up_raw_stdev_avg_raw_chunk(i) = sum(vmon_avg_ramp_up_raw_stdev_avg_raw_chunk_wt(i,:))/sum(vmon_weight_ramp_up_raw_stdev_raw_chunk(i,:));
        lcm1_avg_ramp_up_raw_stdev_stdev_raw_chunk(i) =1/sqrt(sum(lcm1_weight_ramp_up_raw_stdev_raw_chunk(i,:)));
        lcm1_avg_ramp_up_raw_stdev_avg_raw_chunk(i) = sum(lcm1_avg_ramp_up_raw_stdev_avg_raw_chunk_wt(i,:))/sum(lcm1_weight_ramp_up_raw_stdev_raw_chunk(i,:));
end


for i = 1:num_files
    for j = 1:num_down_chunks(i)
        chunk_begin = down_chunk_array(i,2,j) +1;
        chunk_end = down_chunk_array(i,2,j+1);
        vmon_avg_ramp_down_raw_avg_chunk(i,j) =( sum(vmon_avg_ramp_down_raw_avg_wt(i,...
            chunk_begin:chunk_end))/sum(vmon_weight_avg_ramp_down_raw(i,chunk_begin:chunk_end)));
        lcm1_avg_ramp_down_raw_avg_chunk(i,j) =( sum(lcm1_avg_ramp_down_raw_avg_wt(i,...
            chunk_begin:chunk_end))/sum(lcm1_weight_avg_ramp_down_raw(i,chunk_begin:chunk_end)));
    end
end

for i = 1:num_files
    for j = 1:num_up_chunks(i)
        chunk_begin = up_chunk_array(i,2,j) +1;
        chunk_end = up_chunk_array(i,2,j+1);
        for k = chunk_begin:chunk_end
            vmon_avg_ramp_up_raw_stdev_raw_chunk(i,j) =( std(vmon_avg_ramp_up_raw(i,chunk_begin:chunk_end),...
                vmon_weight_avg_ramp_up_raw(i,chunk_begin:chunk_end)));
            lcm1_avg_ramp_up_raw_stdev_raw_chunk(i,j) =( std(lcm1_avg_ramp_up_raw(i,chunk_begin:chunk_end),...
                lcm1_weight_avg_ramp_up_raw(i,chunk_begin:chunk_end)));
        end
    end
end

for i = 1:num_files
    for j = 1:num_up_chunks(i)
        chunk_begin = up_chunk_array(i,2,j) +1;
        chunk_end = up_chunk_array(i,2,j+1);
        vmon_avg_ramp_up_raw_avg_chunk(i,j) =( sum(vmon_avg_ramp_up_raw_avg_wt(i,...
            chunk_begin:chunk_end))/sum(vmon_weight_avg_ramp_up_raw(i,chunk_begin:chunk_end)));
        lcm1_avg_ramp_up_raw_avg_chunk(i,j) =( sum(lcm1_avg_ramp_up_raw_avg_wt(i,...
            chunk_begin:chunk_end))/sum(lcm1_weight_avg_ramp_up_raw(i,chunk_begin:chunk_end)));
    end
end

for i = 1:num_files
    %overall chunk averages and stdevs of avgs
    vmon_avg_ramp_up_avg(i) = vmon_avg_ramp_up_raw_avg(i)*vmon_avg_scale;
    vmon_avg_ramp_down_avg(i) = vmon_avg_ramp_down_raw_avg(i)*vmon_avg_scale;
    vmon_avg_ramp_up_stdev(i) = vmon_avg_ramp_up_raw_stdev(i)*vmon_avg_scale;
    vmon_avg_ramp_down_stdev(i) =  vmon_avg_ramp_down_raw_stdev(i)*vmon_avg_scale;
    lcm1_avg_ramp_up_avg(i) = lcm1_avg_ramp_up_raw_avg(i)*lcm1_avg_scale;
    lcm1_avg_ramp_down_avg(i) = lcm1_avg_ramp_down_raw_avg(i)*lcm1_avg_scale;
    lcm1_avg_ramp_up_stdev(i) = lcm1_avg_ramp_up_raw_stdev(i)*abs(lcm1_avg_scale);
    lcm1_avg_ramp_down_stdev(i) = lcm1_avg_ramp_down_raw_stdev(i)*abs(lcm1_avg_scale);
    
    %overall chunk averages of stdevs and stdevs of stdevs
    vmon_avg_ramp_up_stdev_stdev_chunk(i) = vmon_avg_ramp_up_raw_stdev_stdev_raw_chunk(i)*abs(vmon_avg_scale);
    vmon_avg_ramp_up_stdev_avg_chunk(i) = vmon_avg_ramp_up_raw_stdev_avg_raw_chunk(i)*abs(vmon_avg_scale);
    lcm1_avg_ramp_up_stdev_stdev_chunk(i) = lcm1_avg_ramp_up_raw_stdev_stdev_raw_chunk(i)*abs(lcm1_avg_scale);
    lcm1_avg_ramp_up_stdev_avg_chunk(i) = lcm1_avg_ramp_up_raw_stdev_avg_raw_chunk(i)*abs(lcm1_avg_scale);
    vmon_avg_ramp_down_stdev_stdev_chunk(i) = vmon_avg_ramp_down_raw_stdev_stdev_raw_chunk(i)*abs(vmon_avg_scale);
    vmon_avg_ramp_down_stdev_avg_chunk(i) = vmon_avg_ramp_down_raw_stdev_avg_raw_chunk(i)*abs(vmon_avg_scale);
    lcm1_avg_ramp_down_stdev_stdev_chunk(i) = lcm1_avg_ramp_down_raw_stdev_stdev_raw_chunk(i)*abs(lcm1_avg_scale);
    lcm1_avg_ramp_down_stdev_avg_chunk(i) = lcm1_avg_ramp_down_raw_stdev_avg_raw_chunk(i)*abs(lcm1_avg_scale);
    

end

for i =1:num_files
    for j = 1:num_down_chunks(i)
        vmon_avg_ramp_down_avg_chunk(i,j) = vmon_avg_ramp_down_raw_avg_chunk(i,j)*vmon_avg_scale;
        vmon_avg_ramp_down_stdev_chunk(i,j) = vmon_avg_ramp_down_raw_stdev_raw_chunk(i,j)*vmon_avg_scale;
        lcm1_avg_ramp_down_avg_chunk(i,j) = lcm1_avg_ramp_down_raw_avg_chunk(i,j)*lcm1_avg_scale;
        lcm1_avg_ramp_down_stdev_chunk(i,j) = lcm1_avg_ramp_down_raw_stdev_raw_chunk(i,j)*abs(lcm1_avg_scale);
    end
end

for i =1:num_files
    for j = 1:num_up_chunks(i)
       vmon_avg_ramp_up_avg_chunk(i,j) = vmon_avg_ramp_up_raw_avg_chunk(i,j)*vmon_avg_scale;
       vmon_avg_ramp_up_stdev_chunk(i,j) = vmon_avg_ramp_up_raw_stdev_raw_chunk(i,j)*vmon_avg_scale;
       lcm1_avg_ramp_up_avg_chunk(i,j) = lcm1_avg_ramp_up_raw_avg_chunk(i,j)*lcm1_avg_scale;
       lcm1_avg_ramp_up_stdev_chunk(i,j) = lcm1_avg_ramp_up_raw_stdev_raw_chunk(i,j)*abs(lcm1_avg_scale);
    end
end

for i = 1:num_files
    for j = 1:num_ramp_up_points(i)
        lcm1_avg_ramp_up(i,j) = lcm1_avg_ramp_up_raw(i,j)*lcm1_avg_scale;
        lcm1_weight_avg_ramp_up(i,j) = lcm1_weight_avg_ramp_up_raw(i,j)*lcm1_avg_scale;
        vmon_avg_ramp_up(i,j) = vmon_avg_ramp_up_raw(i,j)*vmon_avg_scale;
        vmon_weight_avg_ramp_up(i,j) = vmon_weight_avg_ramp_up_raw(i,j)*vmon_avg_scale;
    end
end

for i = 1:num_files
    for j = 1:num_ramp_down_points(i)
        lcm1_avg_ramp_down(i,j) = lcm1_avg_ramp_down_raw(i,j)*lcm1_avg_scale;
        lcm1_weight_avg_ramp_down(i,j) = lcm1_weight_avg_ramp_down_raw(i,j)*lcm1_avg_scale;
        vmon_avg_ramp_down(i,j) = vmon_avg_ramp_down_raw(i,j)*vmon_avg_scale;
        vmon_weight_avg_ramp_down(i,j) = vmon_weight_avg_ramp_down_raw(i,j)*vmon_avg_scale;
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

fit_params = zeros(num_files,4);
%a: offset, offset uncertainty
a = zeros(num_files,2);
%b: slope, slope uncertainty
b = zeros(num_files,2);
%calculate residuals and the line
leakage_residual = zeros(num_files,num_rows,2);
fit_line_x = zeros(num_files,2);
fit_line_y = zeros(num_files,2);

for i =1:num_files
    fit_params(i,:) = array_linfit(current_source_avg(i,:),lcm1_avg_raw(i,:),lcm1_avg_stdev_raw(i,:),numpoints(i));
    %convert y values to nA, slope will be overall unitless
    a(i,1) = fit_params(i,1)*lcm1_avg_scale*1e-3;
    a(i,2) = fit_params(i,2)*abs(lcm1_avg_scale)*1e-3;
    b(i,1) = fit_params(i,3)*lcm1_avg_scale*1e-3;
    b(i,2) = fit_params(i,4)*abs(lcm1_avg_scale)*1e-3;
    fit_line_x(i,1) = floor(min(current_source_avg(i,:)));
    fit_line_x(i,2) = floor(max(current_source_avg(i,:)))+1;
    fit_line_y(i,1) = a(i,1)+b(i,1)*fit_line_x(i,1);
    fit_line_y(i,2) = a(i,1)+b(i,1)*fit_line_x(i,2);
    for j = 1:numpoints(i)
        leakage_residual(i,j,1) = lcm1_avg(i,j)*1e-3 - (a(i) + b(i)*current_source_avg(i,j));
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
    red_list_hsv(i,1) = 0;
    red_list_hsv(i,2) = sat;
    red_list_hsv(i,3) = val;
    red_list_rgb(i,:,:) = hsv2rgb(red_list_hsv(i,:,:));
end
cmap = colormap(red_list_rgb);

grnhsv = [0.33 1 1];     % green HSV
grn_list_hsv = zeros(num_files + 1,3);
grn_list_rgb = zeros(num_files + 1,3);
grn_step = round(0.75/(num_files + 1),4);
for i = 1:num_files+1
    sat = 1 + grn_step*(i - num_files - 1);
    val = 1 + grn_step*(i - num_files - 1)/3.;
    grn_list_hsv(i,1) = 0;
    grn_list_hsv(i,2) = sat;
    grn_list_hsv(i,3) = val;
    grn_list_rgb(i,:,:) = hsv2rgb(grn_list_hsv(i,:,:));
end
cmap = colormap(red_list_rgb);



%titles and shit
title_string = 'leakage current comparison plot';
x_label = 'PS voltage (-kV)';
y_label = 'residual (pA)';
legend_titles = [filenames];
xtick_numbers = [ 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75];
ytick_numbers = [-125 -100 -75 -50 -25 0 25 50 75 100 125];
xmin = 0.0;
xmax = 80;
ymin = -15;
ymax = 25;
plot_bounds = [xmin xmax ymin ymax];


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


%%%%%%%%%%%% %grid of leakage current v. time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %   subplot(num_files,1,num_files + 1 - i)
%    plot(time(i,start_point(i):end_point(i)) - time(i,start_point(i)),...
%        lcm1_avg(i,start_point(i):end_point(i))-lcm1_avg_offset(i),...
%        '-','Color', cmap(i+1,:),'MarkerSize', 3, 'LineWidth', 2.0);
% %   axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
% %   text('Position',[xmax ymax],'String',legend_titles(i),...
% %       'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
%    ax.TickDir = 'out'; % make ticks point out



%%%%%%%%%%%%%%%%% grid of power ps voltage v. time %%%%%%%%%%%%%%%%%%%%%%%%
% %    subplot(num_files,1,num_files + 1 - i)
%    plot(time(i,start_offset(i):end_offset(i)) - time(i,start_offset(i)),...
%        (vmon_avg(i,start_offset(i):end_offset(i)) - vmon_avg_offset(i))/.14,...
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

%%%%%%%%%%%%%% ramp data plot of psvoltage v time %%%%%%%%%%%%%%%%%%%%%%%%%
%    subplot(num_files,1,num_files + 1 - i)
%    plot(time_ramp_up(i,:), vmon_avg_ramp_up(i,:) - vmon_avg_offset(i),...
%    'o','Color', cmap(1+i,:),'MarkerSize', 8, 'LineWidth', 2.0);
%    %axis(plot_bounds)
%    pbaspect([1.33 1 1])
%    ax = gca; % current axes
% %    ax.FontSize = 16;
%    ax.TickDir = 'out'; % make ticks point out



% %%%%%%%%%% ramp data plot of mean up chunk voltage v chunk # %%%%%%%
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
% %%%%%%%%% ramp data plot of stdev hi/lo chunk ps voltage v chunk # %%%%%%%
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
% %%%%%%%%%% ramp data plot of mean down chunk leakage current v chunk # %%%%%%%
%     figure4= figure('Units','normalized')
%     plot([ 1:1:num_down_chunks(i) ], lcm1_avg_ramp_down_avg_chunk(i,1:num_down_chunks(i)) - lcm1_avg_offset(i),...
%     'o','Color', 'blue','MarkerSize', 8, 'LineWidth', 2.0);
%     %errorbar([1:1:num_up_chunk_rows(i), lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,'o');
% %    axis(plot_bounds)
%     pbaspect([1.33 1 1])
%     ax = gca; % current axes
%     ax.FontSize = 32;
%     ax.TickDir = 'out'; % make ticks point out
%     title('lcm1 LO ramp means','FontSize',40)
%     xlabel('chunk #','FontSize',32)
%     ylabel('leakage current (pA)','FontSize',32)
%     annotation(figure4,'textbox',...
%     outside_plot,'String',{['avg LO (pA):'],...
%     [sprintf('%.1f',lcm1_avg_ramp_down_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_down_stdev(i))]},...
%     'FontSize',32,'BackgroundColor',[1 1 1]);
% 
% %%%%%%%%%% ramp data plot of mean up chunk leakage current v chunk # %%%%%%%
%     figure5= figure('Units','normalized')
%     plot([1:1:num_up_chunks(i)], lcm1_avg_ramp_up_avg_chunk(i,1:num_up_chunks(i)) - lcm1_avg_offset(i),...
%     'x','Color', 'red','MarkerSize', 8, 'LineWidth', 2.0);
%     %errorbar([[1:1:num_up_chunk_rows(i)], lcm1_avg_ramp_up_avg_chunk, lcm1_avg_ramp_up_stdev_chunk,'o');
% %    axis(plot_bounds)
%     pbaspect([1.33 1 1])
%     ax = gca; % current axes
%     ax.FontSize = 32;
%     ax.TickDir = 'out'; % make ticks point out
%     title('lcm1 HI ramp means','FontSize',40)
%     xlabel('chunk #','FontSize',32)
%     ylabel('leakage current (pA)','FontSize',32)
%     annotation(figure5,'textbox',...
%     outside_plot,'String',{['avg HI (pA):'],...
%     [sprintf('%.1f',lcm1_avg_ramp_up_avg(i) - lcm1_avg_offset(i)) ' \pm ' sprintf('%.1f',lcm1_avg_ramp_up_stdev(i))]},...
%     'FontSize',32,'BackgroundColor',[1 1 1]);
% 
% %%%%%%%%%% ramp data plot of stdev hi/lo chunk leakage current v chunk # %%%%%%%
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
% 
% % %%%%%%% test chunk counting algorithm psvoltage v time %%%%%%%%%%%%%%%%%%%
% %    subplot(num_files,1,num_files + 1 - i)
%    plot(time(i,start_offset(i):end_offset(i)), vmon_avg(i,start_offset(i):end_offset(i)) - vmon_avg_offset(i),...
%    'o','Color', cmap(1+i,:),'MarkerSize', 6, 'LineWidth', 1.0); hold on;
% %    plot(time_ramp_up(i,1:num_ramp_up_points), vmon_avg_ramp_up(i,1:num_ramp_up_points) - vmon_avg_offset(i),...
% %    '^','Color', 'blue','MarkerSize', 7, 'LineWidth', 2.0); hold on;
%    for j = 1:num_up_chunks(i)
%        chunk_begin = up_chunk_array(i,2,j) +1;
%        chunk_end = up_chunk_array(i,2,j+1);
%        plot(time_ramp_up(i,chunk_begin:chunk_end), vmon_avg_ramp_up(i,chunk_begin:chunk_end) - vmon_avg_offset(i),...
%        'x','Color', 'green','MarkerSize', 6, 'LineWidth', 2.0);
%    end
%    for j = 1:num_down_chunks(i)
%        chunk_begin = down_chunk_array(i,2,j) +1;
%        chunk_end = down_chunk_array(i,2,j+1);
%        plot(time_ramp_down(i,chunk_begin:chunk_end), vmon_avg_ramp_down(i,chunk_begin:chunk_end) - vmon_avg_offset(i),...
%       'x','Color', 'blue','MarkerSize', 6, 'LineWidth', 2.0);
%    end
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
% plot(fit_line_x(1,:),...
%      fit_line_y(1,:),...
%     'k-', 'LineWidth', 1.0); hold on; %nA
% errorbar(current_source_avg(1,:),...         
%      lcm1_avg(1,:)*1e-3,lcm1_avg_stdev(1,:)*1e-3,...
%     'o','Color', 'blue','MarkerSize', 10, 'LineWidth', 1.); hold on; %nA
%     %axis(plot_bounds)
%  pbaspect([1.33 1 1])
%  ax = gca; % current axes
%  ax.TickDir = 'out'; % make ticks point out
 
  %%%%%%%% plot leakage sensitivity fit line residual %%%%%%%%%%%%%%%%%%%%%
%convert to pA
figure1 = figure('Units','normalized')
plot(fit_line_x(1,:)*1e3,...
     [0 0]*1e3,...
    'k--', 'LineWidth', 2.0); hold on; %nA
errorbar(current_source_avg(1,:)*1e3,...         
     leakage_residual(1,1:numpoints(1),1)*1e3,lcm1_avg_stdev(1,:),...
    'o','Color', 'blue','MarkerSize', 12, 'LineWidth', 1.5); hold on; %nA
%     %axis(plot_bounds)
 pbaspect([1.33 1 1])
 ax = gca; % current axes
 ax.TickDir = 'out'; % make ticks point out

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

title(title_string,'FontSize',40)
xlabel(x_label,'FontSize',32)
ylabel(y_label,'FontSize',32)
l = legend('show'); 
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
annotation(figure1,'textbox',...
   outside_plot,'String',str,'FontSize',32,'BackgroundColor',[1 1 1]);