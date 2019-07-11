clearvars -except masterlist binsize
close all

RaEDM_leakage_analysis('2018-10-24-164446-hv-1.txt',1); % runs
 
close all;

RaEDM_leakage_analysis('2018-10-25-153857-hv-1.txt',1); % runs but
% optimization quits. Update: increasing sample_rate time fixed this.

% RaEDM_leakage_analysis('2018-10-26-114809-hv-1.txt'); % doesn't run. 
% Note: voltage only goes up to 1 kV. Only ~50 minutes of data, not sure we
% need to look at this one.

% RaEDM_leakage_analysis('2018-10-26-121136-hv-1.txt'); % doesn't run.
% Note: only 30 minutes of data, probably don't need to look at this one.

close all;

RaEDM_leakage_analysis('2018-10-26-142009-hv-1.txt',1); % runs

close all;

RaEDM_leakage_analysis('2018-10-26-162012-hv-1.txt',1); % runs

close all;

RaEDM_leakage_analysis('2018-10-31-162157-hv-1.txt',1); % doesn't run. 
% Update: runs after increasing sample_rate time.

close all;

RaEDM_leakage_analysis('2018-11-01-174934-hv-1.txt',1); % doesn't run. 
% Update: runs after increasing sample_rate time 0.04 --> 0.06. Updated 
% chunk detecting algorithm to be independent of sampling time, so there
% shouldn't be any more patchwork solutions for this issue.

close all;

RaEDM_leakage_analysis('2018-11-02-140707-hv-1.txt',1); % runs

close all;

RaEDM_leakage_analysis('2018-11-02-162509-hv-1.txt',1); % doesn't run.
% Update: runs after updating chunk detecting algorithm.

close all;

RaEDM_leakage_analysis('2018-11-07-161056-hv-1.txt',1); % doesn't run.
% Update: runs after changing chunk detecting algorithm.

close all;

RaEDM_leakage_analysis('2018-11-08-162446-hv-1.txt',1); % runs

% RaEDM_leakage_analysis('2018-11-08-174949-hv-1.txt'); % doesn't run. note:
% this is less than 30 minutes of live data, probably not worth looking at.

% RaEDM_leakage_analysis('2018-11-09-124150-hv-1.txt'); % doesn't run.
% note: network storage space saturated and we lost about 4 hours of data,
% this file's no good.

close all;

RaEDM_leakage_analysis('2018-11-14-163806-hv-1.txt',1); % runs