clearvars -except masterlist binsize
close all

RaEDM_leakage_analysis('2018-12-13-144501-hv-1.txt',4); % runs, but does 
% not fully optimize some chunks. Note: This is set to 2 uA range. 
% stdev_binwidth --> 4 pA

RaEDM_leakage_analysis('2018-12-18-160059-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-05-165127-hv-1.txt',4); % had to change 
% binwidth settings avg_binwidth --> 50 pA stdev_binwidth --> 10 pA, 
% summed_stdev_binwidth --> 2.0 pA

RaEDM_leakage_analysis('2019-02-07-165243-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-12-120123-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-12-140215-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-12-165220-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-14-135705-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-14-155710-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-26-170742-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-28-165856-hv-1.txt',4);

RaEDM_leakage_analysis('2019-03-12-170048-hv-1.txt',4);

RaEDM_leakage_analysis('2019-03-14-165910-hv-1.txt',4);

RaEDM_leakage_analysis('2019-03-19-170320-hv-1.txt',4);

RaEDM_leakage_analysis('2019-03-21-165502-hv-1.txt',4);

RaEDM_leakage_analysis('2019-03-26-171237-hv-1.txt',4);

RaEDM_leakage_analysis('2019-03-28-171324-hv-1.txt',4);

RaEDM_leakage_analysis('2019-04-02-170346-hv-1.txt',4);

RaEDM_leakage_analysis('2019-04-04-165853-hv-1.txt',4);

RaEDM_leakage_analysis('2019-04-09-165437-hv-1.txt',4);

RaEDM_leakage_analysis('2019-04-11-171358-hv-1.txt',4); % doesn't run.
% Update: now runs after fixing sort_state/

RaEDM_leakage_analysis('2019-04-12-162334-hv-1.txt',4);

RaEDM_leakage_analysis('2019-04-16-164342-hv-1.txt',4); % runs, but doesn't
% finish chunk optimization

RaEDM_leakage_analysis('2019-04-18-132346-hv-1.txt',4); % runs, but doesn't
% finish chunk optimization

RaEDM_leakage_analysis('2019-04-26-151416-hv-1.txt',4); % runs, but doesn't
% finish chunk optimization

RaEDM_leakage_analysis('2019-05-01-171913-hv-1-trimmed.txt',4); % runs, doesn't 
% correctly identify chunks. Also had to remove the first 50 minutes of off
% data to get it to run without doing some extensive changing of code.
