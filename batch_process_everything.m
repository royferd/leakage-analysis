% NET

clearvars -except masterlist binsize

RaEDM_leakage_analysis('2019-05-23-161940-hv-1.txt',4);

RaEDM_leakage_analysis('2019-05-31-154720-hv-1.txt',4);

RaEDM_leakage_analysis('2019-06-06-153104-hv-1.txt',4);

RaEDM_leakage_analysis('2019-06-20-150540-hv-1.txt',4);

RaEDM_leakage_analysis('2019-07-18-140330-hv-1.txt',4);

RaEDM_leakage_analysis('2019-07-18-153953-hv-1.txt',4);

RaEDM_leakage_analysis('2019-07-25-144731-hv-1.txt',4);

RaEDM_leakage_analysis('2019-08-01-134407-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-08-01-151610-hv-1.txt',4);

% Nb56

clearvars -except masterlist binsize
close all

RaEDM_leakage_analysis('2017-10-16-180518-hv-1.txt',1); 
 
RaEDM_leakage_analysis('2017-10-17-105238-hv-1.txt',1); 

RaEDM_leakage_analysis('2017-10-17-165207-hv-1.txt',1); 

RaEDM_leakage_analysis('2017-10-18-114452-hv-1.txt',1); 

RaEDM_leakage_analysis('2017-10-18-171909-hv-1.txt',1); 

RaEDM_leakage_analysis('2017-10-19-104525-hv-1.txt',1); 

RaEDM_leakage_analysis('2017-10-19-171143-hv-1.txt',1); 

RaEDM_leakage_analysis('2017-10-20-114432-hv-1.txt',1); 

RaEDM_leakage_analysis('2017-10-23-175202-hv-1.txt',1); 

% Nb23

clearvars -except masterlist binsize
close all

RaEDM_leakage_analysis('2019-11-26-171435-hv-1.txt',4); 
 
RaEDM_leakage_analysis('2019-12-04-171229-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-12-13-154803-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-12-16-135104-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-12-18-161830-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-03-142124-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-06-154742-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-08-154619-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-14-142319-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-15-150026-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-20-154812-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-21-133242-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-23-145403-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-27-142302-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-28-124706-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-01-30-150547-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-03-144707-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-06-145329-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-10-150401-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-11-133322-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-13-144854-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-17-144755-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-18-131540-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-24-145116-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-02-25-131832-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-03-05-145833-hv-1.txt',4); 

RaEDM_leakage_analysis('2020-03-12-144554-hv-1.txt',4); 

% Nb78

clearvars -except masterlist binsize
close all

RaEDM_leakage_analysis('2018-10-24-164446-hv-1.txt',4); 
 
RaEDM_leakage_analysis('2018-10-25-153857-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-10-26-142009-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-10-26-162012-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-10-31-162157-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-11-01-174934-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-11-02-140707-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-11-02-162509-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-11-07-161056-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-11-08-162446-hv-1.txt',4); 

RaEDM_leakage_analysis('2018-11-14-163806-hv-1.txt',4); 

% Ti13. Make sure to do the other two files that require different bin
% width settings manually!

clearvars -except masterlist binsize
close all

% RaEDM_leakage_analysis('2018-12-13-144501-hv-1.txt',4); % runs, but does 
% not fully optimize some chunks. Note: This is set to 2 uA range. 
% stdev_binwidth --> 4 pA

RaEDM_leakage_analysis('2018-12-18-160059-hv-1.txt',4);

RaEDM_leakage_analysis('2019-02-05-165127-hv-1.txt',4); 

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

RaEDM_leakage_analysis('2019-04-11-171358-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-04-12-162334-hv-1.txt',4);

RaEDM_leakage_analysis('2019-04-16-164342-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-04-18-132346-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-04-26-151416-hv-1.txt',4); 

RaEDM_leakage_analysis('2019-05-01-171913-hv-1-trimmed.txt',4); % runs, doesn't 
% correctly identify chunks. Also had to remove the first 50 minutes of off
% data to get it to run without doing some extensive changing of code.