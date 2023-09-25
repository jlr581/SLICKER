addpath("..")
% read in data from files
fileID_target = fopen('test5_input_target.txt','r');
target = fscanf(fileID_target,'%d %f',[2 100])';
fclose(fileID_target);
fileID_timebase = fopen('test5_input_timebase.txt','r');
timebase = fscanf(fileID_timebase,'%d',[1 100])';
fclose(fileID_timebase);
fileID_p1 = fopen('test5_input_proxy1.txt','r');
p1 = fscanf(fileID_p1,'%d %f',[2 80])';
fclose(fileID_p1);
fileID_p2 = fopen('test5_input_proxy2.txt','r');
p2 = fscanf(fileID_p2,'%d %f',[2 80])';
fclose(fileID_p2);

% call SLICKER with default settings except for number of ensemble members, time limit for each ensemble member and the solution tolerance
[time,central_tendency,uncertainty,spread]=SLICKER(target,timebase,p1,p2,'num_ensemble',512,'time_limit',20,'tol',1e-7);

% write the output to file
t=[time,central_tendency,uncertainty,spread]';
fileID = fopen('test5_MATLAB_output.txt','w');
fprintf(fileID,'%12s %12s %12s %8s\n','time base','M-estimator','+/- for_95%_CI',' Qn');
fprintf(fileID,'%12.4e %12.4e %12.4e %12.4e\n',t);
fclose(fileID);


