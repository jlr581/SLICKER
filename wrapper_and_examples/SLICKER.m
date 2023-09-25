function [time, central_tendency, uncertainty, spread] = SLICKER(target,rtime,proxy1,varargin)

% SLICKER - runs Fortran SLICKER reconstruction code "reconstruct"
%
%   Brief description : This is a Matlab wrapper that generates input files to run the Fortran
%                       SLICKER reconstruction code described in Roberts et al. (2023).
%
%   User inputs:
%       target     2d array; 1st column is time, second column is the corresponding data value
%       rtime      1d vector of list of times that you want the reconstruction for
%       proxy1     2d array; 1st column is time, 2nd column is proxy value
%       varargin   up to an additional 9 input proxy time series (proxy2, ..., proxy9) and a list of name-value pairs (listed below) to change defaults
%
%   Input arrays may be different lengths. Data for each MUST be in strictly increasing time order.
%
%   Name-value pairs and default values
%       num_ensemble    Number of ensemble members (default 4096)
%       slick_width     Slick width parameters (default [0.4 1.6])
%       time_limit      Time limit (seconds) for solution of individual 
%                        ensemble member (default 60 seconds)
%       tol             Solution tolerance (default 1e-6)
%       nonlinear       Attempt a non-linear solution (default "n")
%       frac_subset     Fraction of ensemble with best stationarity used in
%                        reported solution (default 0.5)
%       proxy_inversion Test for significantly stronger correlation if 
%                        individual proxies are inverted (default "n" for 
%                        linear solution and "y" for non-linear)
%
%
%   Usage:
%       [time, central_tendency, uncertainty, spread] = SLICKER(target,rtime,proxy1,varargin)
%
%
%   Examples:
%       [time, central_tendency, uncertainty, spread] = SLICKER(target,rtime,proxy1,'num_ensemble',512);
%       [time, central_tendency, uncertainty, spread] = SLICKER(target,rtime,proxy1,proxy2,proxy3,'time_limit',20,'tol',1e-7);
%
%
%   Citations:
%   If you use this Matlab wrapper, please cite
%   Roberts et al. (submitted), 2023.
%
%   Author Info
%   This function was written by Felicity S. McCormack of Monash University,
%   and Jason Roberts of the Australian Antarctic Division
%   Version 2.0: May 2023
%

% Some checks
p=inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validMatrix = @(x) isnumeric(x) && ismatrix(x);
validVector = @(x) isnumeric(x) && isvector(x);
validVectorPosNum = @(x) isnumeric(x) && isvector(x) && (x > 0);
validyn = @(x) ischar(x) && (x=="y" || x=="Y" || x=="n" || x=="N");
addRequired(p,'target',validMatrix);
addRequired(p,'rtime',validVector);
addRequired(p,'proxy1',validMatrix);
addOptional(p,'proxy2',[],validMatrix);
addOptional(p,'proxy3',[],validMatrix);
addOptional(p,'proxy4',[],validMatrix);
addOptional(p,'proxy5',[],validMatrix);
addOptional(p,'proxy6',[],validMatrix);
addOptional(p,'proxy7',[],validMatrix);
addOptional(p,'proxy8',[],validMatrix);
addOptional(p,'proxy9',[],validMatrix);
addOptional(p,'proxy10',[],validMatrix);
addParameter(p,'num_ensemble',4096,validScalarPosNum);
addParameter(p,'slick_width', [0.4 1.6] ,validVectorPosNum);
addParameter(p,'time_limit', 60 , validScalarPosNum);
addParameter(p,'tol', 1e-6 , validScalarPosNum);
addParameter(p,'nonlinear', "n", validyn);
addParameter(p,'frac_subset', 0.5 , validScalarPosNum);
addParameter(p,'proxy_inversion', "e" , validyn);
parse(p,target,rtime,proxy1,varargin{:});

num_ensemble=p.Results.num_ensemble;
slick_width=p.Results.slick_width;
time_limit=p.Results.time_limit;
tol=p.Results.tol;
isnonlinear=lower(p.Results.nonlinear);
frac_subset=p.Results.frac_subset;
if (p.Results.proxy_inversion=="e")
  if (isnonlinear=="y")
    ispinversion="y";
  else
    ispinversion="n";
  end
else
  ispinversion=lower(p.Results.proxy_inversion);
end
proxy{1}=proxy1;
num_proxy_files=1;

% Get proxies 
fields = fieldnames(p.Results);
for i=4:11
	fieldName = fields{i};
	fieldValue = p.Results.(fieldName);

	if ~isempty(fieldValue)
		num_proxy_files = num_proxy_files + 1;
		proxy{num_proxy_files} = p.Results.(fieldName);
	end
end

if (nargin<3)
   help SLICKER
   error('SLICKER error message: bad usage')
end
if size(target, 2) ~= 2
    error('SLICKER error message: reconstruction target data should be Mx2.');
end

for i=1:num_proxy_files
  if size(proxy{i},2) ~= 2
    error('SLICKER error message: reconstruction proxy ',string(i),' should be Mx2.');
  end
end


% Generate user input data

% Filename of reconstruction target data
basename=strcat("SLICKER_",date,"_",datestr(now,'HH_MM_SS_FFF'),"_");
rfilename = strcat(basename,"target.txt");

% Filenames of proxy files
for i = 1:num_proxy_files
    pfilenames{i} = strcat(basename,"proxy",int2str(i),".txt");
end

% Filename of times to be used for reconstruction
tfilename = strcat(basename,"timebase.txt");


% Filename of output file
ofilename = strcat(basename,"output.txt");

% Filename of master input file
mfilename = strcat(basename,"input_keyboard_file.txt");


% Generate user input file
fid = fopen(mfilename, 'w');
fprintf(fid, '%s\n', rfilename);
fprintf(fid, '%d\n', num_proxy_files);
for i = 1:num_proxy_files
   fprintf(fid, '%s\n', pfilenames{i});
end
fprintf(fid, '%s\n', tfilename);
fprintf(fid, '%d\n', num_ensemble);
fprintf(fid, '%f %f\n', slick_width);
fprintf(fid, '%d\n', time_limit);
fprintf(fid, '%e\n', tol);
fprintf(fid, '%s\n', ofilename);
fprintf(fid, '%s\n', isnonlinear);
fprintf(fid, '%f\n', frac_subset);
fprintf(fid, '%s\n', ispinversion);
fclose(fid);


% Write target file
dlmwrite(rfilename, target, 'delimiter', '\t');

% Write reconstruction time file
dlmwrite(tfilename, rtime, 'delimiter', '\t');

% Write proxy files
for i = 1:num_proxy_files
    dlmwrite(pfilenames{i}, proxy{i}, 'delimiter', '\t');
end


% Run Fortran code
SLICK_path=which("SLICKER");
FORTRAN_path=strrep(SLICK_path,"SLICKER.m","reconstruct");
systemcommand = [FORTRAN_path,' < ', mfilename ,' > /dev/null'];
system(systemcommand);

% Read the output file 
fid = fopen(ofilename,'r');
tline=fgetl(fid);
i=1; j=1;
while ischar(tline)
   tline=fgetl(fid);
   if ischar(tline)
      index = cellfun(@(x) ~startsWith(x, '!'), {tline});
      if index
         values=strsplit(tline);
         data(j,:)=str2double({values{find(~cellfun('isempty', values))}});
         j=j+1;
      end
      i=i+1;
   end
end
fclose(fid);

time = data(:,1);
central_tendency = data(:,2);
uncertainty = data(:,3);
spread = data(:,4);

