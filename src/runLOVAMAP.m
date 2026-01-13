clearvars;

% FILE PARAMETERS
filename          = 'particle_assembly.txt';                % input file name
file_path         = './particle_domains/';                  % input file path (include final forward-slash)
excel_path        = './outputs/';                           % file path for Excel output (include final forward-slash)

% INPUT PARAMETERS
%%%%%%%%%%%%%%%%% Only applicable to spherical data (.txt,.csv,.dat file inputs) %%%%%%%%%%%%%%%%%
voxel_size        = 2;                                      % desired mesh size in micrometers   
voxel_range       = [1e7, 1e8];                             % desired resolution range
%%%%%%%%%%%%%%%%%
combine_edge_subs = true;                                   % set to false if you do not want to merge any surface pores

% OUTPUT PARAMETERS
generate_raw_data = true;                                   % set to true to output Excel data
interior_only     = true;                                   % set to false to output both interior and surface pore data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECOMMEND: DO NOT MODIFY BELOW THIS LINE

% Additional input parameters
crop_percent      = 1;          % percentage of the domain to analyze
dip_percent       = 0.8;
hall_cutoff       = 6;          % radius, in um
shell_thickness   = 4;          % in um
num_2D_slices     = 30;

% Creating file path
file_path = [file_path, filename];

% Ensure required folders are on path
srcDir  = fileparts(mfilename('fullpath'));
baseDir = fileparts(srcDir);

% Define required folders
requiredPaths = {
    fullfile(baseDir, 'src'), ...
    fullfile(baseDir, 'external'), ...
    fullfile(baseDir, 'mex-lovamap')
};

% Add to path if not already included
for i = 1:numel(requiredPaths)
    p = requiredPaths{i};
    if exist(p, 'dir') && ~contains(path, p)
        addpath(genpath(p));
        fprintf('✔ Added to path: %s\n', p);
    end
end

% Label output files with date stamp
dateStamp = datestr(now, 'yymmdd-HHMM');

% Run LOVAMAP 
[data, time_log] = LOVAMAP(file_path, voxel_size, voxel_range, crop_percent, dip_percent, ...
    hall_cutoff, shell_thickness, num_2D_slices, combine_edge_subs);

% Output data to Excel file
if generate_raw_data
    % accommodate old naming system
    if strcmp(filename(1:8), 'beadInfo') || strcmp(filename(1:13), 'labeledDomain')
        excel_filename = replace(filename, {'beadInfo_', 'labeledDomain_', '.dat', '.txt', '.json'}, ...
                                           {'stats_', 'stats_', '.xlsx', '.xlsx', '.xlsx'});
    else
        excel_filename = replace(filename, {'.dat', '.txt', '.json'}, ...
                                           {'.xlsx', '.xlsx', '.xlsx'});
        excel_filename = horzcat('stats_', excel_filename);
    end
    write2excel(data, excel_path, excel_filename, filename, interior_only, 0);
end

[tt, lab] = secondsTime(time_log(end).Time);
fprintf('%30s %.5f %s\n', 'Total Time:', tt, lab);
