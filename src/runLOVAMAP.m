clearvars;

% Excel raw data
generate_raw_data = 0;
interior_only     = 1;

% File
% data_filename = 'beadInfo_square100.dat'; % for figs, use 130_5.dat
% data_filename_full = ['../lindsay_data/domain_spheres/square_packing/', data_filename];

input_file = '../lindsay_data/domain_spheres/square_packing/beadInfo_square160.dat';

% Parameters
voxel_size      = 2;
voxel_range     = [1e7, 1e8]; % desired resolution
crop_percent    = 1; % 0.5 100_0.dat for wings gif
hall_cutoff     = 6; % radius, in um
shell_thickness = 4; % in um
num_2D_slices   = 30;

% label output files with date stamp
dateStamp = datestr(now, 'yymmdd-HHMM');

% runtimes file
% runtimes_file = fopen(['../lindsay_data/outputs/runtimes_', dateStamp, '.txt'], 'wt+');
% output the date and time for the run
% printdatetime(runtimes_file);

% output parameters used in the run
fprintf('%30s %s\n\n %29s %.2f\n %29s %.2f\n %29s %.1f\n %29s %i\n %29s %i\n\n', ...
        'Input file:', input_file, ...
        'crop_percent:', crop_percent, ...
        'hall_cutoff:', hall_cutoff);

% Analyze void space
[data, time_log] = LOVAMAP(input_file, voxel_size, voxel_range, crop_percent, hall_cutoff, ...
    shell_thickness, num_2D_slices);

% Output data to Excel file
% if strcmp(generate_raw_data, 'on') || strcmp(generate_raw_data, 'yes') || sum(generate_raw_data == 1) == 1
%     % accommodate old naming system
%     if strcmp(data_filename(1:8), 'beadInfo') || strcmp(data_filename(1:13), 'labeledDomain')
%         excel_filename = replace(data_filename, {'beadInfo_', 'labeledDomain_', '.dat', '.txt', '.json'}, ...
%                                                 {'stats_', 'stats_', '.xlsx', '.xlsx', '.xlsx'});
%     else
%         excel_filename = replace(data_filename, {'.dat', '.txt', '.json'}, ...
%                                                 {'.xlsx', '.xlsx', '.xlsx'});
%         excel_filename = horzcat('stats_', excel_filename);
%     end
%     excel_path = '../lindsay_data/outputs/';

%     write2excel(data, excel_path, excel_filename, data_filename, interior_only, 0);
% end

[tt, lab] = secondsTime(time_log(end).Time);
fprintf('%30s %.5f %s\n', 'Total Time:', tt, lab);
