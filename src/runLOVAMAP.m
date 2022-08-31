clearvars;

% Excel raw data
generate_raw_data = true;
interior_only     = true;

% File
file_name  = 'beadInfo_{100,100}_100_{0,100}_0.dat';
file_path  = ['./particle_domains/', file_name];
excel_path = './lovamap_outputs/';

% Parameters
voxel_size        = 2;
voxel_range       = [1e7, 1e8]; % desired resolution
crop_percent      = 1;          % percentage of the domain to analyze
dip_percent       = 0.8;
hall_cutoff       = 6;          % radius, in um
shell_thickness   = 4;          % in um
num_2D_slices     = 30;
combine_edge_subs = true;

% label output files with date stamp
dateStamp = datestr(now, 'yymmdd-HHMM');

% Analyze void space
[data, time_log] = LOVAMAP(file_path, voxel_size, voxel_range, crop_percent, dip_percent, ...
    hall_cutoff, shell_thickness, num_2D_slices, combine_edge_subs);

% Output data to Excel file
if generate_raw_data
    % accommodate old naming system
    if strcmp(file_name(1:8), 'beadInfo') || strcmp(file_name(1:13), 'labeledDomain')
        excel_filename = replace(file_name, {'beadInfo_', 'labeledDomain_', '.dat', '.txt', '.json'}, ...
                                           {'stats_', 'stats_', '.xlsx', '.xlsx', '.xlsx'});
    else
        excel_filename = replace(file_name, {'.dat', '.txt', '.json'}, ...
                                           {'.xlsx', '.xlsx', '.xlsx'});
        excel_filename = horzcat('stats_', excel_filename);
    end
    write2excel(data, excel_path, excel_filename, file_name, interior_only, 0);
end

[tt, lab] = secondsTime(time_log(end).Time);
fprintf('%30s %.5f %s\n', 'Total Time:', tt, lab);
