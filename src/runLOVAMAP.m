clearvars;

% File parameters
filename          = 'YS_test_1_segment_300.json';
file_path         = ['./domain_labelled/real/', filename];
excel_path        = './outputs/';

% Input parameters (will get this data from json file if .json)
voxel_size        = 2;
voxel_range       = [1e7, 1e8]; % desired resolution
crop_percent      = 1;          % percentage of the domain to analyze
dip_percent       = 0.8;
hall_cutoff       = 6;          % radius, in um
shell_thickness   = 4;          % in um
num_2D_slices     = 30;

% Json required input
combine_edge_subs = false; % typically needed for real scaffold images

% Output parameters
generate_raw_data = true;       % export data to excel
interior_only     = false;       % output interior subunits only

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
