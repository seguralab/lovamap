% Write data outputs to Excel file
%
% mult_domains - 1 if running multiple domains stored in single excel file;
%                0 if not
 
 
function write2excel(data, excel_path, excel_filename, sheet_name, ...
                     interior_only, mult_domains)
 
    warning('off','MATLAB:xlswrite:AddSheet');
 
    num_subs = data.Descriptors.Global.numSubs;
    
    % for subunit data
    which_subslog = true(num_subs, 1);
    if interior_only || strcmp(interior_only, 'yes') || strcmp(interior_only, 'on') || (isnumeric(interior_only) && interior_only == 1)
        % select only interior data.Subunits
        for i = 1 : num_subs
            if data.Subunits{i}.edge
                which_subslog(i) = false;
            end
        end
        % add 'interior' to file name
        excel_filename = replace(excel_filename, {'.'}, {'_interior.'});
        excel_filename_full = [excel_path, excel_filename];
    else
        excel_filename_full = [excel_path, excel_filename];
    end
    if mult_domains ~= 1
        version = 1;
        while isfile(excel_filename_full)
            excel_filename_full = replace(excel_filename_full, ...
                                          {'.xlsx', ['(', num2str(version - 1), ').xlsx']}, ...
                                          ['(', num2str(version), ').xlsx']);
            version = version + 1;
        end
    end
        
    % use domain filename to name sheet
    sheet_name = replace(sheet_name, {'beadInfo_', 'labeledDomain_', '.dat', '.txt', '.json'}, ...
                                     {'', '', '', '', ''});
    sheet_name = sheet_name(1 : min(length(sheet_name), 31)); % Excel has character limit
 
    % FIRST: Populate with GLOBAL data
    global_fields     = fieldnames(data.Descriptors.Global);
    num_global        = length(data.Descriptors.Global.names);
    global_data       = cell(2, num_global);
    global_data(1, :) = data.Descriptors.Global.names;
    for i = 1 : length(global_data)
        % remember first fieldname contains the names of the descriptors
        global_data{2, i} = data.Descriptors.Global.(global_fields{i + 1});
    end
    writecell(global_data, excel_filename_full, 'Sheet', sheet_name);
    
    % SECOND: Populate with INTER-SUBUNIT descriptor data
    intersub_fields = fieldnames(data.Descriptors.NonSubs);
    num_intersubs   = length(data.Descriptors.NonSubs.names);
    
    writecell(data.Descriptors.NonSubs.names', excel_filename_full, 'Sheet', sheet_name, ...
              'Range', 'A3'); 
    for i = 1 : num_intersubs
        % remember first fieldname contains the names of the descriptors
        writematrix(data.Descriptors.NonSubs.(intersub_fields{i + 1})(:), excel_filename_full, 'Sheet', sheet_name, ...
                    'Range', [num2ExcelCol(i), '4']);
    end
    
    % THIRD: Populate with SUBUNIT descriptor data
    sub_fields  = fieldnames(data.Descriptors.Subs);
    num_subcats = length(data.Descriptors.Subs.names);
    
    writecell(data.Descriptors.Subs.names', excel_filename_full, 'Sheet', sheet_name, ...
              'Range', [num2ExcelCol(num_intersubs + 1), '3']); 
    for i = 1 : num_subcats
        % remember first fieldname contains the names of the descriptors
        writematrix(data.Descriptors.Subs.(sub_fields{i + 1})(which_subslog), excel_filename_full, 'Sheet', sheet_name, ...
                    'Range', [num2ExcelCol(num_intersubs + i), '4']);
    end
end
