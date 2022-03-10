% Check whether ridges1D exists between 2 peaks. If not, combine peaks into
% one

% e.g.:
% i          - current ridge1D
% peaks      - peaks.L6
% peaks_key  - peaks.pr1D_key
% peaks_cntr - peaks.pr1D_cntr

function [orig_pks, pot_pks, peaks, peaks_key, peaks_cntr, ...
          ridges1D_ind, ridges1D_key, ridges1D_split] = ...
          combinePks(i, pot_pks, peaks, ridges1D, ...
                     peaks_key, peaks_cntr, voxels)
    
    check_pks = peaks.indices(pot_pks);
    orig_pks  = 1 : peaks.num;
    pkperms   = nchoosek(check_pks, 2);
    no_match  = false(size(pkperms, 1), 1);
    for j = 1 : size(pkperms, 1)
        matches_tot = arrayfun(@(x) fastIntersect(ridges1D.rp_key(x, 1:2), ...
                                                  sort(pkperms(j, :)), 'all'), ...
                               ridges1D.pks2);
        if sum(matches_tot) == 0
            % no ridge exists
            % confirm peaks are physically close, then store to combine
            phys_dist = vecnorm(voxels(pkperms(j, 1), :) - voxels(pkperms(j, 2), :), 2, 2);
            pk_ind1   = binarySearch(peaks.indices, pkperms(j, 1));
            pk_ind2   = binarySearch(peaks.indices, pkperms(j, 2));
            EDT_radii = peaks.EDT(pk_ind1) + peaks.EDT(pk_ind2);
            if phys_dist < EDT_radii
                no_match(j) = true;
            end
        end
    end
    % combine all peaks that do not share ridge1D
    combine_pks = unique(pkperms(no_match, :));
    if ~isempty(combine_pks)
        % keep the peak with the higher EDT and remove the rest
        combine_pks_ind = arrayfun(@(x) binarySearch(peaks.indices, combine_pks(x)), ...
                           1 : length(combine_pks));
        [~, max_ind] = max(peaks.EDT(combine_pks_ind));
        remain_ind   = [1 : max_ind - 1, max_ind + 1 : length(combine_pks)];
        % store an associated ridges1D to later move voxel to ridge
        push2r1D = zeros(length(combine_pks), 1);
        for k = 1 : ridges1D.num % including ridge i
            column_ind = fastIntersect(ridges1D.rp_key(k, 1 : ridges1D.rp_counter(k)), ...
                                       combine_pks(remain_ind), 'indices');
            if ~isempty(column_ind)
                % replace
                if length(column_ind) > 1
                    error('A single ridge1D is flanked by 2 peaks that will be removed.')
                else
                    % store ridge number of voxels to be removed
                    % (exactly which ridge among those that match shouldn't matter)
                    if k ~= i
                        % skip ridge1D i so don't have to worry about
                        % shifting ridge voxels when leaving fnc
                        push2r1D(combine_pks == ridges1D.rp_key(k, column_ind)) = ...
                                 k;
                    end
                    % key
                    ridges1D.rp_key(k, column_ind) = combine_pks(max_ind);
                    % split
                    if fastIntersect(ridges1D.split{1, k}, combine_pks(remain_ind), 'bool')
                        ridges1D.split{1, k}(fastIntersect(ridges1D.split{1, k}, ...
                                             combine_pks(2 : end), 'indices')) = combine_pks(max_ind);
                    elseif fastIntersect(ridges1D.split{4, k}, combine_pks(remain_ind), 'bool')
                        ridges1D.split{4, k}(fastIntersect(ridges1D.split{4, k}, ...
                                             combine_pks(2 : end), 'indices')) = combine_pks(max_ind);
                    else
                        error('Ridge1D was associated with a peak in its key that does not exist in split.')
                    end
                end
            end
        end
        % remove peaks: combine_pks(remain_ind)
        % track changes
        if length(combine_pks) > 1
            for k = remain_ind(:)'
                ind = find(peaks.indices == combine_pks(k));
                peaks.indices(ind)  = [];
                peaks.beads(ind)    = [];
                peaks.numBeads(ind) = [];
                peaks.num           = peaks.num - 1;
                peaks.EDT(ind)      = [];
                % peaks key
                peaks_key(ind, :)   = [];
                peaks_cntr(ind)     = [];
                % tracking changes
                orig_pks(ind)       = [];
                % move old peak indices onto appropriate ridge
                ridges1D.indices{push2r1D(k)} = sort([ridges1D.indices{push2r1D(k)}; ...
                                                      combine_pks(k)]);
            end
        end
        % update check_pks
        rmv_pks = fastIntersect(check_pks, combine_pks(remain_ind), ...
                                'indices');
        pot_pks(rmv_pks) = [];
    end
    ridges1D_ind   = ridges1D.indices;
    ridges1D_key   = ridges1D.rp_key;
    ridges1D_split = ridges1D.split;
end