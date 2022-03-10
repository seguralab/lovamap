% Display time and write to file

function total_time = writeTime(toc_time, total_time, file, description)
    total_time = total_time + toc_time;
    [tt, lab] = secondsTime(total_time);
    fprintf('%30s %.5f %s %20s %.5f %s\n', description, toc_time, 'sec', '', tt, lab);
    fprintf(file, '%30s %.5f %s %20s %.5f %s\n', description, toc_time, 'sec', '', tt, lab);
end