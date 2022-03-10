% Print out date and time in line with runtime outputs and
% Write to file

function printdatetime(runtimes_file)
    c = clock;
    
    % day
    if c(3) < 10
        % hour
        if c(4) == 0
            fprintf('\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5))
            fprintf(runtimes_file, '\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5));
        elseif c(4) < 10
            fprintf('\n%22d-%1d-%4i  %1i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5))
            fprintf(runtimes_file, '\n%22d-%1d-%4i  %1i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5));
        elseif c(4) < 12
            fprintf('\n%22d-%1d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5))
            fprintf(runtimes_file, '\n%22d-%1d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5));
        elseif c(4) == 12
            fprintf('\n%22d-%1d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),c(4),c(5))
            fprintf(runtimes_file, '\n%22d-%1d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),c(4),c(5));
        elseif c(4) < 22
            fprintf('\n%22d-%1d-%4i  %1i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5))
            fprintf(runtimes_file, '\n%22d-%1d-%4i  %1i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5));
        elseif c(4) < 24
            fprintf('\n%22d-%1d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5))
            fprintf(runtimes_file, '\n%22d-%1d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5));
        else
            fprintf('\n%22d-%1d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5))
            fprintf(runtimes_file, '\n%22d-%1d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5));
        end
    else
        % hour
        if c(4) == 0
            fprintf('\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5))
            fprintf(runtimes_file, '\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5));
        elseif c(4) < 10
            fprintf('\n%21d-%2d-%4i  %1i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5))
            fprintf(runtimes_file, '\n%21d-%2d-%4i  %1i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5));
        elseif c(4) < 12
            fprintf('\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5))
            fprintf(runtimes_file, '\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),c(4),c(5));
        elseif c(4) == 12
            fprintf('\n%21d-%2d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),c(4),c(5))
            fprintf(runtimes_file, '\n%21d-%2d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),c(4),c(5));
        elseif c(4) < 22
            fprintf('\n%21d-%2d-%4i  %1i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5))
            fprintf(runtimes_file, '\n%21d-%2d-%4i  %1i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5));
        elseif c(4) < 24
            fprintf('\n%21d-%2d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5))
            fprintf(runtimes_file, '\n%21d-%2d-%4i  %2i:%02ipm\n\n', c(2),c(3),c(1),mod(c(4), 12),c(5));
%         else
%             fprintf('\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5))
%             fprintf(runtimes_file, '\n%21d-%2d-%4i  %2i:%02iam\n\n', c(2),c(3),c(1),12,c(5));
        end
    end
end