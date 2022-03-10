% Peter function

function [idx] = binarySearch(sorted_array, num)
    left = 1;
    right = length(sorted_array);
    flag = 0;
    
    while (left <= right)
        mid = ceil((left + right) / 2);
        
        if (sorted_array(mid) == num)
            idx = mid;
            flag = 1;
            break;
        elseif (sorted_array(mid) > num)
            right = mid - 1;
        else
            left = mid + 1;
        end
    end
    
    if (flag == 0)
        idx = -1;
    end
end