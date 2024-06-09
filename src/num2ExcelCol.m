% Convert numbers to capital letters that correspond to Excel columns

function l = num2ExcelCol(n)
    if n <= 26
        l = char(n + 64);
    else
        first_letter = char(floor((n - 1) / 26) + 64);
        second_letter = char(mod((n - 1), 26) + 65);
        l = [first_letter, second_letter];
    end
end
