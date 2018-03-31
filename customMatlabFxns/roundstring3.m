function numstring = roundstring3(num)
%% function numstring = roundstring3(num)
% converts a decimal number to a string that includes <= 3 decimals after
% zero

if num >= .001
    num = 1000*num;
    num = round(num);
    num = num/1000;
    numstring = num2str(num);
else
    numstring = num2str(num,'%10.1E');
end