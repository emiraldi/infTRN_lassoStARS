function numstring = roundstring2(num)
% converts a decimal number to a string that includes <= 2 decimals after
% zero

num = 100*num;
num = round(num);
num = num/100;
numstring = num2str(num);
