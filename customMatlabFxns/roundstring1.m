function numstring = roundstring2(num)
% converts a decimal number to a string that includes <= 1 decimals after
% zero

num = 10*num;
num = round(num);
num = num/10;
numstring = num2str(num);
