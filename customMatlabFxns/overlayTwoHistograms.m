function overlayTwoHistograms(histbins,y1,y2,color1,color2)
%% overlayTwoHistograms
% Create an overlaid histogram from two frequencies distributed over
% histbins.  Where the histograms overlap, the color will be an average
% of the colors supplied for the two histograms
% Inputs:
%   histbins -- the domain of the histograms, column vector
%   y1 -- the frequencies for histogram 1, column vector
%   y2 -- the frequencies for histogram 2, column vector
%   color1 -- the color for the histogram bars for histogram 1, 3D row
%       vector
%   color2 -- the color for the histogram bars for histogram 2, 3D row
%       vector
% Outputs:
%   plots the overlaid histograms

%% debugging inputs:


%% end 

both = min([y1,y2],[],2);
hold on
bar(histbins,y1,'FaceColor',color1)
bar(histbins,y2,'FaceColor',color2)
bar(histbins,both,'FaceColor', mean([color1;color2]))