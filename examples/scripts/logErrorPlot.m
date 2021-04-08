function [hP,hE] = logErrorPlot(x_data,y_data,e_data,bottom,Pstyle,Estyle)
%
% function h = logErrorPlot(x_data,y_data,e_data,bottom,Pstyle,Estyle)
%
% Plots data and vertical errorbars without endcaps in y-log scale.
%
% inputs:
% x_data: the x data
% y_data: the y data
% e_data: the symmetric vertical errors
% bottom: the minimum y value of plotted error bars (this is needed when
% error is greater than data in log-scale plot)
% Pstyle: a cell array of variable length containing stylings for
% the plot itself (optional)
% Estyle: a cell array of variable length containing stylings for
% the errorbars (optional)
%
% outputs:
% hP: the handle for the data
% hE: the handle for the errorbars
%
% Example:
% >>myErrorbar(1:10,1:10,1:-.1:.1,{'rd','markersize',3})
  
  holdisOn = ishold;
  if nargin<5
    Pstyle = {};
  end;
  if isempty(Pstyle)
    Pstyle = {'b.'};
  end;
  
  hP = semilogy(x_data,y_data,Pstyle{:});
  
  hold on;
  if nargin<6
    Estyle = {};
  end;
  if isempty(Estyle)
    plotColor = get(hP,'color');
    Estyle = {'color',plotColor,'linewidth',2};
  end;
  
  hE = line([x_data;x_data],[y_data+e_data;max(y_data-e_data,bottom)],Estyle{:});
  
  if ~holdisOn
    hold off;
  end; 
