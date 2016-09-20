% Produces a blank figure with everything turned off
%
% Optional argument: vector specifying axis limits (will be passed directly
% to axis() ).

function hf = blankFigure(varargin)

hf = figure; hold on; 
set(gca,'visible', 'off');
set(hf, 'color', [1 1 1]);
if ~isempty(varargin)
  axis(varargin{1});
end
axis square;
