  function c = horzcat(varargin)
%|function c = horzcat(varargin)
%| called for c = [a, b, ...] where any of them is a fatrix

ob = cat(2, varargin{:});
% ob = block_fatrix(varargin, 'type', 'row'); % todo
