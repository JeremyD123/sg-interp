function val = posvars(d, varpos, varargin)
% POSVARS   Shift interpolation parameters to right position
%    Shifts the interpolation parameters to the right position, and 
%    fills field up with any parameters in varargin. (Internal 
%    function).	

if ~isempty(varpos)
	val = cell(1,d+length(varargin));
	for l = 1:d
		val{varpos(l)} = 1;
	end
	m = 1;
	if length(varargin) > 0
		for l = 1:d+length(varargin)
			if isempty(val{l})
				val{l} = varargin{m};
				m = m + 1;
			end
		end
	end
else
	val = varargin;
end
