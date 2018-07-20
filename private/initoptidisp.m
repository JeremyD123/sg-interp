function [isdispiter, iterstr] = initoptidisp(dispopt)
% INITOPTIDISP   Initializes iteration information display
%   (Internal function)

if strcmpi(dispopt, 'iter') 
	disp(' Iteration   Func-count Grad-count     f(x)            Procedure');
	iterstr =' %5.0f        %5.0f     %5.0f     %12.6g         %s';
	isdispiter = true;
else
	iterstr = [];
	isdispiter = false;
end
