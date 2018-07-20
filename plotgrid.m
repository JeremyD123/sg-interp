function varargout = plotgrid(n,d,options)
% PLOTGRID   Plot the sparse grid
%    PLOTGRID(N,D) Plots the sparse grid points of dimension D = 2
%    or D = 3 up to level N.
%    
%    PLOTGRID(N,D,OPTIONS) Plots as above, but with default grid
%    type replaced by the grid type specified in OPTIONS, an
%    argument created with the SPSET function. See SPSET for
%    details.
%
%    H = PLOTGRID(...) Returns a vector of handles to the grid points
%    (useful for changing the look of the plotted grid).  
	
if nargin ~= 2 && nargin ~= 3
	error(['Wrong number of arguments. Please see ''help plotgrid''' ...
				 ' for usage information.']);
end

if nargin < 3, options = []; end

holdstate = ishold;
h = [];

switch d
 case 2
	for k = 0:n
		x = spgrid(k,2,options);
		th = plot(x(:,1),x(:,2),'k.');
		h = [h; th];
		axis equal;
		axis([0 1 0 1]);
		hold on;
	end
 case 3
	for k = 0:n
		x = spgrid(k,3,options);
		th = plot3(x(:,1),x(:,2),x(:,3),'k.');
		h = [h; th];
		axis equal;
		axis([0 1 0 1 0 1]);
		hold on;
	end
 otherwise
	error(['Invalid parameter d = ' num2str(d) '.']);
end

if ~holdstate
	hold off;
end
	
if nargout > 0
  varargout{1} = h;
end