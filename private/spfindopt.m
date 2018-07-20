function [xopt, fval] = spfindopt(z, xbox, isminimize, ismaximize)
% SPFINTOPT   Compute the extrema over the sparse grid points for
% the box XBOX

if nargin < 3;
  isminimize = 1;
  ismaximize = 0;
end

if isfield(z, 'indices')
	sparseIndices = 1;
else
	sparseIndices = 0;
end

if isfield(z, 'selectOutput')
	output = z.selectOutput;
else
	output = 1;
end

xoptmin = []; xoptmax = [];
ymin = []; ymax = [];

if isminimize, ymin = inf; end
if ismaximize, ymax = -inf; end

d = z.d;

n = size(z.vals, 2);

for k = 1:n
	% Get the sparse grid points
	y = [];
	if isfield(z, 'fvals')
		y = z.fvals{output, k};
	end	
	if isfield(z, 'grid')
		sx = z.grid{k};
	else
		if ~sparseIndices
			sx = spgrid(k-1,d,spset('GridType', z.gridType));
		else
			% Make sure that the entire grid is returned!
			% TODO: This should be broken down in the future for very
      % high-dimensional interpolants and many support nodes, since
      % it may require to much memory do store the entire grid in
      % one single array.
			z.indices.currentindex = 1;
			sx = spgrid(z.indices(k,:),[],spset('GridType', z.gridType));
		end
		for l = 1:d
			% Rescale sparse grid to actual range
			sx(:,l) = z.range(l,1) + (z.range(l,2)-z.range(l,1)).*sx(:,l);
		end
	end
	% Crop all points outside of search box
	skip = 0;
	np = size(sx,1);
	id = zeros(np,1);
	nid = 0;
	for k = 1:np
	  crop = 0;
		for l = 1:d
			if sx(k,l) < xbox(l,1) - 10*eps*(z.range(l,2)-z.range(l,1)) ...
			   || sx(k,l) > xbox(l,2) + 10*eps*(z.range(l,2)-z.range(l,1))
			  crop = 1;
				break;
			end
		end
		if crop == 0
			nid = nid + 1;
			id(nid) = k;
		end
	end
  sx = sx(id(1:nid),:);
  
	if ~isempty(sx)
		if ~isempty(y)
			y = y(id(1:nid));
		else
			sxcell = num2cell(sx,1);
			y = spinterp(z, sxcell{:});
		end
	else
		% no points for this index set; continue with next 
		skip = 1;
		break;
	end
	
	if ~skip
		if isminimize
			[ymintemp, id] = min(y);
			if ymintemp < ymin
				ymin = ymintemp;
				xoptmin = sx(id,:);
			end
		end
		if ismaximize
			[ymaxtemp, id] = max(y);
			if ymaxtemp > ymax
				ymax = ymaxtemp;
				xoptmax = sx(id,:);
			end
		end
	end
end

xopt = [xoptmin' xoptmax'];
fval = [ymin ymax];
