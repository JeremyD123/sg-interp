function w = spquadwcc(levelseq)
% SPQUADWCC   Compute quadrature weights, Clenshaw-Curtis grid 
%    W = SPQUADWCC(LEVELSEQ)  Computes the quadrature weights
%    for the given sequence of index sets LEVELSEQ. One row of 
%    column vectow W represents one weight. 
%    (Internal function)

% Get the number of levels
nlevels = uint32(size(levelseq,1));

% Get the dimension
d = uint16(size(levelseq,2));

% Init weights vector
w = [];
	
for kl = 1:nlevels
	npoints = 1;
	wval = 1;
	for k = 1:d
		lval = double(levelseq(kl,k));
		if lval == 0
			np = 1;
			nw = 1;
		elseif lval < 3
			np = 2;
			nw = 4;
		else
			np = 2^(lval-1);
			nw = 2^lval;
		end
		wval = wval * nw;
		npoints = npoints * np;
	end
	w = [w; ones(npoints,1) / wval];
end
