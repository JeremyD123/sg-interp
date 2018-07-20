function w = spquadwm(levelseq)
% SPQUADWM   Compute quadrature weights, Maximum grid 
%    W = SPQUADWM(LEVELSEQ)  Computes the quadrature weights
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
			np = 3;
			nw = [1/4 1/2 1/4];
		else
			np = 2^lval;
			nw = ones(1,np) * 1/(2^(lval+1));
		end
		wval = wval(:) * nw;
		npoints = npoints * np;
	end
	w = [w; wval(:)];
end
