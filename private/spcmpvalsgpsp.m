function ip = spcmpvalsgpsp(d,z,y,seq,fromindex,toindex)
% SPCMPVALSGPSP   Compute surpluses, Gauss-Patterson grid, 
%    sparse indices (internal function)

ip = spcmpvalscbgpsp(d,z,y,seq,fromindex,toindex,true);
