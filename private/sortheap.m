function A = sortheap(A, na, G)
% SORTHEAP  Reorder heap after a new element has been added at end
%    (Internal function)

k = na;
while k > 1
	prev = floor(k/2);
	if G(A(prev)) < G(A(k))
		temp = A(prev);
		A(prev) = A(k);
		A(k) = temp;
	else
		break;
	end
	k = floor(k/2);
end
