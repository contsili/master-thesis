function [pnt, dhk] = tetrahedron642()

[pnt, dhk] = tetrahedron;
[pnt, dhk] = refine(pnt, dhk);
[pnt, dhk] = refine(pnt, dhk);
[pnt, dhk] = refine(pnt, dhk);

pnt = pnt ./ repmat(sqrt(sum(pnt.^2,2)), 1,3);