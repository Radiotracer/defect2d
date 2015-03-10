function [ val ] = actConstraint(p)
% This function computes the activity constraint to suppress noise.
%  Parameters: 
% (a) Center :p(1),p(2)
% (b) Blood pool activity; p(3)
% (c) Background activity; p(4)
% (d) Endocardium radii and corresponding myocaridum thicknesses(#=nRad):
%     [p(5),p(6)] --> [p(5+2*(nRad-1), p(6+2*(nRad-1)]
% (e) Activities of myocardium segments (#=nSeg):
%  p(4+2*nRad+1)-->p(4+2*nRad+nSeg)
global nRad;
global nSeg;
global delta;
p=p(:);
act=p(4+2*nRad+1:4+2*nRad+nSeg);
act=[act;act(1)];
dact=diff(act);
val=sum(clqegyGM(abs(dact),delta));

end

