function [ egy ] = clqegyGR(r,delta )
% The clique energy function by  Green
% Reference:
% Lalush, David S., and Benjamin MW Tsui.
%"Simulation evaluation of Gibbs prior distributions for use in maximum a posteriori SPECT reconstructions." 
% Medical Imaging, IEEE Transactions on 11.2 (1992): 267-275.1

if isscalar(delta) && delta>0 && min(r(:))>=0
    egy=zeros(size(r));
    for k=1:numel(r)
        egy(k)=delta*log(cosh(r(k)/delta));
    end    
else
    error('All entries of r must be non-negative and delta must be positive scalar');
end

end

