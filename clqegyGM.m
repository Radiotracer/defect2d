function [ egy ] = clqegyGM(r,delta )
% The clique energy function by Geman and McClure
% Reference:
% Lalush, David S., and Benjamin MW Tsui.
%"Simulation evaluation of Gibbs prior distributions for use in maximum a posteriori SPECT reconstructions." 
% Medical Imaging, IEEE Transactions on 11.2 (1992): 267-275.1

if isscalar(delta) && delta>0 && min(r(:))>=0
    egy=zeros(size(r));
    for k=1:numel(r)
        if r(k)<delta
            egy(k)=8*delta^2*r(k)^2/(3*sqrt(3)*(r(k)^2+delta^2));
        else
            egy(k)=0.0;
        end
    end    
else
    error('All entries of r must be non-negative and delta must be positive scalar');
end

end

