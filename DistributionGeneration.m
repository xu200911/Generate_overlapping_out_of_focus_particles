function y = DistributionGeneration(x,mu,sigma,Distribution_type)

if nargin<1
    error(message('stats:lognpdf:TooFewInputs'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end
% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
% Negative data would create complex values, potentially creating spurious
% NaNi's in other elements of y.  Map them, and zeros, to the far right
% tail, whose pdf will be zero.
x(x <= 0) = Inf;

try
 switch Distribution_type   
    case 'LogNormal'
    y = exp(-0.5 * ((log(x) - log(mu))./sigma).^2) ./ (x .* sqrt(2*pi) .* sigma);
    case 'Normal'
    y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);            
    case 'RR'
    y = sigma .* (x/(mu)).^(sigma-1) .* exp(-(x/(mu)).^(sigma)) ./ mu;
    case 'Equality' 
%         x_min = min(x);
%         x_max = max(x);  
        y = ones(length(x),1)/length(x);

 end
catch
    error(message('stats:lognpdf:InputSizeMismatch'));
end
