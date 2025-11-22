function xdash = fftdiff(x, k, fs, filtf)

% computes the kth derivative of x. If k is negative, does integral of
% order -k.
% x is to one window of a periodic function
% length(x) is to be an even number
% fs = sampling frequency
% filtf = [flow fhigh]; components of xdash outside this frequency window
%                       are set to zero

x = x(:); % in case it is a row vector, turn it into a column vector
n = length(x);

f = (-n/2:n/2-1)'*(fs/n);
iw = 1i*2*pi*f;

fftx = fftshift(fft(x));
fftxdash = (iw).^(k).*fftx;

% if integrating with non-zero DC comp in x
fftxdash(isnan(fftxdash)|isinf(fftxdash)) = 0; 

% apply filter
if nargin<4 % filtf not given
    filtf = [0 fs/2];
end
fftxdash((abs(f)<filtf(1))|(abs(f)>filtf(2))) = 0;

xdash = ifft(ifftshift(fftxdash));
xdash = real(xdash); % in case there is any imaginary part due to roundoff
