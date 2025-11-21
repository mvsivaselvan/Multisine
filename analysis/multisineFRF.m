function [frf, coh] = multisineFRF(x_, y_, numFFTpts, numFreqExcite)

x = x_ - mean(x_); % input
y = y_ - mean(y_); % output

nwindows = floor(length(x)/numFFTpts);

Sx = repmat(complex(0,0),numFreqExcite,1);
Sy = repmat(complex(0,0),numFreqExcite,1);
Sxx = repmat(complex(0,0),numFreqExcite,1);
Syy = repmat(complex(0,0),numFreqExcite,1);
Sxy = repmat(complex(0,0),numFreqExcite,1);

selct = numFFTpts/2+1+[1:numFreqExcite]; % +1 to ignore DC comp

for n = 1:nwindows
    range = (n-1)*numFFTpts+[1:numFFTpts];
    X = fftshift(fft(x(range)));
    X = X(selct)*(2/numFFTpts);
    Y = fftshift(fft(y(range)));
    Y = Y(selct)*(2/numFFTpts);
    
    Sx = Sx + (X-Sx)/n;
    Sy = Sy + (Y-Sy)/n;
    Sxx = Sxx + (abs(X).^2-Sxx)/n;
    Syy = Syy + (abs(Y).^2-Syy)/n;
    Sxy = Sxy + (conj(X).*Y-Sxy)/n;
end

frf = Sy./Sx;
coh = abs(Sxy).^2./(Sxx.*Syy);
