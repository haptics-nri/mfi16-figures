% DFT321 on Nx3 arrays
% also returns freq. domain data
function [out, freq] = dft321(in, filt)

    % unity-gain first-order high-pass filter
    a = [.02-1 0];
    b = [1 .02-1];

    if nargin > 1 && filt
        A = fft(filtfilt(b, a, in), 2*size(in,1));
    else
        A = fft(bsxfun(@minus, in, mean(in)), 2*size(in,1));
    end
    A_s = sqrt(sum(A .* conj(A), 2));
    theta = angle(sum(A,2));
    
    freq = A_s .* exp(1j * theta);
    out = real(ifft(freq, size(in,1)));
    freq = abs(freq(end-size(in,1):end));

end
