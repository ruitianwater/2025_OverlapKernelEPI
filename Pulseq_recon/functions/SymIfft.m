function result = SymIfft( data, dims )
% Applies "symmetric" fft to selected dimensions of a matrix.
% Aassumes central point in data domain and sets
% the zero point in the Fourier domain at 
% size/2 + 1 for even sizes and (size+1)/2 for odd.
% Selected dimensions are given in the vector dims.

if nargin < 2
    dims = 1:ndims(data);
end

result = data;
for n=dims
    result = fftshift(ifft(ifftshift(result,n),[],n),n);

end

