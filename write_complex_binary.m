function write_complex_binary (filename, z)

% Divide into real and imaginary components
z_real = real(z);
z_imag = imag(z);

% Preallocate the interleaved array
interleaved = zeros(size(z,1)*2, size(z,2));

% Alternate real and imaginary data
interleaved(1,:) = z_real(1,:);
interleaved(2,:) = z_imag(1,:);

% Write the interleaved values
fid = fopen(filename,'w');
fwrite(fid, interleaved, 'double');
fclose(fid);

% for some reason the converter by Youngjune
% can see one more redundant number at the end
% need to be careful on this

