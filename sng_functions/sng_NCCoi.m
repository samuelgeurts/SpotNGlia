function NCC = sng_NCCoi(I1,I2,maskfull)
%% computes normalized correlation coefficient
%{
    I1 = ImageSlice{k4}(:,:,1);
    I2 = ImageSlice{k4}(:,:,3);
%}

[s1,s2,s3] = size(I1);

if nargin < 3
    maskfull = true(s1,s2);
end

maskfull = repmat(maskfull,1,1,s3/size(maskfull,3));    

n = sum(maskfull(:));

I1b = reshape(double(I1(maskfull)),[n/s3,s3]);
I2b = reshape(double(I2(maskfull)),[n/s3,s3]);

mean1 = mean(I1b);
mean2 = mean(I2b);

for k = 1:s3
    I1b(:,k) = I1b(:,k) - mean1(k);
    I2b(:,k) = I2b(:,k) - mean2(k);
end
    
a = sum(I1b(:).*I2b(:));
b = sum(I1b(:).^2);
c = sum(I2b(:).^2);

NCC = a/sqrt(b*c);
end


