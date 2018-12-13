function NCC = sng_NCC(I1,I2)
%% computes normalized correlation coefficient
%{
    I1 = ImageSlice{k4}(:,:,1);
    I2 = ImageSlice{k4}(:,:,3);
%}

I1 = double(I1);
I2 = double(I2);

a = sum(I1(:).*I2(:));

b = sum(I1(:).^2);
c = sum(I2(:).^2);

NCC = a/sqrt(b*c);

end
