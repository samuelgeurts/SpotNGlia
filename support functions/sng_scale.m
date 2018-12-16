function [Img2,xcom2,ycom2] = sng_scale(Img,xcom,ycom,scale)
%resize image ans scale center of mass coordinates

%scale = 1/8;

Img2 = imresize(Img,scale);
xcom2 = xcom * scale;
ycom2 = ycom * scale;

end
