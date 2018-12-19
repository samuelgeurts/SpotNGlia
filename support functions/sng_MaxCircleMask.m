function cm = sng_MaxCircleMask(I40)
   %creates a mask with the largest circle possible in an rectangular
   %box

[iy,ix,iz] = size(I40);
cx = (ix/2)+0.5;
cy = (iy/2)+0.5;
r = min(ix,iy)/2; 
[x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
c_mask=((x.^2+y.^2)<=r^2);
cm = repmat((c_mask),1,1,iz);
end