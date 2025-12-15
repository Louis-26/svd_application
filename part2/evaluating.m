% This code use 256_256_buildings.png as an example
% measure method1
x=[];
y=[];
z=[];
for i = 240:250
    [time,PSNR]=deblur_and_recover_picture(X_name,Tkl,Tkr,l_trunc,r_trunc,svdme);
    x(i-239)=i;
    y(i-239)=PSNR;
    z(i-239)=time;
end

figure; plot(x,y,'-r*')
xlim([240 250])
ylim([30 50])
figure; plot(x,z,'-r*')
xlim([240 250])
zlim([0 0.5])
% measure method2 and method3
%{
x=[];
y=[];
z=[];
for i = 245:255
    [time,PSNR]=deblur_and_recover_picture(X_name,Tkl,Tkr,l_trunc,r_trunc,svdme);
    x(i-244)=i;
    y(i-244)=PSNR;
    z(i-244)=time;
end

figure; plot(x,y,'-r*')
xlim([245 255])
ylim([30 50])
figure; plot(x,z,'-r*')
xlim([245 255])
zlim([0 0.5])
%}