vid = VideoReader("test_videos\2560_1440\street.mp4");
A = zeros(vid.Height*vid.Width,vid.NumFrames);
tic
for frame = 1:vid.NumFrames
     pic = rgb2gray(read(vid,frame));
     A(:,frame) = vec(pic);
end
[U,S,V] = approach2(A);
B = uint8(reshape(S*V(1)*U,size(pic)));
imshow(B)
toc
function [U,S,V] = approach2(A)
C = A'*A;
[~,n] = size(C);
x = zeros(n,1);
x(1) = 1;
threhold  = 0.001;
iter = 0;
eigenvalue = 0;
rate = 1;
while rate > threhold
    iter = iter+1;
    a = C*x;
    xhat = a/norm(a);
    eigenvalue1 = xhat'*C*xhat;
    rate = abs(eigenvalue1- eigenvalue);
    x = xhat;
    eigenvalue = eigenvalue1;
end
S = rate;
V = x;
U = A*V/S;
end
