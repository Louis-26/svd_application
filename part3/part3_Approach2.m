vid = VideoReader("test_videos\2560_1440\street.mp4");
%A = zeros(vid.Height*vid.Width,vid.NumFrames);
tic
%for frame = 1:vid.NumFrames
%     pic = rgb2gray(read(vid,frame));
%     A(:,frame) = vec(pic);
%end
%[U,S,V] = svds(A,1);
%B = uint8(reshape(S*V(1)*U,size(pic)));
imshow(B)
toc