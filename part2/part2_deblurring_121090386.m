disp_picture("256_256_buildings.png",110,160,210,3)
% disp_picture("512_512_stars_01.png",200,300,400)
% disp_picture("640_640_carousel.png",300,400,500)

% use different truncation value
function disp_picture(fileName,trunc_1,trunc_2,trunc_3,mode)
    if(~exist('mode','var'))
        mode = 1; 
    end
    time_list_picture=[];
    PSNR_list_picture=[];
    [time,PSNR]=deblur_and_recover_picture(fileName,trunc_1,1,14,mode);
    time_list_picture=[time_list_picture,time];
    PSNR_list_picture=[PSNR_list_picture,PSNR];
    [time,PSNR]=deblur_and_recover_picture(fileName,trunc_2,2,[5,5,0.6,0.3,0.1],mode);
    time_list_picture=[time_list_picture,time];
    PSNR_list_picture=[PSNR_list_picture,PSNR];
    [time,PSNR]=deblur_and_recover_picture(fileName,trunc_3,3,2,mode);
    time_list_picture=[time_list_picture,time];
    PSNR_list_picture=[PSNR_list_picture,PSNR];
    fprintf("The time for these three truncation numbers is:\n")
    disp(time_list_picture)
    fprintf("The PSNR for these three truncation numbers is:\n")
    disp(PSNR_list_picture)
end

function [T1,T2]=get_T(k1,k2,X_name)
    % get the matrix T
    % T1, T2 mean the blurring kernel for Al and Ar
    % n1,n2,n3 mean the parameter for method 2
    n=get_size_of_X(X_name);
    if k1==1 % meaning we use method 1
        % here k2 denotes the prescription
        a=linspace(1/4.1,1/4.1,n-1);
        d=linspace(2.1/4.1,2.1/4.1,n);
        T=zeros(n);
        for i = 1:(n-1)
            T(i,i)=d(i);
            T(i+1,i)=a(i);
            T(i,i+1)=a(i);
        end
        T(n,n)=d(n);
        T1=T^k2;
        T2=T^k2;
    elseif k1==2 % use method 2, k2=(kl,kr,n1,n2,n3)
        kl=k2(1);% prescription for T1 kernel
        kr=k2(2);% prescription for T2 kernel
        n1=k2(3);% number for the central diagonal 
        n2=k2(4);% number for the diagonal below the central diagonal for T1
        n3=k2(5);% number for the diagonal above the central diagonal for T1
        a=linspace(n2,n2,n-1);
        a1=linspace(n3,n3,n-1);
        d=linspace(n1,n1,n);
        T1=zeros(n);
        T2=zeros(n);
        for i = 1:(n-1)
            T1(i,i)=d(i);
            T2(i,i)=d(i);
            T1(i+1,i)=a(i);
            T2(i+1,i)=a1(i);
            T2(i,i+1)=a(i);
            T1(i,i+1)=a1(i);
        end
        T1(n,n)=d(n);
        T2(n,n)=d(n);
        T1=T1^kl;
        T2=T1^kr;
    elseif k1==3 % use geometric sequence
        % here k2 denotes the ratio of the sequence
        an=(1-1/k2)/(2*(1-(1/k2)^n)-(1-1/k2));
        T=zeros(n);
        for i = 0:(n-1)
            for j=1:(n-i)
                T(j,j+i)=an/k2^i;
                T(j+i,j)=an/k2^i;
            end
        end
        T1=T;
        T2=T;
    end
end

function [time,PSNR]=deblur_and_recover_picture(X_name,trunc,k1,k2,mode) 
    tic % mode: denote the source of svd, 
    % 1: built-in, 2: svd from problem 1 phase II-A, 3: svd from problem 1 phase II-B
    [T1,T2]=get_T(k1,k2,X_name);
    % k1: denote the version of blurring kernel

    % k2: denote the prescription
    % get the blurry kernel A_l and A_r
    A_l=T1;
    A_r=T2;
    n=size(A_l);
    n=n(1);

    % use built-in function
    if mode==1
        fprintf("When we use the svd result in built-in function svd\n")
        % get SVD decomposition of A_l, A_r from built-in svd 
        [U_l,sigma_l,V_l]=svd(A_l);
        [U_r,sigma_r,V_r]=svd(A_r);    
    % use function by phase II-A approach
    elseif mode==2
        fprintf("When we use the svd result in problem 1 phase II-A\n")
        [U_l,sigma_l,V_l]=mySVD1(A_l);
        [U_r,sigma_r,V_r]=mySVD1(A_r);
    % use function by phase II-B appraoch
    else
        fprintf("When we use the svd result in problem 2 phase II-B\n")
        [U_l,sigma_l,V_l]=mySVD2(A_l);
        [U_r,sigma_r,V_r]=mySVD2(A_r);
    end
    [time,PSNR]=get_recover_picture(X_name,A_l,A_r,U_l,sigma_l,V_l,U_r,sigma_r,V_r,trunc);

end

function [time,PSNR]=get_recover_picture(X_name,A_l,A_r,U_l,sigma_l,V_l,U_r,sigma_r,V_r,trunc)
    % print the sigular values of A_l,A_r
    fprintf("The sigular values of A_l are\n")
    n=get_size_of_X(X_name);
    for i=1:n
        disp(sigma_l(i,i))
    end
    
    fprintf("The sigular values of A_r are\n")
    for i=1:n
        disp(sigma_r(i,i))
    end
    
    % determine the value of truncation number
    l_trunc=trunc;
    r_trunc=trunc;
    %compute the pseudoinverse of A_l, A_r
    A_l_pseudoinverse=zeros(n);
    A_r_pseudoinverse=zeros(n);
    for i=1:l_trunc
        A_l_pseudoinverse=A_l_pseudoinverse+V_l(:,i)*U_l(:,i)'/sigma_l(i,i);
    end
    for i=1:r_trunc
        A_r_pseudoinverse=A_r_pseudoinverse+V_r(:,i)*U_r(:,i)'/sigma_r(i,i);
    end

    % In the following operation, the first picture is the orginal clear one
    % the second picture is the blurry one,
    % and the third picture is the recovered one by SVD operation

    % get the matrix regarding the picture
    X=imread(X_name);
    % convert the form of X
    X=double(X)/255;
    % show the original clear image
    figure; imshow(X)
    % get the blurry picture
    B=[];
    if colorful(X_name)
        for i=1:3    
            B(:,:,i)=A_l*X(:,:,i)*A_r;
        end
    else
        B=A_l*X*A_r;
    end
    % show the blurry image
    figure;imshow(B)
    % recover the original image
    X_trunc=[];
    if colorful(X_name)
        for i=1:3    
            X_trunc(:,:,i)=A_l_pseudoinverse*B(:,:,i)*A_r_pseudoinverse;
        end
    else
        X_trunc=A_l_pseudoinverse*B*A_r_pseudoinverse;
    end
    fprintf("When the truncation value is %d\n",trunc)
    % show the recovered image
    figure;imshow(X_trunc)
    % measure the quality of reconstruction by peak-signal-to-noise ratio(PSNR)
    toc
    time=toc;
    PSNR=10*log10(n^2/norm(X_trunc-X,"fro")^2);
    fprintf("The computed PSNR is %f\n",PSNR)
end

function value=get_size_of_X(X_name)
    X=imread(X_name);
    n=size(X);
    value=n(1);
end

function judge=colorful(X_name)
    X=imread(X_name);
    num_of_dim=ndims(X);
    if num_of_dim==2
        judge=false;
    else
        judge=true;
    end
end