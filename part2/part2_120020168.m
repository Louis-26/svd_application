%input picture name,(256_256_buildings.png,640_640_lion.png,512_512_town_02.png)
X_name=input('please input the picture name ','s');
n=get_size_of_X(X_name);
%input truncation number
l_trunc=input('please input left truncation number ');
r_trunc=input('please input right truncation number ');

%input the blurry matrix you want
%blurry="original"(method1) or "shift"(method2) or "geometric"(method3);
blurry=input('please input which blurry matrix you want to use ','s');

%input the svd method you want
%svdme="builtin"(matlab built in svd) or "mySVD1"(phase A) or "mySVD2"(phase B);
svdme=input('please input which svd method you want to use ','s');

%Choose the blurry matrix you want
if blurry=="original"
    %input blurry matrix's power(i.e.k)
    kl=input('please input kl ');
    kr=input('please input kr ');
    [Tkl Tkr]=originalTk(n,kl,kr); %original Tk(i.e.method1)
elseif blurry=="shift"
    %input parameter of blurry matrix A_l and A_r
    %n1 = number in the diagonal(i.e.an), 
    %n2 = number in one subdiagonal(i.e.an+1 for A_l)
    %n3 = number in the other subdiagonal(i.e.an-1 for A_l)
    %n2>=n3  n1+n2+n3=1
    n1=input('please input n1');
    n2=input('please input n2');
    n3=input('please input n3');
    %input blurry matrix's power(i.e.kl,kr)
    kl=input('please input kl');
    kr=input('please input kr');
    [Tkl Tkr]=shiftTk(n,kl,kr,n1,n2,n3); %shifted Tkl and Tkr(i.e.method2)
else
    q=input('please input q');
    [Tkl Tkr]=geometricTk(n,q);
end

[time,PSNR]=deblur_and_recover_picture(X_name,Tkl,Tkr,l_trunc,r_trunc,svdme);
function [time,PSNR]=deblur_and_recover_picture(X_name,Tkl,Tkr,l_trunc,r_trunc,svdme)
    n=get_size_of_X(X_name);
    tic
    % get the blurry kernel A_l and A_r
    A_l=Tkl;
    A_r=Tkr;
    
    % get SVD decomposition of A_l, A_r from built-in svd 
    if svdme=="builtin"
        [U_l,sigma_l,V_l]=svd(A_l);
        [U_r,sigma_r,V_r]=svd(A_r);
        fprintf("When we use the svd result in built-in function svd")
    elseif svdme=="mySVD1"
        [U_l,sigma_l,V_l]=mySVD1(A_l);
        [U_r,sigma_r,V_r]=mySVD1(A_r);
        fprintf("When we use the svd result in phaseA")
    else
        [U_l,sigma_l,V_l]=mySVD2(A_l);
        [U_r,sigma_r,V_r]=mySVD2(A_r);
        fprintf("When we use the svd result in phaseB")
    end
    % save the sigular values of A_l,A_r 
    save sigma_l
    save sigma_r
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
    subplot(1,3,1); imshow(X),title('original picture')
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
    subplot(1,3,2);imshow(B),title('blurry picture')
    % recover the original image
    X_trunc=[];
    if colorful(X_name)
        for i=1:3    
            X_trunc(:,:,i)=A_l_pseudoinverse*B(:,:,i)*A_r_pseudoinverse;
        end
    else
        X_trunc=A_l_pseudoinverse*B*A_r_pseudoinverse;
    end
    fprintf("When the left truncation value is %d\n",l_trunc)
    fprintf("When the right truncation value is %d\n",r_trunc)
    % show the recovered image
    subplot(1,3,3); imshow(X_trunc),title('reconstruction picture')
    toc
    time=toc;
    % compute PSNR and cputime
    XF=X_trunc-X;
    if colorful(X_name)
        num1=trace(XF(:,:,1).'*XF(:,:,1));
        num2=trace(XF(:,:,2).'*XF(:,:,2));
        num3=trace(XF(:,:,3).'*XF(:,:,3));
        num=num1+num2+num3;
    else
        num=trace(XF);
    end
    PSNR=10*log10(n^2/num);
    fprintf("The computed PSNR is %f\n",PSNR)
    fprintf("The total time is %f\n",time)
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
function value=get_size_of_X(X_name)
    X=imread(X_name);
    n=size(X);
    value=n(1);
end
function [Tkl Tkr]=originalTk(n,kl,kr)
    a=linspace(1/4.1,1/4.1,n-1);
    d=linspace(2.1/4.1,2.1/4.1,n);
    T=zeros(n);
    for i = 1:(n-1)
        T(i,i)=d(i);
        T(i+1,i)=a(i);
        T(i,i+1)=a(i);
    end
    T(n,n)=d(n);
    Tkl=T;
    Tkr=T;
    for i=1:(kl-1)
        Tkl=Tkl*T;
    end
    for i=1:(kr-1)
        Tkr=Tkr*T;
    end
end
function [Tkl Tkr]=shiftTk(n,kl,kr,n1,n2,n3)
    a=linspace(n2,n2,n-1);
    a1=linspace(n3,n3,n-1);
    d=linspace(n1,n1,n);
    Tl=zeros(n);
    Tr=zeros(n);
    for i = 1:(n-1)
        Tl(i,i)=d(i);
        Tr(i,i)=d(i);
        Tl(i+1,i)=a(i);
        Tr(i+1,i)=a1(i);
        Tr(i,i+1)=a(i);
        Tl(i,i+1)=a1(i);
    end
    Tl(n,n)=d(n);
    Tr(n,n)=d(n);

    Tkl=Tl;
    Tkr=Tr;
    for i=1:(kl-1)
        Tkl=Tkl*Tl;
    end
    for i=1:(kr-1)
        Tkr=Tkr*Tr;
    end
end
function [Tkl Tkr]=geometricTk(n,q)
    %Geometric sequence
    an=(1-1/q)/(2*(1-(1/q)^n)-(1-1/q));
    T=zeros(n);
    for i = 0:(n-1)
        for j=1:(n-i)
            T(j,j+i)=an/q^i;
            T(j+i,j)=an/q^i;
        end
    end
    Tkl=T;
    Tkr=T;
end
