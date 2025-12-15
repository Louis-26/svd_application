function [U,Sigma,V] = mySVD2(A)
    %A = [1,2,3,4,5;6,7,8,9,0;1,2,3,4,5];
    %A = randn(5,5);
    %A = [1,2,3,4,5;1,2,3,4,5;1,2,3,4,5;1,2,3,4,5;1,2,3,4,5];
    [m, n] = size(A);
    U = eye(m);
    V = eye(n);     
    B = A;
    
   
    for k = 1:(max(m,n)-2)
        if k <= min(m-1, n)
            w = B(k:m, k);
            alpha = -sign(w(1))*norm(w);
            w(1) = w(1) - alpha;
            w = w/norm(w);
            B(k:m, :) = B(k:m, :) - 2*w*(w'*B(k:m, :));
            U(:, k:m) = U(:, k:m) - 2*(U(:, k:m)*w)*w';
        end
        if k <= min(n-1, m)
            x = B(k, (k+1):n);
            w = x;
            w(1) = x(1) + norm(x)*sign(x(1));
            w = w/norm(w);
            B(:, (k+1):n) = B(:, (k+1):n) - 2*(B(:, (k+1):n)*w')*w;
            V(:, (k+1):n) = V(:, (k+1):n) - 2*(V(:, (k+1):n)*w')*w;
        end        
    end
    
    % Phase II QR iteration with Wilkinson shift
    % Phase II B
    k = min(m,n);
    Sigma = zeros(m, n);
    if m>n 
        X = B;
        while (abs(X(1, 2)) + abs(X(k-1, k))) > 1e-12 
            [~, Rk] = qr(X'); 
            Lk = chol(Rk*Rk');
            X = Lk';
        end
        Sigma(1:n,1:n) = X; %The Sigma on here is actually the same value with the later Sigma(line 77).
        
        
    else
        X = B';
        while (abs(X(2, 1)) + abs(X(k,k-1))) > 1e-12 
            [~, Rk] = qr(X');
            Lk = chol(Rk*Rk');
            X = Lk';
        end
        Sigma(1:m,1:m) = X';
        
    end 
    T = B'*B; %tridiagonal matrix 
    Q0 = eye(k); 
    while (abs(T(1, 2)) + abs(T(k-1, k))) > 1e-12 
        [Q, R] = qr(T); 
        T = R*Q;
        Q0 = Q0*Q;
    end
    s = sqrt(diag(T)); 
    [~, Idx]=sort(s, 'descend'); 
    % find V
    V = V*Q0(:, Idx);

    T = B*B';
    Q0 = eye(k);
    while (abs(T(1, 2)) + abs(T(k-1, k))) > 1e-12 
        [Q, R] = qr(T);
        T = R*Q;
        Q0 = Q0*Q;
    end
    s = sqrt(diag(T)); 
    [~, Idx]=sort(s, 'descend');
    % find U
    U = U*Q0(Idx, :);
    Sigma = U'*A*V;
    for i = 1:min(m,n)
        if Sigma(i,i) < 0
            U(:, i) = -U(:, i);
            Sigma(i,i) = -Sigma(i,i);
        end
    end
end
