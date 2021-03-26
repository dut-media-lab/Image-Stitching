function [H, A, C1, C2] = vgg_H_from_x_lin(xs1,xs2,A,W,C1,C2)
    % H = vgg_H_from_x_lin(xs1,xs2)
    %
    % Compute H using linear method (see Hartley & Zisserman Alg 3.2 page 92 in
    %                              1st edition, Alg 4.2 page 109 in 2nd edition). 
    % Point preconditioning is inside the function.
    %
    % The format of the xs [p1 p2 p3 ... pn], where each p is a 2 or 3
    % element column vector.

    [r,c] = size(xs1);

    if (size(xs1) ~= size(xs2))
     error ('Input point sets are different sizes!')
    end

    if (size(xs1,1) == 2)
      xs1 = [xs1 ; ones(1,size(xs1,2))];
      xs2 = [xs2 ; ones(1,size(xs2,2))];
    end


    % condition points
    if nargin == 2
        C1 = vgg_conditioner_from_pts(xs1);
        C2 = vgg_conditioner_from_pts(xs2);
        xs1 = vgg_condition_2d(xs1,C1);
        xs2 = vgg_condition_2d(xs2,C2);
    end

    if nargin == 6
        B = A;
        B(1:2:end,:)=W*A(1:2:end,:);
        B(2:2:end,:)=W*A(2:2:end,:);
        
        % Extract nullspace
        [u,s,v] = svd(B, 0); s = diag(s);
    else
        A = [];
        ooo  = zeros(1,3);
        for k=1:c
          p1 = xs1(:,k);
          p2 = xs2(:,k);
          A = [ A;
            p1'*p2(3) ooo -p1'*p2(1)
            ooo p1'*p2(3) -p1'*p2(2)
           ];
        end
        
        % Extract nullspace
        [u,s,v] = svd(A, 0); s = diag(s);
    end

    nullspace_dimension = sum(s < eps * s(1) * 1e3);
    if nullspace_dimension > 1
      fprintf('Nullspace is a bit roomy...');
    end

    h = v(:,9);

    H = reshape(h,3,3)';

    %decondition
    H = C2\H*C1;

    H = H/H(3,3);
end
