function [u] = FeedbackLinearizationController(X,q_ref,omega_ref,COM_X,COM_Y,COM_Z,Jbx,Jby,Jbz,Jk,Jw,Mb,Mk,Bvb,Bvk,Bvm,rk,rw,g)

    %% Function handlers and variables
    beta = 0;
    Phi = @(q)[q(1) -q(2) -q(3) -q(4);     % for q o p = Phi(q) * p
              q(2) q(1)  -q(4) q(3);
              q(3) q(4)  q(1)  -q(2);
              q(4) -q(3) q(2)  q(1)];           
    Gamma = @(p)[p(1) -p(2) -p(3) -p(4);   % for q o p = Gamma(p) * q
                 p(2) p(1) p(4) -p(3);
                 p(3) -p(4) p(1) p(2);
                 p(4) p(3) -p(2) p(1)];  


    %Weight of the poles for linear controller
    p = 100;
    s = 0.5*p;
    pos = [p p p]; %q2 q3 q4
    speed = [s s s];
    k = [pos speed]; 

    %Quaternion assignment
    q1 = X(3);
    q2 = X(4);
    q3 = X(5);
    q4 = X(6);
    q = [q1 q2 q3 q4]';

    dx = X(7);
    dy = X(8);
    dq1 = X(9);
    dq2 = X(10);
    dq3 = X(11);
    dq4 = X(12);
    dq = [dq1 dq2 dq3 dq4]';


    dchi = reshape(X(7:12), 6, 1);

    % Specify variable sizes by initializing to zero
    q_err = [1;0;0;0];
    dq_ref = [0;0;0;0];
    
    % Quaternion error based on reference input    
    q_err = Phi(q_ref)' * q; % quaternion error in body frame    

    if (q_err(1) < 0)
       q_err = -q_err; % inverse quaternion gives the shorter way to the same rotation
    end

    % Recompute q_ref based on q_err (makes a difference if q_err was negated)    
    q_ref = Gamma(q_err)' * q;
    dq_ref = 1/2 * Phi(q_ref) * [0;omega_ref]; % body angular velocity
    
    % Quaternion derivative error based on reference
    dq_err = Phi(dq_ref)'*q + Phi(q_ref)'*dq;
    % TODO, Consider computing dqref based in dq err, otherwise,why this?
    
    
    M = mass(COM_X,COM_Y,COM_Z,Jbx,Jby,Jbz,Jk,Jw,Mb,Mk,q1,q2,q3,q4,rk,rw);
    C = coriolis(COM_X,COM_Y,COM_Z,Jbx,Jby,Jbz,Jw,Mb,beta,dq1,dq2,dq3,dq4,dx,dy,q1,q2,q3,q4,rk,rw);
    G = gravity(COM_X,COM_Y,COM_Z,Mb,beta,g,q1,q2,q3,q4);
    D = friction(Bvb,Bvk,Bvm,beta,dq1,dq2,dq3,dq4,dx,dy,q1,q2,q3,q4,rk,rw);
    Q = input_forces(q1,q2,q3,q4,rk,rw);

    fq = [zeros(4,2), eye(4)] * (M \ (-C*dchi - G - D)); % M \ b = inv(M)*b
    gq = [zeros(4,2), eye(4)] * (M \ Q);
    
    %Consider reconstruct this with the proper notation from the report
    fc = [dchi;fq];
    gc = [zeros(4,3);gq];

    %Old was fc 6/7/8. Changed to 8 9 10
    D = [fc(8);fc(9);fc(10)];

    E = [gc(6,1), gc(6,2), gc(6,3);
        gc(7,1), gc(7,2), gc(7,3);
        gc(8,1), gc(8,2), gc(8,3)];

    % Old used to be ref
    q1 = q_err(1);
    q2 = q_err(2);
    q3 = q_err(3);
    q4 = q_err(4);
    dq1 = dq_err(1);
    dq2 = dq_err(2);
    dq3 = dq_err(3);
    dq4 = dq_err(4);
    %% Evaluate if it is done with q_err or ref From thesis it was q_err
    e = [q_err;dq_err];

    v = [-k(1)*q2-k(4)*dq2;
        -k(2)*q3-k(5)*dq3;
        -k(3)*q4-k(6)*dq4];
    
    u = E\(v-D);