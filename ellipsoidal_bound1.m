function [ P_ellipsoid ] = ellipsoidal_bound1( system, As, attack_type )
% system: system structure containing F,G,C,K,mu1,R1,mu2,R2
% A: false alarm rate
% attack_ztype = {'none', 'zero alarm', 'hidden'}
t=max(size(system.F));
m=min(size(system.L));
if strcmp(attack_type,'none')
    r=(ncx2inv(1-As,t,0));
    [VR1,ER1]=eig(system.R1);
    A_transfer=VR1'*sqrtm(ER1);
    A1=0.01:0.01:0.99;
    for i=1:max(size(A1))
        clear p Pi sol SolverInfo LMI
        A=system.F;
        B=eye(t);
        R=inv(system.R1)/r;
        a = A1(i);
        p = sdpvar(t,t);
        systemlmi=lmi([a*p-A'*p*A  -A'*p*B; -B'*p*A  (1-a)*R-B'*p*B])
        systemlmi=addlmi(systemlmi,'p>=0');
        h=-logdet(p);
        sol = solvesdp(systemlmi,h);
        SolverInfo = sol.info;
        P(:,:,i) = double(p);
        Pi = inv(P(:,:,i));
        if sol.problem == 1;
            Vol(i) = NaN;
        else sol.problem == 0;
            Vol(i) = det(Pi);
        end
    end
    [Y,k]=nanmin(Vol);
     P_ellipsoid=P(:,:,k);
%     A1=0.01:0.01:0.99;
%     for i=1:max(size(A1))
%         clear p Pi sol SolverInfo LMI
%         A=system.F+system.G*system.K;
%         B=[-system.G*system.K eye(t)];
%         R=0.5*[Pe zeros(t,t);zeros(t,t) inv(system.R1)/r]
%         clear p Pi sol SolverInfo LMI
%         a = A1(i);
%         p = sdpvar(t,t);
%         systemlmi=lmi([[a*p-A'*p*A  -A'*p*B; -B'*p*A  (1-a)*R-B'*p*B]])
%         systemlmi=addlmi(systemlmi,'p>=0');
%         h=-logdet(p);
%         sol = solvesdp(systemlmi,h);
%         SolverInfo = sol.info;
%         P(:,:,i) = double(p);
%         Pi = inv(P(:,:,i));
%         if sol.problem == 1;
%             Vol(i) = NaN;
%         else sol.problem == 0;
%             Vol(i) = det(Pi);
%         end
%     end
%     [Y,k]=nanmin(Vol);
%     P_ellipsoid=P(:,:,k);
    
elseif strcmp(attack_type,'zero alarm only') 
    kappa = system.alpha;
    A1=0.01:0.01:0.99;
    for i=1:max(size(A1))
        clear p Pi sol SolverInfo LMI
        A=system.F;
        B=-system.L*sqrtm(system.Sigma);
        R=eye(t)/kappa;
        a = A1(i);
        p = sdpvar(t,t);
        systemlmi=lmi([a*p-A'*p*A  -A'*p*B; -B'*p*A  (1-a)*R-B'*p*B])
        systemlmi=addlmi(systemlmi,'p>=0');
        h=-logdet(p);
        sol = solvesdp(systemlmi,h);
        SolverInfo = sol.info;
        P(:,:,i) = double(p);
        Pi = inv(P(:,:,i));
        if sol.problem == 1;
            Vol(i) = NaN;
        else sol.problem == 0;
            Vol(i) = det(Pi);
        end
    end
    [Y,k]=nanmin(Vol);
    Pe=P(:,:,k);
        A=system.F+system.G*system.K;
        B=-system.G*system.K;
        R=Pe;
        
    A1=0.01:0.01:0.99;
    for i=1:max(size(A1))
        clear p Pi sol SolverInfo LMI
        a = A1(i);
        p = sdpvar(t,t);
        systemlmi=lmi([[a*p-A'*p*A  -A'*p*B; -B'*p*A  (1-a)*R-B'*p*B]])
        systemlmi=addlmi(systemlmi,'p>=0');
        h=-logdet(p);
        sol = solvesdp(systemlmi,h);
        SolverInfo = sol.info;
        P(:,:,i) = double(p);
        Pi = inv(P(:,:,i));
        if sol.problem == 1;
            Vol(i) = NaN;
        else sol.problem == 0;
            Vol(i) = det(Pi);
        end
    end
    [Y,k]=nanmin(Vol);
    P_ellipsoid=P(:,:,k);
elseif strcmp(attack_type,'zero alarm Total')
        P_noise     = ellipsoidal_bound1( system, As, 'none' );
        P_attack    = ellipsoidal_bound1( system, As, 'zero alarm only' );
        struct(1).p = inv(P_noise);
        struct(2).p = inv(P_attack);
        P_ellipsoid    =  minkowskisumn(struct);
end