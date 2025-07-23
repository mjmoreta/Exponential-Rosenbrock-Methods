% Program that implements the exponential Rosenbrock Euler method avoiding
% the order reduction with p=2 to solve problem
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(t,1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)
% For the spatial discretization we have used the standard second-order
% difference squeme. 
% Here, we use subrutine phipm

% This exponential Rosenbrock Euler method appears in 
% M. Caliari and A. Ostermann, Implementation of exponential Rosenbrock-type 
% integrators, Appl. Num. Math. 59 (2009) 568-581.

% Some initial values, such as h, matrix A_h0 and some more that are used
% several times along the program
% N is such that h/1/N is the grid diameter in [0,1]
N=1000;
A=-2*diag([ones(N-1,1)])+diag (ones (N-2,1),1)+diag (ones (N-2,1),-1);
h=1/N;
h2=h^2;
A=A./(h^2);
A=sparse(A);
x=[h:h:1-h]';
hdiv=1/h2;

% n is such that the time step size is k=1/n.
n=5;

% The program runs for 8 different values of k, from k=1/5, in order to calculate the
% error and the order of the method
for ll=1:8

    k=1/n;

    % U is the approximation to the exact solucion at time t_n. The initial 
    % u(0,x) value is known
    U=zeros(N-1,1);
    for ii=1:N-1
        U(ii)=cos(x(ii));
    end

    t=0;
    
    Phh=zeros(N-1,1);
    vecb=zeros(N-1,4);

    % CPU time calculus starts. 
    tstart=tic;

    % For the local error r=1. For the global one, r=n
    for r=1:n
      
        % J is the Jacobian  
        J=A+2*diag(U);
        J=sparse(J);
          
        % funh is funcion h(t,x)
        for ii=1:N-1
            Phh(ii,1)=funh(t,x(ii));       
        end

        % vecb contains the vectors that are multiplied by the
        % exponential funcions. The first column is multiplied by the 
        % exponential, the second one is multiplied by phi_1(kJ), 
        % the third one by phi_2(kJ) and the last one by phi_3(kJ).

        vecb(:,1)=U;
    
        vecb(:,2)=-U.^2+Phh;
        vecb(1,2)=vecb(1,2)+hdiv*uu(t,0);
        vecb(N-1,2)=vecb(N-1,2)+hdiv*uu(t,1);
    
        vecb(1,3)=funht(t,x(1))+hdiv*uut(t,0);
        for ii=2:N-2
            vecb(ii,3)=funht(t,x(ii));         
        end
        vecb(N-1,3)=funht(t,x(N-1))+hdiv*uut(t,1);

        vecb(1,4)=hdiv*uutt(t,0);
        vecb(N-1,4)=hdiv*uutt(t,1);

        % Aproximation to the exact solution at time t_{n+1}, by using
        % soubrutine phipm
        U=phipm(k,J,vecb,10^(-8),1,10); 

        t=t+k;        
            
    end

    % CPU time calculus finishes
    telapsed=toc(tstart)

    % Sol contains the exact solution at time T
    sol=zeros(N-1,1);

    % Error in the infinite norm
    for ii=1:N-1
        sol(ii)=cos(x(ii)+t);
    end

    err=norm(sol-U,inf);

    % Order. When ll=1, as there is not a previous error, it can't be calculated. 
    % Two consecutive errors are compared
    if ll==1
        err
        err0=err;
    else
        [err log2(err0/err)]
        err0=err;
    end

    % The new value of n is 2*n, and the new value of k_n=k/2. It is
    % calculated at the beginning os the next iteration
    n=2*n;

end
