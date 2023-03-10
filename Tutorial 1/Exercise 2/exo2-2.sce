/******************************************************************************/
/* Scilab script to solve equation from exo2 via finite element method.       */
/* More precisely, the problem to solve is -u"+u=f, u(0)=u(1), u'(0)=u'(1).   */
/* Here, we take f(x)=(4π^2+1)*sin(2πx).                                      */
/******************************************************************************/

/* Parameters of the problem */
pi = 4*atan(1); // Definition of pi
a  = 0;         // Left bound of the domain
b  = 1;         // Right bound of the domain

// Function for f:
function y = f(x)
    y = (4*pi^2+1)*sin(2*pi*x);
endfunction

// Function for the exact solution:
function y = exact_sol(x)
    y = sin(2*pi*x);
endfunction

// Function to solve the problem for a given N:
function [X, Uh, h] = solve_equa(N)
    h = 1/N; // Size of a cell (uniform mesh)

    /* Construction of the matrices */
    v = ones(1, N-1); // Useful to construct Kh and Mh

    // Construction of the stiffness matrix:
    Kh       = 2*eye(N, N) - diag(v, -1) - diag(v, 1);
    Kh(1, N) = -1;
    Kh(N, 1) = -1;
    Kh       = 1/h * Kh;

    // Construction of the mass matrix:
    Mh       = 4*eye(N, N) + diag(v, -1) + diag(v, 1);
    Mh(1, N) = 1;
    Mh(N, 1) = 1;
    Mh       = h/6 * Mh;

    // Construction of the matrix Ah of the system to be solved:
    Ah = Kh + Mh;

    // Construction of the rhs Lh of the system to be solved:
    X   = a:h:b;     // x-coordinate of the unknowns
    pi  = 4*atan(1); // Definition of π
    Xj  = X(1, 1:N);
    Xj1 = X(1, 2:N+1);
    Lh  = (4*pi^2+1)/h * ((sin(2*pi*Xj)-sin(2*pi*Xj1))/(4*pi^2)+h*cos(2*pi*Xj)/(2*pi));
    Int = (4*pi^2+1)/h*((sin(2*pi*Xj1)-sin(2*pi*Xj))/(4*pi^2)-h*cos(2*pi*Xj1)/(2*pi));
    Lh(1, 2:N) = Lh(1, 2:N) + Int(1, 1:N-1);
    Lh(1, 1)   = Lh(1, 1) + Int(1, N);

    // Construction of the rhs Lh of the system to be solved,
    // using the trapezoidal rule:
    //Lh_bis = h*(4*pi^2+1)*sin(2*pi*Xj);

    Uh      = linsolve(Ah, -Lh'); // Solution of the system Ah x = Lh
    Uh(N+1) = Uh(1);
    Uh = Uh'

    //Uh_bis      = linsolve(Ah, -Lh_bis'); // Solution of the system Ah x = Lh_bis
    //Uh_bis(N+1) = Uh_bis(1);
endfunction

/******************************************************************************/
/* Visual comparison of the calculated solution with the exact solution.      */
/******************************************************************************/

N          = 10;            // N = Number of cells (N+1 points)
[X, Uh, h] = solve_equa(N); // Computed solution

// Plot of the solution:
//plot(X, Uh, a:0.01:b, exact_sol(a:0.01:b));
//legend(["$y=u_{h}(x)$", "$y=u(x)$"]);

// Computation of L2 norm of the error:
err = (h*sum(abs(Uh(1:N)-exact_sol(X(1:N)))^2))^0.5;
printf("L2 norm of the error = %f\n", err);

/******************************************************************************/
/* Plot of the error in L2-norm as N grows                                    */
/******************************************************************************/

Liste_N = 10:10:500;
Err     = [];
for N = Liste_N
    [X, Uh, h] = solve_equa(N); // Computed solution
    err        = (h*sum((Uh(1:N)-exact_sol(X(1:N)))^2))^0.5;
    Err        = [Err, err];
end

//plot(log(Liste_N), log(Err), 'k', log(Liste_N), -2*log(Liste_N), '--k');

/******************************************************************************/
/* Plot of the error in H1-seminorm as N grows                                */
/******************************************************************************/

// Function for u':
function y = up(x)
    y = 2*pi*cos(2*pi*x);
endfunction

Liste_N = 10:10:500;
Err     = [];
for N = Liste_N
    [X, Uh, h] = solve_equa(N); // Computed solution

    Up  = up(X(1:N));
    Uhp = (Uh(2:N+1)-Uh(1:N))/h;

    err = (h*sum((Up-Uhp)^2))^0.5;
    Err = [Err, err];
end

//plot(log(Liste_N), log(Err), 'b', log(Liste_N), -log(Liste_N), '--b');
//legend(["L2 error", "Slope = -2", "H1 error", "Slope = -1"]);
//title("Log-log plot of the errors as N grows");
