/******************************************************************************/
/* Scilab script to solve_equa equation from exo1 via finite element method.  */
/* The problem to solve_equa is -u"=f, u'(α)=u'(β)=0, mean value = 0.         */
/* Here, we take f(x)=π^2cos(πx). The exact solution is cos(πx).              */
/******************************************************************************/

/* Parameters of the problem */
pi = 4*atan(1); // Definition of pi
a  = 0;         // Left bound of the domain
b  = 1;         // Right bound of the domain

// Function for f:
function y = f(x)
    y = pi^2*cos(pi*x);
endfunction

// Function for the exact solution:
function y = exact_sol(x)
    y = cos(pi*x);
endfunction

// Function to solve the problem for a given N:
function [X, Uh, h] = solve_equa(N);
    h = 1/N; // Size of a cell (uniform mesh)

    /* Construction of the matrices */
    v = ones(1, N-1); // Useful to construct Kh

    // Construction of the stiffness matrix:
    Kh           = 2*eye(N, N) - diag(v, -1) - diag(v, 1);
    Kh(1, 1)     = 1;
    Kh(N, N)     = 4;
    Kh(N, 1:N-1) = Kh(N, 1:N-1) + 2*v;
    Kh(N, 1)     = 1;
    Kh           = 1/h * Kh;

    // Construction of the matrix Ah of the system to be solve_equad:
    Ah = Kh;

    // Construction of the rhs Lh of the system to be solve_equad:
    X        = a:h:b; // x-coordinate of the unknowns
    Lh       = h*f(X(1, 1:N));
    Lh(1, 1) = 0.5*Lh(1, 1);

    Uh         = linsolve(Ah, -Lh'); // Solution of the system Ah x = Lh
    Uh(N+1, 1) = - (Uh(1, 1)+2*sum(Uh(2:N, 1)));
    Uh = Uh'
endfunction

/******************************************************************************/
/* Visual comparison of the calculated solution with the exact solution.      */
/******************************************************************************/

N          = 10;            // N = Number of cells (N+1 points)
[X, Uh, h] = solve_equa(N); // Computed solution

// Plot of the solution:
//plot(X, Uh, a:0.01:b, exact_sol(a:0.01:b));
//legend(["$y=u_{h}(x)$", "$y=u(x)$"]);

// Computation of the mean value:
L  = Uh(1:N);
R  = Uh(2:N+1);
mv = h/2*sum(L+R);
printf("Mean value of the computed solution = %f\n", mv);

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
    y = -pi*sin(pi*x);
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
