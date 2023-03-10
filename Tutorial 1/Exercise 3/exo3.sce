/******************************************************************************/
/* Scilab script to solve equation from exo3 via finite element method.       */
/* More precisely, the problem to solve is -u"+cu'=0, u(0)=0, u(1)=1.         */
/******************************************************************************/

// Function for the exact solution:
function y = exact_sol(x, c)
    y = (exp(c*x)-1)/(exp(c)-1)
endfunction

/* Parameters of the problem */
a = 0   // Left bound of the domain
b = 1   // Right bound of the domain
c = 8

// Function to solve the problem for a given N:
function [X, Uh, h] = solve_equa(N);
    h = 1/N; // Size of a cell (uniform mesh)

    /* Construction of the matrices */
    v = ones(1, N-2) // Useful to construct Kh and Mh

    // Construction of the stiffness matrix:
    Kh0 = 2*eye(N-1, N-1) - diag(v, -1) - diag(v, 1)
    Kh0 = 1/h * Kh0

    // Construction of the Th matrix:
    Th0 = - diag(v, -1) + diag(v, 1)
    Th0 = c/2 * Th0

    // Construction of the matrix Ah of the system to be solved:
    Ah0 = Kh0 + Th0

    // Construction of the rhs Lh of the system to be solved:
    Lh0      = zeros(N-1)
    Lh0(N-1) = 1/h - c/2

    Uh0        = linsolve(Ah0, -Lh0) // Solution of the system Ah x = Lh
    Uh         = zeros(N+1)
    Uh(1, 2:N) = Uh0
    Uh(N+1)    = 1

    X = a:h:b // x-coordinate of the unknowns
endfunction

/******************************************************************************/
/* Visual comparison of the calculated solution with the exact solution.      */
/******************************************************************************/

N          = 4;            // N = Number of cells (N+1 points)
[X, Uh, h] = solve_equa(N); // Computed solution

// Plot of the solution:
//plot(X, Uh, a:0.01:b, exact_sol(a:0.01:b, c));
//legend(["$y=u_{h}(x)$", "$y=u(x)$"]);

/******************************************************************************/
/* Plot of the error in L2-norm as N grows                                    */
/******************************************************************************/

Liste_N = 10:10:500;
Err     = [];
for N = Liste_N
    [X, Uh, h] = solve_equa(N); // Computed solution
    err        = (h*sum((Uh(1:N)-exact_sol(X(1:N), c))^2))^0.5;
    Err        = [Err, err];
end

//plot(log(Liste_N), log(Err), 'k', log(Liste_N), -2*log(Liste_N), '--k');

/******************************************************************************/
/* Plot of the error in H1-seminorm as N grows                                */
/******************************************************************************/

// Function for u':
function y = up(x, c)
    y = c*exp(c*x)/(exp(c)-1)
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

/******************************************************************************/
/* Now we fix N and observe the influence of the speed c.                     */
/******************************************************************************/

N = 50;

C = 10:30:130;
col = ['b', 'r', 'g', 'm', 'k'];
i = 1;
for c = C
    [X, Uh, h] = solve_equa(N); // Computed solution
    //plot(X, Uh, col(i));
    i = i+1;
end

//legend(["b = 10", "b = 40", "b = 70", "b = 100", "b = 130"]);
//title("N = 50");
