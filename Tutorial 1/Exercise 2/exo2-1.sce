/******************************************************************************/
/* Scilab script to solve equation from exo2 via finite element method.       */
/* More precisely, the problem to solve is -u"+u=f, u(0)=u(1), u'(0)=u'(1).   */
/* Here, we take f≡1. The exact solution is 1.                                */
/******************************************************************************/

/* Parameters of the problem */
a = 0; // Left bound of the domain
b = 1; // Right bound of the domain

// Function to solve the problem for a given N:
function [X, Uh, h] = solve_equa(N);
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

    // Construction of the matrix Ah of the system to be solve_equad:
    Ah = Kh;

    // Construction of the matrix Ah of the system to be solved:
    Ah = Kh + Mh;

    Lh = h * ones(N, 1); // Rhs Lh of the system to be solved

    X       = a:h:b;             // x-coordinate of the unknowns
    Uh      = linsolve(Ah, -Lh); // Solution of the system Ah x = Lh
    Uh(N+1) = Uh(1);
    Uh = Uh'
endfunction

/******************************************************************************/
/* Visual comparison of the calculated solution with the exact solution.      */
/******************************************************************************/

N          = 10;            // N = Number of cells (N+1 points)
[X, Uh, h] = solve_equa(N); // Computed solution

// Plot of the solution:
//plot(X, Uh, X, ones(X));
//legend(["$y=u_{h}(x)$", "$y=u(x)$"]);
//a=gca();
//a.data_bounds=[0, 0; 1, 2];
