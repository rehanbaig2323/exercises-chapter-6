"""A module providing numerical solvers for nonlinear equations."""


class ConvergenceError(Exception):
    """Exception raised if a solver fails to converge."""

    pass


def newton_raphson(f, df, x_0, eps=1.0e-5, max_its=20):
    """Solve a nonlinear equation using Newton-Raphson iteration.

    Solve f==0 using Newton-Raphson iteration.

    Parameters
    ----------
    f : function(x: float) -> float
        The function whose root is being found.
    df : function(x: float) -> float
        The derivative of f.
    x_0 : float
        The initial value of x in the iteration.
    eps : float
        The solver tolerance. Convergence is achieved when abs(f(x)) < eps.
    max_its : int
        The maximum number of iterations to be taken before the solver is taken
        to have failed.

    Returns
    -------
    float
        The approximate root computed using Newton iteration.
    """
    count = 0
    while count <= max_its:
        x = x_0 - (f(x_0) / df(x_0))
        if abs(f(x)) < eps:
            return x
        else:
            count += 1
            x_0 = x
            continue
    raise ConvergenceError("Too many iterations, hasn't converged")


def bisection(f, x_0, x_1, eps=1.0e-5, max_its=20):
    """Solve a nonlinear equation using bisection.

    Solve f==0 using bisection starting with the interval [x_0, x_1]. f(x_0)
    and f(x_1) must differ in sign.

    Parameters
    ----------
    f : function(x: float) -> float
        The function whose root is being found.
    x_0 : float
        The left end of the initial bisection interval.
    x_1 : float
        The right end of the initial bisection interval.
    eps : float
        The solver tolerance. Convergence is achieved when abs(f(x)) < eps.
    max_its : int
        The maximum number of iterations to be taken before the solver is taken
        to have failed.

    Returns
    -------
    float
        The approximate root computed using bisection.
    """
    count = 0
    while count <= max_its:
        if f(x_0) * f(x_1) > 0:
            raise ValueError("Points not on either side of x-axis.")
        x = (x_0 + x_1) / 2
        if abs(f(x)) < eps:
            return x
        else:
            if f(x_0) * f(x) > 0:
                x_0 = x
            else:
                x_1 = x
            count += 1
            continue
    raise ConvergenceError("Too many iterations.")


def solve(f, df, x_0, x_1, eps=1.0e-5, max_its_n=20, max_its_b=20):
    """Solve a nonlinear equation.

    solve f(x) == 0 using Newton-Raphson iteration, falling back to bisection
    if the former fails.

    Parameters
    ----------
    f : function(x: float) -> float
        The function whose root is being found.
    df : function(x: float) -> float
        The derivative of f.
    x_0 : float
        The initial value of x in the Newton-Raphson iteration, and left end of
        the initial bisection interval.
    x_1 : float
        The right end of the initial bisection interval.
    eps : float
        The solver tolerance. Convergence is achieved when abs(f(x)) < eps.
    max_its_n : int
        The maximum number of iterations to be taken before the newton-raphson
        solver is taken to have failed.
    max_its_b : int
        The maximum number of iterations to be taken before the bisection
        solver is taken to have failed.

    Returns
    -------
    float
        The approximate root.
    """
    # Delete these two lines when implementing the method.

# hahashgiajhg
