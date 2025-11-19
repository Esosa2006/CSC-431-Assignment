import math

# ============================
# Utility: Evaluate function
# ============================
def f(x, coeffs):
    """
    Evaluate polynomial: a_n x^n + ... + a1 x + a0
    coeffs = [a_n, ..., a0]
    """
    result = 0
    power = len(coeffs) - 1
    for c in coeffs:
        result += c * (x ** power)
        power -= 1
    return result


def derivative(x, coeffs):
    """Compute derivative of polynomial."""
    result = 0
    power = len(coeffs) - 1
    for c in coeffs[:-1]:
        result += c * power * (x ** (power - 1))
        power -= 1
    return result


# ==================================
# 1. BISECTION METHOD
# ==================================
def bisection(a, b, tol, max_iter, coeffs):
    if f(a, coeffs) * f(b, coeffs) >= 0:
        print("Error: f(a) and f(b) must have opposite signs!")
        return None, None, None

    print("\nITER |     a       b       c       f(c)       Error")
    print("----------------------------------------------------------")

    for i in range(1, max_iter + 1):
        c = (a + b) / 2
        fc = f(c, coeffs)

        error = abs(b - a)

        print(f"{i:>4} | {a:>7.5f} {b:>7.5f} {c:>7.5f} {fc:>10.6f} {error:>10.6f}")

        if fc == 0 or error < tol:
            return c, error, i

        if f(a, coeffs) * fc < 0:
            b = c
        else:
            a = c

    return c, error, max_iter


# ==================================
# 2. REGULA FALSI
# ==================================
def regula_falsi(a, b, tol, max_iter, coeffs):
    print("\nITER |     a       b      x       f(x)       Error")
    print("----------------------------------------------------------")

    for i in range(1, max_iter + 1):
        fa = f(a, coeffs)
        fb = f(b, coeffs)

        x = b - fb * (b - a) / (fb - fa)
        fx = f(x, coeffs)

        error = abs(fx)

        print(f"{i:>4} | {a:>7.5f} {b:>7.5f} {x:>7.5f} {fx:>10.6f} {error:>10.6f}")

        if abs(fx) < tol:
            return x, error, i

        if fa * fx < 0:
            b = x
        else:
            a = x

    return x, error, max_iter


# ==================================
# 3. SECANT METHOD
# ==================================
def secant(x0, x1, tol, max_iter, coeffs):
    print("\nITER |    x0        x1        x2        f(x2)       Error")
    print("-----------------------------------------------------------------")

    for i in range(1, max_iter + 1):
        f0 = f(x0, coeffs)
        f1 = f(x1, coeffs)

        if f1 - f0 == 0:
            print("Error: division by zero")
            return None

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        error = abs(x2 - x1)

        print(f"{i:>4} | {x0:>7.5f} {x1:>7.5f} {x2:>7.5f} {f(x2, coeffs):>10.6f} {error:>10.6f}")

        if error < tol:
            return x2, error, i

        x0, x1 = x1, x2

    return x2, error, max_iter


# ==================================
# 4. NEWTON–RAPHSON
# ==================================
def newton_raphson(x0, tol, max_iter, coeffs):
    print("\nITER |     x        f(x)      Error")
    print("------------------------------------------")

    for i in range(1, max_iter + 1):
        fx = f(x0, coeffs)
        dfx = derivative(x0, coeffs)

        if dfx == 0:
            print("Derivative became zero. Stopping.")
            return None

        x1 = x0 - fx / dfx
        error = abs(x1 - x0)

        print(f"{i:>4} | {x0:>7.5f} {fx:>10.6f} {error:>10.6f}")

        if error < tol:
            return x1, error, i

        x0 = x1

    return x1, error, max_iter


# ==================================
# 5. FIXED POINT ITERATION
# Requires user-defined g(x)
# ==================================
def fixed_point(g, x0, tol, max_iter):
    print("\nITER |     x        g(x)      Error")
    print("------------------------------------------")

    for i in range(1, max_iter + 1):
        x1 = g(x0)
        error = abs(x1 - x0)

        print(f"{i:>4} | {x0:>7.5f} {x1:>10.6f} {error:>10.6f}")

        if error < tol:
            return x1, error, i

        x0 = x1

    return x1, error, max_iter


# ==================================
# 6. MODIFIED SECANT METHOD
# x_new = x - f(x)*δx / (f(x+δx)-f(x))
# ==================================
def modified_secant(x0, delta, tol, max_iter, coeffs):
    print("\nITER |     x       f(x)       Error")
    print("------------------------------------------")

    for i in range(1, max_iter + 1):
        fx = f(x0, coeffs)
        fx_delta = f(x0 + delta * x0, coeffs)

        x1 = x0 - fx * (delta * x0) / (fx_delta - fx)
        error = abs(x1 - x0)

        print(f"{i:>4} | {x0:>7.5f} {fx:>10.6f} {error:>10.6f}")

        if error < tol:
            return x1, error, i

        x0 = x1

    return x1, error, max_iter


# ==================================
# MAIN CLI
# ==================================
def main():
    print("\n===== ROOT FINDING CALCULATOR (ZOF CLI) =====")
    print("Choose a method:")
    print("1. Bisection")
    print("2. Regula Falsi")
    print("3. Secant")
    print("4. Newton–Raphson")
    print("5. Fixed Point Iteration")
    print("6. Modified Secant\n")

    choice = int(input("Enter method number: "))

    # Polynomial coefficients
    coeffs = list(map(float, input("\nEnter polynomial coefficients (highest degree first): ").split()))

    tol = float(input("Enter tolerance: "))
    max_iter = int(input("Enter max iterations: "))

    if choice == 1:
        a = float(input("Enter a: "))
        b = float(input("Enter b: "))
        root, error, it = bisection(a, b, tol, max_iter, coeffs)
        
        if root is None:
            return

    elif choice == 2:
        a = float(input("Enter a: "))
        b = float(input("Enter b: "))
        root, error, it = regula_falsi(a, b, tol, max_iter, coeffs)

    elif choice == 3:
        x0 = float(input("Enter x0: "))
        x1 = float(input("Enter x1: "))
        root, error, it = secant(x0, x1, tol, max_iter, coeffs)

    elif choice == 4:
        x0 = float(input("Enter initial guess x0: "))
        root, error, it = newton_raphson(x0, tol, max_iter, coeffs)

    elif choice == 5:
        print("\nEnter g(x) as a Python function. Example: lambda x: math.cos(x)")
        g_input = input("g(x) = ")
        g = eval(g_input)
        x0 = float(input("Enter initial guess x0: "))
        root, error, it = fixed_point(g, x0, tol, max_iter)

    elif choice == 6:
        x0 = float(input("Enter initial guess x0: "))
        delta = float(input("Enter delta value (e.g., 0.01): "))
        root, error, it = modified_secant(x0, delta, tol, max_iter, coeffs)

    else:
        print("Invalid choice.")
        return

    print("\n===== FINAL RESULT =====")
    print(f"Estimated root: {root}")
    print(f"Final error: {error}")
    print(f"Iterations: {it}")


# Run program
if __name__ == "__main__":
    main()
