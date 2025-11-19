# app.py
from flask import Flask, render_template, request, redirect, url_for, flash
import math, ast, operator
from typing import List, Tuple, Callable, Any

app = Flask(__name__)
app.secret_key = "replace-this-with-a-secret-in-prod"

# ---------- Safe expression evaluator for g(x) ----------
# We convert a user expression like "math.cos(x)" or "x**2 + 1" to a callable safely.
SAFE_MATH = {k: getattr(math, k) for k in dir(math) if not k.startswith("__")}
SAFE_GLOBALS = {"__builtins__": None}
SAFE_GLOBALS.update(SAFE_MATH)

def make_function_from_expr(expr: str) -> Callable[[float], float]:
    """
    Make a function g(x) from a user expression string.
    Uses eval with a restricted globals dictionary (math functions only).
    Allowed variable: x
    """
    # Basic safety: only allow certain AST node types
    allowed_nodes = (
        ast.Expression, ast.BinOp, ast.UnaryOp, ast.Call,
        ast.Name, ast.Load, ast.Constant, ast.Pow, ast.Mult,
        ast.Add, ast.Sub, ast.Div, ast.Mod, ast.USub, ast.UAdd,
        ast.Tuple, ast.List
    )
    parsed = ast.parse(expr, mode="eval")
    for node in ast.walk(parsed):
        if not isinstance(node, allowed_nodes):
            raise ValueError(f"Disallowed expression element: {type(node).__name__}")
    def g(x):
        local = {"x": x}
        return eval(compile(parsed, filename="<expr>", mode="eval"), SAFE_GLOBALS, local)
    # Test it quickly
    _ = g(1.0)
    return g

# ---------- Polynomial f(x) and derivative ----------
def parse_coeffs(text: str) -> List[float]:
    # Accept space- or comma-separated floats
    tokens = text.replace(",", " ").split()
    return [float(t) for t in tokens]

def f_poly(x: float, coeffs: List[float]) -> float:
    res = 0.0
    power = len(coeffs) - 1
    for c in coeffs:
        res += c * (x ** power)
        power -= 1
    return res

def derivative_poly(x: float, coeffs: List[float]) -> float:
    res = 0.0
    power = len(coeffs) - 1
    for c in coeffs[:-1]:
        res += c * power * (x ** (power - 1))
        power -= 1
    return res

# ---------- Root-finding methods (return root, final_error, iters, iter_rows) ----------
def bisection(a: float, b: float, tol: float, max_iter: int, coeffs: List[float]):
    if f_poly(a, coeffs) * f_poly(b, coeffs) >= 0:
        return None, None, 0, [("error", "f(a) and f(b) must have opposite signs")]
    rows = []
    for i in range(1, max_iter+1):
        c = (a + b) / 2.0
        fc = f_poly(c, coeffs)
        error = abs(b - a) / 2.0
        rows.append((i, a, b, c, fc, error))
        if abs(fc) < tol or error < tol:
            return c, error, i, rows
        if f_poly(a, coeffs) * fc < 0:
            b = c
        else:
            a = c
    return c, error, max_iter, rows

def regula_falsi(a: float, b: float, tol: float, max_iter: int, coeffs: List[float]):
    rows = []
    for i in range(1, max_iter+1):
        fa = f_poly(a, coeffs)
        fb = f_poly(b, coeffs)
        if fb - fa == 0:
            return None, None, i-1, [("error", "division by zero in formula")]
        x = b - fb * (b - a) / (fb - fa)
        fx = f_poly(x, coeffs)
        error = abs(fx)
        rows.append((i, a, b, x, fx, error))
        if abs(fx) < tol:
            return x, error, i, rows
        if fa * fx < 0:
            b = x
        else:
            a = x
    return x, error, max_iter, rows

def secant(x0: float, x1: float, tol: float, max_iter: int, coeffs: List[float]):
    rows = []
    for i in range(1, max_iter+1):
        f0 = f_poly(x0, coeffs)
        f1 = f_poly(x1, coeffs)
        denom = (f1 - f0)
        if denom == 0:
            return None, None, i-1, [("error", "division by zero in secant")]
        x2 = x1 - f1 * (x1 - x0) / denom
        error = abs(x2 - x1)
        rows.append((i, x0, x1, x2, f_poly(x2, coeffs), error))
        if error < tol:
            return x2, error, i, rows
        x0, x1 = x1, x2
    return x2, error, max_iter, rows

def newton_raphson(x0: float, tol: float, max_iter: int, coeffs: List[float]):
    rows = []
    for i in range(1, max_iter+1):
        fx = f_poly(x0, coeffs)
        dfx = derivative_poly(x0, coeffs)
        if dfx == 0:
            return None, None, i-1, [("error", "derivative zero")]
        x1 = x0 - fx / dfx
        error = abs(x1 - x0)
        rows.append((i, x0, fx, dfx, x1, error))
        if error < tol:
            return x1, error, i, rows
        x0 = x1
    return x1, error, max_iter, rows

def fixed_point(g_func: Callable[[float], float], x0: float, tol: float, max_iter: int):
    rows = []
    for i in range(1, max_iter+1):
        x1 = g_func(x0)
        error = abs(x1 - x0)
        rows.append((i, x0, x1, error))
        if error < tol:
            return x1, error, i, rows
        x0 = x1
    return x1, error, max_iter, rows

def modified_secant(x0: float, delta: float, tol: float, max_iter: int, coeffs: List[float]):
    rows = []
    for i in range(1, max_iter+1):
        fx = f_poly(x0, coeffs)
        fx_delta = f_poly(x0 + delta * x0, coeffs)
        denom = (fx_delta - fx)
        if denom == 0:
            return None, None, i-1, [("error", "division by zero in modified secant")]
        x1 = x0 - fx * (delta * x0) / denom
        error = abs(x1 - x0)
        rows.append((i, x0, fx, x1, error))
        if error < tol:
            return x1, error, i, rows
        x0 = x1
    return x1, error, max_iter, rows

# ---------- Flask routes ----------
@app.route("/", methods=["GET"])
def index():
    return render_template("index.html", result=None)

@app.route("/solve", methods=["POST"])
def solve():
    try:
        method = int(request.form.get("method"))
        coeffs_text = request.form.get("coeffs", "").strip()
        coeffs = parse_coeffs(coeffs_text) if coeffs_text else []
        tol = float(request.form.get("tol", "1e-6"))
        max_iter = int(request.form.get("max_iter", "50"))

        # method-specific inputs
        a = request.form.get("a", "")
        b = request.form.get("b", "")
        x0 = request.form.get("x0", "")
        x1 = request.form.get("x1", "")
        delta = request.form.get("delta", "")
        g_expr = request.form.get("g_expr", "")

        # run chosen method
        if method == 1:
            a_f = float(a); b_f = float(b)
            root, final_error, iters, rows = bisection(a_f, b_f, tol, max_iter, coeffs)
        elif method == 2:
            a_f = float(a); b_f = float(b)
            root, final_error, iters, rows = regula_falsi(a_f, b_f, tol, max_iter, coeffs)
        elif method == 3:
            x0_f = float(x0); x1_f = float(x1)
            root, final_error, iters, rows = secant(x0_f, x1_f, tol, max_iter, coeffs)
        elif method == 4:
            x0_f = float(x0)
            root, final_error, iters, rows = newton_raphson(x0_f, tol, max_iter, coeffs)
        elif method == 5:
            # build g(x)
            if not g_expr:
                raise ValueError("g(x) expression required for Fixed Point Method")
            g_func = make_function_from_expr(g_expr)
            x0_f = float(x0)
            root, final_error, iters, rows = fixed_point(g_func, x0_f, tol, max_iter)
        elif method == 6:
            x0_f = float(x0)
            delta_f = float(delta)
            root, final_error, iters, rows = modified_secant(x0_f, delta_f, tol, max_iter, coeffs)
        else:
            raise ValueError("Invalid method selected")

        # detect early error returns (format: rows starting with ("error", ...))
        if rows and isinstance(rows[0], tuple) and rows[0][0] == "error":
            flash(f"Method error: {rows[0][1]}", "danger")
            return render_template("index.html", result=None)

        result = {
            "method": method,
            "root": root,
            "final_error": final_error,
            "iterations": iters,
            "rows": rows,
            "coeffs_text": coeffs_text,
            "g_expr": g_expr
        }
        return render_template("index.html", result=result)
    except Exception as e:
        flash(f"Error: {str(e)}", "danger")
        return redirect(url_for("index"))

if __name__ == "__main__":
    # for local debugging
    app.run(debug=True, host="0.0.0.0", port=5000)
