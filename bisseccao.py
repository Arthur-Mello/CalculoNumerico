from flask import Flask, request, render_template_string
import math

app = Flask(__name__)

# ─── Funções de cálculo ───────────────────────────────────────────────────────

def f(x, expr):
    """Avalia a expressão matemática para um dado x."""
    try:
        # Permite usar funções matemáticas comuns na expressão
        env = {
            "x": x,
            "sin": math.sin, "cos": math.cos, "tan": math.tan,
            "exp": math.exp, "log": math.log, "log10": math.log10,
            "sqrt": math.sqrt, "abs": abs, "pi": math.pi, "e": math.e,
            "pow": pow,
        }
        return float(eval(expr, {"__builtins__": {}}, env))
    except Exception:
        return float("nan")


def bisseccao(a, b, tol, max_iter, expr):
    fa = f(a, expr)
    fb = f(b, expr)

    if math.isnan(fa) or math.isnan(fb):
        return {"erro": "Erro ao avaliar a função nos extremos do intervalo."}
    if fa * fb >= 0:
        return {"erro": "O intervalo não satisfaz f(a) * f(b) < 0"}

    iteracoes = []

    for i in range(1, max_iter + 1):
        p = (a + b) / 2
        fp = f(p, expr)
        erro = abs(b - a) / 2

        iteracoes.append({
            "k": i,
            "a": a,
            "b": b,
            "p": p,
            "fp": fp,
            "erro": erro,
        })

        if erro < tol or fp == 0.0:
            return {"raiz": p, "iteracoes": iteracoes}

        if fa * fp < 0:
            b = p
            fb = fp
        else:
            a = p
            fa = fp

    return {"raiz": (a + b) / 2, "iteracoes": iteracoes}


# ─── Template HTML ────────────────────────────────────────────────────────────

TEMPLATE = """
<!DOCTYPE html>
<html lang="pt-br">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Método da Bisseção</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.css" rel="stylesheet">
</head>
<body class="p-4">
  <h2 class="mb-3">Método da Bisseção</h2>

  {% if stage == 'menu' %}
  <form method="post" class="card card-body mb-4" autocomplete="off">
    <input type="hidden" name="stage" value="solve">
    <div class="mb-3">
      <label class="form-label">
        Função
        <span class="text-muted">(ex: <code>x**3 - x - 2</code> ou <code>sin(x) - x/2</code>)</span>
      </label>
      <input type="text" class="form-control" name="funcao"
             value="{{ funcao }}" required>
    </div>
    <div class="row mb-3">
      <div class="col">
        <label class="form-label">Início do intervalo (a)</label>
        <input type="number" step="any" class="form-control" name="a" value="{{ a }}" required>
      </div>
      <div class="col">
        <label class="form-label">Fim do intervalo (b)</label>
        <input type="number" step="any" class="form-control" name="b" value="{{ b }}" required>
      </div>
    </div>
    <div class="row mb-3">
      <div class="col">
        <label class="form-label">Tolerância</label>
        <input type="number" step="any" class="form-control" name="tol" value="{{ tol }}" required>
      </div>
      <div class="col">
        <label class="form-label">Máximo de iterações</label>
        <input type="number" class="form-control" name="maxIter" value="{{ maxIter }}" required>
      </div>
    </div>
    <div class="d-flex gap-2">
      <button type="submit" class="btn btn-success rounded-pill">
        <i class="bi bi-play-fill me-1"></i>Calcular
      </button>
      <a href="/" class="btn btn-secondary rounded-pill">
        <i class="bi bi-arrow-counterclockwise me-1"></i>Reiniciar
      </a>
    </div>
  </form>

  {% elif stage == 'solve' %}
    {% if erro %}
      <div class="alert alert-danger">{{ erro }}</div>
    {% else %}
      <div class="alert alert-success">
        Raiz aproximada: <strong>{{ "%.10g"|format(raiz) }}</strong>
      </div>
      <table class="table table-striped">
        <thead>
          <tr>
            <th>Iteração</th><th>a</th><th>b</th>
            <th>p</th><th>f(p)</th><th>Erro</th>
          </tr>
        </thead>
        <tbody>
          {% for row in iteracoes %}
          <tr>
            <td>{{ row.k }}</td>
            <td>{{ "%.6f"|format(row.a) }}</td>
            <td>{{ "%.6f"|format(row.b) }}</td>
            <td>{{ "%.6f"|format(row.p) }}</td>
            <td>{{ "%.6f"|format(row.fp) }}</td>
            <td>{{ "%.6f"|format(row.erro) }}</td>
          </tr>
          {% endfor %}
        </tbody>
      </table>
    {% endif %}

    <div class="d-flex gap-2 mt-3">
      {# Voltar mantendo os valores preenchidos #}
      <form method="post">
        <input type="hidden" name="stage" value="menu">
        <input type="hidden" name="funcao" value="{{ funcao }}">
        <input type="hidden" name="a" value="{{ a }}">
        <input type="hidden" name="b" value="{{ b }}">
        <input type="hidden" name="tol" value="{{ tol }}">
        <input type="hidden" name="maxIter" value="{{ maxIter }}">
        <button type="submit" class="btn btn-outline-primary">Voltar</button>
      </form>
      <a href="/" class="btn btn-outline-secondary">Reiniciar</a>
    </div>
  {% endif %}

</body>
</html>
"""

# ─── Rota Flask ───────────────────────────────────────────────────────────────

DEFAULTS = {
    "funcao": "x**3 - x - 2",
    "a": "1",
    "b": "2",
    "tol": "1e-6",
    "maxIter": "50",
}

@app.route("/", methods=["GET", "POST"])
def index():
    stage = request.form.get("stage", "menu")

    # Lê inputs (mantém defaults se não enviado)
    funcao  = request.form.get("funcao",  DEFAULTS["funcao"])
    a_str   = request.form.get("a",       DEFAULTS["a"])
    b_str   = request.form.get("b",       DEFAULTS["b"])
    tol_str = request.form.get("tol",     DEFAULTS["tol"])
    max_str = request.form.get("maxIter", DEFAULTS["maxIter"])

    ctx = dict(stage=stage, funcao=funcao,
               a=a_str, b=b_str, tol=tol_str, maxIter=max_str,
               erro=None, raiz=None, iteracoes=[])

    if stage == "solve":
        try:
            a   = float(a_str.replace(",", "."))
            b   = float(b_str.replace(",", "."))
            tol = float(tol_str.replace(",", "."))
            max_iter = int(max_str)
        except ValueError:
            ctx["erro"] = "Valores inválidos nos campos numéricos."
            return render_template_string(TEMPLATE, **ctx)

        resultado = bisseccao(a, b, tol, max_iter, funcao)

        if "erro" in resultado:
            ctx["erro"] = resultado["erro"]
        else:
            ctx["raiz"]     = resultado["raiz"]
            ctx["iteracoes"] = resultado["iteracoes"]

    return render_template_string(TEMPLATE, **ctx)


if __name__ == "__main__":
    app.run(debug=True)