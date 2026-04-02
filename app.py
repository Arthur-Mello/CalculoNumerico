from flask import Flask, request, render_template_string
import math

app = Flask(__name__)

# ════════════════════════════════════════════════════════════════
#  FUNÇÕES DE CÁLCULO
# ════════════════════════════════════════════════════════════════

# ── Bisseção ─────────────────────────────────────────────────────

def f_eval(x, expr):
    try:
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
    fa = f_eval(a, expr)
    fb = f_eval(b, expr)
    if math.isnan(fa) or math.isnan(fb):
        return {"erro": "Erro ao avaliar a função nos extremos do intervalo."}
    if fa * fb >= 0:
        return {"erro": "O intervalo não satisfaz f(a) * f(b) < 0"}
    iteracoes = []
    for i in range(1, max_iter + 1):
        p  = (a + b) / 2
        fp = f_eval(p, expr)
        erro = abs(b - a) / 2
        iteracoes.append({"k": i, "a": a, "b": b, "p": p, "fp": fp, "erro": erro})
        if erro < tol or fp == 0.0:
            return {"raiz": p, "iteracoes": iteracoes}
        if fa * fp < 0:
            b = p; fb = fp
        else:
            a = p; fa = fp
    return {"raiz": (a + b) / 2, "iteracoes": iteracoes}


# ── Eliminação de Gauss ──────────────────────────────────────────

def is_near_zero(v, eps=1e-10):
    return abs(v) < eps


def render_matrix_state(A, b):
    n = len(A)
    html = '<p><strong>Estado da Matriz</strong></p>'
    html += '<div class="table-responsive"><table class="table table-dark table-sm align-middle mt-1 mb-3"><thead><tr>'
    for j in range(n):
        html += f'<th class="text-white">U[·,{j+1}]</th>'
    html += '<th class="text-white">= c</th></tr></thead><tbody>'
    for i in range(n):
        html += '<tr>'
        for j in range(n):
            html += f'<td class="text-white">{A[i][j]:.6f}</td>'
        html += f'<td class="text-white">{b[i]:.6f}</td></tr>'
    html += '</tbody></table></div>'
    return html


def gauss_elimination(A_orig, b_orig, use_P=True):
    n = len(A_orig)
    A = [[float(A_orig[i][j]) for j in range(n)] for i in range(n)]
    b = [float(v) for v in b_orig]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    P = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)] if use_P else None
    steps = []
    for k in range(n):
        p = k
        max_val = abs(A[k][k])
        if use_P:
            for i in range(k + 1, n):
                if abs(A[i][k]) > max_val:
                    max_val = abs(A[i][k]); p = i
            if p != k:
                A[k], A[p] = A[p], A[k]
                b[k], b[p] = b[p], b[k]
                P[k], P[p] = P[p], P[k]
                for col in range(k):
                    L[k][col], L[p][col] = L[p][col], L[k][col]
                steps.append({'type': 'text', 'content': f'Troca L{k+1} ⇄ L{p+1} (pivotamento parcial)'})
        if is_near_zero(A[k][k]):
            steps.append({'type': 'text', 'content': f'Coluna {k}: pivô ~ 0 ⇒ sistema singular.'})
            return None, L, A, P, b, steps
        for i in range(k + 1, n):
            if is_near_zero(A[i][k]):
                continue
            m = A[i][k] / A[k][k]
            L[i][k] = m
            for j in range(k, n):
                A[i][j] -= m * A[k][j]
            b[i] -= m * b[k]
            steps.append({'type': 'text', 'content': f'L{i+1} ← L{i+1} − ({m:.6g})·L{k+1}'})
            steps.append({'type': 'matrix', 'content': render_matrix_state(A, b)})
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i + 1, n))
        if is_near_zero(A[i][i]):
            steps.append({'type': 'text', 'content': 'Pivô zero na retrossubstituição ⇒ singular.'})
            return None, L, A, P, b, steps
        x[i] = (b[i] - s) / A[i][i]
    return x, L, A, P, b, steps


# ── Decomposição LU ──────────────────────────────────────────────

def lu_decomposition(A_orig, b_orig):
    """Decompõe A = L·U (sem pivotamento) e resolve L·U·x = b."""
    n = len(A_orig)
    A = [[float(A_orig[i][j]) for j in range(n)] for i in range(n)]
    b = [float(v) for v in b_orig]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[0.0] * n for _ in range(n)]
    steps = []

    # Fatoração LU (algoritmo de Doolittle)
    for k in range(n):
        # Linha de U
        for j in range(k, n):
            s = sum(L[k][p] * U[p][j] for p in range(k))
            U[k][j] = A[k][j] - s

        if is_near_zero(U[k][k]):
            steps.append(f'Pivô nulo em U[{k+1},{k+1}] ⇒ sistema singular.')
            return None, L, U, steps

        # Coluna de L
        for i in range(k + 1, n):
            s = sum(L[i][p] * U[p][k] for p in range(k))
            L[i][k] = (A[i][k] - s) / U[k][k]

        steps.append(f'Passo {k+1}: coluna {k+1} de L e linha {k+1} de U calculadas.')

    # Substituição direta: L·y = b
    y = [0.0] * n
    for i in range(n):
        s = sum(L[i][j] * y[j] for j in range(i))
        y[i] = (b[i] - s) / L[i][i]

    # Retrossubstituição: U·x = y
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(U[i][j] * x[j] for j in range(i + 1, n))
        if is_near_zero(U[i][i]):
            steps.append('Pivô zero na retrossubstituição ⇒ singular.')
            return None, L, U, steps
        x[i] = (y[i] - s) / U[i][i]

    steps.append('Substituição direta (L·y = b) e retrossubstituição (U·x = y) concluídas.')
    return x, L, U, steps


# ── Jacobi ───────────────────────────────────────────────────────

def criterio_das_linhas(n, A):
    for i in range(n):
        soma = sum(abs(A[i][j]) for j in range(n) if j != i)
        if A[i][i] == 0 or soma >= abs(A[i][i]):
            return False
    return True


def jacobi(n, A, b, x0, tol, max_iter=1000):
    convergencia = criterio_das_linhas(n, A)
    x     = [v for v in x0]
    x_novo = [0.0] * n
    iteracoes = []
    for k in range(max_iter):
        for i in range(n):
            soma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_novo[i] = (b[i] - soma) / A[i][i]
        diferencas = [abs(x_novo[i] - x[i]) for i in range(n)]
        max_dif = max(diferencas)
        max_x   = max(abs(v) for v in x_novo)
        erro    = 0.0 if max_x == 0 else max_dif / max_x
        iteracoes.append({"k": k+1, "x": [round(v, 8) for v in x_novo], "erro": erro})
        x = [v for v in x_novo]
        if erro < tol:
            return {"x": x, "iteracoes": iteracoes, "convergencia": convergencia, "convergiu": True}
    return {"x": x, "iteracoes": iteracoes, "convergencia": convergencia, "convergiu": False}


# ── Gauss-Seidel ─────────────────────────────────────────────────

def criterio_de_sassenfeld(n, A):
    beta = [0.0] * n
    for i in range(n):
        soma = sum(abs(A[i][j]) * (beta[j] if j < i else 1.0) for j in range(n) if j != i)
        if A[i][i] == 0:
            return False
        beta[i] = soma / abs(A[i][i])
    return max(beta) < 1


def gauss_seidel(n, A, b, x0, tol, max_iter=1000):
    convergencia = criterio_de_sassenfeld(n, A)
    x = [v for v in x0]
    iteracoes = []
    for k in range(max_iter):
        x_anterior = [v for v in x]
        for i in range(n):
            soma1 = sum(A[i][j] * x[j]          for j in range(i))
            soma2 = sum(A[i][j] * x_anterior[j] for j in range(i + 1, n))
            x[i]  = (b[i] - soma1 - soma2) / A[i][i]
        diferencas = [abs(x[i] - x_anterior[i]) for i in range(n)]
        max_dif = max(diferencas)
        max_x   = max(abs(v) for v in x)
        erro    = 0.0 if max_x == 0 else max_dif / max_x
        iteracoes.append({"k": k+1, "x": [round(v, 8) for v in x], "erro": erro})
        if erro < tol:
            return {"x": x, "iteracoes": iteracoes, "convergencia": convergencia, "convergiu": True}
    return {"x": x, "iteracoes": iteracoes, "convergencia": convergencia, "convergiu": False}


# ════════════════════════════════════════════════════════════════
#  TEMPLATES HTML
# ════════════════════════════════════════════════════════════════

BASE_HEAD = """
<!doctype html>
<html lang="pt-br">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Cálculo Numérico</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.css" rel="stylesheet">
  <style>
    body { background: #f4f6fb; }
    .navbar-brand { font-weight: 700; letter-spacing: 1px; }
    h1,h2,h3,h4,h5,h6 { color: #1a2340; }
    .card { border: none; border-radius: 1rem; box-shadow: 0 2px 16px rgba(0,0,0,.08); }
    .badge-conv { font-size: .85rem; }
  </style>
</head>


"""

BASE_FOOT = """
</div>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
</body></html>
"""

# ── Menu Principal ───────────────────────────────────────────────

MENU_TEMPLATE = BASE_HEAD + """
<div class="row justify-content-center">
  <div class="col-md-8">
    <h1 class="mb-1">Menu Principal</h1>
    <p class="text-muted mb-4">Selecione o método numérico desejado:</p>
    <div class="row g-3">

      <div class="col-md-6">
        <a href="/bisseccao" class="text-decoration-none">
          <div class="card p-4 h-100 text-center">
            <i class=""></i>
            <h5>Método da Bisseção</h5>
            <p class="text-muted small mb-0"></p>
          </div>
        </a>
      </div>

      <div class="col-md-6">
        <a href="/gauss" class="text-decoration-none">
          <div class="card p-4 h-100 text-center">
            <i class=""></i>
            <h5>Eliminação de Gauss</h5>
            <p class="text-muted small mb-0"></p>
          </div>
        </a>
      </div>

      <div class="col-md-6">
        <a href="/lu" class="text-decoration-none">
          <div class="card p-4 h-100 text-center">
            <i class=""></i>
            <h5>Decomposição LU</h5>
            <p class="text-muted small mb-0"></p>
          </div>
        </a>
      </div>

      <div class="col-md-6">
        <a href="/jacobi" class="text-decoration-none">
          <div class="card p-4 h-100 text-center">
            <i class=""></i>
            <h5>Método de Jacobi</h5>
            <p class="text-muted small mb-0"></p>
          </div>
        </a>
      </div>

      <div class="col-md-6 offset-md-3">
        <a href="/seidel" class="text-decoration-none">
          <div class="card p-4 h-100 text-center">
            <i class=""></i>
            <h5>Gauss-Seidel</h5>
            <p class="text-muted small mb-0"></p>
          </div>
        </a>
      </div>

    </div>
  </div>
</div>
""" + BASE_FOOT

# ── Bisseção ─────────────────────────────────────────────────────

BISSECCAO_TEMPLATE = BASE_HEAD + """
<h1 class="mb-1"><i class=""></i>Método da Bisseção</h1>
<p class="text-muted mb-4">Zeros de funções por divisão de intervalos</p>

{% if stage == 'menu' %}
<div class="card p-4" style="max-width:600px">
  <form method="post" autocomplete="off">
    <input type="hidden" name="stage" value="solve">
    <div class="mb-3">
      <label class="form-label fw-semibold">Função f(x)</label>
      <input type="text" class="form-control" name="funcao" value="{{ funcao }}" required>
      <div class="form-text">Ex: <code>x**3 - x - 2</code> &nbsp;|&nbsp; <code>sin(x) - x/2</code></div>
    </div>
    <div class="row g-3 mb-3">
      <div class="col">
        <label class="form-label">a (início)</label>
        <input type="number" step="any" class="form-control" name="a" value="{{ a }}" required>
      </div>
      <div class="col">
        <label class="form-label">b (fim)</label>
        <input type="number" step="any" class="form-control" name="b" value="{{ b }}" required>
      </div>
    </div>
    <div class="row g-3 mb-4">
      <div class="col">
        <label class="form-label">Tolerância</label>
        <input type="number" step="any" class="form-control" name="tol" value="{{ tol }}" required>
      </div>
      <div class="col">
        <label class="form-label">Máx. de iterações</label>
        <input type="number" class="form-control" name="maxIter" value="{{ maxIter }}" required>
      </div>
    </div>
    <div class="d-flex gap-2">
      <button class="btn btn-success rounded-pill px-4"><i class="bi bi-play-fill me-1"></i>Calcular</button>
      <a href="/" class="btn btn-outline-secondary rounded-pill">Menu</a>
    </div>
  </form>
</div>

{% elif stage == 'solve' %}
  {% if erro %}
    <div class="alert alert-danger">{{ erro }}</div>
  {% else %}
    <div class="alert alert-success fs-5">
      Raiz aproximada: <strong>{{ "%.10g"|format(raiz) }}</strong>
    </div>
    <div class="card p-3">
      <div class="table-responsive">
        <table class="table table-striped table-hover align-middle mb-0">
          <thead class="table-dark">
            <tr><th>k</th><th>a</th><th>b</th><th>p</th><th>f(p)</th><th>Erro</th></tr>
          </thead>
          <tbody>
            {% for r in iteracoes %}
            <tr>
              <td>{{ r.k }}</td>
              <td>{{ "%.8f"|format(r.a) }}</td>
              <td>{{ "%.8f"|format(r.b) }}</td>
              <td>{{ "%.8f"|format(r.p) }}</td>
              <td>{{ "%.8f"|format(r.fp) }}</td>
              <td>{{ "%.2e"|format(r.erro) }}</td>
            </tr>
            {% endfor %}
          </tbody>
        </table>
      </div>
    </div>
  {% endif %}
  <div class="d-flex gap-2 mt-3">
    <form method="post">
      <input type="hidden" name="stage" value="menu">
      <input type="hidden" name="funcao" value="{{ funcao }}">
      <input type="hidden" name="a" value="{{ a }}">
      <input type="hidden" name="b" value="{{ b }}">
      <input type="hidden" name="tol" value="{{ tol }}">
      <input type="hidden" name="maxIter" value="{{ maxIter }}">
      <button class="btn btn-outline-primary rounded-pill">Voltar</button>
    </form>
    <a href="/bisseccao" class="btn btn-outline-secondary rounded-pill">Reiniciar</a>
    <a href="/" class="btn btn-outline-dark rounded-pill">Menu</a>
  </div>
{% endif %}
""" + BASE_FOOT

# ── Gauss ────────────────────────────────────────────────────────

GAUSS_TEMPLATE = BASE_HEAD + """
<h1 class="mb-1"></i>Eliminação de Gauss</h1>
<p class="text-muted mb-4">Método direto {{ '(pivotamento parcial)' if use_P else '(simples)' }}</p>

{% if stage == 'choose_n' %}
<div class="card p-4" style="max-width:500px">
  <form method="post" class="row g-3">
    <input type="hidden" name="stage" value="enter_matrix">
    <div class="col-auto">
      <label class="form-label fw-semibold">Tamanho n × n</label>
      <select name="n" class="form-select rounded-pill">
        {% for i in range(2,11) %}
          <option value="{{ i }}" {{ 'selected' if i==n }}>{{ i }}</option>
        {% endfor %}
      </select>
    </div>
    <div class="col-12">
      <div class="form-check">
        <input class="form-check-input" type="checkbox" name="use_P" value="1" id="use_P" {{ 'checked' if use_P }}>
        <label class="form-check-label" for="use_P">Pivotamento parcial (matriz P)</label>
      </div>
    </div>
    <div class="col-12">
      <button class="btn btn-info rounded-pill"><i class="bi bi-grid me-1"></i>Criar matriz</button>
      <a href="/" class="btn btn-outline-secondary rounded-pill ms-2">Menu</a>
    </div>
  </form>
</div>

{% elif stage == 'enter_matrix' %}
<div class="card p-4">
  <form method="post">
    <input type="hidden" name="stage" value="solve">
    <input type="hidden" name="n" value="{{ n }}">
    <input type="hidden" name="use_P" value="{{ 1 if use_P else 0 }}">
    <p class="mb-2">Informe <strong>A</strong> e <strong>b</strong> (sistema A·x = b):</p>
    <div class="table-responsive mb-3">
      <table class="table table-dark table-striped align-middle">
        <thead><tr>
          {% for j in range(n) %}<th>A[·,{{ j+1 }}]</th>{% endfor %}
          <th>= b</th>
        </tr></thead>
        <tbody>
          {% for i in range(n) %}<tr>
            {% for j in range(n) %}
              <td><input required type="text" class="form-control form-control-sm"
                   name="a_{{ i }}_{{ j }}" value="{{ '1' if i==j else '0' }}"></td>
            {% endfor %}
            <td><input required type="text" class="form-control form-control-sm"
                 name="b_{{ i }}" value="0"></td>
          </tr>{% endfor %}
        </tbody>
      </table>
    </div>
    <div class="d-flex gap-2">
      <button class="btn btn-success rounded-pill"><i class="bi bi-play-fill me-1"></i>Resolver</button>
      <a href="/gauss" class="btn btn-outline-secondary rounded-pill">Reiniciar</a>
      <a href="/" class="btn btn-outline-dark rounded-pill">Menu</a>
    </div>
  </form>
</div>

{% elif stage == 'solve' %}
<div class="row g-4">
  <div class="col-lg-6">
    <div class="card p-3 mb-3">
      <h5><i class="bi bi-list-check me-2"></i>Passos da eliminação</h5>
      {% if steps %}
      <ol>
        {% for s in steps %}
          {% if s.type == 'matrix' %}<div>{{ s.content|safe }}</div>
          {% else %}<li>{{ s.content }}</li>{% endif %}
        {% endfor %}
      </ol>
      {% else %}<p class="text-muted">Nenhuma operação registrada.</p>{% endif %}
    </div>

    {% if use_P and P %}
    <div class="card p-3 mb-3">
      <h5><i class="bi bi-clipboard-data me-2"></i>Matriz de Permutação P</h5>
      <div class="table-responsive">
        <table class="table table-dark table-sm">
          <thead><tr>{% for j in range(n) %}<th>P[·,{{ j+1 }}]</th>{% endfor %}</tr></thead>
          <tbody>{% for i in range(n) %}<tr>
            {% for j in range(n) %}<td>{{ "%.6f"|format(P[i][j]) }}</td>{% endfor %}
          </tr>{% endfor %}</tbody>
        </table>
      </div>
    </div>
    {% endif %}

    <div class="card p-3 mb-3">
      <h5><i class="bi bi-clipboard-data me-2"></i>Matriz L (triangular inferior)</h5>
      <div class="table-responsive">
        <table class="table table-dark table-sm">
          <thead><tr>{% for j in range(n) %}<th>L[·,{{ j+1 }}]</th>{% endfor %}</tr></thead>
          <tbody>{% for i in range(n) %}<tr>
            {% for j in range(n) %}<td>{{ "%.6f"|format(L[i][j]) }}</td>{% endfor %}
          </tr>{% endfor %}</tbody>
        </table>
      </div>
    </div>

    <div class="card p-3">
      <h5><i class="bi bi-clipboard-data me-2"></i>Matriz U e vetor c</h5>
      <div class="table-responsive">
        <table class="table table-dark table-sm">
          <thead><tr>
            {% for j in range(n) %}<th>U[·,{{ j+1 }}]</th>{% endfor %}<th>= c</th>
          </tr></thead>
          <tbody>{% for i in range(n) %}<tr>
            {% for j in range(n) %}<td>{{ "%.6f"|format(U[i][j]) }}</td>{% endfor %}
            <td>{{ "%.6f"|format(c[i]) }}</td>
          </tr>{% endfor %}</tbody>
        </table>
      </div>
    </div>
  </div>

  <div class="col-lg-6">
    </div>
  </div>
</div>
{% endif %}
""" + BASE_FOOT

# ── Decomposição LU ──────────────────────────────────────────────

LU_TEMPLATE = BASE_HEAD + """
<h1 class="mb-1"></i>Decomposição LU</h1>
<p class="text-muted mb-4">Fatoração A = L·U e resolução por substituição (sem pivotamento)</p>

{% if stage == 'choose_n' %}
<div class="card p-4" style="max-width:400px">
  <form method="post" class="row g-3">
    <input type="hidden" name="stage" value="enter_matrix">
    <div class="col-auto">
      <label class="form-label fw-semibold">Tamanho n × n</label>
      <select name="n" class="form-select rounded-pill">
        {% for i in range(2,11) %}
          <option value="{{ i }}" {{ 'selected' if i==n }}>{{ i }}</option>
        {% endfor %}
      </select>
    </div>
    <div class="col-12">
      <button class="btn btn-warning rounded-pill"><i class="bi bi-grid me-1"></i>Criar matriz</button>
      <a href="/" class="btn btn-outline-secondary rounded-pill ms-2">Menu</a>
    </div>
  </form>
</div>

{% elif stage == 'enter_matrix' %}
<div class="card p-4">
  <form method="post">
    <input type="hidden" name="stage" value="solve">
    <input type="hidden" name="n" value="{{ n }}">
    <p class="mb-2">Informe <strong>A</strong> e <strong>b</strong> (sistema A·x = b):</p>
    <div class="table-responsive mb-3">
      <table class="table table-dark table-striped align-middle">
        <thead><tr>
          {% for j in range(n) %}<th>A[·,{{ j+1 }}]</th>{% endfor %}
          <th>= b</th>
        </tr></thead>
        <tbody>
          {% for i in range(n) %}<tr>
            {% for j in range(n) %}
              <td><input required type="text" class="form-control form-control-sm"
                   name="a_{{ i }}_{{ j }}" value="{{ '1' if i==j else '0' }}"></td>
            {% endfor %}
            <td><input required type="text" class="form-control form-control-sm"
                 name="b_{{ i }}" value="0"></td>
          </tr>{% endfor %}
        </tbody>
      </table>
    </div>
    <div class="d-flex gap-2">
      <button class="btn btn-success rounded-pill"><i class="bi bi-play-fill me-1"></i>Resolver</button>
      <a href="/lu" class="btn btn-outline-secondary rounded-pill">Reiniciar</a>
      <a href="/" class="btn btn-outline-dark rounded-pill">Menu</a>
    </div>
  </form>
</div>

{% elif stage == 'solve' %}
<div class="row g-4">
  <div class="col-lg-6">
    <div class="card p-3 mb-3">
      <h5><i class="bi bi-list-check me-2"></i>Passos da fatoração</h5>
      <ol>{% for s in steps %}<li>{{ s }}</li>{% endfor %}</ol>
    </div>
    <div class="card p-3 mb-3">
      <h5>Matriz L (triangular inferior)</h5>
      <div class="table-responsive">
        <table class="table table-dark table-sm">
          <thead><tr>{% for j in range(n) %}<th>L[·,{{ j+1 }}]</th>{% endfor %}</tr></thead>
          <tbody>{% for i in range(n) %}<tr>
            {% for j in range(n) %}<td>{{ "%.6f"|format(L[i][j]) }}</td>{% endfor %}
          </tr>{% endfor %}</tbody>
        </table>
      </div>
    </div>
    <div class="card p-3">
      <h5>Matriz U (triangular superior)</h5>
      <div class="table-responsive">
        <table class="table table-dark table-sm">
          <thead><tr>{% for j in range(n) %}<th>U[·,{{ j+1 }}]</th>{% endfor %}</tr></thead>
          <tbody>{% for i in range(n) %}<tr>
            {% for j in range(n) %}<td>{{ "%.6f"|format(U[i][j]) }}</td>{% endfor %}
          </tr>{% endfor %}</tbody>
        </table>
      </div>
    </div>
  </div>
  </div>
</div>
{% endif %}
""" + BASE_FOOT

# ── Jacobi / Gauss-Seidel (template compartilhado) ───────────────

ITER_TEMPLATE = BASE_HEAD + """
{% if metodo == 'jacobi' %}
  <h1 class="mb-1"></i>Método de Jacobi</h1>
  <p class="text-muted mb-4">Método iterativo — critério das linhas (diagonal dominante)</p>
{% else %}
  <h1 class="mb-1">Gauss-Seidel</h1>
  <p class="text-muted mb-4">Método iterativo — critério de Sassenfeld</p>
{% endif %}

{% if stage == 'choose_n' %}
<div class="card p-4" style="max-width:400px">
  <form method="post" class="row g-3">
    <input type="hidden" name="stage" value="enter_matrix">
    <div class="col-auto">
      <label class="form-label fw-semibold">Tamanho n × n</label>
      <select name="n" class="form-select rounded-pill">
        {% for i in range(2,11) %}
          <option value="{{ i }}" {{ 'selected' if i==n }}>{{ i }}</option>
        {% endfor %}
      </select>
    </div>
    <div class="col-12">
      <button class="btn btn-primary rounded-pill"><i class="bi bi-grid me-1"></i>Criar matriz</button>
      <a href="/" class="btn btn-outline-secondary rounded-pill ms-2">Menu</a>
    </div>
  </form>
</div>

{% elif stage == 'enter_matrix' %}
<div class="card p-4">
  <form method="post">
    <input type="hidden" name="stage" value="solve">
    <input type="hidden" name="n" value="{{ n }}">
    <p class="mb-2">Informe <strong>A</strong>, <strong>b</strong> e o <strong>chute inicial x⁰</strong>:</p>
    <div class="table-responsive mb-3">
      <table class="table table-dark table-striped align-middle">
        <thead><tr>
          {% for j in range(n) %}<th>A[·,{{ j+1 }}]</th>{% endfor %}
          <th>= b</th><th>x⁰</th>
        </tr></thead>
        <tbody>
          {% for i in range(n) %}<tr>
            {% for j in range(n) %}
              <td><input required type="text" class="form-control form-control-sm"
                   name="a_{{ i }}_{{ j }}" value="{{ '1' if i==j else '0' }}"></td>
            {% endfor %}
            <td><input required type="text" class="form-control form-control-sm"
                 name="b_{{ i }}" value="0"></td>
            <td><input required type="text" class="form-control form-control-sm"
                 name="x0_{{ i }}" value="0"></td>
          </tr>{% endfor %}
        </tbody>
      </table>
    </div>
    <div class="row g-3 mb-4">
      <div class="col-auto">
        <label class="form-label">Tolerância</label>
        <input type="number" step="any" class="form-control" name="tol" value="0.05" required>
      </div>
      <div class="col-auto">
        <label class="form-label">Máx. de iterações</label>
        <input type="number" class="form-control" name="maxIter" value="1000" required>
      </div>
    </div>
    <div class="d-flex gap-2">
      <button class="btn btn-success rounded-pill"><i class="bi bi-play-fill me-1"></i>Resolver</button>
      <a href="/{{ metodo }}" class="btn btn-outline-secondary rounded-pill">Reiniciar</a>
      <a href="/" class="btn btn-outline-dark rounded-pill">Menu</a>
    </div>
  </form>
</div>

{% elif stage == 'solve' %}
  {% if erro %}
    <div class="alert alert-danger">{{ erro }}</div>
  {% else %}
    {% if convergencia %}
      <div class="alert alert-success"> Critério de convergência satisfeito.</div>
    {% else %}
      <div class="alert alert-warning"> Critério de convergência NÃO satisfeito — resultado pode divergir.</div>
    {% endif %}

    {% if convergiu %}
      <div class="alert alert-info">Convergência atingida na iteração {{ iteracoes|length }}.</div>
    {% else %}
      <div class="alert alert-warning">Limite de iterações atingido sem convergência.</div>
    {% endif %}

    <div class="row g-4">
      <div class="col-lg-8">
        <div class="card p-3">
          <h5><i class="bi bi-table me-2"></i>Iterações</h5>
          <div class="table-responsive">
            <table class="table table-striped table-hover align-middle mb-0">
              <thead class="table-dark">
                <tr>
                  <th>k</th>
                  {% for i in range(n) %}<th>x{{ i+1 }}</th>{% endfor %}
                  <th>Erro</th>
                </tr>
              </thead>
              <tbody>
                {% for r in iteracoes %}
                <tr>
                  <td>{{ r.k }}</td>
                  {% for v in r.x %}<td>{{ "%.6f"|format(v) }}</td>{% endfor %}
                  <td>{{ "%.2e"|format(r.erro) }}</td>
                </tr>
                {% endfor %}
              </tbody>
            </table>
          </div>
        </div>
      </div>
      
    </div>
  {% endif %}
{% endif %}
""" + BASE_FOOT


# ════════════════════════════════════════════════════════════════
#  ROTAS FLASK
# ════════════════════════════════════════════════════════════════

@app.route("/")
def menu():
    return render_template_string(MENU_TEMPLATE)


# ── Bisseção ─────────────────────────────────────────────────────

@app.route("/bisseccao", methods=["GET", "POST"])
def rota_bisseccao():
    DEFAULTS = {"funcao": "x**3 - x - 2", "a": "1", "b": "2", "tol": "1e-6", "maxIter": "50"}
    stage   = request.form.get("stage", "menu")
    funcao  = request.form.get("funcao",  DEFAULTS["funcao"])
    a_str   = request.form.get("a",       DEFAULTS["a"])
    b_str   = request.form.get("b",       DEFAULTS["b"])
    tol_str = request.form.get("tol",     DEFAULTS["tol"])
    max_str = request.form.get("maxIter", DEFAULTS["maxIter"])
    ctx = dict(stage=stage, funcao=funcao, a=a_str, b=b_str,
               tol=tol_str, maxIter=max_str, erro=None, raiz=None, iteracoes=[])
    if stage == "solve":
        try:
            a = float(a_str.replace(",", ".")); b = float(b_str.replace(",", "."))
            tol = float(tol_str.replace(",", ".")); max_iter = int(max_str)
        except ValueError:
            ctx["erro"] = "Valores inválidos."; return render_template_string(BISSECCAO_TEMPLATE, **ctx)
        res = bisseccao(a, b, tol, max_iter, funcao)
        if "erro" in res:
            ctx["erro"] = res["erro"]
        else:
            ctx["raiz"] = res["raiz"]; ctx["iteracoes"] = res["iteracoes"]
    return render_template_string(BISSECCAO_TEMPLATE, **ctx)


# ── Gauss ────────────────────────────────────────────────────────

@app.route("/gauss", methods=["GET", "POST"])
def rota_gauss():
    stage = request.form.get("stage", "choose_n")
    n     = max(2, min(10, int(request.form.get("n", 3))))
    use_P = request.form.get("use_P", "0") == "1"
    ctx   = dict(stage=stage, n=n, use_P=use_P,
                 steps=[], x=None, L=None, U=None, P=None, c=None)
    if stage == "solve":
        A = [[float(request.form.get(f"a_{i}_{j}", 0)) for j in range(n)] for i in range(n)]
        b = [float(request.form.get(f"b_{i}", 0)) for i in range(n)]
        x, L, U, P, c, steps = gauss_elimination(A, b, use_P=use_P)
        ctx.update(x=x, L=L, U=U, P=P, c=c, steps=steps)
    return render_template_string(GAUSS_TEMPLATE, **ctx)


# ── Decomposição LU ──────────────────────────────────────────────

@app.route("/lu", methods=["GET", "POST"])
def rota_lu():
    stage = request.form.get("stage", "choose_n")
    n     = max(2, min(10, int(request.form.get("n", 3))))
    ctx   = dict(stage=stage, n=n, steps=[], x=None, L=None, U=None)
    if stage == "solve":
        A = [[float(request.form.get(f"a_{i}_{j}", 0)) for j in range(n)] for i in range(n)]
        b = [float(request.form.get(f"b_{i}", 0)) for i in range(n)]
        x, L, U, steps = lu_decomposition(A, b)
        ctx.update(x=x, L=L, U=U, steps=steps)
    return render_template_string(LU_TEMPLATE, **ctx)


# ── Jacobi ───────────────────────────────────────────────────────

@app.route("/jacobi", methods=["GET", "POST"])
def rota_jacobi():
    stage = request.form.get("stage", "choose_n")
    n     = max(2, min(10, int(request.form.get("n", 3))))
    ctx   = dict(metodo="jacobi", stage=stage, n=n,
                 erro=None, x=None, iteracoes=[], convergencia=False, convergiu=False)
    if stage == "solve":
        try:
            A   = [[float(request.form.get(f"a_{i}_{j}", 0)) for j in range(n)] for i in range(n)]
            b   = [float(request.form.get(f"b_{i}", 0)) for i in range(n)]
            x0  = [float(request.form.get(f"x0_{i}", 0)) for i in range(n)]
            tol = float(request.form.get("tol", 0.05))
            mi  = int(request.form.get("maxIter", 1000))
        except ValueError:
            ctx["erro"] = "Valores inválidos."; return render_template_string(ITER_TEMPLATE, **ctx)
        res = jacobi(n, A, b, x0, tol, mi)
        ctx.update(x=res["x"], iteracoes=res["iteracoes"],
                   convergencia=res["convergencia"], convergiu=res["convergiu"])
    return render_template_string(ITER_TEMPLATE, **ctx)


# ── Gauss-Seidel ─────────────────────────────────────────────────

@app.route("/seidel", methods=["GET", "POST"])
def rota_seidel():
    stage = request.form.get("stage", "choose_n")
    n     = max(2, min(10, int(request.form.get("n", 3))))
    ctx   = dict(metodo="seidel", stage=stage, n=n,
                 erro=None, x=None, iteracoes=[], convergencia=False, convergiu=False)
    if stage == "solve":
        try:
            A   = [[float(request.form.get(f"a_{i}_{j}", 0)) for j in range(n)] for i in range(n)]
            b   = [float(request.form.get(f"b_{i}", 0)) for i in range(n)]
            x0  = [float(request.form.get(f"x0_{i}", 0)) for i in range(n)]
            tol = float(request.form.get("tol", 0.05))
            mi  = int(request.form.get("maxIter", 1000))
        except ValueError:
            ctx["erro"] = "Valores inválidos."; return render_template_string(ITER_TEMPLATE, **ctx)
        res = gauss_seidel(n, A, b, x0, tol, mi)
        ctx.update(x=res["x"], iteracoes=res["iteracoes"],
                   convergencia=res["convergencia"], convergiu=res["convergiu"])
    return render_template_string(ITER_TEMPLATE, **ctx)


if __name__ == "__main__":
    app.run(debug=True)