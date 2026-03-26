from flask import Flask, request, render_template_string
import math

app = Flask(__name__)

# ─── Funções de cálculo ───────────────────────────────────────────────────────

def is_near_zero(v, eps=1e-10):
    return abs(v) < eps


def render_matrix_state(A, b):
    n = len(A)
    html = '<p class="text-black"><strong>Estado Atual da Matriz</strong></p>'
    html += '<table class="table table-dark table-sm align-middle mt-2 mb-3"><thead><tr>'
    for j in range(n):
        html += f'<th class="text-white">U[·,{j+1}]</th>'
    html += '<th class="text-white">= c</th></tr></thead><tbody>'
    for i in range(n):
        html += '<tr>'
        for j in range(n):
            html += f'<td class="text-white">{A[i][j]:.6f}</td>'
        html += f'<td class="text-white">{b[i]:.6f}</td>'
        html += '</tr>'
    html += '</tbody></table>'
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
                    max_val = abs(A[i][k])
                    p = i
            if p != k:
                A[k], A[p] = A[p], A[k]
                b[k], b[p] = b[p], b[k]
                P[k], P[p] = P[p], P[k]
                for col in range(k):
                    L[k][col], L[p][col] = L[p][col], L[k][col]
                steps.append({'type': 'text', 'content': f'Troca L{k+1} ⇄ L{p+1} (pivotamento parcial)'})

        if is_near_zero(A[k][k]):
            steps.append({'type': 'text', 'content': f'Coluna {k}: pivô ~ 0 ⇒ sistema singular/indeterminado.'})
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

    # Retrossubstituição
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i + 1, n))
        if is_near_zero(A[i][i]):
            steps.append({'type': 'text', 'content': 'Pivô zero na retrossubstituição ⇒ sistema singular.'})
            return None, L, A, P, b, steps
        x[i] = (b[i] - s) / A[i][i]

    return x, L, A, P, b, steps


# ─── Template HTML ────────────────────────────────────────────────────────────

TEMPLATE = """
<!doctype html>
<html lang="pt-br">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Eliminação de Gauss</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.css" rel="stylesheet">
  <style>
    h1, h2, h3, h4, h5, h6 { color: #000 !important; }
  </style>
</head>
<body>
<div class="container py-4">
  <div class="d-flex align-items-center mb-4">
    <i class="bi bi-calculator-fill fs-2 me-2 text-info"></i>
    <h1 class="h3 m-0">Eliminação de Gauss
      <small class="text-secondary">({{ 'pivotamento parcial' if use_P else 'simples' }})</small>
    </h1>
  </div>

  <div class="card rounded-2xl shadow p-3 mb-4">
    <div class="card-body">

      {# ── ETAPA 1: escolher n ── #}
      {% if stage == 'choose_n' %}
      <form method="post" class="row g-3">
        <input type="hidden" name="stage" value="enter_matrix">
        <div class="col-auto">
          <label class="form-label" for="n">Tamanho do sistema (n × n):</label>
          <select name="n" id="n" class="form-select rounded-pill">
            {% for i in range(2, 11) %}
              <option value="{{ i }}" {{ 'selected' if i == n }}>{{ i }}</option>
            {% endfor %}
          </select>
        </div>
        <div class="col-auto align-self-end">
          <div class="form-check">
            <input class="form-check-input" type="checkbox" value="1" id="use_P" name="use_P" {{ 'checked' if use_P }}>
            <label class="form-check-label" for="use_P">
              Usar matriz de permutação (pivotamento parcial)
            </label>
          </div>
        </div>
        <div class="col-auto align-self-end">
          <button class="btn btn-info rounded-pill">
            <i class="bi bi-grid-3x3-gap me-1"></i>Criar matriz
          </button>
        </div>
      </form>

      {# ── ETAPA 2: preencher matriz ── #}
      {% elif stage == 'enter_matrix' %}
      <form method="post">
        <input type="hidden" name="stage" value="solve">
        <input type="hidden" name="n" value="{{ n }}">
        <input type="hidden" name="use_P" value="{{ 1 if use_P else 0 }}">
        <div class="mb-3">
          <div class="d-flex align-items-center mb-2">
            <span class="me-2">Informe a matriz <strong>A</strong> e o vetor <strong>b</strong> (sistema A·x = b):</span>
            <span class="hint">Use números reais (p. ex.: 2, -3.5, 1e-3).</span>
          </div>
          <div class="table-responsive">
            <table class="table table-dark table-striped align-middle">
              <thead>
                <tr>
                  {% for j in range(n) %}
                    <th class="text-white">A[·,{{ j+1 }}]</th>
                  {% endfor %}
                  <th class="text-white">= b</th>
                </tr>
              </thead>
              <tbody>
                {% for i in range(n) %}
                  <tr>
                    {% for j in range(n) %}
                      <td class="text-white">
                        <input required type="text" class="form-control form-control-sm"
                               name="a_{{ i }}_{{ j }}"
                               value="{{ '1' if i == j else '0' }}">
                      </td>
                    {% endfor %}
                    <td class="text-white">
                      <input required type="text" class="form-control form-control-sm"
                             name="b_{{ i }}" value="0">
                    </td>
                  </tr>
                {% endfor %}
              </tbody>
            </table>
          </div>
        </div>
        <div class="d-flex gap-2">
          <button class="btn btn-success rounded-pill">
            <i class="bi bi-play-fill me-1"></i>Resolver
          </button>
          <a href="/" class="btn btn-secondary rounded-pill">
            <i class="bi bi-arrow-counterclockwise me-1"></i>Reiniciar
          </a>
        </div>
      </form>

      {# ── ETAPA 3: resultado ── #}
      {% elif stage == 'solve' %}
      <div class="row g-3">
        <div class="col-12 col-lg-6">

          {# Passos #}
          <div class="card rounded-2xl mb-3">
            <div class="card-body">
              <h2 class="h5 mb-3 text-black"><i class="bi bi-list-check me-2"></i>Passos da eliminação</h2>
              {% if steps %}
                <ol class="mb-0">
                  {% for s in steps %}
                    {% if s.type == 'matrix' %}
                      <div>{{ s.content | safe }}</div>
                    {% else %}
                      <li>{{ s.content }}</li>
                    {% endif %}
                  {% endfor %}
                </ol>
              {% else %}
                <p class="text-secondary">Nenhuma operação registrada.</p>
              {% endif %}
            </div>
          </div>

          {# Matriz P #}
          {% if use_P and P %}
          <div class="card rounded-2xl mb-3">
            <div class="card-body">
              <h2 class="h5 mb-3 text-black"><i class="bi bi-clipboard-data me-2"></i>Matriz de Permutação P</h2>
              <div class="table-responsive">
                <table class="table table-dark table-sm align-middle">
                  <thead><tr>
                    {% for j in range(n) %}<th class="text-white">P[·,{{ j+1 }}]</th>{% endfor %}
                  </tr></thead>
                  <tbody>
                    {% for i in range(n) %}
                      <tr>{% for j in range(n) %}
                        <td class="text-white">{{ "%.6f"|format(P[i][j]) }}</td>
                      {% endfor %}</tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
          {% endif %}

          {# Matriz L #}
          <div class="card rounded-2xl mb-3">
            <div class="card-body">
              <h2 class="h5 mb-3 text-black"><i class="bi bi-clipboard-data me-2"></i>Matriz triangular inferior L</h2>
              <div class="table-responsive">
                <table class="table table-dark table-sm align-middle">
                  <thead><tr>
                    {% for j in range(n) %}<th class="text-white">L[·,{{ j+1 }}]</th>{% endfor %}
                  </tr></thead>
                  <tbody>
                    {% for i in range(n) %}
                      <tr>{% for j in range(n) %}
                        <td class="text-white">{{ "%.6f"|format(L[i][j]) }}</td>
                      {% endfor %}</tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>
            </div>
          </div>

          {# Matriz U e vetor c #}
          <div class="card rounded-2xl">
            <div class="card-body">
              <h2 class="h5 mb-3 text-black"><i class="bi bi-clipboard-data me-2"></i>Matriz triangular superior U e vetor c</h2>
              <div class="table-responsive">
                <table class="table table-dark table-sm align-middle">
                  <thead><tr>
                    {% for j in range(n) %}<th class="text-white">U[·,{{ j+1 }}]</th>{% endfor %}
                    <th class="text-white">= c</th>
                  </tr></thead>
                  <tbody>
                    {% for i in range(n) %}
                      <tr>
                        {% for j in range(n) %}
                          <td class="text-white">{{ "%.6f"|format(U[i][j]) }}</td>
                        {% endfor %}
                        <td class="text-white">{{ "%.6f"|format(c[i]) }}</td>
                      </tr>
                    {% endfor %}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        </div>

        {# Resultado #}
        <div class="col-12 col-lg-6">
          <div class="card rounded-2xl h-100">
            <div class="card-body d-flex flex-column">
              <h2 class="h5 mb-3 text-black"><i class="bi bi-check2-circle me-2"></i>Resultado</h2>
              {% if x is none %}
                <div class="alert alert-warning rounded-3" role="alert">
                  O sistema parece singular ou indeterminado. Verifique os coeficientes.
                </div>
              {% else %}
                <p class="mb-2 text-black">Solução aproximada (precisão 6 casas decimais):</p>
                <ul class="list-group mb-3">
                  {% for i in range(n) %}
                    <li class="list-group-item d-flex justify-content-between align-items-center">
                      x{{ i+1 }}
                      <span class="badge bg-info rounded-pill">{{ "%.6f"|format(x[i]) }}</span>
                    </li>
                  {% endfor %}
                </ul>
                <div class="mt-auto d-flex gap-2">
                  <form method="post">
                    <input type="hidden" name="stage" value="enter_matrix">
                    <input type="hidden" name="n" value="{{ n }}">
                    <input type="hidden" name="use_P" value="{{ 1 if use_P else 0 }}">
                    <button class="btn btn-secondary rounded-pill">
                      <i class="bi bi-arrow-counterclockwise me-1"></i>Voltar
                    </button>
                  </form>
                  <a href="/" class="btn btn-outline-secondary rounded-pill">
                    <i class="bi bi-x-circle me-1"></i>Reiniciar
                  </a>
                </div>
              {% endif %}
            </div>
          </div>
        </div>
      </div>
      {% endif %}

    </div>
  </div>
</div>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
"""

# ─── Rotas Flask ──────────────────────────────────────────────────────────────

@app.route("/", methods=["GET", "POST"])
def index():
    # Valores padrão
    stage = request.form.get("stage", "choose_n")
    n = int(request.form.get("n", 3))
    n = max(2, min(10, n))
    use_P = request.form.get("use_P", "0") == "1"

    ctx = dict(stage=stage, n=n, use_P=use_P,
               steps=[], x=None, L=None, U=None, P=None, c=None)

    if stage == "solve":
        A = []
        b = []
        for i in range(n):
            row = []
            for j in range(n):
                row.append(float(request.form.get(f"a_{i}_{j}", 0)))
            A.append(row)
            b.append(float(request.form.get(f"b_{i}", 0)))

        x, L, U, P, c, steps = gauss_elimination(A, b, use_P=use_P)
        ctx.update(x=x, L=L, U=U, P=P, c=c, steps=steps)

    return render_template_string(TEMPLATE, **ctx)


if __name__ == "__main__":
    app.run(debug=True)