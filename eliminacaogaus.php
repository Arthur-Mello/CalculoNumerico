<?php

function is_near_zero($v, $eps = 1e-10) { return abs($v) < $eps; }

function gauss_elimination(array $A, array $b, &$steps = []) {
    $n = count($A);
    for ($i = 0; $i < $n; $i++) { $A[$i] = array_map('floatval', $A[$i]); $b[$i] = floatval($b[$i]); }
    for ($k = 0; $k < $n; $k++) {
        $p = $k; $maxVal = abs($A[$k][$k]);
        for ($i = $k + 1; $i < $n; $i++) {
            if (abs($A[$i][$k]) > $maxVal) { $maxVal = abs($A[$i][$k]); $p = $i; }
        }

        if (is_near_zero($maxVal)) {
            $steps[] = "Coluna $k: pivô ~ 0 ⇒ sistema singular/indeterminado.";
            return [null, $A, $b];
        }

        if ($p !== $k) {
            [$A[$k], $A[$p]] = [$A[$p], $A[$k]];
            [$b[$k], $b[$p]] = [$b[$p], $b[$k]];
            $steps[] = sprintf("Troca L%d ⇄ L%d (pivotamento parcial)", $k+1, $p+1);
        }


        for ($i = $k + 1; $i < $n; $i++) {
            if (is_near_zero($A[$i][$k])) continue;
            $m = $A[$i][$k] / $A[$k][$k];
            for ($j = $k; $j < $n; $j++) {
                $A[$i][$j] -= $m * $A[$k][$j];
            }
            $b[$i] -= $m * $b[$k];
            $steps[] = sprintf("L%d ← L%d − (%.6g)·L%d", $i+1, $i+1, $m, $k+1);
        }
    }

    $x = array_fill(0, $n, 0.0);
    for ($i = $n - 1; $i >= 0; $i--) {
        $sum = 0.0;
        for ($j = $i + 1; $j < $n; $j++) { $sum += $A[$i][$j] * $x[$j]; }
        if (is_near_zero($A[$i][$i])) {
            $steps[] = "Encontrado pivô zero na retrossubstituição ⇒ sistema singular.";
            return [null, $A, $b];
        }
        $x[$i] = ($b[$i] - $sum) / $A[$i][$i];
    }

    return [$x, $A, $b];
}

function h($s) { return htmlspecialchars((string)$s, ENT_QUOTES, 'UTF-8'); }
function get_post($key, $default = null) { return $_POST[$key] ?? $default; }

$stage = get_post('stage', 'choose_n');
$n = intval(get_post('n', 3));
if ($n < 2) $n = 2; if ($n > 8) $n = 8;

?>
<!doctype html>
<html lang="pt-br">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Eliminação de Gauss em PHP</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.css" rel="stylesheet">
  <style>
    body { background: #0f172a; color: #e2e8f0; }
    .card { background: #111827; border: 1px solid #1f2937; }
    .form-control, .form-select { background: #0b1220; color: #e5e7eb; border-color: #1f2937; }
    .rounded-2xl { border-radius: 1rem; }
    .matrix-cell { width: 7rem; }
    .hint { color:#93c5fd; }
    a, .btn-link { color:#93c5fd; }
  </style>
</head>
<body>
<div class="container py-4">
  <div class="d-flex align-items-center mb-4">
    <i class="bi bi-calculator-fill fs-2 me-2 text-info"></i>
    <h1 class="h3 m-0">Eliminação de Gauss <small class="text-secondary">(pivotamento parcial)</small></h1>
  </div>

  <div class="card rounded-2xl shadow p-3 mb-4">
    <div class="card-body">
      <?php if ($stage === 'choose_n'): ?>
        <form method="post" class="row g-3">
          <input type="hidden" name="stage" value="enter_matrix">
          <div class="col-auto">
            <label class="form-label" for="n">Tamanho do sistema (n × n):</label>
            <select name="n" id="n" class="form-select rounded-pill">
              <?php for ($i=2; $i<=8; $i++): ?>
                <option value="<?= $i ?>" <?= $i===$n? 'selected':'' ?>><?= $i ?></option>
              <?php endfor; ?>
            </select>
          </div>
          <div class="col-auto align-self-end">
            <button class="btn btn-info rounded-pill"><i class="bi bi-grid-3x3-gap me-1"></i>Criar matriz</button>
          </div>
        </form>
      <?php elseif ($stage === 'enter_matrix'): ?>
        <form method="post">
          <input type="hidden" name="stage" value="solve">
          <input type="hidden" name="n" value="<?= $n ?>">
          <div class="mb-3">
            <div class="d-flex align-items-center mb-2">
              <span class="me-2">Informe a matriz <strong>A</strong> e o vetor <strong>b</strong> (sistema A·x = b):</span>
              <span class="hint">Use números reais (p. ex.: 2, -3.5, 1e-3).</span>
            </div>
            <div class="table-responsive">
              <table class="table table-dark table-striped align-middle">
                <thead>
                  <tr>
                    <?php for ($j=0; $j<$n; $j++): ?>
                      <th>A[·,<?= $j+1 ?>]</th>
                    <?php endfor; ?>
                    <th>= b</th>
                  </tr>
                </thead>
                <tbody>
                  <?php for ($i=0; $i<$n; $i++): ?>
                    <tr>
                      <?php for ($j=0; $j<$n; $j++): ?>
                        <td>
                          <input required type="text" class="form-control form-control-sm matrix-cell" name="a_<?= $i ?>_<?= $j ?>" value="<?= $i===$j? '1': '0' ?>">
                        </td>
                      <?php endfor; ?>
                      <td>
                        <input required type="text" class="form-control form-control-sm matrix-cell" name="b_<?= $i ?>" value="0">
                      </td>
                    </tr>
                  <?php endfor; ?>
                </tbody>
              </table>
            </div>
          </div>
          <div class="d-flex gap-2">
            <button class="btn btn-success rounded-pill"><i class="bi bi-play-fill me-1"></i>Resolver</button>
            <a href="?" class="btn btn-secondary rounded-pill"><i class="bi bi-arrow-counterclockwise me-1"></i>Reiniciar</a>
          </div>
        </form>
      <?php elseif ($stage === 'solve'): ?>
        <?php
          // Coletar A e b do POST
          $A = []; $bvec = [];
          for ($i=0; $i<$n; $i++) {
            $row = [];
            for ($j=0; $j<$n; $j++) {
              $row[] = floatval(get_post("a_{$i}_{$j}", 0));
            }
            $A[] = $row;
            $bvec[] = floatval(get_post("b_{$i}", 0));
          }
          $steps = [];
          [$x, $U, $c] = gauss_elimination($A, $bvec, $steps);
        ?>
        <div class="row g-3">
          <div class="col-12 col-lg-6">
            <div class="card rounded-2xl mb-3">
              <div class="card-body">
                <h2 class="h5 mb-3"><i class="bi bi-list-check me-2"></i>Passos da eliminação</h2>
                <?php if (!empty($steps)): ?>
                  <ol class="mb-0">
                    <?php foreach ($steps as $s): ?>
                      <li><?= h($s) ?></li>
                    <?php endforeach; ?>
                  </ol>
                <?php else: ?>
                  <p class="text-secondary">Nenhuma operação registrada.</p>
                <?php endif; ?>
              </div>
            </div>

            <div class="card rounded-2xl">
              <div class="card-body">
                <h2 class="h5 mb-3"><i class="bi bi-clipboard-data me-2"></i>Matriz triangular superior U e vetor c</h2>
                <div class="table-responsive">
                  <table class="table table-dark table-sm align-middle">
                    <thead>
                      <tr>
                        <?php for ($j=0; $j<$n; $j++): ?><th>U[·,<?= $j+1 ?>]</th><?php endfor; ?><th>= c</th>
                      </tr>
                    </thead>
                    <tbody>
                      <?php for ($i=0; $i<$n; $i++): ?>
                        <tr>
                          <?php for ($j=0; $j<$n; $j++): ?>
                            <td><?= h(number_format($U[$i][$j], 6, '.', '')) ?></td>
                          <?php endfor; ?>
                          <td><?= h(number_format($c[$i], 6, '.', '')) ?></td>
                        </tr>
                      <?php endfor; ?>
                    </tbody>
                  </table>
                </div>
              </div>
            </div>
          </div>

          <div class="col-12 col-lg-6">
            <div class="card rounded-2xl h-100">
              <div class="card-body d-flex flex-column">
                <h2 class="h5 mb-3"><i class="bi bi-check2-circle me-2"></i>Resultado</h2>
                <?php if ($x === null): ?>
                  <div class="alert alert-warning rounded-3" role="alert">
                    O sistema parece singular ou indeterminado. Verifique os coeficientes.
                  </div>
                <?php else: ?>
                  <p class="mb-2">Solução aproximada (precisão 6 casas decimais):</p>
                  <ul class="list-group mb-3">
                    <?php for ($i=0; $i<$n; $i++): ?>
                      <li class="list-group-item d-flex justify-content-between align-items-center">
                        x<?= $i+1 ?>
                        <span class="badge bg-info rounded-pill"><?= h(number_format($x[$i], 6, '.', '')) ?></span>
                      </li>
                    <?php endfor; ?>
                  </ul>
                  <div class="mt-auto d-flex gap-2">
                    <form method="post">
                      <input type="hidden" name="stage" value="enter_matrix">
                      <input type="hidden" name="n" value="<?= $n ?>">
                      <?php for ($i=0; $i<$n; $i++): for ($j=0; $j<$n; $j++): ?>
                        <input type="hidden" name="a_<?= $i ?>_<?= $j ?>" value="<?= h($A[$i][$j]) ?>">
                      <?php endfor; endfor; for ($i=0; $i<$n; $i++): ?>
                        <input type="hidden" name="b_<?= $i ?>" value="<?= h($bvec[$i]) ?>">
                      <?php endfor; ?>
                      <button class="btn btn-outline-info rounded-pill" title="Editar valores"><i class="bi bi-pencil-square me-1"></i>Editar</button>
                    </form>
                    <a class="btn btn-secondary rounded-pill" href="?"><i class="bi bi-arrow-counterclockwise me-1"></i>Novo cálculo</a>
                  </div>
                <?php endif; ?>
              </div>
            </div>
          </div>
        </div>
      <?php endif; ?>
    </div>
  </div>
</div>
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
