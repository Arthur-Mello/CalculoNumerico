<?php
function f($x, $expr) {
    $expr = str_replace('^', '**', $expr);
    try {
        $result = null;
        eval('$result = ' . $expr . ';');
        return $result;
    } catch (Throwable $e) {
        return NAN;
    }
}
function bisseccao($a, $b, $tol, $maxIter, $expr) {
    $fa = f($a, $expr);
    $fb = f($b, $expr);
    if (!is_numeric($fa) || !is_numeric($fb)) {
        return ["erro" => "Erro ao avaliar a função nos extremos do intervalo."];
    }
    if ($fa * $fb >= 0) {
        return ["erro" => "O intervalo não satisfaz f(a) * f(b) < 0"];
    }

    $iteracoes = [];

    for ($i = 1; $i <= $maxIter; $i++) {
        $p = ($a + $b) / 2;
        $fp = f($p, $expr);
        $erro = abs($b - $a) / 2;

        $iteracoes[] = [
            "k" => $i,
            "a" => $a,
            "b" => $b,
            "p" => $p,
            "f(p)" => $fp,
            "erro" => $erro
        ];

        if ($erro < $tol || $fp == 0.0) {
            return ["raiz" => $p, "iteracoes" => $iteracoes];
        }

        if ($fa * $fp < 0) {
            $b = $p;
            $fb = $fp;
        } else {
            $a = $p;
            $fa = $fp;
        }
    }
    return ["raiz" => ($a + $b) / 2, "iteracoes" => $iteracoes];
}

$stage = isset($_REQUEST['stage']) ? $_REQUEST['stage'] : 'menu';

$default_funcao = 'pow($x,3) - $x - 2';
$default_a = '1';
$default_b = '2';
$default_tol = '1e-6';
$default_maxIter = '50';
$input_funcao = isset($_REQUEST['funcao']) ? $_REQUEST['funcao'] : $default_funcao;
$input_a = isset($_REQUEST['a']) ? $_REQUEST['a'] : $default_a;
$input_b = isset($_REQUEST['b']) ? $_REQUEST['b'] : $default_b;
$input_tol = isset($_REQUEST['tol']) ? $_REQUEST['tol'] : $default_tol;
$input_maxIter = isset($_REQUEST['maxIter']) ? $_REQUEST['maxIter'] : $default_maxIter;
$resultado = null;
if ($stage === 'solve') {
    $funcao = trim($_POST['funcao'] ?? $default_funcao);
    $a = floatval(str_replace(',', '.', $_POST['a'] ?? $default_a));
    $b = floatval(str_replace(',', '.', $_POST['b'] ?? $default_b));
    $tol = floatval(str_replace(',', '.', $_POST['tol'] ?? $default_tol));
    $maxIter = intval($_POST['maxIter'] ?? $default_maxIter);
    $resultado = bisseccao($a, $b, $tol, $maxIter, $funcao);
}
?>

<!DOCTYPE html>
<html lang="pt-br">
<head>
    <meta charset="UTF-8">
    <title>Método da Bisseção</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
</head>
<body class="p-4">
    <h2 class="mb-3">Método da Bisseção</h2>

    <?php if ($stage === 'menu'): ?>
        <form method="post" class="card card-body mb-4" autocomplete="off">
            <input type="hidden" name="stage" value="solve">
            <div class="mb-3">
                <label class="form-label">Função <span class="text-muted">(ex: <code>pow($x,3) - $x - 2</code>)</span></label>
                <input type="text" class="form-control" name="funcao" value="<?= htmlspecialchars($input_funcao) ?>" required>
            </div>
            <div class="row mb-3">
                <div class="col">
                    <label class="form-label">Início do intervalo (a)</label>
                    <input type="number" step="any" class="form-control" name="a" value="<?= htmlspecialchars($input_a) ?>" required>
                </div>
                <div class="col">
                    <label class="form-label">Fim do intervalo (b)</label>
                    <input type="number" step="any" class="form-control" name="b" value="<?= htmlspecialchars($input_b) ?>" required>
                </div>
            </div>
            <div class="row mb-3">
                <div class="col">
                    <label class="form-label">Tolerância</label>
                    <input type="number" step="any" class="form-control" name="tol" value="<?= htmlspecialchars($input_tol) ?>" required>
                </div>
                <div class="col">
                    <label class="form-label">Máximo de iterações</label>
                    <input type="number" class="form-control" name="maxIter" value="<?= htmlspecialchars($input_maxIter) ?>" required>
                </div>
            </div>
            <button type="submit" class="btn btn-success rounded-pill"><i class="bi bi-play-fill me-1"></i>Calcular</button>
            <a href="index.php" class="btn btn-secondary rounded-pill ms-2"><i class="bi bi-arrow-counterclockwise me-1"></i>Voltar</a>
        </form>
    <?php elseif ($stage === 'solve'): ?>
        <?php if (isset($resultado["erro"])): ?>
            <div class="alert alert-danger"><?= htmlspecialchars($resultado["erro"]) ?></div>
        <?php else: ?>
            <div class="alert alert-success">
                Raiz aproximada: <strong><?= $resultado["raiz"] ?></strong>
            </div>
            <table class="table table-striped">
                <thead>
                    <tr>
                        <th>Iteração</th>
                        <th>a</th>
                        <th>b</th>
                        <th>p</th>
                        <th>f(p)</th>
                        <th>Erro</th>
                    </tr>
                </thead>
                <tbody>
                    <?php foreach ($resultado["iteracoes"] as $linha): ?>
                        <tr>
                            <td><?= $linha["k"] ?></td>
                            <td><?= number_format($linha["a"], 6) ?></td>
                            <td><?= number_format($linha["b"], 6) ?></td>
                            <td><?= number_format($linha["p"], 6) ?></td>
                            <td><?= number_format($linha["f(p)"], 6) ?></td>
                            <td><?= number_format($linha["erro"], 6) ?></td>
                        </tr>
                    <?php endforeach; ?>
                </tbody>
            </table>
        <?php endif; ?>
        <form method="post" class="d-inline">
            <input type="hidden" name="stage" value="menu">
            <input type="hidden" name="funcao" value="<?= htmlspecialchars($input_funcao) ?>">
            <input type="hidden" name="a" value="<?= htmlspecialchars($input_a) ?>">
            <input type="hidden" name="b" value="<?= htmlspecialchars($input_b) ?>">
            <input type="hidden" name="tol" value="<?= htmlspecialchars($input_tol) ?>">
            <input type="hidden" name="maxIter" value="<?= htmlspecialchars($input_maxIter) ?>">
            <button type="submit" class="btn btn-outline-primary mt-3">Voltar</button>
        </form>
        <form method="post" class="d-inline">
            <input type="hidden" name="stage" value="menu">
            <button type="submit" class="btn btn-outline-secondary mt-3 ms-2">Reiniciar</button>
        </form>
        <a href="index.php" class="btn btn-outline-dark mt-3 ms-2">Voltar ao Menu Principal</a>
    <?php endif; ?>
</body>
</html>