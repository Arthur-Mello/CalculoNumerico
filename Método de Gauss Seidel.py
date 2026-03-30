# Trabalho de Cálculo Numérico 
# Método Iterativo: Gauss-Seidel
# Nomes dos Alunos:
# 1. Arthur Mello
# 2. Arthur Luiz
# 3. Eduardo Bolotari
# 4. Gabriel Fiori

def criterio_de_sassenfeld(n, A):
    """
    Testa a convergência do sistema para o método de Gauss-Seidel 
    usando o critério de Sassenfeld.
    """
    beta = [0.0] * n
    
    for i in range(n):
        soma = 0.0
        for j in range(n):
            if i != j:
                # Se j < i, usamos o beta[j] recém-calculado (como no próprio Gauss-Seidel)
                if j < i:
                    soma += abs(A[i][j]) * beta[j]
                # Se j > i, usamos apenas o valor do coeficiente
                else:
                    soma += abs(A[i][j])
        
        # O elemento da diagonal principal não pode ser nulo
        if A[i][i] == 0:
            return False 
            
        beta[i] = soma / abs(A[i][i])
    
    # Se o maior valor de beta for menor que 1, a convergência é garantida
    if max(beta) < 1:
        return True
    else:
        return False

def gauss_seidel(n, A, b, x0, tol, max_iter=1000):
    """
    Resolve um sistema linear Ax = b utilizando o método iterativo de Gauss-Seidel.
    """
    # 1. Teste de convergência
    print("Realizando teste de convergência (Critério de Sassenfeld)...")
    if not criterio_de_sassenfeld(n, A):
        print("AVISO: O sistema não satisfaz o critério de Sassenfeld. A convergência não é garantida.\n")
    else:
        print("Teste OK! O critério de Sassenfeld foi satisfeito. A convergência é garantida.\n")

    # 2. Inicialização do vetor de incógnitas com o chute inicial
    # Cria uma cópia do chute inicial para não alterar a variável original
    x = [valor for valor in x0]
    
    # Variável para armazenar os valores da iteração anterior
    x_anterior = [0.0] * n

    # 3. Loop de iterações
    for k in range(max_iter):
        # Atualiza o vetor anterior antes de calcular o novo
        for i in range(n):
            x_anterior[i] = x[i]
            
        # Calcula as novas aproximações para cada incógnita
        for i in range(n):
            soma1 = 0.0
            soma2 = 0.0
            
            # Somatório usando os valores de x JÁ atualizados nesta iteração (j < i)
            for j in range(i):
                soma1 += A[i][j] * x[j]
                
            # Somatório usando os valores de x da iteração anterior (j > i)
            for j in range(i + 1, n):
                soma2 += A[i][j] * x_anterior[j]
                
            # Atualiza o x[i] atual
            x[i] = (b[i] - soma1 - soma2) / A[i][i]
            
        # 4. Cálculo do Erro (Erro Relativo)
        # Erro = max(|x_atual - x_anterior|) / max(|x_atual|)
        diferencas = [abs(x[i] - x_anterior[i]) for i in range(n)]
        max_diferenca = max(diferencas)
        
        valores_absolutos_x = [abs(valor) for valor in x]
        max_x = max(valores_absolutos_x)
        
        # Evita divisão por zero se a raiz for exatamente [0, 0, ..., 0]
        if max_x == 0:
            erro = 0
        else:
            erro = max_diferenca / max_x
            
        print(f"Iteração {k+1:02d} | x = {[round(v, 5) for v in x]} | Erro = {erro:.5f}")
        
        # 5. Verifica a condição de parada (precisão atingida)
        if erro < tol:
            print(f"\nConvergência atingida na iteração {k+1}!")
            return x

    # Caso atinja o limite máximo de iterações sem convergir
    print("\nLimite máximo de iterações atingido sem alcançar a precisão desejada.")
    return x

# Bloco Principal para Testar o Algoritmo

if __name__ == "__main__":
    # Dados de entrada
    n = 3
    A = [
        [5, 1, 1],
        [3, 4, 1],
        [3, 3, 6]
    ]
    b = [5, 6, 0]
    
    # Chute inicial (vetor nulo, por exemplo)
    chute_inicial = [0, 0, 0]
    
    # Precisão (Tolerância)
    precisao = 0.05
    
    # Chamada da função
    print("=== RESOLUÇÃO POR GAUSS-SEIDEL ===")
    vetor_x = gauss_seidel(n, A, b, chute_inicial, precisao)
    
    # Saída
    print("\nVetor de incógnitas final (x):")
    for i, val in enumerate(vetor_x):
        print(f"x_{i+1} = {val:.5f}")