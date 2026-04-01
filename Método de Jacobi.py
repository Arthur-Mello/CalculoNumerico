# Trabalho de Cálculo Numérico
# Método Iterativo: Jacobi
# Nomes dos Alunos:
# 1. Arthur Mello
# 2. Arthur Luiz
# 3. Eduardo Bolotari   
# 4. Gabriel Fiori

def criterio_das_linhas(n, A):
    """
    Testa a convergência do sistema para o método de Jacobi 
    verificando se a matriz é estritamente diagonal dominante.
    """
    for i in range(n):
        soma_linha = 0.0
        for j in range(n):
            if i != j:
                soma_linha += abs(A[i][j])
        
        # O elemento da diagonal principal não pode ser nulo
        if A[i][i] == 0:
            return False
            
        # Se a soma dos elementos não-diagonais for maior ou igual
        # ao elemento da diagonal principal, o critério falha
        if soma_linha >= abs(A[i][i]):
            return False
            
    # Se passar por todas as linhas, a convergência é garantida
    return True

def jacobi(n, A, b, x0, tol, max_iter=1000):
    """
    Resolve um sistema linear Ax = b utilizando o método iterativo de Jacobi.
    """
    # 1. Teste de convergência
    print("Realizando teste de convergência (Critério das Linhas)...")
    if not criterio_das_linhas(n, A):
        print("AVISO: A matriz não é estritamente diagonal dominante. A convergência não é garantida.\n")
    else:
        print("Teste OK! O critério das linhas foi satisfeito. A convergência é garantida.\n")

    # 2. Inicialização dos vetores
    # x recebe os valores da iteração atual (inicia com o chute)
    x = [valor for valor in x0]
    # x_novo guardará os resultados calculados para a próxima iteração
    x_novo = [0.0] * n

    # 3. Loop de iterações
    for k in range(max_iter):
        
        # Calcula as novas aproximações para cada incógnita
        for i in range(n):
            soma = 0.0
            
            # Somatório usando APENAS os valores de x da iteração anterior
            for j in range(n):
                if j != i:
                    soma += A[i][j] * x[j]
            
            # Atualiza o valor no vetor auxiliar x_novo
            x_novo[i] = (b[i] - soma) / A[i][i]
            
        # 4. Cálculo do Erro (Erro Relativo)
        # Erro = max(|x_novo - x|) / max(|x_novo|)
        diferencas = [abs(x_novo[i] - x[i]) for i in range(n)]
        max_diferenca = max(diferencas)
        
        valores_absolutos_x_novo = [abs(valor) for valor in x_novo]
        max_x_novo = max(valores_absolutos_x_novo)
        
        # Evita divisão por zero
        if max_x_novo == 0:
            erro = 0
        else:
            erro = max_diferenca / max_x_novo
            
        print(f"Iteração {k+1:02d} | x = {[round(v, 5) for v in x_novo]} | Erro = {erro:.5f}")
        
        # 5. Atualiza o vetor x para a próxima iteração
        for i in range(n):
            x[i] = x_novo[i]
            
        # 6. Verifica a condição de parada (precisão atingida)
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
        [10, 2, 1],
        [1, 5, 1],
        [2, 3, 10]
    ]
    b = [7, -8, 6]
    
    # Chute inicial (vetor nulo)
    chute_inicial = [0, 0, 0]
    
    # Precisão (Tolerância)
    precisao = 0.05
    
    # Chamada da função
    print("=== RESOLUÇÃO PELO MÉTODO DE JACOBI ===")
    vetor_x = jacobi(n, A, b, chute_inicial, precisao)
    
    # Saída
    print("\nVetor de incógnitas final (x):")
    for i, val in enumerate(vetor_x):
        print(f"x_{i+1} = {val:.5f}")