#### Regressao Nao Linear - Capacitor###
# Aplica uma regressao nao linear baseada na equação da descarga do capacitor
# Entrada
# matrix 2xN, com a primeira linha sendo voltagem e a segunda tempo.

regressao_nao_linear_MMQ = function(dados){
  tolerancia = 0.0000001
  iteracoes = 250
  num_dados = ncols(dados);
  
  #chutes iniciais
  voltagem_inicial = 4
  constante_tempo = 1
  
  
  
  #calculo residuais
  for (i in 1:iteracoes){
    matriz_jacob = matrix(0, num_dados, 2)
    residuais = numeric(num_dados)
    
    #calculo residual e jacobiana
    for (j in 1:num_dados){
      tempo = dados[1, j];
      voltagem = dados[2, j];
      
      #predicao voltagem
      voltagem_predita = voltagem_inicial * exp(-tempo/constante_tempo);
      
      # calculo residuo
      residuais[j] = voltagem - voltagem_predita;
    
      
      #derivadas parciais
      der_voltagem_inicial = -exp(-tempo / constante_tempo)
      
      der_constante_tempo = voltagem_inicial * tempo * exp(-tempo/constante_tempo) / (constante_tempo^2)
      
      #construcao da matriz jacobiana
      
      matriz_jacob[j, 1] = der_voltagem_inicial;
      matriz_jacob[j, 2] = der_constante_tempo;
    }
    # Atualização dos parâmetros
    Jt_J = t(matriz_jacob) %*% matriz_jacob
    Jt_r = t(matriz_jacob) %*% residuais
    delta_p = solve(Jt_J) %*% Jt_r
    
    # Atualizar parâmetros
    voltagem_inicial = voltagem_inicial + delta_p[1]
    constante_tempo = constante_tempo + delta_p[2]
    
    # Verificação de convergência
    if (sqrt(sum(delta_p^2)) < tolerancia) {
      break
    }
  }
  
  return(c(voltagem_inicial, constante_tempo))
}