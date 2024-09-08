########### Eliminação Gaussiana ###########
# Entradas
# matriz = uma matriz a ser resolvida
# Saida
# Um vetor de solucao

resolver_gaussiana = function(matriz){
  tempo = proc.time();
  # tamanho da matriz
  numero_colunas = ncol(matriz);
  numero_linhas = nrow(matriz);
  
  # Tranforma em triangular
  for (i in 1:(numero_colunas - 1)){
    # Elimina elementos abaixo do pivo atual
    for (j in (i + 1):numero_linhas){
      if(numero_linhas < j){
        break;
        }
      multiplicador = matriz[j, i]/matriz[i,i];
      
      for (k in i:(numero_colunas)){
        matriz[j, k] = matriz[j, k] - multiplicador * matriz[i, k];  
      }
    }
    
  }
  
  # resolve as incognitas
  vetor_solucao = numeric(numero_linhas);
  
  for (l in numero_linhas:1){
    # resolve para cada variavel a partir da ultima linha
    soma = matriz[l, numero_linhas + 1];
    
    #subtracao de valores conhecidos
    for (m in (l + 1):numero_linhas) {
      if(m > numero_linhas){
        break;
      }
      soma = soma - matriz[l, m] * vetor_solucao[m];
    }
    
    vetor_solucao[l] = soma / matriz[l, l];
    
  }
  
  print(proc.time() - tempo); #retorna a contagem de tempo
  
  return(vetor_solucao);
}