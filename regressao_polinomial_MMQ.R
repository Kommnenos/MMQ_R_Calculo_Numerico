#### Regressao Polinomial ####
# matriz = um conjunto de dados bidimencional 2xN a ser aplicado o metodo de minimos quadrados
# Saidas
# A solucao

regressao_polinomial_MMQ = function(dados, grau){
  
  num_dados = ncol(dados);
  
  # Matriz de Vandermonde
  vandermonde = matrix(0, num_dados, grau + 1)
  for (i in 1:num_dados){
    for (j in 0:grau){
      vandermonde[i, j + 1] = dados[1, i]^j;
    }
  }
  
  transposta = t(vandermonde);
  normal_esquerda = transposta %*% vandermonde;
  normal_direita = transposta %*% dados[2,];
  
  matrix_aumentada = cbind(normal_esquerda, normal_direita);
  coeficientes = resolver_gaussiana(matrix_aumentada);
  
  return(coeficientes);
  
  
}