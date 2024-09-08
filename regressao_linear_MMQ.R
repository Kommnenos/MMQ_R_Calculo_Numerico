##### Regressao Linear ####
# Entradas
# matriz = um conjunto de dados bidimencional 2xN a ser aplicado o metodo de minimos quadrados
# Saidas
# A solucao

regressao_linear_MMQ = function(dados){
  ## verificaçoes de segurança
  if (nrow(dados) != 2 || is.matrix(dados) == FALSE){
    return("Erro, certifique-se que seu dataset seja uma matriz 2xN.");
  }
  num_dados = ncol(dados);
  
  ## somatorios dos dados
  somatorio_x = sum(dados[1,]);
  somatorio_y = sum(dados[2,]);
  somatorio_x_quadrado = sum(dados[1,]^2);
  somatorio_xy = sum(dados[1,] * dados[2,]);
  
  #criação do sistema linear a ser resolvido
  sistema_linear = matrix(c(num_dados, somatorio_x, somatorio_x, somatorio_x_quadrado, somatorio_y, somatorio_xy), 2, 3);
  
  #resolução do sistema
  
  ajuste_linear = resolver_gaussiana(sistema_linear);
  
  return(list(intersecao = ajuste_linear[1], inclinacao = ajuste_linear[2]));
  
  
  
  
}