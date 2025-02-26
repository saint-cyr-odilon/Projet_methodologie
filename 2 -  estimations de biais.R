### Méthode 1: aucune reponderation

# Fonction pour calculer le biais sur Y 
evaluate_raw_bias_Y <- function(sample, population) {
  
  # Moyenne de Y dans la population complète
  mean_Y_pop <- mean(population$Y, na.rm = TRUE)
  
  # Moyenne de Y dans l'échantillon (chez les répondants)
  mean_Y_sample <- mean(sample$Y, na.rm = TRUE)
  
  # Calcul du biais
  biais_Y <- mean_Y_sample - mean_Y_pop
  
  # Retourner un tableau avec les résultats
  return(data.frame(
    Methode = "Sans correction",
    Moyenne_Y = mean_Y_sample,
    Biais = biais_Y
  ))
}

# Méthode 2: calcul du biais avec reponderation basée sur proba de réponse estimée par modelisation
# Probabilité de réponse modelisée par X0, X1 et Y pour les repondants et par X1 pour les non repondants

evaluate_weighted_bias_Y_global <- function(sample, population) {
  
  mean_Y_pop <- mean(population$Y, na.rm = TRUE)
  
  # Vérifier qu'il y a suffisamment de répondants
  if (sum(sample$repondant_final) == 0) {
    stop("Erreur: Pas assez de répondants.")
  }
  
  # Modèle pour estimer la probabilité de réponse
  model_repondants <- glm(repondant_final ~ X0 + X1 + Y, 
                          data = sample %>% filter(repondant_final), 
                          family = binomial, 
                          control = list(maxit = 500))
  
  model_non_repondants <- glm(repondant_final ~ X1, 
                              data = sample, 
                              family = binomial, 
                              control = list(maxit = 500))
  
  # Prédiction des probabilités de réponse
  sample <- sample %>%
    mutate(
      prob_response_repondants = ifelse(repondant_final, 
                                        predict(model_repondants, newdata = sample, type = "response"), NA_real_),
      prob_response_non_repondants = predict(model_non_repondants, newdata = sample, type = "response"),
      
      # Attribution des probabilités selon le statut de réponse
      prob_response = ifelse(!is.na(prob_response_repondants), prob_response_repondants, prob_response_non_repondants),
      
      # Éviter les valeurs extrêmes
      prob_response = pmax(pmin(prob_response, 0.99), 0.01),
      
      # Calcul des poids IPW
      poids_IPW = 1 / prob_response
    )
  
  # Correction du biais avec pondération IPW
  mean_Y_weighted <- weighted.mean(sample$Y, sample$poids_IPW, na.rm = TRUE)
  
  return(data.frame(
    Methode = "IPW global",
    Moyenne_Y = mean_Y_weighted,
    Biais = mean_Y_weighted - mean_Y_pop
  ))
}

# Exemple d'utilisation :
sample <- generate_sample(population, 10000)
biais_raw <- evaluate_raw_bias_Y(sample, population)
biais_weighted_global <- evaluate_weighted_bias_Y_global(sample, population)


# Affichage du biais brut sans correction
print(rbind(biais_raw,
            biais_weighted_global))
