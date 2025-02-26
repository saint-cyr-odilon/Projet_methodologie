# install.packages('dplyr')
# install.packages('ggplot2')
# install.packages('MASS')
# install.packages('GGally')
# install.packages('tidyverse')

library(dplyr)
library(ggplot2)
library(MASS)  # Pour générer les données corrélées
library(GGally)
library(tidyverse)


# Fonction pour générer une population simulée
generate_population <- function(pop_size, cor_U_Y, cor_U_Pi, cor_U_Pt) {
  
  cor_matrix <- matrix(c(
    1.0,  0.3,  0.2,  cor_U_Y, cor_U_Pi, cor_U_Pt,  # U
    0.3,  1.0,  0.4,  0.5,     0.3,     0.3,      # X0
    0.2,  0.4,  1.0,  0.3,     0.35,    0.25,     # X1
    cor_U_Y,  0.5,  0.3,  1.0,  0.5,     0.4,      # Y
    cor_U_Pi, 0.3,  0.35,  0.5,  1.0,     0.3,      # Pi
    cor_U_Pt, 0.3,  0.25,  0.4,  0.3,     1.0),     # Pt
    nrow = 6, byrow = TRUE)
  
  # Vérification et régularisation si nécessaire
  if (any(eigen(cor_matrix)$values <= 0)) {
    cor_matrix <- cor_matrix + diag(nrow(cor_matrix)) * 0.01
  }
  
  # Moyennes et variances des variables
  mu <- c(0, 5, 25, 22, 0, 0)  
  sigma <- c(1, 6, 3, 5, 1, 1)  # X0 a une variance plus élevée, X1 plus faible
  
  # Génération de la population complète
  population <- mvrnorm(n = pop_size, mu = mu, Sigma = cor_matrix * (sigma %o% sigma), empirical = TRUE)
  colnames(population) <- c("U", "X0", "X1", "Y", "Pi_raw", "Pt_raw")
  
  population <- as.data.frame(population) %>%
    mutate(
      id = seq_len(n()),  # Ajout d'un identifiant unique
      Y  = (0.6 * X0 + 0.2 * X1) + 12*U + rnorm(pop_size, mean = 10, sd = 5),  # Y influencé par X0, X1 et U
      Pi = plogis(Pi_raw / 4  + qlogis(0.4)),  # Centrage sur 0.4 + Réduction de variance
      Pt = plogis(Pt_raw / 4  + qlogis(0.6))   # Centrage sur 0.6 + Réduction de variance
    ) %>%
    select(id, U, X0, X1, Y, Pi, Pt)   # Garder l’ordre propre
  return(population)
}

# Fonction pour générer les histogrammes d'une population donnée
generate_histograms <- function(population) {
  histograms <- list(
    ggplot(population, aes(x = U)) +
      geom_histogram(bins=100, fill = "gray", alpha = 0.7, color = "black") +
      labs(title = "Distribution de U", x = "U", y = "Effectif"),
    ggplot(population, aes(x = X0)) +
      geom_histogram(bins=100, fill = "blue", alpha = 0.7, color = "black") +
      labs(title = "Distribution de X0", x = "X0", y = "Effectif"),
    ggplot(population, aes(x = X1)) +
      geom_histogram(bins=100, fill = "green", alpha = 0.7, color = "black") +
      labs(title = "Distribution de X1", x = "X1", y = "Effectif"),
    ggplot(population, aes(x = Y)) +
      geom_histogram(bins=100, fill = "purple", alpha = 0.7, color = "black") +
      labs(title = "Distribution de Y", x = "Y", y = "Effectif"),
    ggplot(population, aes(x = Pi)) +
      geom_histogram(bins=100, fill = "blue", alpha = 0.7, color = "black") +
      labs(title = "Distribution de Pi", x = "Pi", y = "Effectif"),
    ggplot(population, aes(x = Pt)) +
      geom_histogram(bins=100, fill = "red", alpha = 0.7, color = "black") +
      labs(title = "Distribution de Pt", x = "Pt", y = "Effectif")
  )
  
  return(histograms)
}

# Exemple d'utilisation de la fonction
population <- generate_population(10000, 0.8, 0.2, 0.2)
histograms <- generate_histograms(population)

# Histogrammes des distributions des variables
for (hist in histograms) {
  print(hist)
}

# Résumé des distributions et corrélations des variables
ggpairs(population %>% select(-id),
        upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.2, size = 1))) +
  labs(title = "Résumé des distributions et corrélations des variables")


# Fonction de simulation des échantillons (réponse déterministe uniquement)
generate_sample <- function(population, sample_size, coef_seuil_i = 0.45, coef_seuil_t = 0.65) {
  
  # Tirage d'un échantillon
  data_sample <- population %>% sample_n(sample_size)
  
  # Mécanisme de réponse déterministe
  data_sample <- data_sample %>%
    mutate(
      internet_response = Pi > coef_seuil_i,
      telephone_response = (Pt > coef_seuil_t) & !internet_response,
      repondant_final = internet_response | telephone_response,
      
      # Mettre X0 et Y à NA pour les non-répondants
      X1 = ifelse(repondant_final, X1, NA),
      Y = ifelse(repondant_final, Y, NA)
    )
  
  return(data_sample)
}

# Génération des 500 échantillons
simulations <- replicate(500,
                         generate_sample(population,
                                         1000,
                                         coef_seuil_i = 0.425, # Plus ces seuils sont élevés, plus le taux de non-réponse final est élevé
                                         coef_seuil_t = 0.6),
                         simplify = FALSE)

moyennes_taux <- bind_rows(simulations) %>% 
  summarise(
    interroge_i = n(),
    interroge_t = sum(!internet_response),
    repondants_i = sum(internet_response),
    repondants_t = sum(telephone_response)
  ) %>%
  mutate(
    repondants_tot = repondants_i + repondants_t,
    taux_rep_i = repondants_i / interroge_i,
    taux_rep_t = repondants_t / interroge_t,
    taux_rep_tot = repondants_tot / interroge_i
  ) %>%
  select(interroge_i, repondants_i, taux_rep_i,   
         interroge_t, repondants_t, taux_rep_t,   
         repondants_tot, taux_rep_tot)            

# Affichage
print(moyennes_taux)



