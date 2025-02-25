library(dplyr)
library(ggplot2)
library(MASS)  # Pour générer les données corrélées

# Fonction de simulation avec contrôle des corrélations
simulation <- function(pop_size, sample_size, cor_U_Y, cor_U_Pi, cor_U_Pt) {
  
  # Définition de la matrice de corrélation
  cor_matrix <- matrix(c(1.0,  0.0,  0.0,  cor_U_Y, cor_U_Pi, cor_U_Pt,  # U
                         0.0,  1.0,  0.0,  0.0,     0.0,     0.0,      # X0
                         0.0,  0.0,  1.0,  0.0,     0.0,     0.0,      # X1
                         cor_U_Y,  0.0,  0.0,  1.0,  0.0,     0.0,      # Y
                         cor_U_Pi, 0.0,  0.0,  0.0,  1.0,     0.2,      # Pi
                         cor_U_Pt, 0.0,  0.0,  0.0,  0.2,     1.0),     # Pt
                       nrow = 6, byrow = TRUE)
  
  # Moyennes et variances des variables
  mu <- c(0, 0, 0, 0, 0, 0)
  sigma <- c(1, 1, 1, 1, 1, 1)  # Variances unitaires
  
  # Génération des données corrélées
  data_gen <- mvrnorm(n = pop_size, mu = mu, Sigma = cor_matrix * (sigma %o% sigma), empirical = TRUE)
  colnames(data_gen) <- c("U", "X0", "X1", "Y", "Pi_raw", "Pt_raw")
  
  # Stockage dans un data.frame
  data_pop <- as.data.frame(data_gen) %>%
    mutate(
      id = seq_len(n()),  # Ajout d'un identifiant unique
      Y  = 0.5 * X0 + 0.3 * X1 + Y + rnorm(pop_size, mean = 0, sd = 1),  # Y influencé par X0, X1 et U
      Pi = plogis(Pi_raw),  # Transformation en probabilité entre 0 et 1
      Pt = plogis(Pt_raw)   # Transformation en probabilité entre 0 et 1
    ) %>%
    select(id, U, X0, X1, Y, Pi, Pt)  # Garder l’ordre propre
  
  # Tirage d'un échantillon initial
  data_sample <- data_pop %>% sample_n(sample_size)
  
  # Mécanisme de réponse aléatoire
  bilan_aleatoire <- data_sample %>%
    mutate(internet_response = runif(n()) < Pi,
           telephone_response = (runif(n()) < Pt) & !internet_response,
           response = case_when(
             internet_response ~ "1-Internet (Aléatoire)",
             telephone_response ~ "2-Téléphone (Aléatoire)",
             TRUE ~ "3-Non-répondant (Aléatoire)")) %>%
    count(response) %>%
    rename(Mode = response, Effectif = n) %>%
    mutate(Taux = Effectif / sample_size)
  
  # Mécanisme de réponse déterministe avec seuil fixé
  bilan_deterministe <- data_sample %>%
    mutate(internet_response = mean(Pi) * 1.15 < Pi,
           telephone_response = (mean(Pt) * 0.85 < Pt) & !internet_response,
           response = case_when(
             internet_response ~ "1-Internet (Déterministe)",
             telephone_response ~ "2-Téléphone (Déterministe)",
             TRUE ~ "3-Non-répondant (Déterministe)")) %>%
    count(response) %>%
    rename(Mode = response, Effectif = n) %>%
    mutate(Taux = Effectif / sample_size)
  
  return(list(bilan_aleatoire = bilan_aleatoire, bilan_deterministe = bilan_deterministe))
}

# Exécution de 10 simulations et calcul des moyennes des taux de réponse
simulations <- replicate(10, simulation(100000, 10000, cor_U_Y = 0.7, cor_U_Pi = 0.5, cor_U_Pt = 0.3), simplify = FALSE)

moyennes_taux <- bind_rows(
  lapply(simulations, function(sim) bind_rows(sim$bilan_aleatoire, sim$bilan_deterministe))
) %>%
  group_by(Mode) %>%
  summarise(Moyenne_Taux = mean(Taux))

# Affichage des résultats
print(moyennes_taux)
