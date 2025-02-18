# Installation des packages requis
install.packages("sampleSelection", dependencies = TRUE)
install.packages("mice", dependencies = TRUE)
library(sampleSelection)  # Pour le modèle de Heckman
library(mice)             # Pour l'imputation multiple
library(ggplot2)          # Pour la visualisation

# Fixer la graine pour la reproductibilité
set.seed(123)

# Taille de la population simulée
taille_population <- 10000

# Fonction pour générer des données simulées
generate_data <- function(n, influence_U) {
  # Variables auxiliaires (connues pour tous)
  X0 <- rnorm(n, mean = 50, sd = 10)  # Variable auxiliaire connue
  U <- rnorm(n)                       # Variable latente non mesurée
  
  # Probabilités de réponse modifiées par l'influence de U
  Pi <- plogis(-0.1 * X0 + influence_U * U)  # Proba de répondre par Internet
  Pt <- plogis(0.1 * X0 - influence_U * U)  # Proba de répondre par téléphone
  
  # Déterminer la réponse
  respond <- ifelse(runif(n) < Pi, 1, ifelse(runif(n) < Pt, 2, 0))
  
  # Générer X1 uniquement pour les répondants
  X1 <- rep(NA, n)  # Initialisation avec NA
  X1[respond > 0] <- rnorm(sum(respond > 0), mean = 5, sd = 2)
  
  # Générer Y uniquement pour les répondants
  Y <- rep(NA, n)  # Initialisation avec NA
  Y[respond > 0] <- 0.5 * X0[respond > 0] + 2 * X1[respond > 0] + U[respond > 0] + rnorm(sum(respond > 0))
  
  # Retourner les données dans un data.frame
  data.frame(ID = 1:n, X0, X1, Y, U, Pi, Pt, respond)
}

# Fonction pour appliquer et évaluer les méthodes de redressement
test_methods <- function(data) {
  library(dplyr)
  
  # Séparation des répondants et non-répondants
  respondents <- data %>% filter(respond > 0)
  non_respondents <- data %>% filter(respond == 0)
  
  # Vraie moyenne de Y (pour comparaison)
  true_mean_Y <- mean(data$Y, na.rm = TRUE)
  
  # 1. Pondération inverse de la probabilité (IPW)
  weights <- 1 / data$Pi
  weights[data$respond == 0] <- NA  # Poids non définis pour les non-répondants
  ipw_mean <- sum(respondents$Y * weights[respondents$ID], na.rm = TRUE) / 
    sum(weights[respondents$ID], na.rm = TRUE)
  
  # 2. Imputation moyenne conditionnelle
  imputed_Y <- data$Y
  imputed_Y[is.na(imputed_Y)] <- mean(respondents$Y, na.rm = TRUE)
  mean_conditional <- mean(imputed_Y)
  
  # 3. Imputation multiple
  mice_data <- mice(data[, c("X0", "X1", "Y")], m = 5, method = c("norm.predict", "norm.predict", "norm"), maxit = 5, seed = 123)
  imputed <- complete(mice_data, action = "long")
  multiple_imputation_mean <- mean(imputed$Y, na.rm = TRUE)
  
  # 4. Modèle de Heckman
  heckman_model <- selection(
    selection = respond ~ X0,
    outcome = Y ~ X0 + X1,
    data = data %>% filter(!is.na(Y)),
    method = "ml"
  )
  heckman_mean <- mean(predict(heckman_model, part = "outcome", type = "conditional"), na.rm = TRUE)
  
  # Résultats
  list(
    true_mean_Y = true_mean_Y,
    ipw_mean = ipw_mean,
    mean_conditional = mean_conditional,
    multiple_imputation_mean = multiple_imputation_mean,
    heckman_mean = heckman_mean
  )
}

# Fonction pour simuler plusieurs itérations et calculer les statistiques
simulate_iterations <- function(n_iterations, n_population, influence_U) {
  results <- replicate(n_iterations, {
    data <- generate_data(n_population, influence_U)
    test_methods(data)
  }, simplify = FALSE)
  
  # Calcul des biais et variances
  results_df <- do.call(rbind, lapply(results, function(res) {
    data.frame(
      true_mean_Y = res$true_mean_Y,
      ipw_bias = abs(res$ipw_mean - res$true_mean_Y),
      mean_conditional_bias = abs(res$mean_conditional - res$true_mean_Y),
      multiple_imputation_bias = abs(res$multiple_imputation_mean - res$true_mean_Y),
      heckman_bias = abs(res$heckman_mean - res$true_mean_Y)
    )
  }))
  
  # Moyenne et écart-type des biais
  summary_stats <- data.frame(
    method = c("IPW", "Mean Conditional", "Multiple Imputation", "Heckman"),
    mean_bias = c(mean(results_df$ipw_bias), mean(results_df$mean_conditional_bias),
                  mean(results_df$multiple_imputation_bias), mean(results_df$heckman_bias)),
    sd_bias = c(sd(results_df$ipw_bias), sd(results_df$mean_conditional_bias),
                sd(results_df$multiple_imputation_bias), sd(results_df$heckman_bias))
  )
  
  return(list(summary = summary_stats, detailed = results_df))
}

# Simulation pour chaque niveau d'influence de U
influences <- c(0.1, 0.3, 0.6)  # Faible, modérée, forte
n_iterations <- 100  # Nombre d'itérations
n_population <- 10000

# Résultats pour chaque influence de U
all_results <- lapply(influences, function(inf) {
  cat("Simulation en cours pour influence U:", inf, "\n")
  simulate_iterations(n_iterations, n_population, inf)
})

# Visualisation des biais par méthode
plot_biases <- function(all_results, influences) {
  bias_data <- do.call(rbind, lapply(1:length(influences), function(i) {
    result <- all_results[[i]]$detailed
    data.frame(
      influence_U = influences[i],
      method = rep(c("IPW", "Mean Conditional", "Multiple Imputation", "Heckman"), each = nrow(result)),
      bias = c(result$ipw_bias, result$mean_conditional_bias, result$multiple_imputation_bias, result$heckman_bias)
    )
  }))
  
  ggplot(bias_data, aes(x = factor(influence_U), y = bias, fill = method)) +
    geom_boxplot() +
    labs(title = "Distribution des biais par méthode et influence de U", x = "Influence de U", y = "Biais") +
    theme_minimal()
}

# Afficher le graphique des biais
plot_biases(all_results, influences)

# Analyse par sous-groupes
analyze_subgroups <- function(influence_U, subgroup_variable, subgroup_values, n_iterations, n_population) {
  subgroup_results <- lapply(subgroup_values, function(value) {
    results <- replicate(n_iterations, {
      data <- generate_data(n_population, influence_U)
      data <- data[data[[subgroup_variable]] >= value, ]  # Filtrer par sous-groupe
      test_methods(data)
    }, simplify = FALSE)
    
    # Calculer la moyenne des biais pour ce sous-groupe
    do.call(rbind, lapply(results, function(res) {
      data.frame(
        ipw_bias = abs(res$ipw_mean - res$true_mean_Y),
        mean_conditional_bias = abs(res$mean_conditional - res$true_mean_Y),
        multiple_imputation_bias = abs(res$multiple_imputation_mean - res$true_mean_Y),
        heckman_bias = abs(res$heckman_mean - res$true_mean_Y)
      )
    }))
  })
  
  names(subgroup_results) <- paste0("Subgroup_", subgroup_values)
  return(subgroup_results)
}

# Exemple d'analyse par sous-groupe
subgroup_results <- analyze_subgroups(0.3, "X0", c(40, 50, 60), 50, 10000)
