#' Parse path
#'
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
parse_path <- function(path) {
  pp <- str_split(path, "_", simplify = TRUE)
  tr <- paste(pp[2], pp[3], sep = "_")
  sex <- if_else(pp[1] == "Female", "F", "M")
  trait <- case_when(
    tr == "nasion_basion" ~ "n_b",
    tr == "nasion_menton" ~ "n_m",
    tr == "sella_gonion" ~ "s_g",
    tr == "sella_nasion" ~ "s_n",
    tr == "sella_basion" ~ "s_b",
    tr == "nasion_ans" ~ "n_ans",
    tr == "menton_ans" ~ "ans_m",
    tr == "nsb_angle" ~ "nsb_angle",
    tr == "articulare_pogonion" ~ "art_po",
    tr == "mp_angle" ~ "mpa",
    tr == "gonial_angle" ~ "gonial_angle",
    tr == "A_NB" ~ "ANB",
    tr == "ans_pns" ~ "ans_pns",
    tr == "condylion_pogonion" ~ "c_po",
    tr == "condylion_gonion" ~ "c_go",
    tr == "gonion_pogonion" ~ "g_po"
  )
  return(list(trait = trait, sex = sex))
}


#' Fix capitalization for datasets
#'
#' @param s (string): Stored name for growth study set
#'
#' @return
#' @export
#'
#' @examples
cap_set <- function(s) {
  s <- case_when(
    s == "bolton" ~ "Bolton Brush",
    s == "denver" ~ "Denver",
    s == "fels" ~ "Fels",
    s == "iowa" ~ "Iowa",
    s == "michigan" ~ "Michigan",
    s == "oregon" ~ "Oregon"
  )
  return(s)
}


#' Load data for CF growth model
#'
#' @param variable (string): variable to select, will be renamed to y
#' @param which_sex (string): for sex ("F", "M")
#' @param lower_age (numeric): lower age of range to use for model fitting 
#' @param upper_age (numeric): upper age of range to use for model fitting
#' @param drop_1 
#'
#' @return
#' @export
#'
#' @examples
load_CF_data <- function(variable,
                         which_sex,
                         lower_age = 0,
                         upper_age = 100,
                         drop_1 = FALSE) {
  M <- data.table::fread("../Data/master_0805_fixedDENVER_clean.csv") %>% 
    filter(sex == which_sex) %>% 
    filter(age.yrs > lower_age, age.yrs < upper_age) %>% 
    dplyr::select(age.yrs, individual, set, matches(variable)) %>% 
    rename(age = age.yrs,
           Set = set) %>% 
    drop_na()
  
  # Rename variable to y
  names(M)[names(M) == variable] <- "y"
  
  if (drop_1) { # Drop those with 1 obs
    Obs_1 <- M %>% 
      group_by(individual) %>% 
      tally() %>% 
      filter(n == 1L)
    
    M <- M %>% 
      filter(!(individual %in% Obs_1$individual))
  }
  
  # Drop any y < 1 to get rid of 0 values
  M <- M %>% 
    filter(y > 1)
  
  M <- M %>% 
    mutate(individual = factor(individual),
           individual_int = as.integer(individual),
           Set = factor(cap_set(Set)),
           set_int = as.integer(Set)) %>% 
    as.data.frame()
  
  dat <- data.frame(age = M$age,
                    y = M$y,
                    ID = M$individual_int,
                    set_int = M$set_int)
  
  return(list(M = M, dat = dat))
}

# Color scale

cols <- c("Bolton Brush" = "#E69F00",
          "Bolton\nBrush" = "#E69F00",
          "Denver" = "#56B4E9",
          "Iowa" = "#CC79A7",
          "Fels" = "#F0E442",
          "Michigan" = "#0072B2",
          "Oregon" = "#009E73",
          "Population" = "#000000")

clr_scale <- list(
  scale_color_manual(
    values = cols),
  scale_shape_manual(values = c(16, 1),
                     guide = guide_legend(nrow = 2)),
  theme(axis.text = element_text(size = 9),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold")),
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3),
                               reverse = TRUE,
                               nrow = 3),
         shape = guide_legend(override.aes = list(alpha = 1, size = 3))))


#' Calculate first derivative for different models
#'
#' @param P Fitted map2stan model
#'
#' @return tibble with first derivatives
#' @export
#'
#' @examples
derivative <- function(P) {
  out <- tibble(age = numeric(),
                Set = character(),
                `Preece-Baines` = numeric(),
                `Double Logistic` = numeric())
  
  for (ii in unique(P$Set)) {
    Ps <- P %>% filter(Set == ii) %>% dplyr::select(-Set)
    d_age <- Ps %>% pull(age)
    d_age <- d_age[-1]
    
    Ps <- Ps %>% dplyr::select(-age)
    
    Pd <- apply(Ps, 2, diff) %>%
      as_tibble()
    Pd <- bind_cols(age = d_age,
                    Set = rep(ii, length.out = length(d_age)),
                    Pd)
    out <- rbind(out, Pd)
  }
  return(out)
}

#' Find nearest value in a vector
#'
#' @param vec Vector of values to search
#' @param target Target value to find in `vec`
#'
#' @return The index of `vec` nearest to `target`
#' @export
#'
#' @examples
nearest_value <- function(vec, target) {
  which(abs(vec - target) == min(abs(vec - target)))
}

#' Convert precis output to data.frame 
#'
#' @param pp A `rethinking::precis()` object
#'
#' @return A `data.frame`
#' @export
#'
#' @examples
precis_to_df <- function(pp) {
  D <- data.frame(pp@.Data)
  names(D) <- pp@names
  row.names(D) <- pp@row.names
  return(D)
}


#' Postcheck individual case
#'
#' @param fit A `map2stan` fitted object
#' @param case Numeric value for observation
#' @param prob Probability for the credible interval
#'
#' @export
postcheck_case <- function(d, fit, case, k, prob = 0.90) {
  library(ggplot2, quietly = TRUE)
  library(cowplot, quietly = TRUE)
  theme_set(theme_cowplot())
  
  # Informationa about the participant
  pt_ID <- paste("ID:", d$individual[case])
  pt_set <- paste("Set:", d$Set[case])
  pt_age <- paste("Age:", round(d$age[case], 2))
  pt_y <- paste("Measure:", round(d$y[case], 2))
  
  newdata <- data.frame(
    age = fit@data$age[case],
    ID = fit@data$ID[case],
    set_int = fit@data$set_int[case]
  )

  ID_data <- tibble(
    age = fit@data$age[fit@data$ID == fit@data$ID[case]],
    y = fit@data$y[fit@data$ID == fit@data$ID[case]]
  )
  
  sims <- rethinking::sim(fit = fit, n = 1000, data = newdata,
                          refresh = 0)
  mu_mean <- mean(sims)
  mu_PI <- PI(sims, prob = prob)
  
  p <- ggplot() +
    geom_line(data = data.frame(age = fit@data$age,
                                ID = fit@data$ID,
                                y = fit@data$y),
              aes(x = age, y = y, group = ID),
              alpha = 0.05) +
    geom_line(data = ID_data, aes(x = age, y = y),
              color = "red", size = 1) +
    geom_point(data = ID_data, aes(x = age, y = y),
               color = "red", size = 2) +
    geom_point(data = data.frame(x = fit@data$age[case], y = mu_mean),
               aes(x, y),
               color = "blue", size = 2) +
    geom_linerange(data = data.frame(x = fit@data$age[case],
                                     ymin = as.numeric(mu_PI[1]),
                                     ymax = as.numeric(mu_PI[2])),
                   aes(x = x, ymin = ymin, ymax = ymax),
                   color = "blue", size = 1) +
    labs(x = "Age (y)", y = "Distance (mm)",
         title = paste(paste(pt_ID, pt_set, sep = "; "),
                       paste(pt_age, pt_y, sep = "; "),
                       paste("k =", round(k, 2)), sep = "\n"))
  return(p)
}

pagebreak <- function() {
  # https://stackoverflow.com/a/55064070
  if(knitr::is_latex_output())
    return("\\newpage")
  else
    return('<div style="page-break-before: always;" />')
}