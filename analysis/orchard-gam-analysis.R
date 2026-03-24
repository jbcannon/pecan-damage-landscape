library(tidyverse)
library(mgcv)
library(gratia)
library(ggpubr)
library(magrittr)

# Load clean and prepare data
df <- read.csv("../data/pecan-orchard-data.csv")

df$orchard = as.factor(df$orchard)
df = na.omit(df)
#clamp damage to (0,1) exclusive
clamp = function(x, min_val=-Inf, max_val=Inf)  pmax(pmin(x, max_val), min_val)
df$crownDamage = clamp(df$damage/100, 0.001, 0.999)

df$rel_dens = df$rel_dens_1.8ha/ (20^2*pi/10000) # correct density

# Check for multicollinearity especially among size variables
all_vars = c('ht_m', 'predicted_dbh', 'crown_area_m2', 'rel_ht',
         'rel_dens', 'HELENE_MSW')
cor_mat <- cor(df[,all_vars])
diag(cor_mat) = NA
cor_mat[which.max(cor_mat)]
corrplot::corrplot(cor_mat)

# Height, crown area and predcited dbh are are highly correlated
cor(df$ht_m, df$crown_area_m2) #ht_m and crown_area highly correlated (r = 0.86)
cor(df$predicted_dbh, df$crown_area_m2) #dbh and crown_area highly correlated (r = .97)

# Drop these two (see AIC, analysis comparison belwo)
vars = setdiff(all_vars, c('ht_m','predicted_dbh', 'HELENE_MSW'))
df[, vars] %>% cor() %>% {diag(.) <- NA; .} %>% range(na.rm=TRUE)
vars

# Setup variables and terms for GAM
main_effects = c('crown_area_m2', 'rel_dens', 'rel_ht')
random_vars = 'orchard'
main_terms = paste0('s(', main_effects, ", bs='cs')", collapse=' + ')
random_terms = paste0("s(", random_vars, ", bs = 're')", collapse = ' + ')

#Construct final model
f = paste0("crownDamage ~ ", main_terms, ' + ', random_terms)
mod_beta = gam(as.formula(f), data=df, method='REML', family = betar())
summary(mod_beta)

gratia::draw(mod_beta)

# Output gam outputs to csv.
mod = mod_beta
gam_output.s = summary(mod)$s.table
gam_output.s = data.frame(variable = rownames(gam_output.s), gam_output.s, row.names = NULL)
gam_output.s$variable <- gsub("s\\(|\\)", "", gam_output.s$variable)  # remove 's(' and ')'
gam_output.s$variable <- ifelse(grepl("^ti\\(", gam_output.s$variable), 
                                gsub(",", " x ", gsub("^ti\\(|\\)$", "", gam_output.s$variable)),gam_output.s$variable)
gam_output.s = select(gam_output.s, !Ref.df)
gam_output = gam_output.s
write_csv(gam_output, '../figs/gam_output-orchard.csv', col_names = FALSE)
knitr::kable(gam_output.s, digits = c(0,3,3,4))

var_names = c(expression('crown area (m'^2*')'),
              'relative height', 
              expression('density (ha'^{-1}*')'))

#plot interaction effects..
fit_data <- lapply(seq_along(vars), function(i) {
  v <- vars[i]
  term_label = grep(v, rownames(summary(mod)$s.table), value=TRUE)
  term_label = Filter(function(x) {!grepl("^ti", x)}, term_label)
  pval = summary(mod)$s.table[term_label, 'p-value']
  fitdata = smooth_estimates(mod, select = term_label)
  fitdata <- fitdata %>%
      mutate(
        eta_lo = .estimate - 2*.se,
        eta_hi = .estimate + 2*.se,
        mu     = plogis(.estimate),
        mu_lo  = plogis(eta_lo),
        mu_hi  = plogis(eta_hi),
        variable = v,
        value = get(v)
      )
  fitdata = select(fitdata, mu:value)
  dens <- density(mod$model[[v]])
  dens <- density(mod$model[[strsplit(v, ":")[[1]][1]]])
  return(list(fitdata=fitdata, dens=dens,pval=pval))
})

all_plots = lapply(fit_data, function(fit) {
  pval = fit$pval
  dens = fit$dens
  fitdata = fit$fitdata
  v = fitdata$variable[1]
  sig <- pval < 0.05
  dens <- density(mod$model[[v]])
  dens <- density(mod$model[[strsplit(v, ":")[[1]][1]]])
  smooth_range <- range(with(fitdata, c(mu_lo, mu_hi)))
  smooth_range = c(0,0.75)
  y_scaled <- dens$y / max(dens$y) * diff(smooth_range) + smooth_range[1]
  hist_scaled <- data.frame(x=dens$x, y=y_scaled)
  p = pval
  plab = ifelse(p<0.0001,'p < 0.0001', paste0('p = ', formatC(p, format = 'f', digits = 4)))
  ggplot(fitdata, aes(x=value,y=mu))+
    geom_ribbon(
      data = hist_scaled,
      aes(x = x, ymin = smooth_range[1], ymax = y),
      inherit.aes = FALSE,
      fill = 'lightblue2', alpha = 0.6
    ) +
    geom_line(color = ifelse(sig, 'blue', grey(0.2)),
              linetype = ifelse(sig, "solid", "dashed"), 
              linewidth = 1) +
    geom_ribbon(aes(ymin = mu_lo, ymax = mu_hi),
                fill = "grey", alpha = 0.5) +
    theme_bw() +
    labs(x = var_names[vars==v], y = 'Damage severity') +
    annotate("text",
             x = -Inf, y = +Inf,
             label = plab,
             hjust = -.25, vjust = 2,
             size = 4)
})

#plot 
z = ggarrange(plotlist = all_plots, labels = LETTERS, ncol = 3, nrow = 1)
z
ggsave('../figs/orchard-effects.jpg', z, width=8, height=3, dpi=600)

backtransform = bind_rows(lapply(fit_data, function(f) f$fitdata))
backtransform %>% group_by(variable) %>%
  summarize(
    min_sev = min(mu),
    max_sev = max(mu),
    range_sev = diff(range(mu)),
    max_sev_at = value[which.max(mu)],
    min_x = min(value),
    max_x = max(value)
  ) %>% arrange(-range_sev)


