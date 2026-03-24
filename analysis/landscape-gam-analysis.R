# Script to analyze effect of landscape varibles on pecan damage severith
# foloowing hurrican idalia
library(googledrive)
library(mgcv)
library(gratia)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(terra)
library(sf)

## Download and prepare data
#drive_download(file = 'pecan-landscape-data.csv', 
#               path = '../data/pecan-landscape-data.csv',
#               overwrite = TRUE)
landscape = read_csv('../data/pecan-landscape-data.csv')
wind_direction = rast('../data/IDALIA_directions.tif')

# Select only cover in downwind direction and combine to single value
# get max. wind direction to nearest 45 deg
landscape_sf = landscape %>% mutate(geometry = map(.geo, ~ st_as_sfc(.x, GeoJSON = TRUE))) %>%
  unnest(geometry) %>% st_as_sf()
landscape_sf$MaxWindDir = terra::extract(wind_direction, landscape_sf)$direction
landscape_sf$MaxWindDir = round(landscape_sf$MaxWindDir/45)*45
landscape_sf$MaxWindDir = landscape_sf$MaxWindDir %% 360

# Consider only using close by points, e.g., 100 km
idalia = read_sf('../data/Idalia_HurricanePath_WithPoints.shp')
aoi = st_buffer(idalia,dist=100000)

landscape_sf = st_intersection(landscape_sf, aoi)
nrow(landscape_sf)

study_area = landscape_sf %>% st_union() %>% st_convex_hull() %>% st_area()
study_area / (1000*1000)
landscape = st_drop_geometry(landscape_sf)
names(landscape) <- sub("^X\\d+_", "", names(landscape))

# Clean up column names
colnames(landscape) = gsub("^\\d+_","",colnames(landscape))
colnames(landscape) = gsub("^b(\\d+)$", "sand_\\1cm", colnames(landscape))

bbox = st_bbox(landscape_sf)
ggplot() + geom_sf(data = landscape_sf, aes(color=severity)) +
  geom_sf(data=idalia) +  
  #geom_sf(data=aoi, alpha=0.1) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
                                    ylim = c(bbox["ymin"], bbox["ymax"])) +
  scale_color_viridis_c()

#make an output for plotting subset
#write_sf(landscape_sf, '../data/pecan-landscape-subset.shp')

## Calculate exposure index
exposure = function(aspect, wind, slope) {
  aspect_rad <- aspect * pi / 180
  wind_rad <- wind * pi / 180
  angle_diff <- abs(aspect_rad - wind_rad)
  angle_diff <- ifelse(angle_diff > pi, 2*pi - angle_diff, angle_diff)
  exposure <- slope * cos(angle_diff)
}
landscape$exposure = with(landscape, exposure(aspect, MaxWindDir, slope))

# Get the Forest cover matching with the correspondin upwind direction
landscape$upwind = (landscape$MaxWindDir - 180) %% 360

colname_base = sprintf('cover_%03d_', landscape$upwind)
scales = c(1,3,5,10,30,50,100)
names(scales) = sprintf('upwindCover_%02d',scales)
upwind_cover = bind_cols(
  lapply(scales, function(s) {
    cols = paste0(colname_base, sprintf('%02d', s))
    vals = sapply(seq_along(cols), function(i) landscape[[i, cols[i]]])
    return(vals)
    }))

## remove multiple wind measurements and append the upwind Measurems
landscape = cbind(landscape, upwind_cover)
landscape <- landscape[ , !grepl("^cover_[0-9]+_[0-9]+$", names(landscape))]
landscape$upwind = NULL
landscape$tpi_m = landscape$mTPI
landscape$mTPI = NULL
clamp <- function(x, lo, hi) pmax(lo, pmin(x, hi))
landscape$severity = clamp(landscape$severity, 0.01, 0.99)

#gam
# Select correct predictor for each scale
selectVar = function(data, variable) {
  vars = grep(paste0('^', variable), colnames(data), value=TRUE)
  names(vars) = vars
  mods = lapply(vars, function(v) {
    if(length(unique(data[,v])) < 10) return(NULL)
    main_effects = c('idaliaMSW', 'ba2023', v)
    main_terms = paste0("s(", main_effects, ", bs = 'cs')", collapse=' + ')
    intx_terms = paste0("ti(", main_effects[-1], ", ", main_effects[1],  ", bs = c('cs','cs'))", collapse=' + ')
    f <- as.formula(paste0("severity ~ ", main_terms, ' + ', intx_terms))
    cat('fitting', v, '\n')
    gam(f, data=data, family=betar())
  })
  mods = Filter(Negate(is.null), mods) # drop nulls
  aic_vals = suppressWarnings(sapply(mods, AIC))
  aic_vals = sort(aic_vals)
  print(names(aic_vals)[1])
  invisible(aic_vals)
}

# determine at which scale best to measure each multi-scale variable
# And also determine scale sensitivity
scale_vars = c('cover', 'upwindCover', 'tpi', 'sand')
scale_sensitivity = bind_rows(lapply(scale_vars, function(var) {
  scores = selectVar(landscape, var)
  dAIC = scores - min(scores)
  aic_output = data.frame(variable = var, 
                          scale = as.numeric(gsub("\\D+", "", names(scores))),
                          dAIC = dAIC)
  return(aic_output)  
}))
filter(scale_sensitivity, dAIC == 0)
ggplot(scale_sensitivity, aes(x=scale, y=dAIC)) + 
  facet_wrap(~variable, scales = 'free_x') + 
  geom_line() + geom_point() + theme_bw() +
  labs(y = expression(Delta~AIC))
ggsave('../figs/scale-sensitivity.jpg', width=5, height=5,dpi=600)

#results, generally insensitive to sand depth, cover best at large scales, 
# upwind cover generally best at moderate scales, tpi at small scales.
# Update covers to % for easier labeling
landscape <- landscape %>% 
  mutate(across(contains("cover"),
                ~ if(max(.x, na.rm = TRUE) <= 1) .x * 100 else .x
                ))

selected_scales = rownames(filter(scale_sensitivity, dAIC == 0))

# Determine if there is high multi-collinearity
all_predictors = c( 'idaliaMSW', 'ba2023', 
                    grep('sand', selected_scales,value=TRUE),
                    grep('tpi', selected_scales,value=TRUE),
                    grep('^cover',selected_scales,value=TRUE),
                    grep('upwind',selected_scales,value=TRUE),
                     'exposure')

c=cor(landscape[all_predictors])
diag(c)=NA; range(c,na.rm=TRUE)
corrplot::corrplot(c,diag=FALSE)
max(c,na.rm=TRUE)

# Cover 100 and upwindCover 100 have high correlation r = 0.564
with(landscape, cor(cover_100, upwindCover_10))

model_predictors = setdiff(all_predictors, c('cover_100'))

# Create GAM , main effects plus all interactions with windspeed.
main_terms = paste0("s(", model_predictors, ", bs = 'cs')", collapse=' + ')
intx_terms = paste0("ti(", model_predictors[-1], ", ", model_predictors[1],  ", bs = c('cs','cs'))", collapse=' + ')
f <- as.formula(paste0("severity ~ ", main_terms, ' + ', intx_terms))
mod=gam(f, data = landscape, family=betar())
summary(mod)
gratia::draw(mod)

variables <- model_predictors
var_names = c('wind speed (m/s)',
              expression('basal area (m'^2*'/ha)'),
              'sand content (%)',
              'TPI (m)', 'upwind cover (%)', 'exposure') 
gam_output.s = summary(mod)$s.table
gam_output.s = data.frame(variable = rownames(gam_output.s), gam_output.s, row.names = NULL)
gam_output.s$variable <- gsub("s\\(|\\)", "", gam_output.s$variable)  # remove 's(' and ')'
gam_output.s$variable <- ifelse(grepl("^ti\\(", gam_output.s$variable), gsub(",", " x ",
       gsub("^ti\\(|\\)$", "", gam_output.s$variable)),gam_output.s$variable)
gam_output.s = select(gam_output.s, !Ref.df)
gam_output = gam_output.s
write_csv(gam_output, '../figs/gam_output-landscape.csv', col_names = FALSE)
knitr::kable(gam_output.s, digits = c(0,3,3,4))

vars = variables
df = landscape
seq = seq_along(vars)
names(seq) = vars
fit_data <- lapply(seq, function(i) {
  v <- vars[i]
  type = ifelse(v %in% model_predictors, 'main', 'itnx')
  if(type == 'main') {
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
  }
  if(type == 'intx') {
    
  }
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
plots_to_show = setdiff(names(all_plots), '') #drop NS plot
z = ggarrange(plotlist = all_plots[plots_to_show], labels=LETTERS, ncol=3,nrow=2)
z
ggsave('../figs/landscape-effects.jpg', z, width=8, height=6, dpi=600)

itx_ba = gratia::draw(mod, select='ba2023', partial_match=TRUE)
itx_sand = gratia::draw(mod, select='sand', partial_match=TRUE)
z_itx = ggarrange(itx_sand, itx_ba, labels = LETTERS, ncol=1, nrow=2)
ggsave('../figs/landscape-effects-itx.jpg', z_itx, width=8, height=8, dpi=600)

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