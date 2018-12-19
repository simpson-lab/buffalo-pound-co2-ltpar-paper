## Analyse the weekly Carbon data from BP

## packages
library('mgcv')
library('ggplot2')
library('cowplot')
theme_set(theme_bw())
library('viridis')
library('readr')
library('dplyr')
library('ISOweek')
devtools::install_github('gavinsimpson/gratia')
library('gratia')

## read in data
bp <- read_csv('bp-weekly-data.csv', guess_max = 1900)

## fix up the date info
## this uses ISOweek pkg to try to get proper dates by assigning Thursday as the
## day of the week...
bp <- mutate(bp, date = ISOweek2date(paste(year, sprintf('W%02d', week), 4,
                                           sep = '-')))
## ...but even this results in
with(bp, nrow(bp) - length(unique(date)))
## duplicate dates, though at least they are not NA
## Not sure what to do with these, but drop them for now
bp <- bp[with(bp, -which(duplicated(date))), ]

## format the `date` variable to %j, %Y, %M for modelling
meanDate <-  with(bp, mean(as.numeric(date)))
bp <- mutate(bp,
             modYear  = as.numeric(format(date, '%Y')),
             modMonth = as.numeric(format(date, '%m')),
             modDoy   = as.numeric(format(date, '%j')),
             cDate    = (as.numeric(date) - mean(as.numeric(date))) / 1000)

## plot the data
ggplot(bp, aes(x = modDoy, y = pCO2)) +
    geom_line() +
    facet_wrap(~ year)

ggplot(bp, aes(x = date, y = pCO2)) +
    geom_line()

## use more threads...
ctrl <- gam.control(nthreads = 7,
                    newton = list(maxHalf = 200))

## should we refit?
REFIT <- FALSE

## first go at a model
if (REFIT) {
    m1 <- gam(pCO2 ~ s(modYear, k = 35) + s(modDoy, k = 30, bs = 'cc'),
              data = bp, family = tw, method = 'ML', knots = list(modDoy = c(0.5, 366.5)),
              control = ctrl)

    m2 <- gam(pCO2 ~ s(modYear, k = 35) + s(modDoy, k = 30, bs = 'ad', xt = list(bs = 'cc')),
              data = bp, family = tw, method = 'ML', knots = list(modDoy = c(0.5, 366.5)),
              control = ctrl)

    m3 <- gam(pCO2 ~ s(modDoy, factor(modYear), k = 30, bs = 'fs', xt = list(bs = 'cc')),
              data = bp, family = tw, method = 'ML', knots = list(modDoy = c(0.5, 366.5)),
              control = ctrl)

    m4 <- gam(pCO2 ~ s(cDate, k = 35) + s(modDoy, factor(modYear), k = 30, bs = 'fs',
                                          xt = list(bs = 'cc')),
              data = bp, family = tw, method = 'ML', knots = list(modDoy = c(0.5, 366.5)),
              control = ctrl)

    m5 <- gam(pCO2 ~ te(modYear, modDoy, k = c(30, 40), bs = c('tp', 'cc')),
              data = bp, family = tw, method = 'ML', knots = list(modDoy = c(0.5, 366.5)),
              control = ctrl)

    m6 <- gam(pCO2 ~ s(cDate, k = 35) + s(modDoy, by = factor(modYear), id = 1, k = 30,
                                          bs = 'cc'),
              data = bp, family = tw, method = 'ML', knots = list(modDoy = c(0.5, 366.5)),
              control = ctrl)

    m7 <- gam(pCO2 ~ s(modDoy, factor(modYear), k = 30, bs = 'fs', xt = list(bs = 'cr')),
              data = bp, family = tw, method = 'ML',
              knots = list(modDoy = seq(0.5, 366.5, length = 30)),
              control = ctrl)

    m8 <- gam(pCO2 ~ te(cDate, modDoy, bs = c('tp', 'cc'), k = c(40,40)),
              data = bp, family = tw, method = 'ML', knots = list(modDoy = c(0.5, 366.5)),
              control = ctrl)

    ## just try to fit the entire time series as a smooth of time
    m9 <- gam(pCO2 ~ s(cDate, bs = 'tp', k = 400),
              data = bp, family = tw, method = 'ML', control = ctrl)

    ## adaptive smoother version of m9
    m10 <- gam(pCO2 ~ s(cDate, bs = 'ad', k = 250, m = 30),
               data = bp, family = tw, method = 'ML', control = ctrl)

    ## adaptive smoother version of m9 but with a seasonal cycle
    m11 <- gam(pCO2 ~ s(cDate, k = 250) + s(modDoy, bs = 'cc', k = 30),
               data = bp, family = tw, method = 'REML', control = ctrl)

    ## Gaussian process
    m12 <- gam(pCO2 ~ s(cDate, bs = 'gp', m = c(3, r = 200), k = 200),
               data = bp, family = tw, method = 'REML', control = ctrl,
               optimizer = c('outer', 'nlm'))

    ## write out models
    write_rds(m1, 'm1-gam-refitted.rds')
    write_rds(m2, 'm2-gam-refitted.rds')
    write_rds(m3, 'm3-gam-refitted.rds')
    write_rds(m4, 'm4-gam-refitted.rds')
    write_rds(m5, 'm5-gam-refitted.rds')
    write_rds(m6, 'm6-gam-refitted.rds')
    write_rds(m7, 'm7-gam-refitted.rds')
    write_rds(m8, 'm8-gam-refitted.rds')
    write_rds(m9, 'm9-gam-refitted.rds')
    write_rds(m10, 'm10-gam-refitted.rds')
    write_rds(m11, 'm11-gam-refitted.rds')
    write_rds(m12, 'm12-gam-refitted.rds')
} else {
    ## load models
    m1 <- read_rds('m1-gam-refitted.rds')
    m2 <- read_rds('m2-gam-refitted.rds')
    m3 <- read_rds('m3-gam-refitted.rds')
    m4 <- read_rds('m4-gam-refitted.rds')
    m5 <- read_rds('m5-gam-refitted.rds')
    m6 <- read_rds('m6-gam-refitted.rds')
    m7 <- read_rds('m7-gam-refitted.rds')
    m8 <- read_rds('m8-gam-refitted.rds')
    m9 <- read_rds('m9-gam-refitted.rds')
    m10 <- read_rds('m10-gam-refitted.rds')
    m11 <- read_rds('m11-gam-refitted.rds')
    m12 <- read_rds('m12-gam-refitted.rds')
}

## quick base plot for checking
## plot(pCO2 ~ date, data = bp, pch = 16, cex = 0.6)
## lines(predict(m8, type = "response") ~ date, data = bp, col = "red")

## M8 seems the best fitting - it's also the most complex
AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12)
## but need to double check for discontinuities at ends of year

M <- m8

pred <- with(bp, data.frame(date = seq(min(date), max(date), by = 1L)))
pred <- mutate(pred,
               modYear = as.numeric(format(date, '%Y')),
               modDoy  = as.numeric(format(date, '%j')),
               cDate   = (as.numeric(date) - meanDate) / 1000)
predvals <- predict(M, newdata = pred, se.fit = TRUE)
pred <- cbind(pred, as.data.frame(predvals))
pred <- mutate(pred, Fitted = exp(fit),
               Upper = exp(fit + (2 * se.fit)),
               Lower = exp(fit - (2 * se.fit)))

fit_plt <- ggplot(bp, aes(x = date, y = pCO2)) +
    geom_point(size = 0.9) +
    geom_ribbon(data = pred,
                mapping = aes(x = date, ymin = Lower, ymax = Upper),
                alpha = 0.2, fill = 'black', colour = NA, inherit.aes = FALSE) +
    geom_line(data = pred,
              mapping = aes(y = Fitted),
              colour = 'red') +
    labs(x = NULL, y = expression(italic(p)* CO[2]))

fit_plt

ggsave('pco2-fitted-values-m8-gam.pdf', fit_plt, height = 4, width = 16)
ggsave('pco2-fitted-values-m8-gam.png', fit_plt, height = 3, width = 12)

doy_plt <- ggplot(pred, aes(x = modDoy, y = Fitted, group = modYear, colour = modYear)) +
    geom_line() +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.95, name = 'Year') +
    labs(x = 'Day of year', y = expression(italic(p)* CO[2]))

doy_plt

ggsave('pco2-fitted-values-by-day-of-year.pdf', plot = doy_plt, height = 5, width = 9)

facetted_fit <- ggplot(pred, aes(x = modDoy, y = Fitted, group = modYear)) +
    geom_point(data = bp, mapping = aes(y = pCO2), size = 0.6) +
    geom_line(col = 'red') + facet_wrap(~ modYear) +
    labs(x = 'Day of year', y = expression(italic(p)* CO[2]))

facetted_fit

ggsave('pco2-facetted-fitted-values.pdf', plot = facetted_fit, height = 6, width = 9)

ice_dates <- unique(bp[, c('year','iceoffdoy','iceondoy')])
names(ice_dates)[1] <- 'modYear'

facetted_fit <- ggplot(pred, aes(x = modDoy, y = Fitted, group = modYear)) +
    geom_vline(data = ice_dates, mapping = aes(xintercept = iceoffdoy)) +
    geom_vline(data = ice_dates, mapping = aes(xintercept = iceondoy)) +
    geom_point(data = bp, mapping = aes(y = pCO2), size = 0.6) +
    geom_line(col = 'red') + facet_wrap(~ modYear) +
    labs(x = 'Day of year', y = expression(italic(p)* CO[2]))

facetted_fit

ggsave('pco2-facetted-fitted-values-with-ice-on-off.pdf',
       plot = facetted_fit, height = 6, width = 9)

## pred <- with(bp, expand.grid(date = seq(min(date), max(date), by = 1L),
##                              modDoy  = seq(50, 300, by = 50L)))
## pred <- mutate(pred,
##                modYear = as.numeric(format(date, '%Y')),
##                cDate   = (as.numeric(date) - meanDate) / 1000)

## predvals <- predict(M, newdata = pred, se.fit = TRUE)

## pred <- cbind(pred, as.data.frame(predvals))
## pred <- mutate(pred, Fitted = exp(fit),
##                Upper  = exp(fit + (2 * se.fit)),
##                Lower  = exp(fit - (2 * se.fit)),
##                fmodDoy = factor(modDoy, levels = seq(50, 300, by = 50L)))


pred <- with(bp, data.frame(date = seq(min(date), max(date), by = 1L)))
pred <- mutate(pred,
               modYear = as.numeric(format(date, '%Y')),
               modDoy  = as.numeric(format(date, '%j')),
               cDate   = (as.numeric(date) - meanDate) / 1000)

predvals <- predict(M, newdata = pred, se.fit = TRUE)

pred <- cbind(pred, as.data.frame(predvals))
pred <- mutate(pred, Fitted = exp(fit),
               Upper = exp(fit + (2 * se.fit)),
               Lower = exp(fit - (2 * se.fit)))

ggplot(subset(pred, modDoy %in% seq(70, 120, by = 10L)),
       aes(x = date, y = Fitted)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, x = date),
                alpha = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    facet_grid(modDoy ~ ., scales = 'free_y') +
    labs(y = expression(italic(p) * CO[2]), x = NULL)


## Kerri's question 1b

## generate a sample from the posterior of the model
min_doy <- 50
max_doy <- 160
newd <- with(bp, expand.grid(year = seq(min(year), max(year), by = 1L),
                             modDoy  = 50:160))
newd <- mutate(newd, date = as.Date(paste(year, modDoy, sep = '-'), '%Y-%j'))
newd <- newd[with(newd, order(year, modDoy)), ]
newd <- mutate(newd, cDate = (as.numeric(date) - meanDate) / 1000)

sims <- simulate(m8, nsim = 1000L, seed = 123, newdata = newd,
                 unconditional = TRUE)

colnames(sims) <- paste0('sim', seq_len(ncol(sims)))
aug_sims <- cbind(year = newd[['year']], doy = newd[['modDoy']], sims)

spl <- split(as.data.frame(aug_sims), list(year = factor(aug_sims[, 'year'])))

pco2_drop <- function(mat, error_value = NA, inv_link) {
    pco2_diff <- function(sim) {
        min_pco2 <- min(sim, na.rm = TRUE)
        max_pco2 <- max(sim, na.rm = TRUE)
        min_day <- which(sim == min_pco2)
        max_day <- which(sim == max_pco2)
        if (length(min_day) >1L) {
            min_day <- max(min_day)
        }
        if (length(max_day) >1L) {
            max_day <- max(max_day)
        }
        pco2_diff <- if (max_day > min_day) {
            error_value
        } else {
            max_pco2 - min_pco2
        }
        data.frame(pCO2Change = pco2_diff,
                   minDoY     = min_day,
                   maxDoY     = max_day)
    }
    imat <- inv_link(mat[, -c(1:2)])
    out <- lapply(seq_len(ncol(imat)), function(i, mat) pco2_diff(mat[, i]), mat = imat)
    out <- do.call('rbind', out)
    out
}

tmp <- lapply(spl, pco2_drop, inv_link = family(M)$linkinv)

pco2_change <- do.call('rbind', lapply(tmp, `[[`, 'pCO2Change'))
pco2_mean_change <- rowMeans(pco2_change, na.rm = TRUE)
pco2_interval <- as.data.frame(t(apply(pco2_change, 1L, quantile, probs = c(0.025, 0.975), na.rm = TRUE)))
names(pco2_interval) <- c('Lower', 'Upper')
pco2_change_summary <- cbind(year = as.numeric(names(spl)), mean = pco2_mean_change, pco2_interval)

pco2_change_plt <- ggplot(pco2_change_summary,
       aes(x = year, y = mean)) +
    geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
    labs(y = expression('Change in' ~ italic(p)*CO[2] ~ DoY[50-160]),
         x = NULL)

pco2_change_plt

ggsave('change-in-pco2-doys-50-160.pdf', pco2_change_plt, height = 5, width = 9)

doy_min_change <- do.call('rbind', lapply(tmp, `[[`, 'minDoY'))
doy_min_change <- doy_min_change + (min_doy - 1)
doy_min_mean_change <- rowMeans(doy_min_change, na.rm = TRUE)
doy_min_interval <- as.data.frame(t(apply(doy_min_change, 1L, quantile, probs = c(0.025, 0.975), na.rm = TRUE)))
names(doy_min_interval) <- c('Lower', 'Upper')
doy_min_change_summary <- cbind(year = as.numeric(names(spl)), mean = doy_min_mean_change, doy_min_interval)

doy_min_change_plt <- ggplot(doy_min_change_summary,
       aes(x = year, y = mean)) +
    geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
    labs(y = 'Day of year', x = NULL,
         title = expression(Minimum ~ italic(p)*CO[2]))

doy_min_change_plt

doy_max_change <- do.call('rbind', lapply(tmp, `[[`, 'maxDoY'))
doy_max_change <- doy_max_change + (min_doy - 1)
doy_max_mean_change <- rowMeans(doy_max_change, na.rm = TRUE)
doy_max_interval <- as.data.frame(t(apply(doy_max_change, 1L, quantile, probs = c(0.025, 0.975), na.rm = TRUE)))
names(doy_max_interval) <- c('Lower', 'Upper')
doy_max_change_summary <- cbind(year = as.numeric(names(spl)), mean = doy_max_mean_change, doy_max_interval)

doy_max_change_plt <- ggplot(doy_max_change_summary,
       aes(x = year, y = mean)) +
    geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
    labs(y = 'Day of year', x = NULL,
         title = expression(Peak ~ italic(p)*CO[2]))

doy_max_change_plt

doy_plts <- plot_grid(doy_max_change_plt, doy_min_change_plt, ncol = 1)
doy_plts

ggsave('day-of-year-of-peak-and-min-pco2.pdf', doy_plts, height = 7, width = 9)

write_csv(doy_min_change_summary, 'doy-min-pco2-summary.csv')
write_csv(doy_max_change_summary, 'doy-peak-pco2-summary.csv')
write_csv(pco2_change_summary, 'change-in-pco2-summary.csv')
