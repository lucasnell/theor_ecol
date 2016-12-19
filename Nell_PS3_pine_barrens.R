#########################################################################
# Problem Set 3 for Pine Barrens, with questions 6, 7, and 8
#########################################################################

# load needed libraries
library(ape)
library(phylolm)
library(vegan)
library(lme4)
library(pez)
library(lmerTest)



# Use this version fo the data set that formats the species names differently
d <- read.csv(paste0('./box/ZooEnt_540_2016/Data set analyses/',
                     'Pine barrens data 28Nov16/veg_aggregated_2012_24Nov16.csv'), 
              header = T)
d$site <- as.factor(d$site)
# summary(d)

w <- reshape(d, direction = "wide", timevar = "sp", v.names = "count", idvar = "site")

# I'm going to reduce the number of species.
nlevels(d$sp)
min.occur <- 5
min.mean.count <- .5
for (i in levels(d$sp)) {
	if (sum(d$count[d$sp == i] > 0) < min.occur | mean(d$count[d$sp == i]) < min.mean.count) 
		d$count[d$sp == i] <- NA
}

d <- d[!is.na(d$count), ]
d$sp <- droplevels(d$sp)
d$site <- droplevels(d$site)
nlevels(d$sp)

# add density variable
d$den <- asin((d$count/50)^0.5)

##  Import data "envi.csv"
e <- read.csv(file = paste0('./box/ZooEnt_540_2016/Data set analyses/',
                             'Pine barrens data 28Nov16/envi.csv'), 
                      header = T)
e <- e[e$date == 2012,]
# summary(e)

# Add a subset of environmental variables to d; you can add more variables if you want
d$tmin <- e$tmin[match(d$site, e$site)]
d$H <- e$H.prop[match(d$site, e$site)]
d$Mn <- e$Mn.mg_kg[match(d$site, e$site)]
d$Shade <- e$shade[match(d$site, e$site)]
d$Precip <- e$precip[match(d$site, e$site)]

# Standardize environmental variables to have mean = 0 and variance = 1 so that the coefficients can be compared. This also makes analyses more stable (less likely to have convergence problems).

for (i in 5:dim(d)[2]) d[, i] <- (d[, i] - mean(d[, i]))/sd(d[, i])
# summary(d)

############################################################################
# Previous analyses for species-specific variation in responses to 
# environmental variables 
############################################################################

envir.names <- c("tmin", "H", "Mn", "Shade", "Precip")

mlm.pvals <- NULL

z <- lmer(den ~ 1 + (1 | sp) + tmin + (0 + tmin | sp) + H + (0 + H | sp) + Mn + (0 + Mn | sp) + Shade + (0 + Shade | sp) + Precip + (0 + Precip | sp) + (1 | site), data = d, REML = F)
summary(z)

z0 <- lmer(den ~ 1 + (1 | sp) + tmin + H + (0 + H | sp) + Mn + (0 + Mn | sp) + Shade + (0 + Shade | sp) + Precip + (0 + Precip | sp) + (1 | site), data = d, REML = F)
mlm.pvals <- c(mlm.pvals, anova(z, z0)[2,8])

z0 <- lmer(den ~ 1 + (1 | sp) + tmin + (0 + tmin | sp) + H + Mn + (0 + Mn | sp) + Shade + (0 + Shade | sp) + Precip + (0 + Precip | sp) + (1 | site), data = d, REML = F)
mlm.pvals <- c(mlm.pvals, anova(z, z0)[2,8])

z0 <- lmer(den ~ 1 + (1 | sp) + tmin + (0 + tmin | sp) + H + (0 + H | sp) + Mn + Shade + (0 + Shade | sp) + Precip + (0 + Precip | sp) + (1 | site), data = d, REML = F)
mlm.pvals <- c(mlm.pvals, anova(z, z0)[2,8])

z0 <- lmer(den ~ 1 + (1 | sp) + tmin + (0 + tmin | sp) + H + (0 + H | sp) + Mn + (0 + Mn | sp) + Shade + Precip + (0 + Precip | sp) + (1 | site), data = d, REML = F)
mlm.pvals <- c(mlm.pvals, anova(z, z0)[2,8])

z0 <- lmer(den ~ 1 + (1 | sp) + tmin + (0 + tmin | sp) + H + (0 + H | sp) + Mn + (0 + Mn | sp) + Shade + (0 + Shade | sp) + Precip + (1 | site), data = d, REML = F)
mlm.pvals <- c(mlm.pvals, anova(z, z0)[2,8])

names(mlm.pvals) <- envir.names

mlm.pvals

############################################################################
# Input phylogeny.
############################################################################
phy <- read.tree(paste0('./box/ZooEnt_540_2016/Data set analyses/',
                        'Pine barrens data 28Nov16/pine_barrens_phylo_24Nov16.tre'))

# Remove species from phylogeny that do not occur in the selected species
phy <- drop.tip(phy, tip = phy$tip.label[!is.element(phy$tip.label, as.character(levels(d$sp)))])
length(phy$tip.label)

# Remove species from selected list if they don't occur in the phylogeny
d <- d[is.element(d$sp, phy$tip.label), ]
d$sp <- droplevels(d$sp)

# check the species contained in the phylogeny and data set
cbind(sort(phy$tip.label), as.character(levels(d$sp)))

# # plot phylogeny
# plot(phy)

############################################################################
# Input trait data
############################################################################

trait <- read.csv(paste0('./box/ZooEnt_540_2016/Data set analyses/',
                         'Pine barrens data 28Nov16/trait_24Nov16.csv'), 
                         header = T)

# Drop species with SLA = NA
trait <- trait[!is.na(trait$SLA),]
trait$sp <- droplevels(trait$sp)

# Remove species from trait list that do not occur in the selected species
x <- trait[is.element(trait$sp, levels(d$sp)), ]
x$sp <- droplevels(x$sp)
rownames(x) <- x$sp

# Remove species from trait list that do not occur in the phylogeny
x <- trait[is.element(trait$sp, phy$tip.label), ]
x$sp <- droplevels(x$sp)

# Remove species from the site data set that do not occur in the trait list
dim(d)[1]/30
d <- d[is.element(d$sp, x$sp), ]
d$sp <- droplevels(d$sp)
dim(d)[1]/30

# Remove species from phylogeny that do not occur in the trait list
length(phy$tip.label)
phy <- drop.tip(
    phy, 
    tip = phy$tip.label[!is.element(phy$tip.label, as.character(levels(x$sp)))])
length(phy$tip.label)

# add traits to data set d
d$LDMC <- x$LDMC[match(x$sp, d$sp)]
d$VH <- x$Vegetative.height[match(x$sp, d$sp)]
d$LT <- x$Leaf.thickness[match(x$sp, d$sp)]
d$SLA <- x$SLA[match(x$sp, d$sp)]

# standardize continuous trait variables
for(i in 10:13) d[,i] <- (d[,i] - mean(d[,i], na.rm = T))/sd(d[,i], na.rm = T)
# summary(d)
# tibble::as_data_frame(d)

############################################################################
# Question 6: Investigate the four traits to determine whether any of them
# could be responsible for differences in species relative abundances among communities.
# (Hint: The analyses are similar to those investigating whether environmental
# variables are responsible for differences in relative abundances.)
############################################################################

trait_names <- c('LDMC', 'VH', 'LT', 'SLA')

# all_mod_fun <- function(var_strings, extra_re = NULL) {
#     one_mod_len <- function(i, vars) {
#         var_mat <- t(combn(vars, i))
#         apply(var_mat, 1, function(x) {
#             paste(c('(1 | sp)', as.character(x), 
#                     paste0('(0 + ', as.character(x), ' | sp)'), extra_re), 
#                   collapse = ' + ')
#         })
#     }
#     out_list <- lapply(1:length(var_strings), one_mod_len, vars = var_strings)
#     return(c(out_list, recursive = TRUE))
# }

# all_mods <- all_mod_fun(trait_names, extra_re = '(1 | site)')


remove_re_forms <- function(var_strings, extra_re = NULL, gr_col = 'site') {
    
    full_mod <- paste(c(paste0('(1 | ', gr_col, ')'), var_strings, 
                        paste0('(0 + ', var_strings, ' | ', gr_col, ')'), extra_re), 
                      collapse = ' + ')
    
    mod0_mat <- t(combn(var_strings, length(var_strings) - 1))
    
    all_mods <- c(full_mod,
                  apply(mod0_mat, 1, function(x) {
                      paste(c(paste0('(1 | ', gr_col, ')'), var_strings, 
                              paste0('(0 + ', as.character(x), ' | ', gr_col, 
                                     ')'), extra_re), 
                            collapse = ' + ')
                  }))
    return(all_mods)
}

comp_re <- remove_re_forms(trait_names, extra_re = '(1 | sp)')

z <- lmer(as.formula(paste('den ~ 1 +', comp_re[1])), data = d, REML = F)
mlm_pvals <- sapply(2:length(comp_re), 
                    function(i){
                        form_i <- as.formula(paste('den ~ 1 +', comp_re[i]))
                        z0 <- lmer(form_i, data = d, REML = F)
                        anova(z, z0)[2,8]
                    })
names(mlm_pvals) <- c('SLA', 'LT', 'VH', 'LDMC')
mlm_pvals





############################################################################
# Question 7: How could phylogenetic signal affect your results from the
# analyses above? Use communityPGLMM (with family = 'gaussian') to test
# whether phylogenetic signal changes the results. I've set up some of the
# needed matrices to help.
############################################################################

# Set up random effects for communityPLMM
nspp <- nlevels(d$sp)
nsite <- nlevels(d$site)

# Make covariance matrix needed for LMMs
Vphy <- vcv(phy)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/det(Vphy)^(1/nspp)

# # random effect for site
# re.site <- list(1, site = d$site, diag(nsite))
# re.site.SLA <- list(d$SLA, site = d$site, diag(nsite))
# 
# 
# re.1 <- list(1, sp = d$sp, covar = diag(nspp))
# re.1.phy <- list(1, sp = d$sp, covar = Vphy)


dat <- d
dat$freq <- ifelse(d$count > 0, 1, 0)

# random effect for site
re.site <- list(1, site = dat$site, covar = diag(nsite))

# random intercept with species independent and with phylogenetic covariances
re.1 <- list(1, site = dat$sp, covar = diag(nspp))
re.1.phy <- list(1, site = dat$sp, covar = Vphy)



randos <- append(lapply(trait_names, 
                        function(s){list(dat[,s], sp = dat$sp, covar = diag(nspp))}), 
                 lapply(trait_names, 
                        function(s){list(dat[,s], sp = dat$sp, covar = Vphy)}))
randos <- append(randos, list(re.site, re.1, re.1.phy))

length(randos)



# # Below *should* work, but wasn't run
# z <- communityPGLMM(as.formula(paste0(c('den ~ ', paste0(trait_names, collapse = ' + ')),
#                                       collapse = '')), 
#                     data = dat,
#                     sp = dat$sp, site = dat$site,
#                     random.effects = randos,
#                     REML = TRUE, verbose = FALSE, s2.init = 0.1)
# summary(z)
# 
# z0_mods <- lapply(5:8,
#                   function(i){
#                       z0 <- communityPGLMM(as.formula(paste0(c('den ~ ', paste0(
#                           trait_names, collapse = ' + ')), collapse = '')), 
#                                            data = dat,
#                                            sp = dat$sp, site = dat$site,
#                                            random.effects = randos[-i],
#                                            REML = TRUE, verbose = FALSE, s2.init = 0.1)
#                   })
# 
# pglmm_pvals <- sapply(z0_mods, function(z0){
#     # Divide by two bc variances can't be negative, so you have to do a 1-tailed test
#     pchisq(2*(z$logLik - z0$logLik), df = 1, lower.tail = FALSE) / 2
# })
# pglmm_pvals


# Remove random effects (e.g., randos[[5]]:randos[[8]]) to see if phylogenetic differences
# are significant



############################################################################
# Question 8: Read Li & Ives (manuscript) that discusses how to test whether
# specific traits are associated with specific environmental variables in
# determining the abundances of species among communities. Confining your 
# attention to three traits:

# SLA
# LT
# VH

# and three environmental variables

# Shade
# Precip
# tmin

# test whether there is an association (interaction) between any trait and any
# environmental variable. First ignore phylogeny and then, if the results are
# significant, include phylogeny.
############################################################################


# What happens when you account for correlated data?
#   Correct type I error control
#   Efficiency higher (variance lower)
#   Increased "true" power
#   (Bias reduced)
#   (Consistent estimator)
#   () for those not always affected












