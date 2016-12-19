#########################################################################
# Problem Set 2: Pine Barrens
#########################################################################

# load needed libraries
library(ape)
library(vegan)
library(lme4)
library(pez)

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
                            'Pine barrens data 28Nov16/envi.csv'), header = T)
e <- e[e$date == 2012,]
# summary(e)

# Add a subset of environmental variables to d; you can add more variables if you want
d$tmin <- e$tmin[match(d$site, e$site)]
d$H <- e$H.prop[match(d$site, e$site)]
d$Mn <- e$Mn.mg_kg[match(d$site, e$site)]
d$Shade <- e$shade[match(d$site, e$site)]
d$Precip <- e$precip[match(d$site, e$site)]

# Standardize environmental variables to have mean = 0 and variance = 1 so that the 
# coefficients can be compared. This also makes analyses more stable (less likely to 
# have convergence problems).

for (i in 5:dim(d)[2]) d[, i] <- (d[, i] - mean(d[, i]))/sd(d[, i])
# summary(d)

############################################################################
# Previous analyses for species-specific variation in responses to 
# environmental variables 
############################################################################

envir.names <- c("tmin", "H", "Mn", "Shade", "Precip")

make_form <- function(vars, sp_name = 'sp', prefix = 'den ~ 1 + (1 | sp)') {
    randos <- paste(paste0('(0 + ', vars, ' | ', sp_name, ')'), collapse = ' + ')
    fixed <- paste(vars, collapse = ' + ')
    string <- paste(prefix, fixed, randos, sep = ' + ')
    return(as.formula(string))
}
# make_form(envir.names)

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
# Question 4: Read Ives & Helmus (2011) at least through Methods: Model II.
# Then perform an analysis similar to the one above to investigate
# whether phylogenetically closely related species have the same responses (slopes)
# to tmin, Mn, Shade, and Precip (don't worry about H). How do you interpret the
# results? I have given you code below that uploads the phylogeny and
# matches it to the plant species data set. The question is essentially to construct
# model II in Ives & Helmus (2011). You will need "communityPGLMM()" in the {pez}
# library. I'm not going to set up the model in R for you, so you will have
# to check out the documentation. (I'm curious to see how easy the documentation
# is to follow.) Finally, this will take some time (>1 hour) to run. This is
# mainly because hardcore numerical stuff in communityPGLMM() is written in R
# rather than C++ (which runs 100x faster). You can use the "verbose = T"
# option if you want to follow the numerical progress.
############################################################################

# Input phylogeny.
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

# plot phylogeny
plot(phy)



# communityPGLMM(den ~ 1 + tmin + Mn + Shade + Precip, data = d, 
#                sp = d$species, site = d$site, random.effects = list(sp, ))

sapply(d, class)


dat <- d
dat$freq <- ifelse(d$count > 0, 1, 0)
nspp <- length(unique(d$sp))
nsite <- length(unique(d$site))
Vphy <- vcv(phy)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy / (det(Vphy) ^ (1/nspp))

re.site <- list(1, site = dat$site, covar = diag(nsite))

re.1 <- list(1, site = dat$sp, covar = diag(nspp))
re.1.phy <- list(1, site = dat$sp, covar = Vphy)



randos <- lapply(c('tmin', 'Mn', 'Shade', 'Precip'), 
                 function(s){
                     list(dat[,s], sp = dat$sp, covar = diag(nspp))
                 })

randos <- append(randos, 
                 lapply(c('tmin', 'Mn', 'Shade', 'Precip'), 
                           function(s){list(dat[,s], sp = dat$sp, covar = Vphy)}))

randos <- append(randos, list(re.site, re.1, re.1.phy))

length(randos)



z <- communityPGLMM(den ~ tmin + Mn + Shade + Precip, data = dat,
                    # family = "binomial",
                    sp = dat$sp, site = dat$site, 
                    random.effects = randos, 
                    REML = TRUE, verbose = FALSE, s2.init = 0.1)

head(dat)

summary(z)


# Remove random effects (e.g., tm.phy) to see if phylogenetic differences are significant


############################################################################
# Question 5: For a model that contains only Precip (to make it run faster),
# test whether there is phylogenetic signal in the response of species to Precip
# when you use presence/absence of species among sites as the dependent variable. 
# Do this with both a gaussian and a binomial model. Do
# you get the same results with both? Which do you believe more? 
############################################################################
