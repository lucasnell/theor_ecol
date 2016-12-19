

# Custom ggplot2 theme
theme_lan <- function(base_size = 10, base_family = 'Helvetica') {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold', size = 11),
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid = element_blank()
        )
}

z_trans <- function(x) {(x - mean(x)) / sd(x)}
min_log <- function(x) {
    nz <- min(x[x>0])
    y <- log(x + nz)
    return(y)
}
anti_min_log <- function(y, x_orig){
    nz <- min(x_orig[x_orig>0])
    x <- exp(y) - nz
    return(x)
}
asin_sqrt <- function(x){
    y <- asin(sqrt(x))
    return(y)
}
anti_asin_sqrt <- function(y){
    x <- sin(y)^2
    return(x)
}

