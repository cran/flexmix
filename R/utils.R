#
#  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
#  $Id: utils.R 4044 2008-07-26 12:30:25Z leisch $
#

list2object = function(from, to){
    n = names(from)
    s = slotNames(to)
    p = pmatch(n, s)
    if(any(is.na(p)))
        stop(paste("\nInvalid slot name(s) for class",
                   to, ":", paste(n[is.na(p)], collapse=" ")))
    names(from) = s[p]
    do.call("new", c(from, Class=to))
}

printIter = function(iter, logLik, label="Log-likelihood")
    cat(formatC(iter, width=4),
        label, ":", formatC(logLik, width=12, format="f"),"\n")
    

## library(colorspace)
## dput(x[c(1,3,5,7,2,4,6,8)])

## x = hcl(seq(0, 360*7/8, length.out = 8), c=30)
LightColors <- c("#F9C3CD", "#D0D4A8", "#9DDDD5", "#D1CCF5",
                 "#EDCAB2", "#AFDCB8", "#ACD7ED", "#EFC4E8")
    
## x = hcl(seq(0, 360*7/8, length.out = 8), c=100, l=65)
FullColors <- c("#FF648A", "#96A100", "#00BCA3", "#9885FF",
                "#DC8400", "#00B430", "#00AEEF", "#F45BE1")

