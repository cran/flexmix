#
#  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
#  $Id: z.R 3913 2008-03-13 15:13:55Z gruen $
#

###**********************************************************
## Backward compatibility

## component model driver
FLXglm <- FLXMRglm
FLXglmFix <- FLXMRglmfix
FLXmclust <- FLXMCmvnorm
FLXbclust <- FLXMCmvbinary

## concomitant model driver
FLXmultinom <- FLXPmultinom
FLXconstant <- FLXPconstant
