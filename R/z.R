#
#  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
#  $Id: z.R 4666 2011-02-23 15:52:35Z gruen $
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
