#
#  FlexMix: Flexible mixture modelling in R
#  Copyright (C) 2004 Friedrich Leisch
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

ExNPreg = function(n)
{
    if(n %% 2 != 0) stop("n must be even")
    
    x <- runif(2*n, 0, 10)

    data.frame(x=x,
               yn=c(5*x[1:n], 40-(x[(n+1):(2*n)]-5)^2)+3*rnorm(n),
               yp=rpois(2*n, exp(c(2-0.2*x[1:n], 1+0.1*x[(n+1):(2*n)]))),
               class = rep(1:2, c(n,n)),
               id1 = factor(rep(1:n, rep(2, n))),
               id2 = factor(rep(1:(n/2), rep(4, n/2))))
}


    
ExNclus = function(n=100)
{
    if(n %% 2 != 0) stop("n must be even")

    require("mvtnorm")
    
    rbind(rmvnorm(n, mean=rep(0,2)),
          rmvnorm(n, mean=c(8,0), sigma=diag(1:2)),
          rmvnorm(1.5*n, mean=c(-2,6), sigma=diag(2:1)),
          rmvnorm(2*n, mean=c(4,4), sigma=matrix(c(1,.9,.9,1), 2)))
}

    

    



    
