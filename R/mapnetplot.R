## library(rgdal)
## library(rgeos)
## library(geosphere)


## mapnet.dir <- "~/Documents/networks/code/mapnetplot/"
## load(file.path(mapnet.dir,"manhattan-shapefile.RData"))
## load(file.path(mapnet.dir,"nyc-shapefile.RData"))
## load(file.path(mapnet.dir,"chicago-shapefile.RData"))

#####  EDA Network Plots
## load("chicago-shapefile.RData")

## col.1 <- adjustcolor("black", alpha=0.2)
## col.2 <- adjustcolor("red", alpha=0.6)
## col.3 <- adjustcolor("yellow",alpha=0.8)
## edge.pal <- colorRampPalette(c(col.1, col.2,col.3), alpha = TRUE,bias=1)


convert.to.edgelist <- function(AA){
    nn <- nrow(AA)
    el <- array(NA,c(nn*nn,3))
    rn <- 0
    for(ii in 1:nn){
        for(jj in 1:nn){
            rn <- rn + 1
            el[rn,] <- c(ii,jj,AA[ii,jj])
        }
    }
    return(el)
}

sh.base.plot <- function(shapefile,leg.pal=edge.pal,leg.count=1000,
                         include.legend=TRUE,
                         xlim=c(0,1), ylim=c(0,1),leg.lab="", ...){

    par(bg=grey.colors(100)[60],fg="blue")

    plot(0,0,xaxt="n",yaxt="n",xlab="",ylab="",frame=FALSE,
         xlim=xlim,ylim=ylim,...)
    rgeos::plot(shapefile,add=TRUE,border="white",lwd=2)

    if(include.legend){
        leg.col <- edge.pal(leg.count)
        leg.seq <- seq(from=xlim[1] + .01, to=xlim[2]-.01,length.out=leg.count)
        points(leg.seq,rep(ylim[1] - 0.0,leg.count),col=leg.col,pch=15)
        mtext(text=leg.lab,side=1,cex=1.75)
    }
}


chicago.base.plot <- function(shapefile=NULL,leg.pal=edge.pal,
                              leg.count=1000,include.legend=TRUE,
                              xlim=c(-87.96,-87.5),ylim=c(41.64,42.03),
                              trip.max=50,trip.min=1,...){

    if(is.null(shapefile)){
        shapefile <- chicago.sh
    }

    leg.text <- paste0(trip.min,"                  Trip Count                  ",trip.max)
    sh.base.plot(shapefile=shapefile,leg.pal=leg.pal,leg.count=leg.count,include.legend=include.legend,
                 xlim=xlim,ylim=ylim,leg.lab=leg.text, ...)


}

#' @export
map.net.plot <- function(net,pal=edge.pal,centroids,
                         el.max=max(net),pt.max,pt.min,
                         edge.size=1/25,pt.size=1.5){

    ##  Setting up Color Palette
    edge.col <- pal(100)
    pt.ind <- (rowSums(net) + colSums(net))/2

    if(missing(el.max)) el.max <- max(net)
    if(missing(pt.max)) pt.max <- max(pt.ind)
    if(missing(pt.min)) pt.min <- max(pt.ind) / 3

    ## edge.size <- 5
    ## point.size <- 3

    ##  EL Version
    el <- convert.to.edgelist(net)
    short.el <- el[el[,3] > 0,]
    short.el <- short.el[order(short.el[,3]),]

    for(xx in 1:nrow(short.el)){
        ii <- short.el[xx,1]
        jj <- short.el[xx,2]
        node1 <- centroids[ii,]
        node2 <- centroids[jj,]

        if(net[ii,jj] > 0){
            arc <- gcIntermediate(p1=c(as.numeric(node1[1]),
                                       as.numeric(node1[2])),
                                  p2=c(as.numeric(node2[1]),
                                       as.numeric(node2[2])),
                                  n=1000,addStartEnd=TRUE)
            edge.ind <- round(100 * (short.el[xx,3] / el.max))

            lines(arc,col=edge.col[edge.ind],lwd=edge.ind*edge.size)
        }
    }
    points(centroids, pch=19,col="orange",
           cex=pt.size*(pt.ind+pt.min)/(pt.max + pt.min))

}

#' @export
chicago.full.plot <- function(net,pal=edge.pal,
                             xlim=c(-87.96,-87.5),ylim=c(41.64,42.03),
                             centroids=chicago.centroids,
                             el.max=max(net),pt.max,pt.min,trip.max=max(net),
                             ...){

    chicago.base.plot(leg.pal=pal,xlim=xlim,ylim=ylim,trip.max=el.max)

    if(all(rownames(net) == chicago.short.labs))
        centroids <- chicago.short.centroids
    map.net.plot(net=net,pal=pal,centroids=centroids,
                 el.max=el.max,pt.max=pt.max,pt.min=pt.min, ...)
}

#' @export
chicago.downtown.plot <- function(net,pal=edge.pal,
                             xlim=c(-87.85,-87.5),ylim=c(41.83,41.95),
                             centroids=chicago.centroids,
                             el.max=max(net),pt.max,pt.min,
                             ...){

    chicago.base.plot(leg.pal=pal,xlim=xlim,ylim=ylim,trip.max=el.max)

    if(all(rownames(net) == chicago.short.labs))
        centroids <- chicago.short.centroids

    map.net.plot(net=net,pal=pal,centroids=centroids,
                 el.max=el.max,pt.max=pt.max,pt.min=pt.min, ...)
}


man.base.plot <- function(shapefile=NULL,leg.pal=edge.pal,leg.count=1000,
                          include.legend=TRUE,
                          xlim=c(-74.05,-73.91),ylim=c(40.68,40.88),
                          trip.max=50,trip.min=1,...){

    if(is.null(shapefile)){
        shapefile <- man.sh
    }

    if(include.legend)
        leg.text <- paste0(trip.min,"                Trip Count                ",trip.max)
    sh.base.plot(shapefile=shapefile,leg.pal=leg.pal,leg.count=leg.count,
                 include.legend=include.legend,
                 xlim=xlim,ylim=ylim,leg.lab=leg.text, ...)

}

#' @export
man.net.plot <- function(net,pal=edge.pal,
                         xlim=c(-74.05,-73.91),ylim=c(40.68,40.88),
                         centroids=man.ct.nodes[,4:3],
                         el.max=max(net),pt.max,pt.min,trip.max=max(net),
                         include.legend=TRUE,
                         ...){

    man.base.plot(leg.pal=pal,xlim=xlim,ylim=ylim,trip.max=el.max,
                  include.legend=include.legend)

    map.net.plot(net=net,pal=pal,centroids=centroids,
                 el.max=el.max,pt.max=pt.max,pt.min=pt.min, ...)
}
