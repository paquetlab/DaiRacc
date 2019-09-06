## Eric R. Paquet (eric.paquet@fsaa.ulaval.ca)
## Copyright : Université Laval
##

## Define the night here:
MAX.NIGHT = 20
MIN.NIGHT = 8

MatchFarmCow <- function(str){
    province <- substr(str,1,2) 
    ## stopifnot(province %in% c("QC","AB"))
    str <- gsub(sprintf("^%s",province),"",str)
    str <- gsub("^_","",str)
    str <- strsplit(str,"_")[[1]]
    str[1] <- as.character(as.numeric(str[1])) #farm id
    str[2] <- as.character(as.numeric(str[2])) #animal id
    
    list(province=province,
         farm=str[1],
         animal=str[2])
}

## Read CSV files exported from the HoboWare software
ReadHoboCSV <- function(filename,add20=FALSE){
    ## message(filename)
    ## First line contains sample name
    ## sample.name <- rev(strsplit(read.csv(filename,nrows=1,header=F)," ")[[1]])[1]
    sample.name <- as.character(read.csv(filename,nrows=1,header=F))
    data <- read.csv(filename,skip=1,header=T,check.names=F)
    ## add 20 to the year. Eg. 01/01/11 -> 01/01/2011
    if (add20){
        data[,2] <- gsub("/11 ","/2011 ",data[,2])
    }
    cur.time <- as.POSIXlt(data[,2], format = "%m/%d/%Y %I:%M:%S %p")
    cur.night <- as.numeric(lapply(cur.time,function(x){x$hour}))
    cur.night <- cur.night > MAX.NIGHT | cur.night < MIN.NIGHT
    list(name=sample.name,
         filename=filename,
         data=data,
         time=cur.time,
         night=cur.night,
         info=MatchFarmCow(filename))
}

PlotOneTrace <- function(trace,add.lowess=F){
    plot(trace$time,trace$data[,3],pch=".",ylab="Y acceleration",xlab="Time",ylim=c(-3,1))
    to.p <- trace$time[trace$night]
    points(to.p,rep(1,length(to.p)))
    if (add.lowess){
        lines(lowess(trace$time,trace$data[,3],f=1/500),lwd=1,col="red")
    }
}

CntSeq <- function(string){
    cs.to.ret <- c()
    cur.s <- ""
    cur.cnt <- 0
    
    for (si in string){
        if (si != cur.s){
            if (cur.s != ""){
                cs.to.ret <- rbind(cs.to.ret,c(cur.s,cur.cnt))
            }
            cur.s <- si
            cur.cnt <- 1
        } else {
            cur.cnt <- cur.cnt + 1
        }
    }

    if (cur.s != ""){
        cs.to.ret <- rbind(cs.to.ret,c(cur.s,cur.cnt))
    }
    cs.to.ret
}

HMMMatch <- function(accele.d,
                     prior=c(0.5,0.5),
                     transition=cbind(c(0.70,0.3),c(0.3,0.7)) ## base on Marianne's annotation it seems to be the closest approximate
                     ){
    ## 1 == down
    ## 2 == up
    ## Prior
    ## Up (50%), Down (50%)

    ## Emission
    ## P(Up | d) ~ N(d|-1,1)
    ## P(Down | d) ~ N(d|0,1)

    ## Transition
    ##      Up Down
    ## Up   0.9 0.1
    ## Down 0.1 0.9
    T1 <- matrix(NA,nrow=2,ncol=length(accele.d)) # Hold probabilities
    T2 <- matrix(NA,nrow=2,ncol=length(accele.d)) # hold idx most probable state

    prior <- log(prior)
    A <- log(transition) ## Transition
    Bf <- function(x){log(dnorm(x,mean=c(0,-1),sd=c(0.5,0.5)))} ## Emission function
    
    T1[,1] <- prior

    for (si in 2:length(accele.d)){
        B <- Bf(accele.d[si])
        T.A <- cbind(T1[,si-1],T1[,si-1]) + A
        T2[,si] <- apply(T.A,2,which.max)
        T1[,si] <- apply(T.A,2,max) + B
    }

    decoding <- c(which.max(T1[,ncol(T1)]))
    for (i in ncol(T2):2){
        decoding <- c(decoding,T2[decoding[length(decoding)],i])
    }
    decoding <- rev(decoding)
    list(T1=T1,T2=T2,decoding=decoding)
}

