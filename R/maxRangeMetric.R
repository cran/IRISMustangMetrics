##
##    Metric to calculate maximum sample range for a Stream of seismic data
##

maxRangeMetric <- function(st, window=300, increment=150) {

   snclq <- st@traces[[1]]@id
   starttime <- st@requestedStarttime
   endtime <- st@requestedEndtime
   x <- mergeTraces(st)@traces[[1]]

   n_samp <- window * x@stats@sampling_rate
   n_incr <- increment * x@stats@sampling_rate

   range <- seismicRoll::roll_range(x@data,n_samp,n_incr, align="left")
   maxRange <- max(range,na.rm=T)

   m1 <- new("GeneralValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="max_range", elementNames=c("value"), elementValues=maxRange)
   return(c(m1))
}
