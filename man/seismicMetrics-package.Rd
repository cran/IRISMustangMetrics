\docType{package}

\name{seismicMetrics-package}
\alias{seismicMetrics-package}
\alias{seismicMetrics}

\title{Utilities for calculating seismic metrics from IRIS DMC data}

\description{
This package provides S4 classes and functions for calculating metrics from seismological 
data available from the IRIS Data Management Center (DMC) (\url{http://www.iris.edu/dms/nodes/dmc/}).
This package is part of the MUSTANG project and is intended \strong{for DMC internal use only}.
}

\section{Introduction}{
  The \pkg{seismicMetrics} package depends upon the \pkg{seismic} package which  defines
  new S4 classes and methods for manipulating seismic data.  Please see the "seismic-intro" 
  vignette for introductory examples on using \pkg{seismic}.
}

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}

\references{
IRIS DMC web services: \url{http://service.iris.edu/}
}

\section{History}{

  version 1.0.0 -- First Public Release
  
}

%% \seealso{ }

\examples{
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient", debug=TRUE)

# Get the seismic data
starttime <- as.POSIXct("2010-02-27 06:45:00",tz="GMT")
endtime <- as.POSIXct("2010-02-27 07:45:00",tz="GMT")
st <- getDataselect(iris,"IU","ANMO","00","BHZ",starttime,endtime)

# Apply a metric and show the results
metricList <- basicStatsMetric(st)
dummy <- lapply(metricList, show)
}

%% Optional keywords
