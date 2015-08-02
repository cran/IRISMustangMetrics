\name{basicStatsMetric}

\alias{basicStatsMetric}


\title{
Min, median, mean, and rms max of a signal
}
\description{
The basicStatsMetric() function calculates the \code{min, median, mean, max} and \code{rms} for
the incoming seismic signal.
}

\usage{
basicStatsMetric(st)
}

\arguments{
  \item{st}{ a \code{Stream} object containing a seismic signal}
}

\details{
This metric merely applies the \code{min, median, mean} and \code{max} methods of \code{Stream} objects
to the \code{st} parameter to calculate the following metrics:
\itemize{
  \item \code{sample_min}
  \item \code{sample_median}
  \item \code{sample_mean}
  \item \code{sample_max}
  \item \code{sample_rms}
}

Any error messages generated in the process will pass through untrapped.
}

\value{
A list of \code{SingleValueMetric} objects is returned.
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}
\note{
%%  ~~further notes~~
See the \code{seismic} package for documentation on \code{Stream} objects and the \code{getDataselect} method.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{ }

\examples{
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient")

# Get the waveform
starttime <- as.POSIXct("2012-01-24 00:00:00", tz="GMT")
endtime <- as.POSIXct("2012-01-25 00:00:00", tz="GMT")
st <- getDataselect(iris,"AK","PIN","","VEA",starttime,endtime)

# Calculate some metrics and show the results
metricList <- basicStatsMetric(st)
dummy <- lapply(metricList, show)
}

\keyword{metrics}
