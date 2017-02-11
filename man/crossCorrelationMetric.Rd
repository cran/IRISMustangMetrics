\name{crossCorrelationMetric}

\alias{crossCorrelationMetric}

\title{Correlation between channels}

\description{
The crossCorrelationMetric() function calculates the maximum absolute correlation (\code{polarity_check}) 
and lag at maximum correlation (\code{timing_drift}) associated with two streams of seismic data.
}

\usage{
crossCorrelationMetric(st1, st2, maxLagSecs=10, filter=signal::butter(2,0.2))
}

\arguments{
  \item{st1}{a \code{Stream} object containing a seismic signal}
  \item{st2}{a \code{Stream} object containing a seismic signal}
  \item{maxLagSecs}{maximum number of seconds of lag to use}
  \item{filter}{a \pkg{signal} package filter to be applied before cross-correlating}
}

\details{
%% max preformatted length ---------------------------------------------------
Details of the algorithm are as follows:
\itemize{
  \item{Both signals are demeaned and detrended}
  \item{If one signal has a higher sampling rate, it is decimated to the lower sampling rate using an IIR filter if it is a multiple of the lower sample rate. See (\code{signal::decimate}).}
  \item{Both signals are filtered, by default with a Butterworth low pass filter. See (\code{signal::filter}).}
  \item{Signals are cross-correlated using the \code{stats::ccf()} function.}
}

The maximum absolute correlation is saved as \code{polarity_check} while the lag at peak correlation is saved as \code{timing_drift}.

\bold{Note:} For cross-correlation, seismic signals must not have any gaps -- they must be contained in a single \code{Trace} object.
}

\value{
A list with two \code{SingleValueMetric} objects is returned.
The metric names are \code{polarity_check} and \code{timing_drift}.
}

\note{
\code{timing_drift} is deprecated and will be removed in the future. 
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan@mazamascience.com} (R code),
Mary Templeton \email{met@iris.washington.edu} (algorithm)
}

\note{
The metricList generated for this two-channel metric will have an additional \code{sncl2} attribute identifying the SNCL in \code{st2}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{ }

\examples{
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient")

# Get the same signal, shifted by 3 seconds
starttime <- as.POSIXct("2013-11-12 07:09:45",tz="GMT")
endtime <- starttime + 600
st1 <- getSNCL(iris,"NM.SLM.00.BHZ",starttime,endtime)
st2 <- getSNCL(iris,"NM.SLM.00.BHZ",starttime+3,endtime+3)

# Cross-correlate
crossCorrelationMetric(st1,st2)
}

\keyword{metrics}
