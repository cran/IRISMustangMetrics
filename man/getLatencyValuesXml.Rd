\name{getLatencyValuesXml}

\alias{getLatencyValuesXml}
\alias{getLatencyValuesXml,IrisClient,character,character,character,character,character-method}
\alias{getLatencyValuesXml,IrisClient,character,character,character,character,missing-method}

\title{ Retrieve instantaneous data latency XML from the MUSTANG latency servlet.  }

\description{
The \code{getLatencyValuesXml} method makes a request of the IRIS DMC \code{latency} 
servlet and returns the unprocessed data latency XML for submission to the MUSTANG BSS (Backend Storage System).
}

\usage{
getLatencyValuesXml(obj, network, station, location, channel, url)
}

\arguments{
  \item{obj}{ \code{IrisClient} object }
  \item{network}{ character string with the two letter seismic network code }
  \item{station}{ character string with the station code }
  \item{location}{ character string with the location code }
  \item{channel}{ character string with the three letter channel code }
  \item{url}{ optional url of the latency web service }
}

\details{
Error returns from the webservice will stop evaluation and throw an error message.
}

\value{
An XML string is returned.
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}

\note{
%%  ~~further notes~~
This method is called by \code{latencyMetric.Rscript}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{IrisClient-class}}
}

%% \examples{ }

\keyword{ webservices }
