% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcEATFruitvegRatio.R
\name{calcEATFruitvegRatio}
\alias{calcEATFruitvegRatio}
\title{calcEATFruitvegRatio}
\usage{
calcEATFruitvegRatio(populationweight = "PopulationPast")
}
\arguments{
\item{populationweight}{datasource of populationweight: FAO can be selected in order to better meet exact values.
Normal datasource is PopulationPast}
}
\value{
List of magpie objects with results on country level, weight on country level, unit and description.
}
\description{
Calculates the share of fruits and vegetables in the calorie supply from the others MAgPIE commodity for the past.
Information on the calorie supply from fruits and vegetables is relevant in the context of dietary recommendations,
e.g. as proposed by the EAT-Lancet Commission on healthy diets from sustainable food systems.
}
\examples{
\dontrun{
calcOutput("EATFruitvegRatio")
}
}
\seealso{
\code{\link[madrat]{calcOutput}}, \code{\link{calcEATLancetTargets}}, \code{\link[mrfaocore]{calcFAOharmonized}},
\code{\link{calcEATLancetDiets}}
}
\author{
Isabelle Weindl, Felicitas Beier
}
