#' A simulated short version of the VA Dental Longitudinal Study
#'
#' A dataset containing the longitudinal health and dental records of 250 veteran men. This data was simulated to mimic the VA Dental Longitudinal Study described in Mitani et al.
#'
#'
#' @format A data frame with 20044 rows and 11 variables:
#' \describe{
#'   \item{subject}{subject ID}
#'   \item{tooth}{tooth ID}
#'   \item{visit}{visit number, max is 6}
#'   \item{toothstat}{indicator for whether or not the tooth is present at a specific visit}
#'   \item{baseage}{baseline age for each subject}
#'   \item{basenumteeth}{number of teeth present at baseline for each subject}
#'   \item{basesmoking}{indicator for smoking at baseline for each subject}
#'   \item{basemets}{indicator of presence or absence of metabolic syndrome for each subject}
#'   \item{maxcal5mm}{maximum clinical attachment loss for each tooth, 1: maxcal > 5mm, 0: maxcal <= 5mm}
#'   \item{prevmaxcal5mm}{maxcal5mm at the previous visit}
#'   \item{edu}{level of education for each subject, 1: has college degree, 0: does not have college degree}
#' }
"dental"
