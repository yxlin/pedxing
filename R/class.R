##' An S4 class 
##'
##' @slot choice_RT choice_RT
##' @slot activation activation 
##' @slot inhibition inhibition
##' @slot leakage leakage
##' @slot dx dx
##' @slot counter counter
##' @export
setClass("LCA", slots = c(
           choice_RT  = "array",
           activation = "array",
           inhibition = "array",
           leakage    = "array",
           dx         = "array",
           counter    = "array"))


##' An S4 class for the Input-only model
##'
##' @slot choice_RT choice_RT
##' @slot activation activation 
##' @slot dx dx
##' @slot counter counter
##' @export
setClass("Input", slots = c(
  choice_RT  = "array",
  activation = "array",
  dx         = "array",
  counter    = "numeric"))

##' An S4 class for the Leakage model
##'
##' @slot choice_RT choice_RT
##' @slot activation activation 
##' @slot dx dx
##' @slot leak leakage
##' @slot counter counter
##' @slot samples simulation samples
##' @export
setClass("Leak", slots = c(
  choice_RT  = "array",
  activation = "array",
  dx         = "array",
  leak       = "array",
  counter    = "numeric"))

##' An S4 class for the Inhibition model
##'
##' @slot choice_RT choice_RT
##' @slot activation activation 
##' @slot dx dx
##' @slot leak leakage
##' @slot counter counter
##' @slot samples simulation samples
##' @export
setClass("Inhibition", slots = c(
  choice_RT  = "array",
  activation = "array",
  dx         = "array",
  leak       = "array",
  inhibition     = "array",
  inhibition_rec = "array",
  counter    = "numeric"))
