#' #' @title Functions used only internally 
#' #' @keywords internal
#' #' @NoRd
#' #' Remove negative values
#' .noneg <- function(x){
#'    res <- ifelse(x < 0, NA, x)
#'    return(res)
#' } 
#' 
#' #' Support Substitution Model  
#' .sub.data <- function (mcomp, posub, min, subvar, iv) {
#'   
#'   nd <- NULL
#'   newd <- vector("list")
#'   
#'   for (j in seq_len(min)) {
#'     sub <- posub * j
#'     for (k in seq_len(nrow(sub))) {
#'       newcomp <- mcomp + sub[k, ]
#'       names(newcomp) <- paste0(names(posub))
#'       nd[[k]] <- cbind(mcomp, newcomp, sub[k, ][[i]])
#'       
#'       }
#'     newd[[j]] <- do.call(rbind, nd)
#'     }
#'   newd <- as.data.table(do.call(rbind, newd))
#'   
#'   # Add more useful information to newd
#'   colnames(newd)[ncol(newd)] <- "MinSubstituted"
#'   newd[, Substitute := rep(subvar, length.out = nrow(newd))]
#'   newd$Predictor <- iv
#'   
#'   # Remove impossible reallocation that result in negative values - TODO
#'   cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")
#'   
#'   newd[, (cols) := lapply(.SD, .noneg), .SDcols = cols]
#'   newd <- newd[complete.cases(newd), ]
#'   
#'   return(newd)
#' }
#' 
#' .posub.data <- function(substitute, i) {
#'   posub <- copy(substitute)
#'   posub <- as.data.table(posub)
#'   posub <- posub[(get(i) != 0)]
#'   posub <- posub[order(-rank(get(i)))]
#'   
#'   posub
#' }
#' 
#' .get.mcomp <- function(data) {
#'   
#'   b <- data$CompIlr$BetweenComp
#'   mcomp <- mean(b, robust = TRUE)
#'   mcomp <- clo(mcomp, total = 1440)
#'   mcomp <- as.data.table(t(mcomp))
#'   colnames(mcomp) <- colnames(b)
#' 
#'   mcomp
#' }
