#' Clean points via nearest neighbors
#'
#' This function XXX
#'
#' @param good_points XX
#' @param bad_points XX
#' @param k XX
#'
#' @export
#' @examples
#'
#' ## i am an example
#'
nearestNeighborCleanup = function(good_points,bad_points,k=21){
  pr <- class::knn(train=good_points[,c("longitude","latitude")],
                   test=bad_points[,c("longitude","latitude")],
                   cl=good_points[,c("assigned")],
                   k=k,
                   use.all=T)
  bad_points$KNN = pr
  return(bad_points)
}
