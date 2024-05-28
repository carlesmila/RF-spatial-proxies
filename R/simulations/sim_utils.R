#' Sample simulation
#' @details
#' Simulates a series of sampling points for simulation problem 1.
#' @param nsamples Integer. Number of samples to simulate.
#' @param dsamples Character. Spatial distribution of the samples. 5 are
#' possible: sregular, wregular, random, wclust, sclust.
#' @param sarea sf/sfc polygon where samples will be simulated.
sim1_samples <- function(nsamples, dsamples, sarea){

  if(dsamples=="sregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 2)
  }else if(dsamples=="wregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 5)
  }else if(dsamples=="random"){
    simpoints <- st_sample(sarea, nsamples)
  }else if(dsamples=="wclust"){
    simpoints <- clustered_sample(sarea, nsamples, 25, 8)
  }else if(dsamples=="sclust"){
    simpoints <- clustered_sample(sarea, nsamples, 10, 6)
  }

  simpoints <- st_sf(geometry=simpoints)
  simpoints
}

#' Evaluate quietly
#' @details
#' Function to hide undesired prints. Copied from:
#' https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#' @param x A function to evaluate.
#' @examples
#' quiet(print("Nothing is shown"))
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' Create a raster stack from a point sf object.
#' @details
#' Create a raster stack from a point sf object by indicating the relevant
#' columns to rasterize as well as their names.
#' @param sf_points A point sf object to rasterize.
#' @param cols_indx Integer vector with the columns of the sf point object to
#' rasterize.
#' @param layers_names Character vector with the names of the layers to be
#' applied to the stack. It must have the same dimensions of cols_indx.
rasterise_stack <- function(sf_points, cols_indx, layers_names){

  # Make sure of equal length of indices and names
  if(length(cols_indx)!=length(layers_names)){
    stop("Colummns indeces and layer names must have the same length.")
  }

  # Create stack from first element and name it
  res_stack <- rasterFromXYZ(
    cbind(st_coordinates(sf_points),
          as.matrix(as.data.frame(sf_points)[,cols_indx[1]], ncol=1)))
  names(res_stack) <- layers_names[1]

  # If there are more elements to be stacked, proceed
  if(length(cols_indx)>1){
    for(i in 2:length(cols_indx)){
      cindx <- cols_indx[i]
      temprast <- rasterFromXYZ(cbind(st_coordinates(sf_points),
                                      as.matrix(as.data.frame(sf_points)[,cindx], ncol=1)))
      names(temprast) <- layers_names[i]
      res_stack <- stack(res_stack, temprast)
    }
  }
  return(res_stack)
}


#' Simulates regular samples jittered by an amount of noise.
#' @details
#' Simulates regular samples jittered by an amount of noise ~ U(-amount, amount).
#' @param sarea sf/sfc polygon where samples will be simulated.
#' @param nsamples Integer. Number of samples to simulate.
#' @param amount Numeric. Amount of jitter to apply.
jitterreg_sample <- function(sarea, nsamples, amount){

  # Simulate regular points, jitter
  res <- st_sample(sarea, nsamples, type="regular")
  res <- as.data.frame(st_coordinates(res))
  res$X2 <- res$X + runif(nrow(res), -amount, amount)
  res$Y2 <- res$Y + runif(nrow(res), -amount, amount)

  # Ensure they fall within the sampling window, if not try again until they do
  res_sf <- st_as_sf(res, coords=c("X2", "Y2"), crs=st_crs(sarea))
  interF <- !st_intersects(res_sf, sarea, sparse = FALSE)
  while(any(interF)){
    res$X2[interF] <- res$X[interF] + runif(sum(interF), -amount, amount)
    res$Y2[interF] <- res$Y[interF] + runif(sum(interF), -amount, amount)
    res_sf <- st_as_sf(res, coords=c("X2", "Y2"), crs=st_crs(sarea))
    interF <- !st_intersects(res_sf, sarea, sparse = FALSE)
  }

  # Convert to geometries
  res$X <- NULL
  res$Y <- NULL
  res <- st_as_sf(res, coords=c("X2", "Y2"), crs=st_crs(sarea))
  res
}

#' Simulates clustered samples.
#' @details
#' Simulates clustered samples by simulating a number of randomly sampled
#' parents, and then randomly simulate children within a buffer of the parents.
#' @param sarea sf/sfc polygon where samples will be simulated.
#' @param nsamples Integer. Number of samples to simulate.
#' @param nparents Integer. Number of parents to simulate.
#' @param radius Numeric. Radius of the buffer for children simulation.
clustered_sample <- function(sarea, nsamples, nparents, radius){

  # Number of offspring per parent
  nchildren <- round((nsamples-nparents)/nparents, 0)

  # Simulate parents
  parents <- st_sf(geometry=st_sample(sarea, nparents, type="random"))
  res <- parents

  # Simulate offspring
  for(i in 1:nrow(parents)){

    # Generate buffer and cut parts outside of the area of study
    buf <- st_buffer(parents[i,], dist=radius)
    buf <- st_intersection(buf, sarea)

    # Simulate children
    children <- st_sf(geometry=st_sample(buf, nchildren, type="random"))
    res <- rbind(res, children)
  }

  return(res)
}
