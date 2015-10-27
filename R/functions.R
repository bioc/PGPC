##' Function used to segment a single image
##' 
##' \code{segmentXman} is called by the function \code{segmentAllSpots}
##' to segment the single image of one spot.
##' 
##' @param x A \code{imageHTS} object.
##' @param uname  A character vector, containing the well names to segment. See
##' getUnames for details.
##' @param p List of parameters which are read from the file specified in 
##' \code{segmentationPar} of the \code{segmentWells} function.
##' @param access A character string indicating how to access the data. Valid 
##' values are \code{'local'}, \code{'server'} and \code{'cache'}, the default.
##' See \code{fileHTS} for details.
##' @param spot An single integer, indicating the spot number to segment. If it
##' is specified, the default \code{1} will automatically be used.
##' @return Returns a \code{list} with the following items
##' cal: calibrated RGB image of the different channels
##' nseg: nuclear segmentation mask
##' cseg: cell segmentation mask
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{segmentAllSpots},\link{getUnames}}
##' 
##' @examples
##' 
##' ## see also section 2.1 "Image processing on cluster"
##' localPath <- tempdir()
##' serverURL <- system.file("extdata", package = "PGPC")
##' 
##' x <- parseImageConf(file.path("conf", "imageconf.txt"), 
##'                     localPath=localPath, 
##'                     serverURL=serverURL) 
##'                    
##' well <- "045-01-C23"
##' 
##' ## get segmentation parameter
##' p <- readHTS(x, type = "file", 
##'              filename = file.path("conf", "segmentationpar.txt"), 
##'              format = "dcf")
##'             
##' segmentation <- segmentXman(x, well, p, access="cache", spot=1)
##' 
segmentXman <- function(x, uname, p, access, spot=NULL)
{
  msg <- function(...) cat(..., "\n")

  if (is.null(spot)){
    msg("No spot defined, using spot 1")
    spot = 1
  } 
  fnuc = fileHTS(x, "source", uname = uname, channel = 1, 
                 spot = spot, access = access)
  fcyt = fileHTS(x, "source", uname = uname, channel = 2, 
                 spot = spot, access = access)
  fe = file.exists(c(fnuc, fcyt))
  if (any(!fe))
    stop("Missing source files: cannot open ", fnuc)
  
  cat(unlist(p), "\n")

  pchar <- c(names(p)[grep("seg.", names(p))], "")
  pnum = which(!(names(p) %in% pchar))
  p[pnum] = lapply(p[pnum], as.numeric)
 
  scale.nuc = p[["scale.nuc"]]
  scale.cyt = p[["scale.cyt"]]
  thresh.nuc.peak =  p[["thresh.nuc.peak"]]
  thresh.nuc.mask =  p[["thresh.nuc.mask"]]
  thresh.cyt =  p[["thresh.cyt"]]
  thresh.cyt.global =  p[["thresh.cyt.global"]]

  msg("read images")
  if (is.null(scale.nuc) || length(scale.nuc) == 0) 
    scale.nuc = 16    
  if (is.null(scale.cyt) || length(scale.cyt) == 0)
    scale.cyt = 16
  
  inuc = readImage(fnuc)
  icyt = readImage(fcyt)
  fe = c(length(inuc), length(icyt)) != 0
  if (any(!fe)) 
    stop("Missing source files: cannot open ", fnuc, fcyt)

  if(FALSE){
    ## mask duplicated area
    borderMask = array(1, dim = dim(inuc))
    
    if(spot == 2){
      msg("crop spot2")
      borderMask[1:(p$border.crop),] = 0
    } 
    if(spot == 3){
      msg("crop spot3")
      borderMask[,1:(p$border.crop)] = 0
    }
    if(spot == 4){
      msg("crop spot4")
      borderMask[1:(p$border.crop),] = 0
      borderMask[,1:(p$border.crop)] = 0
    }
  }
  
  ## crop image
  crop = c((p$border.crop+1):(dim(inuc)[1] - p$border.crop))
  inuc = inuc[crop, crop]
  icyt = icyt[crop, crop]		

  msg("smooth image")
  inucSmooth = gblur(inuc,sigma=1)

  msg("computing nseg")
  seed = bwlabel(opening(thresh(inucSmooth, 10, 10, thresh.nuc.peak),
    kern=makeBrush(1,shape="disc")))
  nseg = propagate(inucSmooth, 
    seed, 
    mask=fillHull(thresh(inucSmooth, 20, 20, thresh.nuc.mask)))
  
  nf <- as.data.frame(computeFeatures.shape(nseg))
  rmindexNseg = which(nf$s.area > p$nuc.area.max | nf$s.area < p$nuc.area.min)
  nseg = fillHull(rmObjects(nseg, rmindexNseg))

  msg("segmenting cells using Voronoi tesselation")
  actinMask = thresh(icyt, 100, 100, thresh.cyt) | icyt > thresh.cyt.global
  cseg = propagate(icyt, 
    nseg,
    lambda=1.0e-10, 
    mask = actinMask)
  cseg = fillHull(cseg)

  hf <- as.data.frame(computeFeatures.shape(cseg))
  rmindexCseg = which(hf$s.area > p$cyt.area.max)
  nseg = fillHull(rmObjects(nseg, rmindexCseg))
  cseg = fillHull(rmObjects(cseg, rmindexCseg))
  res = list(cal = rgbImage(green = inuc* scale.nuc,  
                            red = icyt * scale.cyt), 
             nseg = nseg, 
             cseg = cseg)
  
  msg("segmentation OK")
  res
}

##' Function used to segment all spots for a well
##' 
##' \code{segmentAllSpots} is called by the \code{imageHTS} function 
##' \code{segmentWells} to segment the images of the screen
##' 
##' @param x A \code{imageHTS} object.
##' @param uname  A character vector, containing the well names to segment. See 
##' getUnames for details.
##' @param p List of parameters which are read from the file specified in 
##' \code{segmentationPar} of the \code{segmentWells} function.
##' @param access A character string indicating how to access the data. Valid
##' values are \code{'local'}, \code{'server'} and \code{'cache'}, the default.
##' See \code{fileHTS} for details.
##' @return Returns a \code{list} with the following items
##' cal: calibrated RGB image of the different channels
##' nseg: nuclear segmentation mask
##' cseg: cell segmentation mask
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{imageHTS}}, \code{\link{segmentWells}}, 
##' \code{\link{fileHTS}}
##' 
##' @examples
##' 
##' ## see section 2.1 Image processing on cluster for a working example
##' localPath = tempdir()
##' serverURL = system.file("extdata", package = "PGPC")
##' 
##' imageConfFile = file.path("conf", "imageconf.txt")
##' 
##' ## circumvent memory problem on 32bit windows by segementing only spot 1.
##' if(.Platform$OS.type == "windows" & R.Version()$arch == "i386")
##'     imageConfFile = file.path("conf", "imageconf_windows32.txt")
##' 
##' x = parseImageConf(imageConfFile, 
##'                    localPath=localPath, 
##'                    serverURL=serverURL)
##'
##'                    
##' well = "045-01-C23"
##' 
##' ## get segmentation parameter
##' p = readHTS(x, type = "file", 
##'             filename = file.path("conf", "segmentationpar.txt"), 
##'             format = "dcf")
##'             
##' segmentation = segmentAllSpots(x, well, p, access="cache")
##' 
##' 
segmentAllSpots <- function(x, uname, p, access){
  spots = as.numeric(getImageConf(x)$SpotNames)
  resAll = lapply(spots, segmentXman, x=x, uname=uname, p=p, access=access)
  
  cal = EBImage::combine(lapply(resAll, function(z) z[["cal"]])) 
  nseg = EBImage::combine(lapply(resAll, function(z) z[["nseg"]]))
  cseg = EBImage::combine(lapply(resAll, function(z) z[["cseg"]]))

  res = list(cal = cal, nseg = nseg, cseg=cseg)
  cat(" nbcells=", countObjects(nseg), sep = "")              
  res
}

startProcessing <- function(text) {
  time1 = Sys.time()
  cat(text, " ... ")
  if (nchar(text) < 25) {
    for (i in (nchar(text)+1):25) { cat(" ") }
  }
  return(time1)
}

endProcessing <- function(time1) {
  time2 = Sys.time()
  cat("finished. Time: ",format(difftime(time2, time1, tz = "",
                                         units = c("auto"))),".\n")
  return(invisible(NULL))
}


## calculate an estimate of local cell density using the radius of the xth 
## nearest neighbor
localDensityNN <- function(data, k=c(10, 15, 20), 
                           imgSize, 
                           offset=0, 
                           edgeCorrection="none") {

  ## use quadtree to search for NN
  tree = createTree(data)
  ## create matrix for the local density data
  NN = matrix(0,nrow=nrow(data), ncol=length(k))
  ## split calculation (could now work in one go!)
  M = ceiling(nrow(data) / 1000)
  I = sort(rep(1:M,length.out=nrow(data)))
  for (m in seq_len(M)) {
    ## looks up the max(k)th nearest neighbors, selects the ones specified in k 
    NN[I == m] = knnLookup(tree=tree, 
                           newx=data[I == m,1], 
                           newy=data[I == m,2], 
                           k=max(k))[,k]
  }
  ## calculate local cell density for k NN  
  tmpDen = list()
  for (s in seq_along(k)) {
    SqrtEuclDist = (data[,1] - data[NN[,s],1])^2 + 
        (data[,2] - data[NN[,s],2])^2
    if(edgeCorrection == "circular"){
      ## correct edge effects
      distToCenter = sqrt((data[,1] - imgSize[1]/2)^2 + 
                              (data[,2] - imgSize[2]/2)^2)
      maskRadius = imgSize[2]/2 - offset
      atEdge = sqrt(SqrtEuclDist) + distToCenter > maskRadius      
      SqrtEuclDist[atEdge] = 1/pi*
          intersectionAreaCircular(sqrt(SqrtEuclDist)[atEdge],
                                   rep(maskRadius, sum(atEdge)),
                                   distToCenter[atEdge])
    }
    if(edgeCorrection == "square"){
      ## correct edge effects
      distToCenter = apply(abs(data - imgSize/2), 1, max)
      mask = imgSize/2 - offset
      atEdge = sqrt(SqrtEuclDist) + distToCenter > mask      
      SqrtEuclDist[atEdge] = 
        1/pi*intersectionAreaSquared(sqrt(SqrtEuclDist)[atEdge],
                                     abs(distToCenter[atEdge]-mask))
    }
    ## calculate local density by deviding 
    tmpDen[[s]] = k[s] / (pi*SqrtEuclDist)
    names(tmpDen)[s] = paste("lcd.", k[s], "NN", sep="")
  }
  den = do.call(cbind, tmpDen) ## maybe keep as matrix
  den
}


intersectionAreaCircular = function(r, R, d){ 
  ## from http://mathworld.wolfram.com/Circle-CircleIntersection.html
  r^2*acos((d^2+r^2-R^2)/(2*d*r)) + R^2*acos((d^2+R^2-r^2)/(2*d*R)) -
    .5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))  
}


intersectionAreaSquared = function(r, h){ 
  ## from http://mathworld.wolfram.com/CircularSegment.html
  pi*r^2 - (r^2*acos((r - h)/r) - (r - h)*sqrt(2*r*h-h^2))
}

##' Function used to segment all spots for a well
##' 
##' \code{getFeaturesAllSpots} is called by the function \code{extractFeatures}
##' to segment the images of the screen
##' 
##' @param cal Calibrated RGB image matrix.
##' @param seg List of the nuclear and cell segmentation masks.
##' @param p List of parameters which are read from the file specified in 
##' \code{featurePar} of the \code{extractFeatures} function.
##' @return Returns a \code{data.frame} with the extracted features.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{imageHTS}}, \code{\link{extractFeatures}}
##'  
##' @examples
##' 
##' ## see section 2.1 Image processing on cluster for a working example
##' localPath <- tempdir()
##' serverURL <- system.file("extdata", package = "PGPC")
##' 
##' imageConfFile = file.path("conf", "imageconf.txt")
##' 
##' ## circumvent memory problem on 32bit windows by segementing only spot 1.
##' if(.Platform$OS.type == "windows" & R.Version()$arch == "i386")
##'     imageConfFile = file.path("conf", "imageconf_windows32.txt")
##' 
##' x = parseImageConf(imageConfFile, 
##'                    localPath=localPath, 
##'                    serverURL=serverURL)
##'                    
##' well <- "045-01-C23"
##' 
##' ## get segmentation parameter
##' p <- readHTS(x, type = "file", 
##'              filename = file.path("conf", "segmentationpar.txt"), 
##'              format = "dcf")
##'             
##' segmentation <- segmentAllSpots(x, well, p, access="cache")
##' 
##' ## get feature parameter
##' pf <- readHTS(x, type = "file", 
##'               filename = file.path("conf", "featurepar.txt"), 
##'               format = "dcf")
##'             
##' ftrs <- getFeaturesAllSpots(cal=segmentation$cal, 
##'                             seg=list(nseg=segmentation$nseg,
##'                                      cseg=segmentation$cseg), 
##'                             pf)
##' 
getFeaturesAllSpots = function(cal, seg, p){
##  browser(skipCalls=1)
  ## get number of spots
  nf = numberOfFrames(cal, "render")
  spots = 1:nf
  ## select red and blue channel
  nucCh = 2
  cytCh = 1
  nseg = seg$nseg
  cseg = seg$cseg

  if(length(dim(nseg)) < 3){
    dim(nseg) = c(dim(nseg), 1)
    dim(cal) = c(dim(cal), 1)
  }

  knn = as.numeric(p[["knn"]])
  lcdCols = paste("lcd.", knn, "NN", sep="")
  edgeCorrection = p[["edge.correction"]]

  time = startProcessing("Feature Extraction")

  myExpandRef <- function(ref, refnames) ref
  
  ftrsAll = lapply(spots, function(spot){
    print(spot)
    tmpFtrs = cbind(                  
      computeFeatures(x=nseg[,,spot], 
                      ref=cal[,,nucCh,spot], 
                      xname = "nseg", 
                      refnames="dna",
                      expandRef = myExpandRef),
      computeFeatures(x=cseg[,,spot], 
                      ref=cal[,,cytCh,spot], 
                      xname = "cseg", 
                      refnames=c("act"),
                      expandRef = myExpandRef),
      computeFeatures(x=cseg[,,spot], 
                      methods.noref=NULL,
                      ref=(cal[,,nucCh,spot] - mean(cal[,,nucCh,spot]))/
                        sd(as.vector(cal[,,nucCh,spot])) * 
                          (cal[,,cytCh,spot] - mean(cal[,,cytCh,spot]))/
                        sd(as.vector(cal[,,cytCh,spot])),
                      xname="cseg",
                      refnames="dnaact",
                      expandRef = myExpandRef)
    )                 
    if (is.null(tmpFtrs)){
      pFtrs = rbind(
        computeFeatures(x=nseg[,,spot], 
                        ref=cal[,,nucCh,spot], 
                        xname = "nseg", 
                        refnames="dna",
                        properties=TRUE,
                        expandRef = myExpandRef),
        computeFeatures(x=cseg[,,spot], 
                        ref=cal[,,cytCh,spot], 
                        xname = "cseg", 
                        refnames=c("act"), 
                        properties=TRUE,
                        expandRef = myExpandRef),   
        computeFeatures(x=cseg[,,spot],, 
                        methods.noref=NULL, ## do not recompute moments and 
                        ## shape for the cseg mask
                        ref=(cal[,,nucCh,spot] - mean(cal[,,nucCh,spot]))/
                          sd(as.vector(cal[,,nucCh,spot])) * 
                            (cal[,,cytCh,spot] - mean(cal[,,cytCh,spot]))/
                          sd(as.vector(cal[,,cytCh,spot])),
                        xname="cseg",
                        refnames="dnaact", 
                        properties=TRUE,
                        expandRef = myExpandRef)   
      )  
      eFtrs = matrix(nrow=0, ncol=nrow(pFtrs)+2+length(lcdCols))
      colnames(eFtrs) = c("spot", "id", pFtrs[,"name"], lcdCols)
      return(eFtrs)
    } else if (nrow(tmpFtrs) > max(knn))  {
      ## add local cell density
      cat("Calcualting lcd for spot", spot, "\n")
      lcd <- localDensityNN(tmpFtrs[, c("nseg.0.m.cx", "nseg.0.m.cy")], 
                            k=knn, 
                            imgSize=dim(nseg)[1], 
                            edgeCorrection=edgeCorrection)
      return(cbind(spot=spot, id=1:nrow(tmpFtrs), tmpFtrs, lcd))
    } else {
      cat("Cannot calculate lcd for spot", spot, "\n")
      lcd <- matrix(NA, 
                    nrow=nrow(tmpFtrs),
                    ncol=length(knn))
      colnames(lcd) = lcdCols
      return(cbind(spot=spot, id=1:nrow(tmpFtrs), tmpFtrs, lcd))
    }
  })  
  ftrs = do.call(rbind, ftrsAll)
  
  ## remove quantile features calculated over pixels
  pixelQuantileFtr = grep(".q0", colnames(ftrs))
  ftrs = ftrs[,-pixelQuantileFtr]
  
  endProcessing(time)
  ftrs
}

##' Function to merge profiles of extracted features from parallel
##' processing 
##' 
##' Merges all feature profiles of each well and saves the result in the
##' specified file.
##' 
##' @param x  A \code{imageHTS} object.
##' @param profilename pattern of profile file names.
##' @param output File name to save the merged profiles.
##' @param folder folder name in which the profile files are stored and in
##'        which the the result is saved.
##' @param access Access parameter passed to \code{fileHTS}
##' @return None, writes the merged profiles into the specified file on disk.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{imageHTS}}, \code{\link{extractFeatures}},
##' \code{\link{fileHTS}}
##'  
##' @examples
##' 
##' ## see section 2.1 Image processing on cluster for usage
##' 
mergeProfiles = function (x, profilename="profiles", output="profiles.tab", 
                          folder="data", access="cache"){
  
  dataPath = fileHTS(x, type="file", filename=folder, access=access)
  profiles = dir(dataPath, pattern=profilename)
  
  if(file.exists(file.path(dataPath, output)))
     profiles = profiles[-grep(output, profiles)]
  
  allProfiles = lapply(file.path(folder, profiles), function (filename){
                      readHTS(x, type="file", filename=filename, format="tab", 
                              access=access)
                    })
  
  profiles = do.call(rbind, allProfiles)
  profiles = profiles[order(profiles$uname),]
  
  ff = fileHTS(x, type = "file", filename = file.path("data", output), 
            createPath = TRUE, access = "local")
  write.table(profiles, file = ff, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
}


##' Function to summarize the extracted features per cell for each well
##' 
##' The function extends the \code{imageHTS} function \code{summarizeWells}. It
##' calculates summary statistics over all cells for each well.
##' In particular the trimmed mean (trim = 0.1) and sd is calculated for each
##' extracted feature. Additionally the 1%, 5%, 95% and 99% quantiles are 
##' calculated for features that are not the standard deviation, median absolut
##' deviation or Halralick statistics calculated over each cell.
##' 
##' 
##' @param x  A \code{imageHTS} object.
##' @param uname A character vector, containing the well names that will be
##' summarized.
##' @param featurePar File containing the feature parameters used for
##' summarizing the wells.
##' @param profileFilename File name to save the summarized features.
##' @param access Access parameter passed to fileHTS       
##' @return None, writes the summarized well profiles into the specified files 
##' on disk.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{imageHTS}}, \code{\link{summarizeWells}}
##'  
##' @examples
##' 
##' ## see section 2.1 Image processing on cluster for a working example
##' localPath = tempdir()
##' serverURL = system.file("extdata", package = "PGPC")
##' 
##' imageConfFile = file.path("conf", "imageconf.txt")
##' 
##' ## circumvent memory problem on 32bit windows by segementing only spot 1.
##' if(.Platform$OS.type == "windows" & R.Version()$arch == "i386")
##'     imageConfFile = file.path("conf", "imageconf_windows32.txt")
##' 
##' x = parseImageConf(imageConfFile, 
##'                    localPath=localPath, 
##'                    serverURL=serverURL)
##'                    
##' well = "045-01-C23"
##' 
##' segmentWells(x, well, file.path("conf", "segmentationpar.txt"))
##' PGPC:::extractFeaturesWithParameter(x, 
##'                                     well,
##'                                     file.path("conf", "featurepar.txt"))
##' 
##' summary <- summarizeWellsExtended(x,
##'                                   well, 
##'                                   file.path("conf", "featurepar.txt"))
##' 
summarizeWellsExtended <- function (x, 
                                    uname, 
                                    featurePar, 
                                    profileFilename =
                                        file.path("data", "profiles.tab"), 
                                    access = "cache") 
{
    p = readHTS(x, type = "file", filename = featurePar, access = access, 
        format = "dcf")
    profiles = as.list(rep(NA, length(uname)))
    names(profiles) = uname
    for (i in 1:length(uname)) {
        u = uname[i]
        cat(u, ": ", sep = "")
        ftrs = try(readHTS(x, "ftrs", uname = u, access = access))
        if (class(ftrs) != "try-error") {
            ftrs$spot = NULL
            ftrs$id = NULL
            n = nrow(ftrs)
            meanftrs = try(apply(ftrs, 2, mean, trim = 0.1, na.rm=TRUE),
                           silent = TRUE)
            if (class(meanftrs) == "try-error") 
                prof = NA
            else {
                names(meanftrs) = paste(names(meanftrs), ".mean", 
                  sep = "")

                sdftrs = apply(ftrs, 2, sd, na.rm=TRUE)
                names(sdftrs) = paste(names(sdftrs), ".sd", 
                  sep = "")

                ## calculate quantiles for all features except summary features
                ftrnames = names(ftrs)
                haralickFtr = grepl(".h.", ftrnames)
                sdFtr = grepl(".sd", ftrnames)
                madFtr = grepl(".mad", ftrnames)
                toRemove <- haralickFtr | sdFtr | madFtr
                
                probes = c(.01,.05, .95, .99)
                qtftrs = as.vector(apply(ftrs[,!toRemove], 2, quantile, 
                                         probs=probes, na.rm=TRUE))
                names(qtftrs) = paste(rep(names(ftrs[,!toRemove]), 
                                          each=length(probes)), 
                                      ".qt.",
                                      probes, 
                  sep = "")

                if (!is.null(p$cell.classes)) {
                  cfrac = rep(0, length(p$cell.classes))
                  names(cfrac) = p$cell.classes
                  f = fileHTS(x, "clabels", uname = u, access = access)
                  if (file.exists(f)) {
                    clabels = readHTS(x, "clabels", uname = u, 
                      access = access)$label
                    z = table(clabels)
                    cfrac[names(z)] = z/n
                  }
                  else {
                    msg = paste("cannot find the file that contains",
                                " class labels for well=", 
                                u, "\n", sep = "")
                    msg = paste(msg, "  maybe predictCellLabels has",
                                "not be called ?\n", 
                                sep = "")
                    msg = paste(msg, "  maybe the field 'cell.classes' in the",
                                "feature parameters file should be empty ?\n", 
                                sep = "")
                    stop(msg)
                  }
                  prof = c(n = n, meanftrs, sdftrs, qtftrs, cfrac)
                }
                else prof = c(n = n, meanftrs, sdftrs, qtftrs)
            }
            profiles[[i]] = prof
            cat("OK\n")
        }
        else cat("NA\n")
    }
    profiles = do.call(rbind, profiles)
    profiles = data.frame(uname = uname, profiles, stringsAsFactors = FALSE)
    rownames(profiles) = NULL
    if (all(is.na(profiles[, 1]))) 
        stop("no cell features found, no profiles generated.")
    else {
        ff = fileHTS(x, type = "file", filename = profileFilename, 
            createPath = TRUE, access = "local")
        write.table(profiles, file = ff, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
    }
    invisible(profiles)
}

##' Function to calculate chemical genetic interactions
##' 
##' To detect chemcial genetic interactions, the data of each feature is 
##' modeled using a multiplicative model and robust L1 regression to estimate
##' the effects of the cell line and drug treatment using the \code{medpolish}
##' function. In this iterative approach row and column median values are
##' subtracted alternately until the proportional change of the absolute
##' residuals falls below a defined threshold. The final row and column values
##' describe the drug and cell line effect respectively. The residuals 
##' represent the interaction terms. This process is done for each replicate
##' and each feature individually.
##' To detect significant interactions the values of replicates are used to
##' perform a moderated t-test against the null hypothesis t = 0, using the 
##' implementation in the Bioconductor package \code{limma}. p-values are
##' adjusted for multiple testing by controlling for the false discovery rate
##' using the method of Benjamini & Hochberg.
##'  
##' @param d \code{list} containing a \code{array} with four dimensions of 
##' features \code{D}, a \code{list} with annotation. The annotation 
##' \code{list} needs to contain a character \code{vector} \code{ftr} with the
##' feature names and a character \code{vector} \code{drug$Content} defining
##' whether the data comes from "sample" or "other" wells.
##' @param ftrs Parameter to select certain features for the caluclation of 
##' chemical genetic interactions.
##' @param samplesOnly If set to \code{TRUE} only values of "sample" wells will
##' be used for the calculation of chemical genetic interactions.
##' @param scaleByLine If set to \code{TRUE} the interaction terms for each
##' cell line will be scaled by the median absolute deviation of interaction
##' terms for the individual cell line and replicate.
##' @param ... Additional parameters passed to \code{medpolish}
##' @return A list with the annotation \code{anno}, raw data of the selected
##' features \code{D}, the chemical genetic interaction results \code{res},
##' the estimated drug and cell line effect \code{effect} and 
##' calculated p-values and multiple testing adjusted p-values.
##' @author Felix A. Klein, \email{felix.klein@@embl.de}
##' @seealso \code{\link{medpolish}}, \code{\link{interactions}},
##' \code{\link{lmFit}}, \code{\link{eBayes}}
##'  
##' @examples
##' 
##' data(interactions)
##' x <- getInteractions(interactions)
##' 
getInteractions = function(d, 
                           ftrs = NULL, 
                           samplesOnly=FALSE, 
                           scaleByLine=FALSE,
                           ...){

  anno = d$anno$drug$Content

  ## dimension format: 1: drug, 2: cell line, 3: replicate, 4: feature
  D = d$D

  if(samplesOnly){
    samples = anno == "sample" ## annotation is identical for every cellline
  } else {
    samples = TRUE
  }
  
  #########################################
  ## MEDPOLISH for specified features or all
  #########################################

  ftrs = if(!is.null(ftrs)) 
    match(ftrs, d$anno$ftr) else seq_along(d$anno$ftr)

  X = array(NA, dim = c(dim(D)[1:3], length(ftrs)))
  P = array(NA, dim = c(dim(D)[1:2], length(ftrs), 3))
  lineEff = array(NA, dim = c(dim(D)[2:3], length(ftrs)))
  drugEff = array(NA, dim = c(dim(D)[c(1,3)], length(ftrs)))
  
  for (ftr in ftrs){

    cat(d$anno$ftr[ftr], ":\n")

    ## check for inf: sum(is.infinite(D[,,2]))
    M1 = medpolish(D[samples,,1,ftr], na.rm=TRUE, ...)
    M2 = medpolish(D[samples,,2,ftr], na.rm=TRUE, ...)

    res1 = M1$residuals
    res2 = M2$residuals

    if(scaleByLine){
      res1 = apply(res1, 2, function(x) {x = x/mad(x, na.rm=TRUE)})
      res2 = apply(res2, 2, function(x) {x = x/mad(x, na.rm=TRUE)})
    }

    X[samples,,1,ftr==ftrs] = res1
    X[samples,,2,ftr==ftrs] = res2

    lineEff[,1,ftr==ftrs] = M1$col
    lineEff[,2,ftr==ftrs] = M2$col

    drugEff[samples,1,ftr==ftrs] = M1$row
    drugEff[samples,2,ftr==ftrs] = M2$row

    ## limma
    tmp = X[,,,ftr==ftrs]
    dim(tmp) = c(prod(dim(X)[1:2]), dim(X)[3])
    fit = lmFit(tmp)
    efit = eBayes(fit)

    P[samples,,ftr==ftrs, 1] = efit$p.value
    P[samples,,ftr==ftrs, 2] = p.adjust(efit$p.value,method="BH")
    P[samples,,ftr==ftrs, 3] = cor(tmp[,1], tmp[,2],
                                   use="pairwise.complete.obs")
    cat("OK\n")
  }
  
  effect = list(drug=drugEff, line=lineEff)

  if(!is.null(ftrs)) d$anno$ftr = d$anno$ftr[ftrs]

  list(anno=d$anno, D=d$D[,,,ftrs], res=X, effect=effect, pVal=P)       
}

selectPart = function(x, i, n) {
  a = max(ceiling(1+(i-1)*length(x)/n), 0)
  b = max(ceiling(i*length(x)/n), 0)
  z = x[a:b]
  z[!is.na(z)]
}

## function from imageHTS 
## p is passed to the evaluating function
extractFeaturesWithParameter <- function (x, 
                                          uname, 
                                          featurePar,
                                          access = "cache") 
{
  p = readHTS(x, type = "file", filename = featurePar, access = access, 
              format = "dcf")
  efMethod = p$extractfeatures.method
  for (u in uname) {
    cat(u, ": ", sep = "")
    cal = try(readHTS(x, "cal", uname = u, access = access))
    seg = try(readHTS(x, "seg", uname = u, access = access))
    if (class(cal) != "try-error" & class(seg) != "try-error") {
      ftrs = eval(call(efMethod, cal, seg, p))
      write.seg.tab = TRUE
      if (!is.null(p$feature.write.seg.tab)) 
        write.seg.tab = as.logical(p$feature.write.seg.tab)
      if (write.seg.tab) {
        fftrs = fileHTS(x, "ftrs", uname = u, createPath = TRUE, 
                        access = "local")
        write.table(ftrs, file = fftrs, sep = "\t", quote = FALSE, 
                    row.names = FALSE, col.names = TRUE)
      }
      cat(" OK\n")
    }
    else cat("NA\n")
  }
}

## extends the thumbnails written in the segmentWells function from imageHTS 
## may be include this as an option there
writeSegmentationThumbnail = function(x, 
                                      uname, 
                                      segmentationPar, 
                                      access = "local"){
  if (length(uname) > 1) {
    for (u in uname) {
      z = try(writeSegmentationThumbnail(x, 
                                         uname = u,
                                         segmentationPar = segmentationPar, 
                                         access = access), 
              silent = TRUE)
      if (class(z) == "try-error") 
        cat(z, " KO\n")
    }
    return(invisible())
  }
  p = readHTS(x, type = "file", filename = segmentationPar, 
              access = access, format = "dcf")
  cal = try(readHTS(x, "cal", uname = uname, access = access))
  seg = try(readHTS(x, "seg", uname = uname, access = access))
  if (class(cal) != "try-error" & class(seg) != "try-error") {  
    hseg = highlightSegmentation(cal, seg$nseg, seg$cseg)
    ff = fileHTS(x, "viewthumb", uname = uname, 
                 createPath = TRUE, access = "local")
    ff = gsub("_thumb.jpeg", "_seg_thumb.jpeg", ff)
    frame = getFrame(hseg, 1, type = "render")
    imageHTS:::writeThumbnail(x, uname = uname, p = p, input.image = frame, 
                              output.filename=ff)
  }  
}


## plot residuals vs. main effect and compare replicates
plotResiduals <- function(d, ftr="n", ...){
  ftr = match(ftr, d$anno$ftr)
  if(is.na(ftr)) stop("unknown feature selected.")

  forPlot = d$res[,,,ftr]
  dim(forPlot) = c(prod(dim(forPlot)[1:2]), dim(forPlot)[3])

  plot(forPlot[,2] ~ forPlot[,1], 
       pch=".", 
       col = as.factor(rep(d$anno$line$name, dim(d$res)[1])),
       ...)
  abline(a=0,b=1)
}

## plots histograms of pvalues and pAdjusteds
plotPvalues <- function(d, ftr="n", ...){
  ftr = match(ftr, d$anno$ftr)
  if(is.na(ftr)) stop("unknown feature selected.")
  hist(d$pVal[,,ftr,1], xlab="p-value", main=NULL, ...)
  hist(d$pVal[,,ftr,2], xlab="adjusted p-value", main=NULL, breaks=20, ...)
}

## plot heatmap with the given cols, ids as lables and filters
plotHeatmap <- function(d, ftr="n", pAdjustedThresh=NULL, ...){
  ## select feature
  ftr = match(ftr, d$anno$ftr)
  if(is.na(ftr)) stop("unknown feature selected.")
  forHeatmap = d$res[,,,ftr]
  dim(forHeatmap) = c(prod(dim(forHeatmap)[1]), prod(dim(forHeatmap)[2:3]))

  colnames(forHeatmap) <- c(paste(rownames(d$anno$line), "_r1", sep=""),
                            paste(rownames(d$anno$line), "_r2", sep=""))
  rownames(forHeatmap) = d$anno$drug$GeneID

  ## order by colnames for dendogram
  forHeatmap = forHeatmap[, order(colnames(forHeatmap))]

  c("PAR001", "PAR007", "104-007", "02-004", "02-006", "104-004", 
  "104-008", "02-030", "02-031", "104-009", "104-001", "02-008")
  pathwayNumber = c(5, 9, 3, 1, 8, 10, 12, 2, 4, 6, 7, 11)
  pathwayOrder <- paste(rep(pathwayNumber, each=2), c("r1", "r2"), sep="_")

  forHeatmap <- forHeatmap[,pathwayOrder]

  if(!is.null(pAdjustedThresh))
     filter = apply(d$pVal[,,ftr,2] < pAdjustedThresh, 1, any)
  else
     filter=TRUE
     
  heatmap.2(forHeatmap[filter,],
            symbreaks=TRUE,
            col=colorRampPalette(c("yellow", "black", "cornflowerblue"))(64),
            las=2,
            trace="none",
            ...)

  meanVal = apply(d$res[,,,ftr], c(1,2), mean)
  colnames(meanVal) <- rownames(d$anno$line)
  rownames(meanVal) = d$anno$drug$GeneID
  heatmap.2(meanVal[filter,as.character(pathwayNumber)],
            symbreaks=TRUE,
            col=colorRampPalette(c("yellow", "black", "cornflowerblue"))(64),
            las=2,
            trace="none",
            ...)
}

plotEffects <- function(d, ftr="n", ...){
  ftr = match(ftr, d$anno$ftr)
  if(is.na(ftr)) stop("unknown feature selected.")

  ## get data
  drug = d$effect$drug[,,ftr]
  line = d$effect$line[,,ftr]
  anno = d$anno$drug$Content

  plot(line[,1], 
       main = d$anno$ftr[ftr], 
       ylim=c(min(line, na.rm=TRUE) - .1, max(line+1, na.rm=TRUE)),
       xlab="cell line",
       ylab="cell line effect")
  points(line[,2], col=2)

  plot(line[,1],
       line[,2],
       main = paste(d$anno$ftr[ftr], "cell line effect"),
       ylim=c(min(line, na.rm=TRUE) - .1, max(line+1, na.rm=TRUE)),
       xlab="rep 1",
       ylab="rep 2")
  abline(a=0,b=1) 

  plot(seq_len(dim(drug)[1])[anno=="sample"], 
       drug[anno=="sample",1],  
       xlim = c(0, dim(drug)[1] + 1),
       ylim = c(min(drug, na.rm=TRUE) - .1, max(drug, na.rm=TRUE) + .1), 
       main = d$anno$ftr[ftr], 
       xlab="drug",
       ylab="drug effect")
  points(seq_len(dim(drug)[1])[! anno %in% c("sample", "empty")], 
         drug[! anno %in% c("sample", "empty"),1], col=2)

  plot(seq_len(dim(drug)[1])[anno=="sample"], 
       drug[anno=="sample",2],  
       xlim = c(0, dim(drug)[1] + 1),
       ylim = c(min(drug, na.rm=TRUE) - .1, max(drug, na.rm=TRUE) + .1), 
       main = d$anno$ftr[ftr], 
       xlab="drug",
       ylab="drug effect")
  points(seq_len(dim(drug)[1])[! anno %in% c("sample", "empty")], 
         drug[! anno %in% c("sample", "empty"),2], 
         col=2)

  plot(drug[,1],
       drug[,2],
       main = paste(d$anno$ftr[ftr], "drug effect"),
       xlab="rep 1",
       ylab="rep 2")
  abline(a=0,b=1)
  
}

plotResVsEffect <- function(d, ftr="n", ...){
  ftr = match(ftr, d$anno$ftr)
  if(is.na(ftr)) stop("unknown feature selected.")

  ## get data
  drug = d$effect$drug[,,ftr]
  res = d$res[,,,ftr]

  plot(drug[,1], 
       res[,1,1], 
       pch=20, 
       xlim = c(min(drug, na.rm=TRUE) - .1, max(drug, na.rm=TRUE) + .1),
       ylim = c(min(res, na.rm=TRUE) - .1, max(res, na.rm=TRUE) + .1), 
       main = paste(d$anno$ftr[ftr], "rep 1"))
  for(i in 2:dim(res)[2]){
    points(drug[,1], res[,i,1], pch=20, col=i)
  }

  plot(drug[,2], 
       res[,1,2], 
       pch=20, 
       xlim = c(min(drug, na.rm=TRUE) - .1, max(drug, na.rm=TRUE) + .1),
       ylim = c(min(res, na.rm=TRUE) - .1, max(res, na.rm=TRUE) + .1),  
       main = paste(d$anno$ftr[ftr], "rep 2"))
  for(i in 2:dim(res)[2]){
    points(drug[,2], res[,i,2], pch=20, col=i)
  }
}

plotResVsValue <- function(d, ftr="n", ...){
  ftr = match(ftr, d$anno$ftr)
  if(is.na(ftr)) stop("unknown feature selected.")

  ## get and fromat data
  resForPlot = d$res[,,,ftr]
  dim(resForPlot) = c(prod(dim(resForPlot)[1:2]), dim(resForPlot)[3])
  valueForPlot = d$D[,,,ftr]
  dim(valueForPlot) = c(prod(dim(valueForPlot)[1:2]), dim(valueForPlot)[3])
  
  plot(valueForPlot[,1], 
       resForPlot[,1],
       xlab="normalized value",
       ylab="interaction score",
       main=paste(d$anno$ftr[ftr], "rep 1"))

  plot(valueForPlot[,2], 
       resForPlot[,2],
       xlab="normalized value",
       ylab="interaction score",
       main=paste(d$anno$ftr[ftr], "rep 2"))
}


plotAll <- function(d, 
                    type, 
                    localPath=NULL, 
                    prefix=NULL, 
                    pAdjustedThresh=NULL,
                    ...){
switch(type,
       "residuals" = {
         for(ftr in d$anno$ftr){
           plotf <- file.path(localPath, 
                              "plots", 
                              paste(prefix, 
                                    ftr, 
                                    "_residuals_%02d.png",
                                    sep=""))

           png(plotf, 
               width=1024, 
               height=1024, 
               pointsize=40)

           plotResiduals(d, ftr=ftr, ...)

           dev.off()
         }
       },
       "heatmap" = {
         for(ftr in d$anno$ftr){
           plotf <- file.path(localPath, 
                              "plots", 
                              paste(prefix, 
                                    ftr,
                                    "_heatmap_%02d.png",
                                    sep=""))

           png(plotf, 
               width=1024, 
               height=3/2*1024, 
               pointsize=25)

           plotHeatmap(d, ftr=ftr, ...)
           if(!is.null(pAdjustedThresh)){
             if(any(d$pVal[,,match(ftr, d$anno$ftr),2] < pAdjustedThresh)) 
               print(plotHeatmap(d, ftr=ftr, pAdjustedThresh=pAdjustedThresh))
           }
           dev.off()

           plotf <- file.path(localPath, 
                              "plots", 
                              paste(prefix, 
                                    ftr,
                                    "_heatmap_ordered_%02d.png",
                                    sep=""))
           png(plotf, 
               width=1024, 
               height=3/2*1024, 
               pointsize=25)

           plotHeatmap(d, ftr=ftr, Colv=FALSE, ...)
           if(!is.null(pAdjustedThresh)){
             if(any(d$pVal[,,match(ftr, d$anno$ftr),2] < pAdjustedThresh)) 
               print(plotHeatmap(d, 
                                 ftr=ftr,
                                 Colv=FALSE,
                                 pAdjustedThresh=pAdjustedThresh))
           }
           dev.off()
         }
       },
       "pval" = {
         for(ftr in d$anno$ftr){
           plotf <- file.path(localPath, 
                              "plots", 
                              paste(prefix, 
                                    ftr,
                                    "_pval_%02d.png",
                                    sep=""))

           
           png(plotf, 
               width=1024, 
               height=1024, 
               pointsize=40)

           plotPvalues(d, ftr=ftr, ...)

           dev.off()
         }
       },
       "effects" = {
         for(ftr in d$anno$ftr){
           plotf <- file.path(localPath, 
                              "plots", 
                              paste(prefix, 
                                    ftr,
                                    "_effects_%02d.png",
                                    sep=""))

           
           png(plotf, 
               width=1024, 
               height=1024, 
               pointsize=40)

           plotEffects(d, ftr=ftr, ...)

           dev.off()
         }
       },
       "resVsEffect" = {
         for(ftr in d$anno$ftr){
           plotf <- file.path(localPath, 
                              "plots", 
                              paste(prefix, 
                                    ftr,
                                    "_interactionScoreVsEffect_%02d.png",
                                    sep=""))

           
           png(plotf, 
               width=1024, 
               height=1024, 
               pointsize=40)

           plotResVsEffect(d, ftr=ftr, ...)

           dev.off()
         }
       },
       "resVsValue" = {
         for(ftr in d$anno$ftr){
           plotf <- file.path(localPath, 
                              "plots", 
                              paste(prefix, 
                                    ftr,
                                    "_interactionScoreVsValue_%02d.png",
                                    sep=""))
 
           
           png(plotf, 
               width=1024, 
               height=1024, 
               pointsize=40)

           plotResVsValue(d, ftr=ftr, ...)

           dev.off()
         }
       },
       "all" = {
         plotAll(d, localPath=localPath, 
                 type="residuals", prefix, ...)
         plotAll(d, localPath=localPath,  
                 type="heatmap", prefix, pAdjustedThresh, ...)
         plotAll(d, localPath=localPath, 
                 type="pval", prefix, ...)
         plotAll(d, localPath=localPath, 
                 type="effects", prefix, ...)
         plotAll(d, localPath=localPath,
                 type="resVsEffect", prefix, ...)
         plotAll(d, localPath=localPath,
                 type="resVsValue", prefix, ...)         
       })
}

topX = function(d, ftr, x, filterCellno, minFDR = 0.2){

  cat(ftr, ":")
  ftr = match(ftr, d$anno$ftr)
  pAdjusted = d$pVal[,,ftr,2]
  top = 
    pAdjusted <= pAdjusted[!filterCellno][order(pAdjusted[!filterCellno])][x] & 
    pAdjusted < minFDR
#  top = c(top,top)
#  dim(top) = dim(res$res)[1:3]

  if ( sum(top) > 0){
  r1 = d$res[,,1,ftr][which(top)]
  r2 = d$res[,,2,ftr][which(top)]

  unames <- data.frame(well = rep(d$anno$drug$Well, 
                                  length(d$anno$line$startPlate)), 
                       plate= rep(as.numeric(gsub("P", 
                                                  "", 
                                                  d$anno$drug$PlateName)),
                                  length(d$anno$line$startPlate)),
                       startPlate = rep(d$anno$line$startPlate, 
                                        each=nrow(d$anno$drug)))
  unames$uname = paste(sprintf("%03d", 
                               unames$plate + unames$startPlate - 1),
                       unames$well, 
                       sep="-")

  unames = as.matrix(unames$uname)
  dim(unames) = dim(top) 

  uname = unames[which(top)]
  
  ## add annotation  
  selTreatment = which(top) %% dim(top)[1]
  selTreatment[selTreatment == 0] = dim(top)[1]
  drug = d$anno$drug[selTreatment,]
  line = d$anno$line[which(top) %/% dim(top)[1]+1,]

  topHits = data.frame(ftr = d$anno$ftr[ftr],
    uname = gsub("-", "-01-", uname),
    compoundID = drug$compoundID,
    r1,
    r2,
    pAdjusted = pAdjusted[which(top)])

  topHits = cbind(topHits, line, drug[, -match("compoundID", names(drug))])
  topHits = topHits[order(topHits$pAdjusted),]
  rownames(topHits) = 1:nrow(topHits)

  cat(" OK\n")

  } else {
    topHits = NULL
}
  topHits
}

##-------------------------------------------------------------------------
## selectByStability
## select features itteratively based on the correlation of residuals after lm 
## modelling using the previously selected features
##-------------------------------------------------------------------------
selectByStability <- function(subsample,
                               preselect,
                               Rdim = 40, 
                               verbose = TRUE) {

  n = 1
  I = phenotype = 4
  m = 3

  D = subsample$D
  phenotype = subsample$phenotype
  
  Sel = match(preselect, phenotype)
  d = dim(D)
  print(dim(D))
  dimnames(D) = list(NULL,NULL,phenotype)
  correlation = rep(NA_real_, Rdim)
  correlationAll = list()
  ratioPositive = rep(NA_real_, Rdim)
  selected = rep("", Rdim)
  DSel = D[,1,]
  DSel[] = NA_real_
  Dall = D
  for (k in 1:Rdim) {
#     D = D[,,-Sel]
#     D2 = apply(DSel,c(1,3),mean,na.rm=TRUE)
    if (k > 1) {
      for (i in 2:dim(D)[3]) {
        cat("k=",k," i=",i,"\r")
        model = lm(as.vector(D[,1,i]) ~ DSel[,1:(k-1),drop=FALSE]+0)
        D[,1,i] = model$residuals
        model = lm(as.vector(D[,2,i]) ~ DSel[,1:(k-1),drop=FALSE]+0)
        D[,2,i] = model$residuals
      }
    }
    C = apply(D, 3, function(x) {
      a = x[,1]
      b = x[,2]
      I = which(is.finite(a) & is.finite(b))
      cor(a[I],b[I]) 
      } )
    if (k > length(preselect)) {
      I = names(C)[which.max(C)]
    } else {
      I = preselect[k]
    }
    correlation[k] = C[I]
    ratioPositive[k] = sum(C > 0, na.rm=TRUE) / length(C)
    selected[k] = I
    cat("k=",k," selected = ",selected[k], " cor = ", correlation[k], 
        " r = ", ratioPositive[k], "\n")
    correlationAll[[k]] = C
    DSel[,k] = apply(D[,,I,drop=FALSE],1,mean,na.rm=TRUE)
    D = D[,,dimnames(D)[[3]] != I,drop=FALSE]
  }

  res = list(selected = selected, correlation = correlation, 
             ratioPositive = ratioPositive, correlationAll = correlationAll)
  res
}


##-------------------------------------------------------------------------
## toRaster
## x a numeric matrix, returns a matrix of same shape with colour strings
##-------------------------------------------------------------------------
toRaster = function(x, cuts, col) {
  cux =cut(x,cuts,include.lowest = TRUE, labels=FALSE)
  rv = x
  rv[] = col[cux]
  return(rv)
}

##--------------------------------------------------
## transform 3D array x into matrix
##--------------------------------------------------
toMatrix = function(x) {
  dim(x) = c(dim(x)[1], prod(dim(x)[-1]))
  return(x)
}


##--------------------------------------------------------------
## determine 'hclust ordering' of one dimension i of array x
##--------------------------------------------------------------
orderDim = function(x, i) {
   px = toMatrix(aperm(x, c(i, setdiff(seq(along=dim(x)), i))))
   theCorr = cor(t(px), use = "pairwise.complete.obs")
   theCorr[!is.finite(theCorr)] = 0
   hc = hclust(as.dist(trsf(theCorr)))
   return(hc$order)
}

##--------------------------------------------------------------
## transform a correlation coefficient (in [-1,1]) to a distance
##--------------------------------------------------------------
trsf = function(x) {
    exp(-x)-exp(-1)
}

##----------------------------------------------------------
## This function plots dim(x)[3] heatmaps of x[,,i].
## x is a 3D array (dimensions already sorted for plotting)
##----------------------------------------------------------
myHeatmap = function(x, cuts, col,
  fontsize=18, colnames=TRUE, rownames=FALSE) {
    stopifnot(is.array(x), length(dim(x))==3)
    use = 0.99

    rx = toRaster(x, cuts=cuts, col=col)
    pushViewport(viewport(xscale = c(if(rownames) -1 else 0.5, dim(x)[3]+0.5),
                          yscale = c(0, if(colnames) 1/0.8 else 1),
                          y=unit(use/2, "npc"), height=unit(use, "npc"),
                          x=unit(use/2, "npc"), width=unit(use, "npc")))
    colnm = dimnames(x)[[3]]
    if(colnames) {
        stopifnot(!is.null(colnm))
        colnm = strsplit(colnm, split="|", fixed=TRUE)
    }
    rownm = dimnames(x)[[1]]
    if(rownames)
        stopifnot(!is.null(rownm))

    for(i3 in seq_len(dim(x)[3])) {
        grid.raster(rx[,,i3], x=unit(i3,  "native"), width =unit(0.9, "native"),
                              y=unit(0.5, "native"), height=unit(1,   "native"),
                              interpolate=FALSE)

        if(colnames) {
          nLines = length(colnm[[i3]])
          for(k in seq_len(nLines))
            grid.text(paste(" ", colnm[[i3]][k], sep=""),
                      x=unit(i3, "native") + 
                        unit( (if(nLines>1) {(k-1)/(nLines-1)} else 0) -0.5, 
                              "char"),
                      y=unit( 1, "native"),
                      hjust = 0,
                      vjust = 0.5, rot=90,
                      gp=gpar(fontsize=fontsize))
        }
        if(rownames) {
           grid.text(paste(rownm, " ", sep=""),
                     x=unit(0.5, "native"),
                     y=unit(seq(along=rownm)/length(rownm), "native"),
                     hjust = 1,
                     vjust = 1,
                     gp=gpar(fontsize=fontsize))
        }
      }

    popViewport()

}


##----------------------------------------------------------
## This function plots dim(x)[3] heatmaps of x[,,i] for specific i
## x is a 3D array (dimensions already sorted for plotting)
##----------------------------------------------------------
myHeatmapSingle = function(x, cuts, col, i,
  fontsize=10, colnames=TRUE, rownames=TRUE) {
    stopifnot(is.array(x), length(dim(x))==3)
    use = 0.99

    rx = toRaster(x, cuts=cuts, col=col)
    pushViewport(viewport(xscale = c(if(rownames) -1 else 0.5, length(i)+0.5),
                          yscale = c(0, if(colnames) 1/0.8 else 1),
                          y=unit(use/2, "npc"), height=unit(use, "npc"),
                          x=unit(use/2, "npc"), width=unit(use, "npc")))
    colnm = dimnames(x)[[3]]
    if(colnames) {
        stopifnot(!is.null(colnm))
        colnm = strsplit(colnm, split="|", fixed=TRUE)
    }
    rownm = dimnames(x)[[1]]
    if(rownames)
        stopifnot(!is.null(rownm))

    for(i3 in seq_len(i)) {
        grid.raster(rx[,,i[i3]],
                    x=unit(i3,  "native"),
                    width = unit(0.9, "native"),
                    y=unit(0.5, "native"),
                    height=unit(1,   "native"),
                    interpolate=FALSE)

        if(colnames) {
          nLines = length(colnm[[i[i3]]])
          for(k in seq_len(nLines))
            grid.text(paste(" ", colnm[[i[i3]]][k], sep=""),
                      x=unit(i3, "native") +
                        unit( (if(nLines>1) {(k-1)/(nLines-1)} else 0) -0.5, 
                              "char"),
                      y=unit( 1, "native"),
                      hjust = 0,
                      vjust = 0.5, rot=90,
                      gp=gpar(fontsize=fontsize))
        }
        if(rownames) {
           grid.text(paste(rownm, " ", sep=""),
                     x=unit(0.5, "native"),
                     y=unit(seq(along=rownm)/length(rownm), "native"),
                     hjust = 1,
                     vjust = 1,
                     gp=gpar(fontsize=5))
        }
      }

    popViewport()

}


## functions for annotating feature names with human readable names
hrNames = function (names){
  humanReadableNames = 
    c("n" = "Cell number", 
      "cseg.act.m.majoraxis.mean" = "Cell major axis", 
      "cseg.act.h.cor.s1.mean" = "Actin Haralick texture (1)", 
      "nseg.0.m.majoraxis.mean" = "Nuclear major axis", 
      "cseg.dnaact.m.eccentricity.sd" = "Nuclear-Actin eccentricity SD", 
      "cseg.act.h.f12.s2.sd" = "Actin Haralick texture SD (1)", 
      "nseg.dna.h.var.s2.mean" = "Nuclear Haralick texture (1)", 
      "nseg.0.m.eccentricity.mean" = "Nuclear eccentricity", 
      "cseg.0.s.radius.min.qt.0.05" = "5% quantile of cell radius",
      "cseg.dnaact.b.mean.qt.0.05" = "5% quantile of Nuclear-Actin intensity", 
      "cseg.act.m.eccentricity.mean" = "Actin eccentricity",
      "nseg.dna.m.eccentricity.sd" = "Nuclear eccentricity SD",
      "nseg.dna.h.idm.s1.sd" = "Nuclear Haralick texture SD (1)",
      "cseg.dnaact.h.f13.s1.mean" = "Nuclear-Actin Haralick texture (1)",
      "cseg.act.h.asm.s2.mean" = "Actin Haralick texture (2)",
      "cseg.dnaact.h.den.s2.sd" = "Nuclear-Actin Haralick texture SD (1)",
      "nseg.dna.h.cor.s2.sd" = "Nuclear Haralick texture SD (2)",
      "cseg.dnaact.b.mad.mean" = "Nuclear-Actin intensity MAD",
      "cseg.act.h.idm.s2.sd" = "Actin Haralick texture SD (2)",
      "nseg.0.s.radius.max.qt.0.05" = "5% quantile of Nuclear radius",
      "nseg.0.s.radius.min.qt.0.99" = "99% quantile of Nuclear radius",
      "lcd.10NN.sd" = "Local cell density SD",
      "cseg.act.m.majoraxis.qt.0.99" = "99% quantile of cell major axis",
      "nseg.dna.h.var.s2.sd" = "Nuclear Haralick texture SD (3)",
      "nseg.0.s.area.mean" = "Nuclear area",
      "cseg.dnaact.m.majoraxis.mean" = "Nuclear-Actin major axis",
      "cseg.act.h.cor.s2.sd" = "Actin Haralick texture SD (3)",
      "nseg.dna.h.cor.s2.mean" = "Nuclear Haralick texture (2)",
      "cseg.dnaact.b.mean.qt.0.01" = "1% quantile of Nuclear-Actin intensity",
      "cseg.dnaact.h.asm.s1.mean" = "Nuclear-Actin Haralick texture (2)",
      "nseg.0.s.area.qt.0.95" = "95% quantile of nuclear area",
      "cseg.dnaact.m.majoraxis.qt.0.95" = 
          "95% quantile of Nuclear-Actin major axis",
      "cseg.dnaact.h.f12.s1.sd" = "Nuclear-Actin Haralick texture SD (2)",
      "cseg.dnaact.h.con.s1.sd" = "Nuclear-Actin Haralick texture SD (3)",
      "cseg.dnaact.m.eccentricity.qt.0.99" = 
          "99% quantile of Nuclear-Actin eccentricity",
      "cseg.0.s.radius.min.qt.0.99" = "99% quantile of cell radius",
      "nseg.dna.b.mean.qt.0.01" = "1% quantile of nuclear intensity",
      "cseg.act.m.eccentricity.sd" = "Actin eccentricity SD",
      "cseg.dnaact.h.sav.s2.sd" = "Nuclear-Actin Haralick texture SD (4)",
      "cseg.act.m.majoraxis.qt.0.01" = "1% quantile of cell major axis")
  res = names
  I = match(names, names(humanReadableNames))
  res[!is.na(I)] = humanReadableNames[I[!is.na(I)]]
  res
}


## functions for annotating feature names with human readable feature classes
hrClass <- function(classes){
  humanClass = c(`n` = "cell number",
                 `nseg.dna.h.cor.s2.sd` = "nuclear texture",
                 `nseg.dna.h.idm.s1.sd` = "nuclear texture",
                 `nseg.dna.h.var.s2.mean` = "nuclear texture",
                 `nseg.dna.m.eccentricity.sd` = "nuclear shape",
                 `nseg.0.m.majoraxis.mean` = "nuclear shape",
                 `nseg.0.m.eccentricity.mean` = "nuclear shape",
                 `nseg.0.s.radius.max.qt.0.05` = "nuclear shape",
                 `cseg.act.m.majoraxis.mean` = "cell shape",
                 `cseg.act.m.eccentricity.mean` = "cell shape",
                 `cseg.dnaact.m.eccentricity.sd` = "cell shape",
                 `cseg.0.s.radius.min.qt.0.05` = "cell shape",
                 `cseg.act.h.idm.s2.sd` = "cellular texture",
                 `cseg.dnaact.h.f13.s1.mean` = "cellular texture",
                 `cseg.act.h.cor.s1.mean` = "cellular texture",
                 `cseg.dnaact.b.mean.qt.0.05` = "cellular texture",
                 `cseg.dnaact.h.den.s2.sd` = "cellular texture",
                 `cseg.dnaact.b.mad.mean` = "cellular texture",
                 `cseg.act.h.asm.s2.mean` = "cellular texture",
                 `cseg.act.h.f12.s2.sd` = "cellular texture")
  res = classes
  I = match(classes, names(humanClass))
  res[!is.na(I)] = humanClass[I[!is.na(I)]]
  res
}


## plot single values drug combinations
plotSingleValues <- function(df, pthresh){
  levels(df$color) <- gsub("white", "black", levels(df$color))
  p = ggplot(data=df, aes(x=concdrug1.micro_M., 
                          y= value, 
                          color=identifier, 
                          shape = variable))
  p = p + #geom_point(size=3) + 
    scale_x_log10() + 
    scale_y_continuous(limits = c(0.45,1.1)) +#scale_y_log10() +
    theme_bw()  +
    scale_colour_manual(values = levels(df$color), 
                        limits=levels(df$identifier)) +
    labs(title = gsub("_", " ", gsub("2013-10-22platelist_", "", experiment))) +
    geom_jitter(position = position_jitter(width = .04, height=0), 
                size=2.5,
                na.rm=TRUE)
  if(any(df$pvalue < pthresh)) p = p +
    geom_vline(xintercept = df$concdrug1.micro_M.[df$pvalue < pthresh], 
               size=12, 
               alpha=.1, 
               colour="black")
  
  print(p)
}

## summary plot of drug combinations
plotSummary <- function(summary, concentrations, pthresh, ...){
  levels(summary$color) <- gsub("white", "black", levels(summary$color))
  if(!missing(concentrations)) 
    summary = subset(summary, concdrug1.micro_M. %in% concentrations)
  
  theme_new = theme_new + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p = ggplot(data=summary, 
             aes(x=concdrug1.micro_M., 
                 y= mean,
                 ymin=mean-2*sem, 
                 ymax=mean+2*sem,
                 color=identifier,
                 fill=identifier))
  
  p = p + geom_point(size=3) + 
    geom_line() +
    scale_x_log10(...) +
    theme_new +
    scale_y_continuous(limits = c(0.35,1.1)) +
    scale_colour_manual(values = levels(summary$color),
                        limits=levels(summary$identifier)) + 
    geom_errorbar(aes(ymin=mean-2*sem, 
                      ymax=mean+2*sem),
                  width=.065)  +
    labs(title = gsub("_", " ", gsub("2013-10-22platelist_", "", experiment)))
  if(any(summary$pvalue < pthresh)){
    p = p + 
      geom_vline(xintercept = 
                     summary$concdrug1.micro_M.[summary$pvalue < pthresh],
                 size=12, 
                 alpha=.1, 
                 colour="black")
  }
  print(p)  
}  



## barplot of drug combinations
plotSummaryBarplot <- function(summary, concentrations, pthresh){
  if(!missing(concentrations)) 
    summary = subset(summary, concdrug1.micro_M. %in% concentrations)
  
  theme_new = theme_new + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  summary$max = 1.05
  p = ggplot(data=summary, 
             aes(x=factor(concdrug1.micro_M.), 
                 y= mean,
                 ymin=mean-2*sem, 
                 ymax=mean+2*sem,
                 fill=identifier))
  dodge <- position_dodge()
  
  p = p + geom_bar(stat="identity", position=dodge, colour="black") +
    theme_new + 
    scale_y_continuous(limits = c(-0.05,1.1)) +
    scale_fill_manual(values = levels(summary$color), 
                      limits=levels(summary$identifier)) + 
    scale_colour_manual(values = levels(summary$color), 
                        limits=levels(summary$identifier)) + 
    geom_errorbar(position=dodge, color="black")  +
    labs(title = gsub("_", " ", gsub("2013-10-22platelist_", "", experiment)))+
    xlab("drug concentration")
  if(any(summary$pvalue < pthresh)){
    p = p + geom_point(aes(x=factor(concdrug1.micro_M.), 
                           y=max),
                       data = subset(summary, 
                                     pvalue < pthresh),
                       pch=8,
                       col=1,
                       show_guide = FALSE)
  }  
  print(p)    
}    


## barplot of drug combination inhibition values
plotSummaryBarplotInhibition <- function(summary, 
                                         concentrations, 
                                         pthresh, 
                                         ylimits=c(-0.05,.5)){
  if(!missing(concentrations)) 
    summary = subset(summary, concdrug1.micro_M. %in% concentrations)
  
  theme_new = theme_new + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid = element_blank())
  
  summary$inhibition = 1 - summary$mean
  summary$max = max(summary$inhibition+summary$sem)
  p = ggplot(data=summary, 
             aes(x=factor(concdrug1.micro_M.), 
                 y= inhibition,
                 ymin=inhibition-sem, 
                 ymax=inhibition+sem,
                 fill=identifier))
  dodge <- position_dodge()
  
  p = p + geom_bar(stat="identity", position=dodge, colour="black") +
    theme_new + 
    scale_y_continuous(limits = ylimits) +
    scale_fill_manual(values = levels(summary$color), 
                      limits=levels(summary$identifier)) + 
    scale_colour_manual(values = levels(summary$color),
                        limits=levels(summary$identifier)) + 
    geom_errorbar(position=dodge, color="black")  +
    labs(title = gsub("_", " ", gsub("2013-10-22platelist_", "", experiment))) +
    xlab("concentrations (\u00B5M)") +
    ylab("rel. inhibition of cell growth")
  
  if(any(summary$pvalue < pthresh)){
    p = p + geom_point(aes(x=factor(concdrug1.micro_M.), 
                           y=max),
                       data = subset(summary, 
                                     pvalue < pthresh),
                       pch=8,
                       col=1,
                       show_guide = FALSE)
  }  
  print(p) 
}


## calculate distance and annotate
getDist <- function(M, drugAnno){
  template = dimnames(M)[[1]]
  dim(M) = c(dim(M)[1], prod(dim(M)[-1]))
  
  theCorr = cor(t(M), use = "pairwise.complete.obs")
  theCorr[!is.finite(theCorr)] = 0
  
  dist = 1-theCorr
  
  drugName <- drugAnno$Name[match(template, 
                                  drugAnno$GeneID)]
  
  Selectivity <- drugAnno$Selectivity_updated[match(template, 
                                                    drugAnno$GeneID)]
  
  dimnames(dist) = list(paste(template, drugName), 
                        ifelse(is.na(Selectivity), 
                               template, 
                               paste(template, Selectivity)))
  
  invisible(dist)
}

## calculate corrlation and annotate
getCorr <- function(M, drugAnno){
  template = dimnames(M)[[1]]
  dim(M) = c(dim(M)[1], prod(dim(M)[-1]))
  
  theCorr = cor(t(M), use = "pairwise.complete.obs")
  theCorr[!is.finite(theCorr)] = 0
  
  drugName <- drugAnno$Name[match(template, 
                                  drugAnno$GeneID)]
  
  Selectivity <- drugAnno$Selectivity_updated[match(template, 
                                                    drugAnno$GeneID)]
  
  dimnames(theCorr) = list(template, 
                           ifelse(is.na(Selectivity), 
                                  template, 
                                  paste(template, Selectivity)))
  
  invisible(theCorr)
}

checkConsistency <- function(objectName){
  env = new.env()
  data(list=objectName, package='PGPC', envir=env)
  if(!identical(get(objectName, envir=.GlobalEnv), get(objectName, envir=env)))
    warning(paste(sprintf("The stored and calculated version of '%s'", 
                          objectName),
                  "are not exactly identical.\n",
                  "This might influence the downstream results."))
  stopifnot(isTRUE(all.equal(get(objectName, envir=.GlobalEnv), 
                             get(objectName, envir=env))))
}
