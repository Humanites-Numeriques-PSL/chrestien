# Load and clean a text



# Functions for stylometric analysis

## Distance measures

### Cosine distance

cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

### MinMax metric

MinMax =
  function(x){
    myDist = matrix(nrow = ncol(x),ncol = ncol(x), dimnames = list(colnames(x),colnames(x)))
    for(i in 1:ncol(x)){
      for(j in 1:ncol(x)){
        min = sum(apply(cbind(x[,i],x[,j]), 1, min))
        max = sum(apply(cbind(x[,i],x[,j]), 1, max))
        resultat = 1 - (min / max)
        myDist[i,j] = resultat
      }
    }
    return(myDist)
  }


### Normalisations

relativeFreqs = function(x){
  # Relative frequencies
  for(i in 1:ncol(x)){
    x[,i] = x[,i]/sum(x[,i])
  }
  return(x)
}
# Z-scores
ZTransf = function(x){
  for(i in 1:nrow(x)){
    x[i,] = ( x[i,] - mean(x[i,]) )  / sd(x[i,])
  }
  return(x)
}

normalisations = function(x){
  # Z-transformation  
  x = ZTransf(x)
  
  # Vector length normalisation
  for(i in 1:ncol(x)){
    x[,i] = x[,i] / sqrt(sum(x[,i]^2))
  }
  return(x)
}

### Feature selection (Moisl 2011)

selection = function(x, z = 1.96){
  
  # Conversion to probabilities
  probs = relativeFreqs(x)
  
  # Prepare output
  results = matrix(nrow = nrow(probs), ncol = 4, dimnames = list(rownames(probs), c('freq', 'mean prob', 'sample size necessary', 'passes')))
  results = as.data.frame(results)
  results[,1] = rowSums(x)
  results[,2] = rowMeans(probs)
  
  for (i in 1:nrow(probs)){
    var = probs[i,]
    # hist(probs[i,])
    # Calculates mirror-image to compensate for non normality
    mirror = ( max(var) + min(var) ) - var
    var = c(var,mirror)
    # e as function of sigma
    e = 2 * sd(var)
    results[i,3] = mean(var) * (1 - mean(var)) * (z / e )^2 
  }
  
  # And now, mark as false all the rows that would necessit bigger sample than available
  results[,4] = results[,3] <= min(colSums(x))
  
  return(results)
}

### SOM

library("kohonen")

somCluster = function(x, gridx = 10, gridy = 10){
  maSOM = som(x, grid=somgrid(gridx, gridy,"hexagonal"))
  #plot(maSOM, type="mapping", labels=rownames(maSOM$data))
  results = cluster::agnes(t(as.data.frame(maSOM$codes)), metric  = "manhattan", method = "ward")
  return(results)
}


### Robustness checks

# Gives the cluster purity with reference to alledged authors, and an adjusted Rand index in comparison with the analysis showed in fig. 1
robustnessChecks = function(data, refCAH, k = "10"){
  # Get classes from the reference CAH
  refCAHClasses = cutree(refCAH, k = k)
  # Prepare results
  results = matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("N", "CPAuteurs", "CPReference")))
  for (i in list(
    seq(0, 1, 0.001),
    seq(0, 1, 0.01),
    seq(0, 1, 0.1),
    seq(0, 1, 0.25),
    seq(0, 1, 0.5),
    seq(0, 0.5, 0.25),
    seq(0, 1, 1)) ) {
    # First, get the cutoffs: first 1000-quantile, first 100-quantile, first decile, all
    selec = quantile(rowSums(data), probs = i)
    selec = selec[length(selec) - 1]
    
    myData = data[rowSums(data) >= selec, , drop = FALSE]
    myData = normalisations(myData)
    myCAH = cluster::agnes(t(myData), metric = "manhattan", method="ward")
    
    # Classes as per alledged author
    expected = sub("_.*", "", rownames(myCAH$data))
    # Cluster purity
    classes = cutree(myCAH, k = k)
    
    N = nrow(myData)
    purity = NMF::purity(as.factor(classes), expected)
    #NMF::entropy(as.factor(classes), expected)
    purityRef = NMF::purity(as.factor(classes), as.factor(refCAHClasses))
    #Rand = mclust::adjustedRandIndex(classes, refCAHClasses)
    
    MF = paste(round(100 - as.numeric(sub("%", "", names(selec))), digits = 2), "%", sep = "")
    
    #localRes = matrix(data = c(N, purity, purityRef, Rand), nrow = 1, ncol = 4, dimnames = list(MF, NULL))
    localRes = matrix(data = c(N, purity, purityRef), nrow = 1, ncol = 3, dimnames = list(MF, NULL))
    results = rbind(results, localRes)
  }
  return(results)
}

# Gives the cluster purity with reference to alledged authors, and an adjusted Rand index in comparison with the analysis showed in fig. 1
compareHC = function(cahList, k = "10"){
  cahList = lapply(cahList, as.hclust)
  # Prepare results
  results = matrix(ncol = length(cahList), nrow = length(cahList), dimnames = list(labels(cahList), labels(cahList)))
  for (i in 1:length(cahList)){
    for (j in 1:length(cahList)){
      classes1 = cutree(cahList[[i]], k = k)
      classes2 = cutree(cahList[[j]], k = k)
      results[i, j] = purityRef = NMF::purity(as.factor(classes1), as.factor(classes2))
      #Rand = mclust::adjustedRandIndex(classes, refCAHClasses)
    }
  }
  
  return(results)
}

### Gives, for each HC, three indicators
# - Purity with regard to Meyer Leg. A/B/C and Wauchier
# - Quality of clustering with HC
# - Stability: ARI with regard to a reference HC
benchmarkHC = function(refCAH, cahList, k = 10){
  # Prepare results
  results = matrix(ncol = 4, nrow = length(cahList), 
                   dimnames = list(labels(cahList), c("N", "AC", "CPMeyer", "CPREF")))
  for (i in 1:length(cahList)){
    N = ncol(cahList[[i]]$data)
    # 1. Purity with regard to Wauchier/non-Wauchier
    classes1 = cutree(as.hclust(cahList[[i]]), k = k)
    expected = classes1
    # Now we set 1 for Wauchier, 2 for NOT Wauchier
    
    expected[grep("Leg.A", rownames(cahList[[i]]$data))] = "Leg-A"
    expected[grep("Leg.B", rownames(cahList[[i]]$data))] = "Leg-B"
    expected[grep("Leg.C", rownames(cahList[[i]]$data))] = "Leg-C"
    expected[grep("_Wau_", rownames(cahList[[i]]$data))] = "Wauchier"
    CPMeyer = NMF::purity(as.factor(classes1), as.factor(expected))
    
    # 2. Quality of clustering with HC
    AC = cahList[[i]]$ac
    
    # 3. Stability: CP with regard to a reference HC
    
    classes2 = cutree(as.hclust(refCAH), k = k)
    CPRef = NMF::purity(as.factor(classes1), as.factor(classes2))
    
    results[i, ] = c(N, AC, CPMeyer, CPRef)
    
  }
  return(results)
}

volatility = function(cahList, k = 9){
  textsLabels = cahList[[1]]$order.lab
  results = matrix(nrow = length(textsLabels), ncol = 1, 
                   dimnames = list(textsLabels, "V_i"))
  
  classes = lapply(cahList, as.hclust)
  classes = lapply(classes, cutree, k = k)
  nclasses = length(classes)
  for (i in 1:length(textsLabels)){
    # Find all X with their freqs. Start with as much row as labels, will remove
    # unused ones later
    X = matrix(ncol = 1, nrow = length(textsLabels), dimnames = list(textsLabels, 'Freq'), data = 0)
    for (j in 1:length(classes)){
      myMembers = labels(classes[[j]][ classes[[j]] ==  classes[[j]][textsLabels[i]]])
      X[myMembers, ] = X[myMembers, ] + 1
    }
    # Remove 0s
    X = X[X > 0, ]
    # Compute index
    V_i = sum( (X - (nclasses - X) )/ nclasses ) / length(X)
    results[i, ] = V_i
  }
  return(results)
}

###########################################
### Compare results with other datasets ###
###########################################

#### First, a function to replicate all other analyses on an alternate dataSet

## subfunction to perform analysis

readData = function(toKeep, path){
  # toKeep: texts to keep for analysis
  # path: path to data
  # return: data…
  data = read.csv(path, header = TRUE, row.names = 1)
  data = t(data)
  data = data[, toKeep]
  data = data[rowSums(data) > 0, ]
  return(data)
}

readAnnotData = function(toKeep, path){
  # toKeep: texts to keep for analysis
  # path: path to data
  # return: data…
  data = read.csv(path, header = TRUE, row.names = 1, sep = ";")
  data = data[, -1]
  colnames(data) = gsub("^X", "", colnames(data))
  colnames(data) = gsub(".decolumnized", "", colnames(data))
  colnames(data) = gsub("Leg.", "Leg-", colnames(data))
  data = data[, toKeep]
  data = data[rowSums(data) > 0, ]
  data = as.matrix(data)
  return(data)
}

performAnalysis = function(data){
  # data: well, … data.
  # returns: a HC (agnes object)
  # Selection based on Moisl 2011
  select = selection(data, z = 1.645)
  select = select[,4]
  # Normalisations
  data = relativeFreqs(data)
  data = data[select,]
  data = normalisations(data)
  myCAH = cluster::agnes(t(data), metric = "manhattan", method="ward")
  return(myCAH)
}

## and the one to replicate

replicateAnalysis = function(toKeep, path_raw_char3grams, 
                             path_expanded_words, path_pos3_gr, path_lemmas, 
                             functionWords, functionLemmas){
  # toKeep: list of texts to keep for analysis
  # path_raw_char3grams: path to raw_char3grams file;
  # path_expanded_words, path_pos3_gr, path_lemmas: you get the idea…
  # functionWords: list of functionWords
  # functionLemmas: well, in fact, a list of functionLemmas (really)
  
  data3grams = readData(toKeep, path_raw_char3grams)
  CAHRaw3gr = performAnalysis(data3grams)
  dataWords = readData(toKeep, path_expanded_words)
  CAHForms  = performAnalysis(dataWords)
  dataAffs = countAffixes(dataWords)
  CAHAffs = performAnalysis(dataAffs)
  
  # Now, function words, a bit more annoying (should put it in readData function, yeah, right)
  dataFW = relativeFreqs(dataWords)
  dataFW = dataFW[functionWords,]
  dataFWsave = dataFW
  dataFW = normalisations(dataFW)
  CAHFW = cluster::agnes(t(dataFW), metric = "manhattan", method="ward")
  
  # back to normal. Or not.
  dataPOS = readAnnotData(toKeep, path_pos3_gr)
  CAHPOS3gr = performAnalysis(dataPOS)
  dataLemmas = readAnnotData(toKeep, path_lemmas)
  CAHLemmas = performAnalysis(dataLemmas)
  
  # And now back to functionLemmas
  dataFL = relativeFreqs(dataLemmas)
  dataFL = dataFL[functionLemmas,]
  dataFL = normalisations(dataFL)
  CAHFL = cluster::agnes(t(dataFL), metric = "manhattan", method="ward")
  
  # And NOOOOOW: Affixes + POS 3-gr + Function words (unnorm)
  
  # Select relevant data
  dataList = list(dataAffs, dataPOS)
  results = matrix(ncol = ncol(dataAffs), nrow = 0, dimnames = list(NULL, colnames(dataAffs)))
  
  for (i in 1:length(dataList)){
    select = selection(dataList[[i]], z = 1.645)
    select = select[,4]
    # Normalisations
    dataList[[i]] = relativeFreqs(dataList[[i]])
    results = rbind(results, dataList[[i]][select,])
  }
  results = rbind(results, dataFWsave)
  rm(dataList)
  dataGlob = normalisations(results)
  CAHGlob2 = cluster::agnes(t(dataGlob), metric = "manhattan", method="ward")
  
  # AND NOW: lemmas and words
  # Select relevant data
  dataList = list(dataLemmas, dataWords)
  results = matrix(ncol = ncol(dataLemmas), nrow = 0, dimnames = list(NULL, colnames(dataLemmas)))
  for (i in 1:length(dataList)){
    select = selection(dataList[[i]], z = 1.645)
    select = select[,4]
    # Normalisations
    dataList[[i]] = relativeFreqs(dataList[[i]])
    results = rbind(results, dataList[[i]][select,])
  }
  rm(dataList)
  dataWordsLemmas = normalisations(results)
  CAHWordsLemmas = cluster::agnes(t(dataWordsLemmas), metric = "manhattan", method="ward")
  
  cahList = list(raw3grams = CAHRaw3gr, Affs = CAHAffs, FunctWords = CAHFW, FunctLemm = CAHFL, POS3gr = CAHPOS3gr, FWPOSandAffs = CAHGlob2, Forms = CAHForms,  Lemmas = CAHLemmas, WordsLemmas = CAHWordsLemmas)
  return(cahList)
}

### And now for something entirely different
### Compare one to one the analyses

compareReplications = function(refCahList, replicatedCahList, k = 5){
  results = matrix(nrow = length(names(refCahList)), ncol = 1,  
                   dimnames = list(names(refCahList), NULL))
  for(i in 1:nrow(results)){
    classes1 = cutree(as.hclust(refCahList[[i]]), k = k)
    classes2 = cutree(as.hclust(replicatedCahList[[i]]), k = k)
    results[i, 1] = NMF::purity(as.factor(classes1), as.factor(classes2))
  }
  return(results)
}


### Plots and layout

#### ACP

library(ggfortify)
pcaPlot = function(x, main ="Plot"){
  dat = x
  auts = vector(length=ncol(dat))
  auts[startsWith(colnames(dat), "CORNEILLEP_")] = 'CORNEILLEP'
  auts[startsWith(colnames(dat), "CORNEILLET_")] = 'CORNEILLET'
  auts[startsWith(colnames(dat), "MOLIERE_")] = 'MOLIERE'
  auts[startsWith(colnames(dat), "ROTROU_")] = 'ROTROU'
  auts[startsWith(colnames(dat), "SCARRON_")] = 'SCARRON'
  auts[startsWith(colnames(dat), "OUVILLE_")] = 'OUVILLE'
  # rename texts
  colnames(dat) = sub('^[^_]+_', '', colnames(dat))
  colnames(dat) = substring(colnames(dat), 1, 10)
  dat = dat[rowSums(dat) != 0,]
  
  ggplot2::autoplot(prcomp(t(dat), scale. = TRUE, center = TRUE), label = TRUE, data = cbind.data.frame(t(dat), auts),  colour='auts', main = main, label.show.legend=FALSE)
}

#### HC

cahPlot = function(x, main="Plot", xlab = paste(ncol(x$data), "features"), k = 6){
  x$order.lab = sub("CORNEILLEP","CP", x$order.lab)
  x$order.lab = sub("CORNEILLET","CT", x$order.lab)
  x$order.lab = sub("MOLIERE","M", x$order.lab)
  x$order.lab = sub("OUVILLE","O", x$order.lab)
  x$order.lab = sub("ROTROU","R", x$order.lab)
  x$order.lab = sub("SCARRON","S", x$order.lab)
  x$order.lab = sub("BOISSY","B", x$order.lab)
  x$order.lab = sub("DANCOURT","DA", x$order.lab)
  x$order.lab = sub("DUFRESNY","DU", x$order.lab)
  x$order.lab = sub("NIVELLE","N", x$order.lab)
  x$order.lab = sub("REGNARD","R", x$order.lab)
  x$order.lab = sub("VOLTAIRE","V", x$order.lab)
  x$order.lab = sub("BOURSAULT","B", x$order.lab)
  x$order.lab = sub("CHEVALIER","C", x$order.lab)
  x$order.lab = sub("DONNEAUDEVISE","DDV", x$order.lab)
  x$order.lab = sub("DORIMOND","DOR", x$order.lab)
  x$order.lab = sub("GILLET","G", x$order.lab)
  x$order.lab = sub("LAFONTAINE","LF", x$order.lab)
  x$order.lab = sub("QUINAULT","Q", x$order.lab)
  #TODO: modif temporaire pour lisibilité
  x$order.lab = substring(x$order.lab, 1, 10)
  plot(x, main=main, xlab=xlab, which.plots = 2)
  myCAH2 = as.hclust(x)
  # Cut in 6 groups
  rect.hclust(myCAH2, k = k, border = 2:5)
}


cahPlotCol = function(x, main="Plot", xlab = paste(ncol(x$data), "features"), k = 3, lth = 7, lrect = -12, cex = 0.6, ylab = "height"){
  # Redefining labels and  Coloring them
  # Redefining labels
  x$order.lab = substring(x$order.lab, 1, 25)
  # Coloring them
  labels = vector(length = length(x$order.lab))
  labels[grep("GerbM_", x$order.lab)] = "darkgreen"
  labels[grep("ChrTr_", x$order.lab)] = "red"
  labels[grep("WauDen_", x$order.lab)] = "pink"
  labels[grep("JnRen_", x$order.lab)] = "orange"
  labels[grep("ANO_", x$order.lab)] = "grey"
  labels[grep("ChrX_", x$order.lab)] = "lightblue"
  xlab = paste(xlab, "|| Agglomerative coefficient = ", round(x$ac, digits = 2))
  
  factoextra::fviz_dend(x, k = k, 
                        k_colors = rep("black", k), 
                        color_labels_by_k = FALSE, 
                        rect = TRUE, 
                        labels_track_height = lth, 
                        label_cols = labels, 
                        cex = cex,
                        lower_rect = lrect,
                        main = main, xlab = xlab, ylab = ylab) + theme(plot.margin = margin(5,15,5,5))
}


#### Boxplots and descriptive statistics

myDescPlot = function(x, type = "boxplot",  main = "", ylab = "freq", xlab = "", withOuville = FALSE){
  names = c('CP', 'CT', 'M','R', 'S')
  if(withOuville){
    names = c(names, 'O')
  }
  CORNEILLEP = x[,grepl('CORNEILLEP', colnames(x))]
  CORNEILLET = x[,grepl('CORNEILLET', colnames(x))]
  MOLIERE = x[,grepl('MOLIERE', colnames(x))]
  ROTROU = x[,grepl('ROTROU', colnames(x))]
  SCARRON = x[,grepl('SCARRON', colnames(x))]
  if(withOuville){
    OUVILLE = x[,grepl('OUVILLE', colnames(x))]
  }
  if('counts' %in% type){
    return(list(CORNEILLEP, CORNEILLET, MOLIERE, ROTROU, SCARRON, OUVILLE))
  }
  if('boxplot' %in% type){ 
    #boxplot
    if(withOuville){
      boxplot(list(CORNEILLEP, CORNEILLET, MOLIERE, ROTROU, SCARRON, OUVILLE), names=names, main=main,ylab=ylab) 
    }
    else{
      boxplot(list(CORNEILLEP, CORNEILLET, MOLIERE, ROTROU, SCARRON), names=names, main=main,ylab=ylab) 
    }
  }
  if('violinplot' %in% type){ 
    #violinplot
    data = cbind(as.data.frame(t(x)), sub("_.*$", "", colnames(x)))
    colnames(data)[2] = "author"
    
    levels(data[,2]) = sub("CORNEILLEP","CP", levels(data[,2]))
    levels(data[,2]) = sub("CORNEILLET","CT", levels(data[,2]))
    levels(data[,2]) = sub("MOLIERE","M", levels(data[,2]))
    levels(data[,2]) = sub("OUVILLE","O", levels(data[,2]))
    levels(data[,2]) = sub("ROTROU","R", levels(data[,2]))
    levels(data[,2]) = sub("SCARRON","S", levels(data[,2]))
    levels(data[,2]) = sub("BOISSY","B", levels(data[,2]))
    levels(data[,2]) = sub("DANCOURT","DA", levels(data[,2]))
    levels(data[,2]) = sub("DUFRESNY","DU", levels(data[,2]))
    levels(data[,2]) = sub("NIVELLE","N", levels(data[,2]))
    levels(data[,2]) = sub("REGNARD","R", levels(data[,2]))
    levels(data[,2]) = sub("VOLTAIRE","V", levels(data[,2]))
    levels(data[,2]) = sub("BOURSAULT","B", levels(data[,2]))
    levels(data[,2]) = sub("CHEVALIER","C", levels(data[,2]))
    levels(data[,2]) = sub("DONNEAUDEVISE","DDV", levels(data[,2]))
    levels(data[,2]) = sub("DORIMOND","DOR", levels(data[,2]))
    levels(data[,2]) = sub("GILLET","G", levels(data[,2]))
    levels(data[,2]) = sub("LAFONTAINE","LF", levels(data[,2]))
    levels(data[,2]) = sub("QUINAULT","Q", levels(data[,2]))
    
    violinplot <- ggplot(data, aes_(x = quote(author), y = as.name(colnames(data)[1]))) +
      ggtitle(main) +
      ylab(ylab) +
      xlab(xlab) +
      geom_violin() + 
      geom_boxplot(width=0.1) +
      theme(axis.text.x = element_text(size = rel(0.7)))
    
    return(violinplot)
    
  }
  if('barplot' %in% type){ 
    data = cbind(as.data.frame(t(x)), sub("_.*$", "", colnames(x)))
    colnames(data)[2] = "author"
    barplot = ggplot(data, aes_(x = quote(author), y = as.name(colnames(data)[1]))) +
      ggtitle(main) +
      ylab(ylab) +
      xlab("") +
      geom_col() # equivalent to geom_bar(stat=identity)
    
    return(barplot)
  }
}

classesDesc = function(x, y, k = "10"){
  # Classes description
  classes = cutree(x, k = k)
  #Add classes to data frame
  dataClassif = t(y)
  dataClassif = cbind(as.data.frame(dataClassif), as.factor(classes))
  colnames(dataClassif[length(dataClassif)])[] = "Classes"
  dataClassif[length(dataClassif)]
  # Desc
  classDesc = FactoMineR::catdes(dataClassif, num.var = length(dataClassif))
  return(classDesc)
}

specifPlot = function(freq_abs, myCAH, k = 5, nfeats = 5){
  # freq_abs: absolute frequencies original data
  # myCAH : the CAH to cut
  # k : the number of classes
  # nfeats: number of positive and negative feats to 
  # TODO: Mieux, utiliser seuil de banalité
  classes = cutree(myCAH, k = k)
  freq_abs_class = matrix(nrow = nrow(freq_abs), ncol = length(unique(classes)), dimnames = list(rownames(freq_abs), unique(classes)))
  for(i in 1:length(unique(classes))){
    # sum the values for each member of the class
    freq_abs_class[, i] = rowSums(freq_abs[, classes == i])
  }
  specs = textometry::specificities(freq_abs_class)
  
  plots = list()
  for(i in 1:ncol(specs)){
    # et maintenant, on peut regarder les spécificités de la classe 1 
    # positives ou négatives
    values = c(head(sort(specs[, i], decreasing = TRUE), n = nfeats), head(sort(specs[, i]), n = nfeats))
    colors = vector(mode="logical", length = length(values))
    colors[values>0] = TRUE
    colors[values<=0] = FALSE
    df = cbind(labels(values), as.data.frame(values), as.data.frame(colors))
    colnames(df)[1] = "labels"
    plots[[i]] = 
      ggplot(df, aes(x = reorder(labels, values), y = values, fill = colors)) +
      geom_col(position = "identity", colour = "black", size = 0.25) +
      scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = FALSE) +
      ggtitle(paste("Specificities for class ", i)) +
      ylab("Specif.") +
      xlab("feats")
    
  }
  myfun <- get("grid.arrange", asNamespace("gridExtra"))
  do.call("myfun", c(plots, ncol=2))
} 


### Data manipulation

aggrAndClean = function(x){
  # aggregate
  aggr = aggregate(x~rownames(x),FUN = sum)
  # cleaning
  x = as.matrix(aggr[,-1])
  rownames(x) = aggr[,1]
  return(x)
}

countAffixes = function(x){
  # Prefixes
  prefs = x
  # Remove words shorter than 3+1 chars
  prefs = prefs[stringr::str_length(rownames(prefs)) > 3,]
  # Extract the first three chars as new rownames
  rownames(prefs) = paste("$", substr(rownames(prefs), 1, 3), sep = "")
  prefs = aggrAndClean(prefs)
  
  # Space prefixes
  spPrefs = x 
  rownames(spPrefs) = paste("_", substr(rownames(spPrefs), 1, 2), sep = "")
  spPrefs = aggrAndClean(spPrefs)
  
  # Suffixes
  sufs = x
  sufs = sufs[stringr::str_length(rownames(sufs)) > 3,]
  rownames(sufs) = paste(stringr::str_sub(rownames(sufs), -3), "^", sep = "")
  sufs = aggrAndClean(sufs)
  
  # Space suffixes
  spSufs = x
  rownames(spSufs) = paste(stringr::str_sub(rownames(spSufs), -2), "_", sep = "")
  spSufs = aggrAndClean(spSufs)
  
  results = rbind(prefs, spPrefs, sufs, spSufs)
  
  return(results)
}
