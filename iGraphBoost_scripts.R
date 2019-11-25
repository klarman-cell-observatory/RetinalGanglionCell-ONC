XgboostAssign = function(train_object, test_object,train_counts=NULL, var.genes = NULL, id_train="typeID", test.label=NULL, do.scale=TRUE, newID = "xgboost", plot.validation = FALSE){
  
  Idents(train_object) <- id_train
  
  # Find variable genes
  if (is.null(var.genes)){
    var.genes0 = GammaPoissonFit(train_counts[,colnames(train_object)],  num.sd = 0.7, x.high.cutoff = 1, x.low.cutoff = 0.001, do.text = FALSE)
    var.genes = ClusterSpecificGenes(train_object, ident_use = id_train, genes_test = var.genes0)
    var.genes = intersect(var.genes, rownames(test_object@assays$RNA))
  }
  
  
  # XGBoost Step
  bst_atlas= XGBoost_train(train_object, var.genes = var.genes, do.scale=do.scale, plot.validation=FALSE)
  thresh.use=0.7
  test_Data <- as.matrix(test_object@assays$RNA[var.genes,])
  if (do.scale) test_Data <- t(scale(t(test_Data)))
  test_Data <- xgb.DMatrix(data = t(test_Data), label = test.label )
  test_pred_bst <- predict(bst_atlas$bst_model, newdata = test_Data)
  numberOfClasses = length(levels(Idents(train_object)))
  test_pred_bst <- matrix(test_pred_bst, nrow= numberOfClasses,
                          ncol=length(test_pred_bst)/numberOfClasses)
  margins_bst = apply(test_pred_bst,2,max)
  transfer.label_bst = apply(test_pred_bst,2,function(x){ 
    if (which.max(x) <= 40){
      if (max(x) > 0.8){ which.max(x) } else {46}
    } else {
      if (max(x) > 0.7){ which.max(x)} else {46}
    }
    
  })
  names(transfer.label_bst) = rownames(test_object@meta.data)
  test_object@meta.data[,newID] = factor(transfer.label_bst[rownames(test_object@meta.data)])
  
  # Graph voting
  Idents(test_object) = newID
  return(test_object)
  
}

XGBoost_train = function(train_object, var.genes = NULL, do.scale=FALSE, max.cells.per.ident = 1000, train.frac = 0.9, weights = TRUE, plot.validation=FALSE){
  
  # Subsetting and Scaling
  train_Data = as.matrix(train_object@assays$RNA[var.genes,])
  scale.mean = rowMeans(train_Data)
  scale.var = apply(train_Data, 1, sd)
  train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.var))
  
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using either ", train.frac*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training, whichever is smaller"))
  
  for (cc in c(1:length(levels(Idents(train_object))))){
    i = levels(Idents(train_object))[cc]
    cells.in.clust = WhichCells(train_object,idents=i);
    
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(cc-1,length(train.temp))); validation.label = c(validation.label, rep(cc-1, length(validation.temp)));
  }
  print(table(training.label))
  train_matrix <- xgb.DMatrix(data = t(train_Data[,training.set]), label=training.label)
  validation_matrix <- xgb.DMatrix(data = t(train_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  print(numberOfClasses)
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,"max_depth"=6, subsample = 0.6)
  nround    <- 200 # number of XGBoost rounds
  if (weights==TRUE){
    freqs = table(training.label)
    weights1 = c()
    minf = min(freqs)
    maxf = max(freqs)
    for (i in 1:length(freqs)){
      f = freqs[i]
      weights1 = c(weights1, rep(0.6 + 0.4*(f-minf)/(maxf-minf), f))
    }
      bst_model <- xgb.train(params = xgb_params,
                             data = train_matrix,
                             nrounds = nround,weight = weights1)
    } else {
      bst_model <- xgb.train(params = xgb_params,
                             data = train_matrix,
                             nrounds = nround)
    }
 
  
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  validation_prediction <- matrix(validation_pred, nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  valid_predlabels=apply(validation_prediction,2,which.max)-1
  
  # Confusion matrix
  if (plot.validation){
    A = table(validation.label, valid_predlabels)
    colnames(A) = levels(Idents(train_object)); rownames(A) = colnames(A)
    plotConfusionMatrix(A, order="Row", xlab.use = "True", ylab.use = "Predicted")
  }
 
  
  to.return = list()
  to.return$bst_model = bst_model
  to.return$scale_mean = scale.mean
  to.return$scale_var = scale.var
  
  return(to.return)
  
}

HarmonizeXgboostRGC = function(test_object, ID1 = "xgboost01", ID2 = "xgboost02", newID="xgboost"){
  
  cell.names = colnames(test_object)
  label1 = test_object@meta.data[,ID1]; names(label1) = cell.names
  label2 = test_object@meta.data[,ID2]; names(label2) = cell.names
  
  new_label = label1
  for (i in c(1:45)){
    cells.use = union(cell.names[label1 == i & label2 == i], 
                      cell.names[label1 == i & label2 == 46])
    cells.use = union(cells.use, cell.names[label2 == i & label1 == 46] )
    new_label[cells.use] = i
  }
  
  test_object@meta.data[,newID] = new_label
  Idents(test_object) = newID
  return(test_object)
}


GraphAssignment = function(object, k=60, knn=15, ref.cells = NULL, target.cells = NULL, null.ident = 46, reduction.use = "tsne"){
  cells.ident0 = as.character(Idents(object)); names(cells.ident0) = names(Idents(object))
  cells.ident = cells.ident0
  cells.ident.factor = factor(cells.ident, levels = levels(Idents(object)))
  
  if (reduction.use == "tsne") k = 2
  
  l=0
  iter=1
  while (l==0){
    print(paste0("Iteration = ", iter))
    if (is.null(target.cells)){
      cells.to.assign = names(cells.ident.factor)[cells.ident.factor == null.ident]
      cells.unassigned = cells.to.assign
    } else {
      cells.unassigned = target.cells
      cells.to.assign = target.cells
    }
    
    if (is.null(ref.cells)){
      ref.cells = setdiff(rownames(object@reductions[[reduction.use]]@cell.embeddings), cells.unassigned)
    }
    
    if (is.null(target.cells)){ 
      frac.unassigned = length(cells.unassigned)/(length(rownames(object@reductions[[reduction.use]]@cell.embeddings)))
      print(paste0(round(frac.unassigned*100), " percent cells are unassigned"))
      if (frac.unassigned < 0.01){
        return(cells.ident.factor)
      }
    }
    
    Xquery = object@reductions[[reduction.use]]@cell.embeddings[cells.to.assign,1:k]
    X = object@reductions[[reduction.use]]@cell.embeddings[ref.cells,1:k]
  
    # Forward
    nearest = get.knnx(X,Xquery,k=knn+1, algorithm = "kd_tree" )
    new.assignments = apply(nearest$nn.index[,1:knn],1,function(x){
      Nvotes = table(cells.ident.factor[ref.cells][x]);
      if (max(Nvotes) > 0.6*knn){
        return(which.max(Nvotes))
      } else {
        return(null.ident)
      }
    })
    names(new.assignments) = cells.unassigned
    
    # # Reciprocal
    # # Mutual nearest neighbor
    # new.assignments0 = new.assignments
    # nearest_inv = get.knnx(Xquery,X,k=knn+1, algorithm = "kd_tree" )
    # new.assignments=sapply(1:length(new.assignments0),function(x){
    #   y= rowSums(nearest_inv$nn.index[,1:knn] == x) > 0
    #   y=which(y)
    #   if (length(y) >= 10){
    #     z=table(cells.ident.factor[ref.cells][y])
    #     if (max(z) > 0.5*sum(z)){
    #       return(which.max(z))
    #     } else {
    #       return(null.ident)
    #     }
    #   } else {
    #     return(null.ident)
    #   }
    # })
    # names(new.assignments) = cells.unassigned
    # new.assignments[new.assignments != new.assignments0] = nullident
    
  
    cells.ident = as.character(cells.ident); names(cells.ident) = names(Idents(object))
    cells.ident[cells.to.assign] = levels(Idents(object))[new.assignments]
    cells.ident.factor = factor(cells.ident, levels = levels(Idents(object)))
    if (!is.null(target.cells)){
      cells.ident.factor = cells.ident.factor[target.cells]
    }
   
    
    iter=iter+1
    frac.unassigned.new = table(cells.ident.factor)[null.ident] / length(cells.ident)
    percentage.reduction = (frac.unassigned - frac.unassigned.new)  / frac.unassigned
    ref.cells = NULL
    if (percentage.reduction < 0.01){
      l = 1
    }
  }
  
  return(cells.ident.factor)
  
}


findMatchToClust = function(object, clust.to.correct = "45_M1",id.now="GraphBoost", cells.ident = NULL, k.use=10, min.mnn=0.4){
  
  cells.in.clust = which.cells(object, clust.to.correct)
  unassigned.cells = object@cell.names[object@data.info[,id.now] == "Unassigned"]
  other.cells = object@cell.names[object@data.info[,id.now] != "Unassigned" & object@data.info[,"atlas_labels_orig"] != "Cr2w"]
  
  # Find mutual nearest neighbors
  nearest.cells = get.knnx(object@pca.rot[unassigned.cells,],object@pca.rot[cells.in.clust,],k=k.use, algorithm = "kd_tree" )
  nearest.un = get.knnx(object@pca.rot[other.cells,],object@pca.rot[unassigned.cells,],k=k.use, algorithm = "kd_tree" )
  
  # Find unassigned cells with
  cells.index = match(cells.in.clust, other.cells)
  cells.nn = unassigned.cells[apply(nearest.un$nn.index,1,function(x){ sum(cells.index %in% x) > 0} )]
  cells.un.index = match(cells.nn,unassigned.cells)
  cells.use = cells.nn[sapply(cells.un.index,function(x){sum(nearest.cells$nn.index %in% x) >= min.mnn*length(cells.in.clust)})]
  return(cells.use)
}



GraphAssignmentWithinClust = function(object, clust, k=30, knn=10, test.label=NULL){
  cells.use = rownames(object@data.info)[object@data.info$m == clust]
  
  cells.ident0 = as.character(object@ident[cells.use]); names(cells.ident0) = names(object@ident[cells.use])
  cells.ident = cells.ident0
  cells.ident.factor = factor(cells.ident, levels = levels(object@ident))
  
  l=0
  iter=1
  print(iter)
  cells.to.assign = names(cells.ident)[cells.ident == "Unassigned"]
  cells.unassigned = cells.to.assign
  ref.cells = setdiff(cells.use, cells.unassigned)
  ref.cells = ref.cells[object@data.info[ref.cells,"atlas_labels_orig"] != "Cr2w"]
  
  if (length(cells.to.assign) < 3 | length(ref.cells) < 20) return(cells.ident)
  
  frac.unassigned = length(cells.unassigned)/length(cells.use)
  print(paste0(round(frac.unassigned*100), " percent cells are Unassigned"))
  Xquery = object@pca.rot[cells.to.assign,1:k]
  X = object@pca.rot[ref.cells,1:k]
  nearest = get.knnx(X,Xquery,k=knn+1, algorithm = "kd_tree" )
  new.assignments = apply(nearest$nn.index[,1:knn],1,function(x){
    Nvotes = table(cells.ident.factor[ref.cells][x]);
    if (max(Nvotes) > 0.7*knn){
      return(which.max(Nvotes))
    } else {
      return(47)
    }
  })
  
  cells.ident = as.character(cells.ident); names(cells.ident) = names(object@ident[cells.use])
  cells.ident[cells.to.assign] = levels(object@ident)[new.assignments]
  cells.ident.factor = factor(cells.ident, levels = levels(object@ident))
  
  #A_bst = table(test.label,cells.ident)
  #A_bst = A_bst[,levels(object@ident)]
  #row.order = c(18,10,3,5,4,45,2,16,21,13,22,6,11,20,31,7,41,24,12,9,15,17,26,14,29,37,36,32,34,28)
  #row.order = c(row.order,setdiff(1:47, row.order))
  #col.order = c(1:3,4,29,5,18,6,7,8,9,10,24,11,12,21,13:16,35,17,20,19,22,25,23,30,26,32,27,31,33,34,45,28,39:42)
  #col.order = c(col.order,setdiff(1:47, col.order))
  #p2 = plotConfusionMatrix(A_bst[row.order,col.order],order=NULL,ylab.use = "DR", xlab.use="Atlas RGC", x.lab.rot = 45,col.low="white",col.high="darkblue") + ggtitle("XG Boost")
  
  frac.unassigned.new = table(cells.ident.factor)[47] / length(cells.ident)
  percentage.reduction = (frac.unassigned - frac.unassigned.new)  / frac.unassigned
  if (percentage.reduction < 0.12)
  
  return(cells.ident)
  
}


# Plotting the Confusion matrix

plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, col.low="blue", col.high="red", max.size=5, ylab.use="Known", xlab.use="Predicted", order=NULL, x.lab.rot=TRUE, plot.return=TRUE, max.perc=100, title.use = "Validation test"){
  
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  X[is.na(X)] = 0
  if (max(X) > 100){
    X=X/100
  }
  
  orig.rownames = rownames(X)
  orig.colnames = colnames(X)
  if (!is.null(order)){
    if (order == "Row"){  
      factor.levels = c()
      for (i1 in colnames(X)){
        if (max(X[,i1]) < 50) next
        ind.sort = rownames(X)[order(X[,i1], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
    
    if (order == "Col") {
      factor.levels = c()
      for (i1 in rownames(X)){
        if (max(X[i1,]) < 50) next
        ind.sort = rownames(X)[order(X[i1,], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(t(X)), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
  } else {
    factor.levels = rownames(t(X))
  }
  
  factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  #X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  #X$Predicted = factor(X$Predicted, levels = rev(factor.levels))
  
  if (!is.null(order)){
    if (order == "Row"){ 
      X$Known = factor(X$Known, levels=rev(factor.levels));
      X$Predicted = factor(X$Predicted, levels = orig.colnames)
      
    }
    if (order == "Col"){
      X$Predicted = factor(X$Predicted, levels = factor.levels);
      X$Known = factor(X$Known, levels=rev(orig.rownames));
    }
  } else {
    X$Known = factor(X$Known, levels=rev(unique(X$Known)));
    X$Predicted = factor(X$Predicted, levels=unique(X$Predicted));
  }
  
  #print(sum(is.na(X$Known)))
  
  
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low =col.low,   high = col.high, limits=c(0, 100 ))+scale_size(range = c(1, max.size), limits = c(0,max.perc))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic")) + ggtitle(title.use) 
  
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  print(p)
  
  if (plot.return) return(p)
}
