
globalMahalanobisDistances <- function(Table1, coverVariance, minObjects, SType = "vehicle control",  url){
  
    require(tidyverse)
    require(tictoc)
  
    ############## 1. grab all well-level data
    
    Table2 = Table1 %>% filter(n_cells_keep>=minObjects & rel_cell_count>50 & stype!="null") 
    
    ############## 2. Calculate the Eigen features from the well-level data (~ 1 min)
    
    iDataCols = which(grepl("^f_", colnames(Table2)))
    
    PCA = prcomp(Table2[,iDataCols], center=F, scale.=F)
    RotationMatrix = PCA$rotation
    
    CumProportion = cumsum(PCA$sdev^2)/sum(PCA$sdev^2)

    ##############  3. Find the inverse of the covariance matrix
    
    ## 3.a) Model the data to get the covariance matrix
    PC = length(which(CumProportion<coverVariance))+1
    Model = lm(PCA$x[,1:PC] ~ 0+Table2$trt_name) 
    
    ## 3.b) 
    Cov = estVar(Model)
    
    ## 3.c) inverse
    invCov = solve(Cov)
    
    ##############  4. Calculate Mahalanobis for each well on a per plate basis
    
    ## 4.a) transform data from table1 to the eigenfeatures
    iDataCols     = which(grepl("^f_", colnames(Table1)))
    iMetadataCols = which(!grepl("^f_", colnames(Table1)))
    
    data = as.matrix(Table1[,iDataCols]) %*% PCA$rotation[,1:PC]
    transfTable1 = cbind(Table1[,iMetadataCols], data)
    rm(data)
    
    iEigenfeatureCol=length(iMetadataCols)+(1:PC)
    
    PlateID = transfTable1$plate_id[1]
    
    Result= NULL
    for (PlateID in unique(transfTable1$plate_id)){
      print(PlateID)
      
      Subset = transfTable1 %>% filter(plate_id==PlateID)
    
      ctrMean = Subset %>% filter(stype==SType) %>% 
        summarise_at(.vars=colnames(Subset)[iEigenfeatureCol], .funs="mean") 
      
      Delta = Subset
      Delta[iEigenfeatureCol] = sweep(as.matrix(Subset[iEigenfeatureCol]), 2, as.matrix(ctrMean), "-")
      
      D = apply(Delta[iEigenfeatureCol], 1, function(x) sqrt((x %*% invCov %*% x)))
      
      Distance = Subset %>% select(1:length(iMetadataCols)) %>% mutate(d = round(D,3))
      
      Result = bind_rows(Result, Distance)
      
      rm(Subset, ctrMean, Delta, Distance)
      
    }#for each plate
    
    ##############  5. Write results in a new collection
    htpp_global_mah <- mongo(collection="htpp_global_mah", url=url, verbose=getOption("verbose"))
        
    ## 5.a) remove all data that was in this collection
    htpp_global_mah$drop()
    
    ## 5.b) write each row as a document in the collection
    print(dim(Result)[1])
    for(i in 1:dim(Result)[1]){
          newDocument = Result[i,]
          htpp_global_mah$insert(newDocument, auto_unbox=TRUE)
    }#for each row
    
    print(htpp_global_mah$count() )
 
    Output = list(CumProportion = CumProportion,
                  RotationMatrix = RotationMatrix,
                  invCov = invCov)
    return(Output)
}#end function