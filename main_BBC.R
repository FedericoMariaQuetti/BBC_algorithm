require(caret)
require(data.table)
require(ModelMetrics)
require(stats)
require(crank)
require(MCMCpack)
require(MASS)
require(cluster)
require(mvtnorm)
require(mltools)
require(crank)
require(plotly)
require(ggplot2)
library(readxl)

require(factoextra)

generate_ds <- function(prp = c(3,3,3),size=33,means = matrix(0,3,2),sg = 1) {
  
  dim = ncol(means)
  k = length(prp)
  
  df <- data.frame()
  
  for (i in 1:k) {
    set.seed(i)
    df_k <- as.data.frame(mvrnorm(n = size*prp[i], 
                                  mu = means[i,],  
                                  Sigma = diag(dim)*sg))
    df_k$cluster <- i
    df <- rbind(df, df_k)
  }
  
  df_normalized <- df
  # df_normalized <- df_normalized[,-ncol(df_normalized)]
  # 
  # for (j in 1:nrow(df)) {
  # 
  #   df_normalized[j,] <- df_normalized[j,]/norm(df_normalized[j,], type="2")
  # 
  # }
  # 
  # df_normalized$cluster <- df$cluster
  
  cls <- c()
  for (dmn in 1:dim) {cls <- append(cls, paste("X", dmn, sep = "_"))}
  colnames(df_normalized) <- c(cls,"cluster")
  df_normalized$index <- 1:size
  
  df_normalized
  
}

clust <- function(df, k=3){
  
  df <- as.data.frame(df)
  
  # km <- kmeans(subset(df, select = -c(cluster)),k, iter.max = 15)
  
  km <- kmeans(subset(df, select = -c(cluster, index)),k)
  
  df$cluster <- km$cluster
  
  df
  
}

mypermute <- function(dataf){
  
  vec <- dataf$cluster
  
  df <- numeric()
  
  val <- unique(vec)
  
  perm_tot <- permute(val)
  
  for (i in 1:nrow(perm_tot)) {
    
    perm <- perm_tot[i,]
    
    v <- vec
    
    for (j in 1:ncol(perm_tot)) {
      
      v[which(vec == val[j])] <- perm[j]
      
    }
    
    df <- cbind(df,v)
    
  }
  
  as.data.frame(df)
  
}

clu_labels_adj <- function(df1,df2){
  
  # questa funzione è ottimizzabilissima, definisce robe in più
  
  temp <- df2
  temp$original_cluster <- 0
  
  for (i in 1:nrow(temp)) {
    temp$original_cluster[i] <- df1$cluster[which(df1$index == temp$index[i])]
  }
  
  clu_ind <- mypermute(temp)
  
  cnt <- 0
  
  ind <- 1
  
  for (l in 1:ncol(clu_ind)) {
    if (length(which((temp$original_cluster-clu_ind[[l]]) == 0)) > cnt) {
      cnt <- length(which((temp$original_cluster-clu_ind[[l]]) == 0))
      ind <- l
    }
    
  }
  
  df2$cluster <- clu_ind[[ind]]
  
  df2
  
}

cluster_partecipation <- function(df1,df2,numclu=3){ 
  
  df <- subset(df1, select = -c(cluster))
  
  dim <- ncol(df)
  
  for (nc in 1:numclu) {
    
    df$zero <- 0
    
    colnames(df)[nc+dim] <- paste("partecipation_to_cluster_",nc,sep = "")
    
    for (i in 1:nrow(df)) {
      
      df[i,nc+dim] <- length(which(df2$index == i & df2$cluster == nc))
      
    }
    
  }
  
  df
  
}

get_res3 <- function(cl_df,ds_og,dim=3,nclu,PATH="."){
  
  vdim <- 1:dim
  res <- data.table()
  partecipations <- subset(cl_df, select = -copy)[,-vdim]
  res_part <- data.frame()
  
  for (i in 1:length(unique(cl_df$index))) {
    res <- rbind(res,cl_df[i,vdim])
    res_part <- rbind(res_part,(colSums(partecipations[which(partecipations$index == i),-1])))
    colnames(res_part) = colnames(partecipations)[-1]
  }
  
  totals <- rowSums(res_part)
  
  res_part <- res_part / totals 
  
  # res_part$index <- 1:nrow(res_part)
  
  write.csv2(res_part, paste(PATH, "/results_nclu", nclu, ".csv", sep = ""), row.names = FALSE)
  
  res_part
  
}

############# functions for prior part

getcenters <- function(dtst) {
  
  dim <-  ncol(dtst)-2
  
  cls <- length(unique(dtst$cluster))
  
  centers <- data.frame()
  
  for (i in 1:cls) {
    
    centers <- rbind(centers,as.data.frame(colMeans(dtst[which(dtst$cluster==i),]))[1:dim,])
    
  }
  
  centers
  
}

getcovariances <- function(dtst) {
  
  dim <- ncol(dtst) - 2
  
  cls <- length(unique(dtst$cluster))
  
  covariances <- list()
  
  for (i in 1:cls) {
    
    # Store the covariance matrix in the list
    covariances[[i]] <- cov(cl0[which(dtst$cluster == i),1:dim])
    
  }
  
  covariances
}

pbbs_mixture <- function(ds_clusters, w=0.25, mu, sigma) {
  
  # manca parametro m da inserire, per ora preso = n
  
  dts <- as.data.table(subset(ds_clusters, select = -c(cluster)))
  
  # Calculate the number of clusters
  num_clusters <- nrow(mu)
  
  # Calculate the proportions of each cluster
  cluster_proportions <- sapply(1:num_clusters, function(i) {
    length(which(ds_clusters$cluster == i)) / nrow(ds_clusters)
  })
  
  n <- nrow(dts)
  k <- n * w / (1 - w)
  m <- n
  
  weights <- rdirichlet(1, rep((n + k) / m, m))
  
  xstar <- data.table()
  index <- numeric()
  
  for (j in 1:m) {
    soglia <- runif(1, 0, 1)
    
    if (soglia > w) {
      ind <- sample.int(n, size = 1)
      xstar <- rbind(xstar, dts[ind, ])
    } else {
      soglia2 <- runif(1, 0, 1)
      cumulative_proportion <- 0
      
      for (i in 1:num_clusters) {
        cumulative_proportion <- cumulative_proportion + cluster_proportions[i]
        if (soglia2 < cumulative_proportion) {
          media <- mu[i,]
          sigma_cl <- sigma[[i]]
          break
        }
      }
      
      v_0 <- data.table(t(mvrnorm(1, as.numeric(media), sigma_cl)))
      v_0$index <- 0
      colnames(v_0) <- colnames(dts)
      xstar <- rbind(xstar, v_0)
    }
  }
  
  # xstar
  
  resample_cumulative <- cumsum(weights)
  rnd <- runif(m)
  
  ds <- data.table()
  
  for (i in 1:m) {
    
    if(rnd[i]<min(resample_cumulative)){
      ds <- rbind(ds,xstar[1,])
    }
    else {
      ds <- rbind(ds,xstar[max(which(resample_cumulative <= rnd[i]))+1,])
    }
    
  }
  
  ds
  
}

clust2 <- function(df, k=3){
  
  df <- as.data.table(df)
  
  km <- kmeans(subset(df, select = -index),k, iter.max = 15)
  
  # km <- pam(subset(df, select = -index),k)
  
  df$cluster <- km$cluster
  
  df
  
}

clust3 <- function(df, cnt, k=3){
  
  df <- as.data.table(df)
  
  km <- kmeans(subset(df, select = -index), cnt, k, iter.max = 15)
  
  # km <- pam(subset(df, select = -index),k)
  
  df$cluster <- km$cluster
  
  df
  
}

bigdist <- function(dt1,dt2) {
  
  dis <- 0
  
  dt1_ <- subset(dt1, select = -c(index, cluster))
  dt2_ <- subset(dt2, select = -c(index, cluster))
  
  for (i in 1:max(dt1$cluster)) {
    dis = dis + norm(as.matrix(colMeans(dt1_[which(dt1$cluster == i),])) - as.matrix(colMeans(dt2_[which(dt2$cluster == i),])), type = "2")
  }
  
  sqrt(dis)
  
}

relabel <- function(dt1, dt2) {
  
  dt1 <- as.data.frame(dt1)
  dt2 <- as.data.frame(dt2)
  
  ind <- mypermute(dt2)
  
  d <- Inf
  counter <- 1
  
  for (i in 1:ncol(ind)) {
    dt2$cluster <- as.data.frame(ind)[,i]
    
    dst <- bigdist(dt1,dt2)
    
    if (dst <= d) 
    {d <- dst
    counter <- i}
    
  }
  
  dt2$cluster <- as.data.frame(ind)[,counter]
  dt2
}

fin_res <- function(my_data){
  # Initialize an empty dataframe to store the result
  result <- data.frame(index = unique(my_data$index))
  additional_columns <- matrix(0, nrow = nrow(result), ncol = max(my_data$cluster))
  colnames(additional_columns) <- paste0("partecipation_to_cluster_", 1:max(my_data$cluster))
  
  # Concatenate additional columns to my_data
  result <- cbind(result, additional_columns)
  
  # Loop through each unique value of "index"
  for (i in 1:nrow(result)) {
    # Subset the data for the current index value
    subset_data <- my_data[my_data$index == result$index[i], ]
    
    # Count the occurrences of each cluster label
    counts <- table(subset_data$cluster)
    
    counts <- counts[order(as.integer(names(counts)))]
    
    result[i, paste0("partecipation_to_cluster_", names(counts))] <- counts
    
    
  }
  
  # Replace NA values with 0
  result[is.na(result)] <- 0
  
  result <- result[order(result$index),]
  
  result_fin <- result[-1,]
  result_fin_temp <- result_fin[,-1]
  
  
  for (i in 1:nrow(result_fin)) {
    result_fin$cluster[i] <- which(result_fin_temp[i,] == max(result_fin_temp[i,]))[1]
  }
  
  result_fin
  
}

compute_entropy <- function(row) {
  entropy <- -sum(row * log(row + (row == 0)))
  if (is.nan(entropy)) {
    entropy <- 0  # Assign entropy 0 if a component is 1
  }
  return(entropy)
}

compute_matrix_element_relativeentropy <- function(row) {
  n <- length(row)
  result <- numeric(n * (n - 1) / 2)
  index <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      ej <- row[i]
      ek <- row[j]
      ejk <- ej / (ej + ek)
      ekj <- ek / (ej + ek)
      element <- -(ejk * log2(ejk) + ekj * log2(ekj))
      # if (is.nan(element)) {
      #   element <- 0  # Set NaNs to 0
      # }
      result[index] <- element
      index <- index + 1
    }
  }
  return(result)
}

compute_matrix_element_mergingentropy <- function(row) {
  n <- length(row)
  result <- numeric(n * (n - 1) / 2)
  index <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      ej <- row[i]
      ek <- row[j]
      ejk <- ej / (ej + ek)
      ekj <- ek / (ej + ek)
      element <- (ej + ek) * log2(ej + ek) - ej * log2(ej) - ek * log2(ek)
      # if (is.nan(element)) {
      #   element <- 0  # Set NaNs to 0
      # }
      result[index] <- element
      index <- index + 1
    }
  }
  return(result)
}

compute_relative_entropy <- function(ds){
  n <- ncol(ds)
  
  
  result_df <- data.frame(matrix(NA, nrow = nrow(ds), ncol = n*(n-1)/2))
  
  # Iterate over each row of ds
  for (i in 1:nrow(ds)) {
    # Compute the matrix element relative entropy for the current row
    result <- compute_matrix_element_relativeentropy(ds[i,])
    
    # Store the result as a row in result_df
    result_df[i,] <- result
  }
  
  # Convert the result_df to a data frame
  result_df <- as.data.frame(result_df)
  
  
  col_names <- c()
  index <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      col_names[index] <- paste0(i, "_", j)
      index <- index + 1
    }
  }
  
  colnames(result_df) <- col_names
  result_df[is.na(result_df)] <- 0
  return(result_df)
  
}

compute_merging_entropy <- function(ds){
  n <- ncol(ds)
  
  result_df <- data.frame(matrix(NA, nrow = nrow(ds), ncol = n*(n-1)/2))
  
  # Iterate over each row of ds
  for (i in 1:nrow(ds)) {
    # Compute the matrix element relative entropy for the current row
    result <- compute_matrix_element_mergingentropy(ds[i,])
    
    # Store the result as a row in result_df
    result_df[i,] <- result
  }
  
  # Convert the result_df to a data frame
  result_df <- as.data.frame(result_df)
  
  
  
  col_names <- c()
  index <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      col_names[index] <- paste0(i, "_", j)
      index <- index + 1
    }
  }
  
  colnames(result_df) <- col_names
  result_df[is.na(result_df)] <- 0
  return(result_df)
  
}

visual <- function(temp_vis = final){
  
  num_clusters <- length(unique(temp_vis$cluster))
  
  # Dynamically assign cluster names and colors
  cluster_names <- paste0("Cluster", 1:num_clusters)
  cluster_colors <- rainbow(num_clusters)
  
  # Map cluster numbers to cluster names
  for (i in 1:num_clusters) {
    temp_vis$cluster[which(temp_vis$cluster == i)] <- cluster_names[i]
  }
  
  # Convert cluster column to factor
  temp_vis$cluster <- as.factor(temp_vis$cluster)
  # temp_vis$cluster[which(temp_vis$cluster == 1)] <- 'Cluster1'
  # temp_vis$cluster[which(temp_vis$cluster == 2)] <- 'Cluster2'
  # temp_vis$cluster[which(temp_vis$cluster == 3)] <- 'Cluster3'
  # temp_vis$cluster <- as.factor(temp_vis$cluster)
  
  fig <- plot_ly(temp_vis, x = ~X_1, y = ~X_2, color = ~cluster, colors = c('#BF382A', '#0C4B8E','#FFFF00'))
  fig <- fig %>% add_markers()
  # fig <- fig %>% layout(scene = list(xaxis = list(title = 'Cluster1'), yaxis = list(title = 'Cluster2')))
  
  fig
  
}

initialize_df <- function(nclu) {
  # Number of rows and columns
  num_rows <- nclu-1
  num_cols <- nclu * (nclu - 1) / 2
  
  # Initialize a matrix with NA values
  mat <- matrix(NA, nrow = num_rows, ncol = num_cols)
  
  # Convert the matrix to a data frame
  df <- as.data.frame(mat)
  
  # Optionally, set row and column names
  rownames(df) <- 2:nclu
  # col_names <- c()
  # index <- 1
  # for (i in 1:(nclu - 1)) {
  #   for (j in (i + 1):nclu) {
  #     col_names[index] <- paste0(i, "_", j)
  #     index <- index + 1
  #   }
  # }
  # 
  # colnames(df) <- col_names
  
  return(df)
}

########## main

# main <- function(){

setwd("C:/Users/Federico/Desktop/BayesClustering/final")
PATH <- "."

ncl <- 3
ds1 <- generate_ds(prp = c(1,1,1)*3,means = rbind(c(1,0),c(0,sqrt(3)),c(-1,0)))
sample_size <- nrow(ds1)
dimens <- ncol(ds1)-2

### dimension of the observations; the dataset is comprised of the vector 
### and a column for the cluster index

ds1$index <- 1:sample_size

write.csv2(ds1, "./dataset.csv", row.names = FALSE)

### clustering of the original dataset

###### different choices of K number of clusters
###### to be used for partitioning

for (ncl in 2:6) {

  print(ncl)
  
  cl0 <- clust(ds1,ncl)
  
  cl0 <- clu_labels_adj(ds1,cl0)
  
  write.csv2(cl0, paste("./dataset_clustered_K", ncl, ".csv", sep = ""), row.names = FALSE)
  
  ds1 <- read.csv2("./dataset.csv")
  dim <-  ncol(ds1)-2

  cl0 <- clust2(subset(ds1,select=-cluster),ncl)
  cl0$index <- 1:nrow(cl0)
  
  if(max(unique(ds1$cluster)) == ncl)
  {cl0 <- relabel(ds1,cl0)}
  
  centers <- getcenters(cl0)
  
  sig <- getcovariances(cl0)
  
  for(wg in 0.5){
    # tuning parameter of the convex combination of prior and posterior
    print(wg)
    vec <- c(1,1.5,3,4.5,6)
    for (l in vec){
      # tuning parameter of the informativeness of the prior, s
      print(l)
      
      res_final <- data.frame()
      
      for (i in 1:100) {
        
        print(i)
        
        bsdts_mixture <- pbbs_mixture(cl0, w=wg, mu = centers, sigma = sig)
        
        cldts <- clust3(bsdts_mixture, centers, ncl)
        
        final <- cldts
        
        final$copy <- i
        
        res_final <- rbind(res_final, final)
        
      }
      
      # memberships results
      
      df_fin <- fin_res(res_final)
      
      res_df_fin <- subset(df_fin, select = -c(index, cluster))
      
      res_df_fin <- res_df_fin/rowSums(res_df_fin)
      
      write.csv2(res_df_fin, paste("./results_pbbs_nclu", ncl, "_var", l, ".csv", sep = ""), row.names = FALSE)
      
      
    }
    
  }
  
}



########## evaluations

print("Computing results")

nclu_max <- 6

for (l in vec){
  print(l)
  
  mn_en <- data.frame()
  re_en <- initialize_df(nclu_max)
  re_en_pbbs <- initialize_df(nclu_max)
  
  for (nclu in 2:nclu_max) {
    
    print(nclu)
    
    ds1 <- read.csv2(paste("./results_pbbs_nclu", nclu, "_var", l, ".csv", sep = ""))
    
    entropy_values <- apply(ds1, 1, compute_entropy)
    m_pbbs <- mean(entropy_values)
    
    tmp <- ds1
    tmp$entropy <- entropy_values
    hist(entropy_values, breaks = 50, main = paste("Entropies, nclu = ", nclu))
    print(paste0("Mean of Shannon entropies pbbs: ", mean(entropy_values)))
    
    relative_entropy_pbbs <- compute_relative_entropy(ds1)
    
    
    mn <- data.frame(mean_entropy = nclu, mean_entropy_pbbs = m_pbbs)
    rownames(mn) <- nclu
    
    mn_en <- rbind(mn_en, mn)
    
    r_mn_pbbs <- colMeans(relative_entropy_pbbs)
    
    for (i in 1:length(r_mn_pbbs)) {
      
      re_en_pbbs[nclu-1,i] <- r_mn_pbbs[i]
      
    }
    
  }
  
  # nome brutto, è mean entropy dopo pbbs
  write.csv2(mn_en, paste0("./mean_entropy_var", l, ".csv"), row.names = FALSE)
  # write.csv2(mn_en, paste0("./mean_entropy_var", l, ".csv"), row.names = FALSE)
  write.csv2(re_en_pbbs, paste0("./rel_entropy_pbbs_var", l, ".csv"), row.names = FALSE)
  
}