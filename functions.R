pca <- function(x, space = c("rows", "columns"),
                center = TRUE, scale = FALSE) {
  space <- match.arg(space)
  if (space == "columns") {
    x <- t(x)
  }
  x <- t(scale(t(x), center = center, scale = scale))
  x <- x / sqrt(nrow(x) - 1)
  s <- svd(x)
  loading <- s$u
  colnames(loading) <- paste0("Loading", 1:ncol(loading))
  rownames(loading) <- rownames(x)
  pc <- diag(s$d) %*% t(s$v)
  rownames(pc) <- paste0("PC", 1:nrow(pc))
  colnames(pc) <- colnames(x)
  pve <- s$d^2 / sum(s$d^2)
  if (space == "columns") {
    pc <- t(pc)
    loading <- t(loading)
  }
  return(list(pc = pc, loading = loading, pve = pve))
}

make_filtered_data = function(counts, metadata){
  print(colnames(metadata))
  filtered_metadata <- metadata
  filtered_counts <- as_data_frame(counts)
  for (column in names(feature_vec)) {
    if (column %in% colnames(filtered_metadata)){
      control_names_tb = table(filtered_metadata[,column])
      control_names_tb = control_names_tb[names(control_names_tb) %in% feature_vec[[column]]]
      print(control_names_tb)
      control_name = names(control_names_tb)[control_names_tb == max(control_names_tb)]
      filtered_counts <- filtered_counts[,filtered_metadata[,column] == control_name]
      filtered_metadata <- filtered_metadata[filtered_metadata[,column] == control_name,]
    }
  }
  return(list(counts = filtered_counts, 
              metadata = filtered_metadata))
}

get_control_cols <- function(metadata){
  names(metadata)[!(names(metadata) %in% c(columns_to_ignore,"lib.size","norm.factors"))]
}

remove_large_factors = function(metadata){

  cols_to_control = get_control_cols(metadata)
  for(col in cols_to_control){
    n_coef = length(unique(metadata[[col]])) - 1
    print(paste(col, n_coef))
    if(n_coef > n_samples/5 | n_coef <= 0){
      print(paste("Droping column", col, "from metadata because it has ", n_coef, "levels."))
      metadata = metadata[,names(metadata) != col]
      cols_to_control = cols_to_control[cols_to_control!=col]
     }
  }
  return(metadata)
}

remove_redundant_features <- function(metadata){  
    # Removing redundant features - see https://stackoverflow.com/questions/38222318/how-to-remove-duplicate-columns-content-in-data-table-r

  features_pair <- combn(names(metadata), 2, simplify = F) # list all column pairs
  to_remove <- c() # init a vector to store duplicates
  for (pair in features_pair) { # put the pairs for testing into temp objects
    f1 <- pair[1]
    f2 <- pair[2]
    if (!(f1 %in% to_remove) & !(f2 %in% to_remove)) {
      if (all(as.numeric(as.factor(metadata[[f1]])) == as.numeric(as.factor(metadata[[f2]])))) {
        cat(f1, "and", f2, "are equal.\n")
        to_remove <- c(to_remove, f2) # build the list of duplicates
      }
    }
  }

  metadata <- metadata[,!(names(metadata) %in% as.list(to_remove))]
}

pca_plot = function(resids){
  pcar = pca(resids)
  results = t(pcar$pc)
  mean = data.frame(matrix(colMeans(results), 
                    ncol = dim(results)[2]))
  colnames(mean) = colnames(results)
  vars = data.frame(matrix(colVars(results), 
                           ncol = dim(results)[2]))
  colnames(vars) = colnames(results)
  
  plot <-  ggplot(results, aes(PC1, PC2)) + 
    geom_point() + 
    stat_ellipse(level = 0.99) + 
    stat_ellipse(level = 0.99, type = "norm", linetype = 2) + 
    geom_point(data= mean, color = "tomato3", size = 7) +
    coord_fixed(ratio=1) +
    labs(x = paste0("PC1 (", round(pcar$pve[1]*100, 2), "%)"),
         y = paste0("PC2 (", round(pcar$pve[2]*100, 2), "%)")) +
    theme_cowplot()
  plot
}


scree_plot <- function(resids){
  var_explained <- pca(resids)$pve
  scree_on_resids  <- qplot(c(1:10), var_explained[1:10]) +
      geom_line() +
      xlab("Principal Component") +
      ylab("Variance Explained") +
      ggtitle("Scree Plot") +
      ylim(0, 1) + 
      coord_fixed(ratio=1) +
      theme_cowplot()
}

vector_cor = function(x, y) x%*%y/sqrt(x%*%x * y%*%y)
