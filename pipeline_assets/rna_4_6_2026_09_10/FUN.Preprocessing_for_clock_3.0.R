preprocessing_for_clock <- function(exprs_data,gene_list,log_transform=F,
                                    show_plots=T,map_to_orthologs=F,map_from="Homo sapiens",
                                    scaling=T,yugene=T,differential=F,control_group=NULL,
                                    map_to="Mus musculus",map_type="Entrez",
                                    function_path="/Users/Alex/Documents/Science/PhD/Functions/"){
  
  #1. Map to orthologs (if needed)
  source(paste0(function_path,"FUN.Map_to_orthologs.R"))
  
  if (log_transform==T){
    exprs_data <- log10(exprs_data+1)
  }
  
  if (map_to_orthologs==T){
    print("Mapping to orthologs...")
    source(paste0(function_path,"FUN.Map_to_orthologs.R"))
    mapping_table <- Map_to_orthologs(rownames(exprs_data),
                                species_from = map_from,species_to = map_to,
                                type=map_type)
    exprs_data <- exprs_data[!is.na(mapping_table[colnames(exprs_data),]$To_Entrez),]
    rownames(exprs_data) <- as.character(mapping_table[rownames(exprs_data),]$To_Entrez)
  }
  
  #2. Scaling and YuGene transformation
  if (scaling==T){
    print("Scaling...")
    exprs_scale <- scale(exprs_data)

    if (show_plots==T){
      plot(density(as.numeric(exprs_scale[,1]),na.rm = T))
      for (i in 1:ncol(exprs_scale)){
        lines(density(as.numeric(exprs_scale[,i]),na.rm=T),col=i)
      }
    }
  }
  
  if (yugene==T){
    print("Yugene transformation...")
    library(YuGene)
    exprs_yugene <- apply(exprs_scale,2,function(x){x-min(x,na.rm=T)})
    # fill na with 0
    exprs_yugene[is.na(exprs_yugene)] <- 0
    exprs_yugene = YuGene(exprs_yugene)

    if (show_plots==T){
      plot(density(as.numeric(exprs_yugene[,1]),na.rm = T))
      for (i in 1:ncol(exprs_yugene)){
        lines(density(as.numeric(exprs_yugene[,i]),na.rm=T),col=i)
      }
    }
    
  }
  
  #3. Create input table
  print("Creating input table...")
  common_genes <- intersect(gene_list,rownames(exprs_data))
  extra_genes <- setdiff(gene_list,rownames(exprs_data))
  print("Total number of genes in the clock:")
  print(length(gene_list))
  print("Number of genes presented in the data:")
  print(length(common_genes))

  output_list <- list()
  if (scaling==T){
    exprs_scale = t(exprs_scale)
    exprs_scale_final = matrix(nrow=nrow(exprs_scale),ncol=length(gene_list))
    rownames(exprs_scale_final) = rownames(exprs_scale)
    colnames(exprs_scale_final) = gene_list
    
    exprs_scale_final[,common_genes] <- exprs_scale[,common_genes]
    if (differential==T){
      if (is.null(control_group)){
        exprs_scale_control = apply(exprs_scale_final,2,median)
      }else{
        exprs_scale_control = apply(exprs_scale_final[control_group,,drop=F],2,median)
      }
      print(length(exprs_scale_control))
      exprs_scale_diff = apply(exprs_scale_final,1,function(x){x-exprs_scale_control})
      exprs_scale_diff = as.data.frame(t(exprs_scale_diff))
      output_list[["Scaled_differential_data"]] <- exprs_scale_diff
    }
    exprs_scale_final = as.data.frame(exprs_scale_final)
    output_list[["Scaled_data"]] <- exprs_scale_final
  }
  if (yugene==T){
    exprs_yugene = t(exprs_yugene)
    exprs_yugene_final = matrix(nrow=nrow(exprs_yugene),ncol=length(gene_list))
    rownames(exprs_yugene_final) = rownames(exprs_yugene)
    colnames(exprs_yugene_final) = gene_list
    exprs_yugene_final[,common_genes] <- exprs_yugene[,common_genes]
    if (differential==T){
      if (is.null(control_group)){
        exprs_yugene_control = apply(exprs_yugene_final,2,median)
      }else{
        exprs_yugene_control = apply(exprs_yugene_final[control_group,,drop=F],2,median)
      }
      print(length(exprs_yugene_control))
      exprs_yugene_diff = apply(exprs_yugene_final,1,function(x){x-exprs_yugene_control})
      exprs_yugene_diff = as.data.frame(t(exprs_yugene_diff))
      output_list[["YuGene_differential_data"]] <- exprs_yugene_diff
    }
    
    exprs_yugene_final = as.data.frame(exprs_yugene_final)
    output_list[["YuGene_data"]] <- exprs_yugene_final
  }
  
  return(output_list)
}
