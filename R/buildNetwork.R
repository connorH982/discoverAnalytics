library(magrittr)
library(dplyr)
library(tidyr)
library(data.table)
library(discover)
library(mygene)


getOddsRatio<-function(res,events,type = c("mut","co"),as.vec=TRUE){
  oddsRatios<-rep(NA,nrow(res))
  for(i in seq(nrow(res))){
    cont <- events %>% t %>% as.data.table %>% dplyr::select(unlist(res[i,1]),unlist(res[i,2])) %>% table
    coOccur<-(cont[1,1]+cont[2,2])/sum(cont)
    excl<-(cont[1,2]+cont[2,1])/sum(cont)
    if(type=="mut"){oddsRatios[i]<-abs(log(coOccur/excl))}else{oddsRatios[i]<-abs(log(abs(excl/coOccur)))}
  }
  if(as.vec==FALSE){res %>% mutate(oddsRatios=oddsRatios)}else{oddsRatios}
}

#' Make a dictionary of gene ontology terms, given a list of genes.
#'
#' Uses myGene.io to build a dictionary of gene ontology terms. It relies on myGene's query many function to build this list.
#' The difference here is this function also allows for forward and reverse lookups (either go terms associated with a gene, or genes assoicated with a go term).
#' This is largely to improve performance within the application. This function is also necessary to run group-wise tests defined by gene ontology terms.
#' @param mtdt Can either be an indicator matrix of the type used by the discover test (see ?\code{\link{buildDirectory}}) or a list of genes.
#' @import magrittr dplyr data.table mygene
#' @return This will return a nested list. The list can be searched by either genes (searching associated gene ontology terms) or by gene ontology term (in this case searching for associated genes).
#' @export
makeDictionary <- function(events){
  if(is.matrix(events)){
    search.terms<-rownames(events)
  }else if(all(names(events) == c("events","bg"))){
    search.terms<-rownames(events$events)
  }else{
    search.terms<-events
  }


  # Get a list of go terms by gene using an API
  search.terms <- search.terms %>% gsub(pattern = "\\(.*\\)",replacement = "") %>% trimws
  dict         <- queryMany(search.terms, scopes='symbol', fields=c('entrezgene', 'go'), species='human',return.as = "DataFrame")
  dict$notfound[is.na(dict$notfound)]<-FALSE

  # Create a dictionary that allows us to search for the go terms of genes
  dictionary<-list()
  for(g in search.terms){

    ind<-which(dict$query==g)
    if(any(dict$notfound[ind]!=TRUE)){
      if(length(ind)>1){ind=ind[1]}
      BP<-dict$go.BP[[ind]]$term
      MF<-dict$go.MF[[ind]]$term
      CC<-dict$go.CC[[ind]]$term
      g<-list(BP=BP,MF=MF,CC=CC)
    }else{
      BP<-"Not Found"
      MF<-"Not Found"
      CC<-"Not Found"
      g<-list(BP=BP,MF=MF,CC=CC)
    }
    dictionary<-append(dictionary,values = list(g))

  }
  names(dictionary)<-search.terms

  # Now allow the dictionary to be searched in reverse as well

  types<-c("bp","mf","cc")


  suppressWarnings({


    stacks<-list()
    for(type in types){

      stacks<-append(stacks,
                     list({
                       dictionary %>% sapply(FUN = `[[`, toupper(type)) %>% stack %>% as.data.table
                     }))

    }
    names(stacks)<-types

    vals<-list()
    for(type in types){

      vals<-append(vals,
                   list({
                     stacks[[type]] %>% dplyr::select(values) %>% unique %>% unlist(use.names = FALSE)
                   }))

    }
    names(vals)<-types



  })


  dict.reverse<-list()
  for(type in types){
    a <- list()
    for(go in vals[[type]]){
      a <- append(a,
                  list({
                    stacks[[type]] %>% subset(values==go) %>% .$ind %>% unique %>% as.character
                  }))

    }
    names(a)<-vals[[type]]
    dict.reverse<-append(dict.reverse,list(a))
  }
  names(dict.reverse)<-types

 append(list(gene = dictionary),dict.reverse)
}


lookUp<-function(dictionary){
  test<-lapply(dictionary$gene,FUN=function(g){
    BP<-g$BP[1]
    MF<-g$MF[1]
    CC<-g$CC[1]
    if(is.null(BP)) BP <- "Not Found"
    if(is.null(MF)) MF <- "Not Found"
    if(is.null(CC)) CC <- "Not Found"
    data.table(BP,MF,CC)
}) %>% rbindlist
}

indicatorWithEffects<-function(mutations){
  samples_by_genes <- mutations[,.(gene=unique(gene)),by="ID"]
  samples_by_genes$one <- 1
  mtdt <- spread(samples_by_genes, gene, one, fill=F)[,-1] %>% t

  mtdt
}


runDiscover<-function(events, q.threshold = 0.01){
  mtdt         <- events$events
  result.mutex <- pairwise.discover.test(events,alternative = "less")
  result.co    <- pairwise.discover.test(events,alternative = "greater")


  res.mutex<-as.data.frame(result.mutex, q.threshold = q.threshold)
  res.co<-as.data.frame(result.co, q.threshold = q.threshold)

  res.combined<-rbind.data.frame(res.mutex,res.co) %>% as.data.table %>% mutate(Effect=c(rep("red",nrow(res.mutex)),rep("green",nrow(res.co)))) %>% as.data.table
  res.combined %<>%
    mutate(odds.ratio=getOddsRatio(.,mtdt,"mut",as.vec = TRUE))
  return(res.combined)
}

pathway.relabel<-function(genes,data.ref,go){
  a<-lapply(genes,FUN=function(g){
    a<-data.ref[gene==g] %>% dplyr::select(go)
  }) %>% unlist
}


#' Prepares a new network visualization
#'
#' This is the main function of the package, meant to prepare mutation data to be visualized by the discoverInteractiveAnalytics application.
#'
#' @param mtdt An indicator matrix, of the same type used in the DISCOVER test. The dimensions should be n genes (rows) by m columns (samples), with cells being 0 (non-mutated) or 1 (mutated).
#' @param path_to_file Where to save the results. This file can be uploaded to the associated application.
#' @param network.name The name of the network that will be displayed within the application.
#' @param q.threshold The FDR cutoff used by the discover test, defaults to 0.01.
#' @param mutation_effects Supply a list of mutations, with a column "gene" and "functionalEffect" to include mutation effects in the visualization. This is not necessary for the majority of functionality, therefore this function will run if it is not supplied.
#' @import magrittr dplyr data.table mygene discover
#' @export
buildAppFile<-function(mtdt, path_to_file, q.threshold = 0.01, mutation_effects = NULL, res.combined = NULL){

  if(!grepl(path_to_file,pattern = ".rds$")){
    path_to_file<-paste(path_to_file,".rds", sep = "")
  }

  # Generate probability background matrix and run discover test
  events<-discover.matrix(mtdt)


  cat("\nRunning discover test. This may take several minutes...\n")
  if(is.null(res.combined)){
    res.combined<-runDiscover(events = events,q.threshold = q.threshold)
  }



  # Build dictionary
  cat("\nBuilding dictionary of gene ontology terms...\n")
  dictionary <- makeDictionary(events)




  # Used to get nodes of the graph (allowing for disconnected nodes) and fill colors
  cat("\nCreating network elements...\n")
  Group<-lookUp(dictionary = dictionary)
  data.ref<-cbind(rownames(events),Group) %>% as.data.table
  names(data.ref)<-c("gene","BF","MF","CC")



  # Finish building edge list
  res.combined$Effect<-as.factor(res.combined$Effect)

  res.combined %<>% mutate(key.column = seq(nrow(res.combined)),
                           BF1        = pathway.relabel(gene1,data.ref,"BF"),
                           BF2        = pathway.relabel(gene2,data.ref,"BF"),
                           MF1        = pathway.relabel(gene1,data.ref,"MF"),
                           MF2        = pathway.relabel(gene2,data.ref,"MF"),
                           CC1        = pathway.relabel(gene1,data.ref,"CC"),
                           CC2        = pathway.relabel(gene2,data.ref,"CC")) %>%
    as.data.table %>% setkey(key.column) %>%
    as.data.table

  # Save necessary files
  if(!is.null(mutation_effects)){
    effct<-mutation_effects %>% dplyr::select(gene,functionalEffect)
  }



  # Take all elements and combine them into a single RDS file

  if(!exists("effct")){
    main    <- list(dictionary,res.combined,events,data.ref)
    names(main) <- c("dict","edges","events","nodes")

  }else{
    main    <- list(dictionary,res.combined,events,data.ref)
    names(main) <- c("dict","edges","events","nodes")
  }
  print(path_to_file)
  print(length(main))
  saveRDS(object = main,file = path_to_file)

  cat("\nFinished building file.\n")
  return(TRUE)
}


#' Run DISCOVER group-wise tests on gene sets defined by gene ontology terms
#'
#' This implements the discover group wise tests of exclusivity, coverage and impurity \link{"https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1114-x"}.
#' Automatically implements a the Benjamini-Hochberg multiple testing correction.
#' @param events Can either be a discover matrix that results from running the discover.pairwise.test() or an indicator matrix. If it is an indicator matrix, the function will run the discover test, but this will significantly increase execution time.
#' @param dictionary If a dictionary object (\code{\link{makeDictionary}}) is not supplied this function will automatically generate one, though this increases execution time.
#' @param go.type Character. Can be "bp","cc" or "mf" for biological process, cellular component, or molecular function respectively. Defaults to "bp"
#' @param group.test.type The type of group test to be employed: "exclusivity", "coverage" or "impurity"
#' @import magrittr dplyr data.table discover
#' @export
groupTestByGeneOntology<-function(events,dictionary = NULL, go.type = "bp", group.test.type = "exclusivity"){
  if(!(group.test.type %in% c("impurity","coverage","exclusivity"))) stop("Available group tests are “impurity”, “coverage”, or “exclusivity” with a default of exclusivity.")


  exclusions<-lapply(names(dictionary[[go.type]]),FUN = function(b){

    genes.in.bp<-dictionary[[go.type]][[b]][dictionary[[go.type]][[b]] %in% rownames(events$events)]
    if(length(genes.in.bp) > 1){
      result<-data.table(process  = b,
                         genes                = paste0(genes.in.bp,collapse = ", "),
                         n.genes              = length(genes.in.bp),
                         exclusion.p       = groupwise.discover.test(events[genes.in.bp, ],method = group.test.type)[["p.value"]])
    }else{
      result<-NULL
    }
    result


  }) %>% do.call(what = "rbind") %>% as.data.table

  exclusions[,BH := stats::p.adjust(exclusion.p,method = "BH")]
  exclusions %<>% arrange(BH)  %>% as.data.table

  exclusions
}





#' Test for Enrichment of particular gene ontology terms
#'
#'
testEnrichment<-function(res.mutex, dictionary, gene_ontology_domain = c("bp","cc","mf"), resamples = 1000){

  if(gene_ontology_domain == c("bp","cc","mf")){

    warning("No gene ontology domain supplied, defaulting to using molecular function")

  }

  # Mirror results to get a pairs to get a symmetric matrix
  res.mutex %<>% mutate(effect.size = -log(p.value)) %>% as.data.table
  res <- res.mutex %>% dplyr::select(gene1,gene2,effect.size)
  res$gene1<-as.character(res$gene1)
  res$gene2<-as.character(res$gene2)
  m <- res
  m$gene1<-res$gene2
  m$gene2<-res$gene1
  res %<>% rbind(m) %>% as.data.table

  # Have levels be identical
  res$gene1<-factor(res$gene1,levels = {
    res$gene1 %>% unique %>% sort
  })
  res$gene2<-factor(res$gene2,levels = {
    res$gene2 %>% unique %>% sort
  })

  # Create matrix
  mat<-spread(res, key = gene2,value = effect.size,fill = 0) %>% setkey("gene1")
  bp<-dictionary$bp %>% lapply(as.data.table) %>% rbindlist(idcol = "bp")
  names(bp)[2]<-"gene"


  get.score<-function(b,mat,dictionary){
    b %<>% unlist
    genes.in.bp<-b[b %in% names(mat)]

    d<-mat[gene1 %in% genes.in.bp,genes.in.bp,with = FALSE]
    sum(d)
  }



  scores.observed<-lapply(names(dictionary$bp),FUN = function(b){

    genes.in.bp<-dictionary$bp[[b]][dictionary$bp[[b]] %in% names(mat)]
    if(length(genes.in.bp)!=0){
      d<-mat[which(mat$gene1 %in% dictionary$bp[[b]]),genes.in.bp,with = FALSE]

      result<-data.table(bp = b,obs = sum(d))
    }else{
      result<-NULL
    }
    result

  }) %>% rbindlist %>% arrange(bp) %>% as.data.table





  #scores.observed$run<-"observed"

  reruns<-lapply(seq(resamples), FUN = function(g){
    d<-mat
    # Reshuffle gene names
    d$gene1      <-  sample(d$gene1,length(d$gene1))
    names(d)[-1] <-  as.character(d$gene1)

    #names(d)[-1] <- sample(names(d)[-1],length(names(d)[-1]))
    #d$gene1      <- sample(d$gene1,length(d$gene1))


    # Score sampled data
    scores.sampled<-lapply(names(dictionary$bp),FUN = function(b){

      genes.in.bp<-dictionary$bp[[b]][dictionary$bp[[b]] %in% names(mat)]
      if(length(genes.in.bp)!=0){
        a<-d[which(d$gene1 %in% dictionary$bp[[b]]),genes.in.bp,with = FALSE]

        result<-data.table(bp = b,score = sum(a))
      }else{
        result<-NULL
      }
      result
    }) %>% rbindlist %>% arrange(bp) %>% as.data.table

    # Return result
    scores.sampled

  })

  #names(reruns)<-seq(10)
  reruns %<>% purrr::reduce(left_join, by = "bp")
  names(reruns)[-1]<-paste("score",seq(1000),sep = '')
  reruns %<>% as.data.table
  #reruns %<>% rbindlist(idcol = "bp")
  #sample.results<-spread(data = reruns,key = bp,value = score)
  #reruns %<>% rbind(scores.observed) %>% as.data.table

  all.data <- merge(scores.observed,reruns, by = "bp") %>% as.data.table
  dt<-all.data %>% split(f = seq(nrow(all.data)))

  processed     <- dt %>% lapply(FUN = function(f){
    process     <- unlist(f)[1]
    obs         <- unlist(f)[2] %>% as.numeric
    scores      <- unlist(f)[-c(1:2)] %>% as.numeric

    p.value     <- pnorm(
      q                = abs(obs),
      mean             = mean(scores),
      sd               = sd(scores))
    effect.size <- qnorm(
      p                = p.value,
      mean             = mean(scores),
      sd               = sd(scores))
    data.table(process      = process,
               obs          = obs,
               sampled.mean = mean(scores),
               sampled.sd   = sd(scores),
               p.value      = p.value)

  }) %>% rbindlist %>% arrange(p.value) %>% as.data.table


  processed %<>% mutate(process      = factor(process,unique(process)),
                        BH.corrected = stats::p.adjust(processed$p.value,
                                                       method = "BH")) %>% as.data.table
  return(processed)

}



