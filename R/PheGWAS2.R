## Function to add gene for repective rsid, if genes are not provided by user
addgene <- function(gwasmulti){

  ensembl <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")

  human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

  gwasmulti$gene <- "NA"
  for(i in rownames(gwasmulti)){
    xx <- getBM(attributes=c( "ensembl_gene_stable_id"),filters="snp_filter", values=gwasmulti[i, "SNP"],mart=ensembl, uniqueRows=TRUE)
    if (nrow(xx) == 0){
      gwasmulti[i,]$gene <- "NA"
    }else{
      gene <- getBM(attributes = c("hgnc_symbol"),
                    filters = "ensembl_gene_id", values = unlist(xx), mart = human)
      gwasmulti[i,]$gene <- if (nrow(gene) == 0) "NA" else  toString( unlist(gene))
    }}
  print("Finished mapping")
  return(gwasmulti)
}
######################################### ######################################### #########################################
############## FAST PROCESSPHEGWAS - THIS IS USED FOR PROCESISNG THE PHENOTYPE BEFORE THE LANDSCAPE FUNCTION
######################################### ######################################### #########################################
#' Prepare the dataframe to pass to landscape function
#'
#' @import tidyverse
#' @import tidyr
#' @importFrom data.table setDT
#' @param phenos Vector of names of dataframes that need to do PheGWAS on. Arrange the dataframe in the order how the the phenotypes should align in y axis
#' @param LDblock If want to pass a custom LDblock file for division of BP groups (applicable only for chromosomal level)
#' @return A processed dataframe to pass to PheGWAS landscape function
#' @details Make sure there are no duplicate rsid's in any of the dataframe, If there aremake sure to resolve it before passing it to this function.
#' @author George Gittu
#' @examples
#' phenos <- c("HDL", "LDL", "TRIGS", "TOTALCHOLESTROL")
#' ## y is ready to be passed to function landscape
#' y <- fastprocessphegwas(phenos)
#' @export
fastprocessphegwas <- function(phenos,LDblock= FALSE,LDpop= "eur"){
  if(is.vector(phenos) & !is.list(phenos)) {
    list.df_pre = mget(phenos,envir = .GlobalEnv)
    for (df in 1:length(list.df_pre)){
      list.df_pre[[df]] <- list.df_pre[[df]] %>%  rename(BETA = beta, SE = se)
    }
  }else{
    phenos <- names(list.df_pre)
  }
  action3 = . %>% mutate(P.value = as.double(P.value))
  list.df_pre <- lapply(list.df_pre, action3)
  if(LDblock){
    print("Processing for LDBlocks passed by the user")
    bpp = read.table(LDblock,header=TRUE,sep = "\t")
  }else if(LDpop == "asn"){
    print("Processing for AFR LDBlocks")
    bpp = bppasn
  }else if(LDpop == "afr"){
    print("Processing for ASN LDBlocks")
      bpp = bppasn
  }else {
      bpp = bppeur
  }
  colnames(bpp)[1] <- "CHR"
  bpp$CHR <-as.integer(gsub("chr","",bpp$CHR))
  bpp$start <- as.integer(bpp$start)
  bpp$stop <- as.integer(bpp$stop)
  bpp <- na.omit(bpp)
  ##getting all the minimum from each label
  labelfunction <- function(xx,bpp){
    setDT(xx)[, label :=
                setDT(bpp)[setDT(xx), on=.(CHR, start<=BP, stop>BP), paste(CHR, x.start, sep="_")]
              ]
  }
  #this list.df shoudnt be changed
  list.df_pre <- lapply(list.df_pre, labelfunction,bpp =bpp)

  # THis is addeed to the enviornment as we need this in nthe landscape functionn
  assign('list.df_pre',list.df_pre,envir=.GlobalEnv)
  # Following processing the dataframe for the labels
  action = . %>% group_by(label) %>%  slice(which.min(P.value))
  list.df <- lapply(list.df_pre, action)
  list.df <- Map(function(x, nm) {i1 <- !(names(x) %in% c('label','CHR'))
  names(x)[i1] <- paste0(names(x)[i1], "_", nm)
  x
  }, list.df, names(list.df))
  out <- Reduce(function(...) merge(..., by = c('CHR','label'), all = TRUE), list.df)
  names(out)
  for(i in 1:length(phenos))
  {
    a <- grep(paste0("^", phenos[i], "$", collapse = "|"), sapply(strsplit(names( out ), split="\\_"), tail, 1L))
    selection <- colnames(out)[a]
    hju <- paste0(phenos[i])
    out <- unite_(out, hju, selection, sep = "and", remove = TRUE)
  }
  out %>% gather(color,Entire_Val,-CHR, -label) -> gwasmulti.meltF

  if(length(grep("gene", colnames(list.df_pre[[1]]))) == 0){
  gwasmulti.melt <- gwasmulti.meltF %>% separate(Entire_Val, c("BP", "SNP","A1","A2","BETA","SE","P"), "and")
  ## Getting the dataframe that can be used
  d <- data.frame(CHR = as.numeric(gwasmulti.melt$CHR), BP = as.numeric(gwasmulti.melt$BP),A1 =gwasmulti.melt$A1,A2 =gwasmulti.melt$A2,SNP =gwasmulti.melt$SNP,
                  P = as.numeric(gwasmulti.melt$P),BETA = as.numeric(gwasmulti.melt$BETA),SE = as.numeric(gwasmulti.melt$SE),
                  PHENO = gwasmulti.melt$color,label = gwasmulti.melt$label)

  }else{
  gwasmulti.melt <- gwasmulti.meltF %>% separate(Entire_Val, c("BP", "SNP","A1","A2","BETA","SE","P","gene"), "and")
  ## Getting the dataframe that can be used
  d <- data.frame(CHR = as.numeric(gwasmulti.melt$CHR), BP = as.numeric(gwasmulti.melt$BP),A1 =gwasmulti.melt$A1,A2 =gwasmulti.melt$A2,SNP =gwasmulti.melt$SNP,
                  P = as.numeric(gwasmulti.melt$P),BETA = as.numeric(gwasmulti.melt$BETA),SE = as.numeric(gwasmulti.melt$SE),gene = gwasmulti.melt$gene,
                  PHENO = gwasmulti.melt$color,label = gwasmulti.melt$label)
  }
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P) ))
  d <- d[order(d$CHR, d$BP), ]
  d$logp <- round(-log10(d$P),3)
  d$BETA <- round(d$BETA,3)
  d$logp[is.infinite(d$logp)] <- max(d$logp[!is.infinite(d$logp)],
                                     na.rm = TRUE) + 10
  # This dataframe is what used for the landscape function.
  d
}


######################################### ######################################### #########################################
############## FAST PROCESSLANDSCAPE - THIS IS USED FOR VIWEING THE LANDSCAPE
######################################### ######################################### #########################################
#' Interactive 3-D association landscape for many phenotypes
#'
#' @import tidyverse
#' @import tidyr
#' @import plotly
#' @import reshape2
#' @importFrom biomaRt useMart getBM
#' @importFrom httr GET stop_for_status content_type content
#' @importFrom utils read.table
#' @import xml2
#' @import jsonlite
#' @importFrom stringr str_trim
#' @param d DataFrame output from processphegwas
#' @param sliceval Integer to indicate value of -log10(p) to do the sectionalcut. Usually value > -log10 6 is considered to be significant
#' @param chromosome Integer to indicate the chromosome number thats interested, If not given entire chromosome is given
#' @param geneview This checks for the common genes across the section
#' @param LDblock If want to pass a custom LDblock file for division of BP groups (applicable only for chromosomal level)
#' @param calculateLD This shoudld be set to true if the calcualte LD logic needed to be added to the plot
#' @param pop The population to select for calculation the LD (default GB)
#' @param R2 The value to set to calculate LD
#' @param D THe value to set to calcualte LD
#' @param mutualLD Calcualte the mutual LD SNP between the phenotypes
#' @param levelsdown Used to find the independant signals
#' @param phenos Vector of names of dataframes that need to do PheGWAS on. Arrange the dataframe in the order how the the phenotypes should align in y axis
#' chromosome view the max peak is selected
#' @author George Gittu
#' @examples
# \dontrun{
#' phenos <- c("HDL", "LDL", "TRIGS", "TOTALCHOLESTROL")
#' ## pass the dataframe from the processphegwas
#' y <- fastprocessphegwas(phenos)
#'
#'
#' 3D landscape visualization of all the phenotypes across the base pair positions(above a threshold of -log10 (p) 6)
#' landscapefast(y,sliceval = 10,phenos =phenos)
#'
#' 3D landscape visualization of chromosome number 19 (above a threshold of -log10 (p) 10)
#' landscapefast(y,sliceval = 7.5,chromosome = 19,phenos =phenos)
#'
#' 3D landscape visualization of chromosome number 19, gene view active  (above a threshold of -log10 (p) 10)
#' landscape(y,sliceval = 7.5,chromosome = 19, geneview = TRUE,phenos =phenos)
#'
#' 3D visualization with LD block (for european population) passing externally, parameter to pass LD and also calculate the mutualLD block
#' landscapefast(y, sliceval = 30, chromosome = 19,calculateLD= TRUE,mutualLD = TRUE,phenos =phenos)
#}
#' @export
landscapefast <- function(d,sliceval = 7,chromosome = FALSE,pop = "GBR",R2 = 0.75,
                          D = 0.75,calculateLD = FALSE, mutualLD = FALSE,phenos,geneview = FALSE,levelsdown = 0){
  if(chromosome == FALSE){
    print("Processing for the entire chromosome")
    # this is for the MAtrix to process the entire chromosome
    gwasmultifull <- d %>%
      ### dont want to show the phenotypes thats doesnt have anything above a certain threshold
      group_by(PHENO) %>%
      filter(any(logp > sliceval)) %>%
      ungroup() %>%
      group_by(PHENO,CHR) %>%
      slice(which.max(logp))
    if(length(grep("gene", colnames(d))) == 0){
      print("Applying BioMArt module for matching gene to rsid")
      gwasmultifull <- addgene(gwasmultifull)
    }
    gwasmultifull <- gwasmultifull[!is.na(gwasmultifull$CHR) & !is.na(gwasmultifull$BP),]
    ## Selecting only chr 1 to 22 ignoring X and Y, phenos is the order that we wannt the traits to be ordered
    gwas_surface_full <- acast(gwasmultifull, gwasmultifull$PHENO ~ gwasmultifull$CHR, value.var = "logp")[phenos,1:22]
    ## To remove the cases which dont have values fro all the snps
    gwas_surface_full <- gwas_surface_full[complete.cases(gwas_surface_full), ]
    gwas_surface_prime_full_use <- gwas_surface_full_use <- gwas_surface_full
    upperlimit = max(gwas_surface_full_use)
    gwas_surface_prime_full_use[,]= upperlimit + 2
    gwas_surface_copy_full_use <- gwas_surface_full_use
    # preparing labels
    for (i in 1:nrow(gwas_surface_copy_full_use)) {
      for (j in 1:ncol(gwas_surface_copy_full_use))
        gwas_surface_copy_full_use[i, j] <-
          paste("Phenotype = ", rownames(gwas_surface_copy_full_use)[i], '\n',"Chromosome = ", colnames(gwas_surface_copy_full_use)[j],'\n',
                "BP =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$BP,"\n",
                "A1 =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$A1,' ',
                "A2 =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$A2,'\n',
                "Effect Size =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$BETA,' ',
                "SE =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$SE,'\n',
                "-log10 (P-value) = ", gwas_surface_copy_full_use[i, j],'\n',
                "SNPID =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$SNP,'\n',
                "Gene =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$gene)
    }

    p <- plot_ly(x=colnames(gwas_surface_full_use),y= rownames(gwas_surface_full_use),z = gwas_surface_full_use) %>%
      add_surface(cmin = sliceval,cmax = upperlimit,opacity=.95,surfacecolor = gwas_surface_full_use,text = gwas_surface_copy_full_use,hoverinfo = "text") %>%
      layout(title = 'PheGWAS',scene = list(camera = list(eye = list(x=2, y=2, z=2)),yaxis = list(title = "Phenotypes", tickmode = "array",nticks = 8),xaxis =list(title= "Chromosomes",autotick = F, dtick = 1),
                                            zaxis =list(title= "-log 10(P)",range=c(sliceval,upperlimit+2)),aspectmode = "manual",aspectratio = list( x = 2, y = .8, z = .8))) %>%
      add_trace(x=colnames(gwas_surface_prime_full_use),y= rownames(gwas_surface_prime_full_use),z = gwas_surface_prime_full_use, type = "surface",surfacecolor = gwas_surface_full_use,showscale = FALSE,
                text = gwas_surface_copy_full_use,hoverinfo = "text",cmin = sliceval,cmax = upperlimit) %>% colorbar(title = "-log10 (P-value)")
    p
  }else{
    print(paste0("Processing for chromosome ",chromosome))
    whileivar = 0
    while (whileivar <= levelsdown) {
      dchrom1initial <- d[d$CHR==chromosome,] %>% rename(lab = label)
      gwasmulti <- dchrom1initial %>%
        group_by(lab) %>%
        filter(any(logp > sliceval)) %>%
        ungroup() %>%
        group_by(lab, PHENO) %>%
        slice(which.max(logp))

      if(length(grep("gene", colnames(d))) == 0){
        print("Applying BioMArt module for matching gene to rsid")
        gwasmulti <- addgene(gwasmulti)
      }

      ###rewriting calculating ld logic
      if (calculateLD){
        gwasmulti$snpsinld <- "NA"
        gwasmulti$snpsinldcount <-0
        gwasmulti$snpsinldup <- "NA"
        z<-NULL
        for (i in 1:nrow(gwasmulti)) {
          if(gwasmulti[i,]$logp > sliceval){
            stringcreate  <- sprintf("/1000GENOMES:phase_3:%s?",pop)
            server <- "http://grch37.rest.ensembl.org"
            ext <- paste0("/ld/human/",gwasmulti[i,]$SNP,stringcreate)
            r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
            tryCatch({
              stop_for_status(r)
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            lddf <- fromJSON(toJSON(content(r)))
            ldlist <- unlist(lddf$variation2)
            if(!is.null(ldlist)){
              ldlist <- unlist(lddf[lddf$r2 >= R2 | lddf$d_prime >= D,]$variation2)
              var <- eval(parse(text=paste0("list.df_pre$",gwasmulti[i,]$PHENO)))
              hhh <-  var[var$label == gwasmulti[i,]$lab & var$P.value <= as.numeric(10^-sliceval)]
              snpsinld <- intersect(hhh$rsid,ldlist)
              snpsinldup <- length(snpsinld) / nrow(hhh)
              gwasmulti[i,]$snpsinld <- paste0(toString(snpsinld))
              gwasmulti[i,]$snpsinldcount <- length(snpsinld)
              gwasmulti[i,]$snpsinldup <- paste0(snpsinldup)
            }
          }
        }
        if(mutualLD){
          print("Calculating the mutually shared SNP's")
          gwasmulti$common <- "NA"
          for(i in unique(gwasmulti$lab)){
            for (j in 1:length(phenos)) {
              for (k in phenos[-j]) {
                uuu <- strsplit(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == phenos[j], ]$snpsinld,",")
                vfg <- strsplit(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == k, ]$snpsinld,",")
                vfg <- lapply(vfg, str_trim)
                uuu <- lapply(uuu, str_trim)
                unld <- intersect(append(unlist(uuu),toString(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == phenos[j], ]$SNP)), append(unlist(vfg), toString(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == k, ]$SNP)))
                if(length(unld)== 0){
                  unld <- NA
                }
                z <- c(z,paste(phenos[j],"with",k,":" ,paste0(toString(unld)), collapse = "," ))
              }
              gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == phenos[j],]$common <- paste(toString(z))
              z <- NULL
            }
          }
        }
      }
      if (levelsdown != 0 & calculateLD) {
        for (ipheno in unique(gwasmulti$PHENO)) {
          bgt <- strsplit(gwasmulti[gwasmulti$PHENO ==
                                      ipheno, ]$snpsinld, ",")
          xxxx <- lapply(bgt, str_trim)
          ppp <- union(unlist(xxxx), gwasmulti[gwasmulti$PHENO ==
                                                 ipheno, ]$SNP)
          xxlist <- as.list(ppp)
          list.df_pre[[ipheno]] <- list.df_pre[[ipheno]][!list.df_pre[[ipheno]]$rsid %in% xxlist, ]
        }
        assign('list.df_pre',list.df_pre,envir=.GlobalEnv)
        d <- fastprocessphegwas(list.df_pre)
      }
      whileivar = whileivar + 1
    }
    gwas_surface <- acast(gwasmulti, gwasmulti$PHENO ~ gwasmulti$lab, value.var = "logp")[phenos,]
    gwas_surface_use <- gwas_surface
    gwas_surface_copy_use <- gwas_surface
    gwas_surface_copy_gene_use <- gwas_surface_use
    colnames(gwas_surface) <- gsub(paste0(chromosome,"_"), "", colnames(gwas_surface), fixed = TRUE)
    yy <- order(as.numeric(colnames(gwas_surface)))
    gwas_surface <- gwas_surface[,yy]
    gwas_surface_copy_use <- gwas_surface_copy_use[,yy]
    gwas_surface_use <- gwas_surface_use[,yy]
    gwas_surface_copy_gene_use <- gwas_surface_copy_gene_use[,yy]

    for (i in 1:nrow(gwas_surface_copy_use)) {
      for (j in 1:ncol(gwas_surface_copy_use))
        if(calculateLD){
          gwas_surface_copy_use[i, j] <-
            paste("Phenotype: ", rownames(gwas_surface_copy_use)[i], ' ',"kbp: ", colnames(gwas_surface_use)[j],' ',
                  "Effect Size: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_use)[j],]$BETA,'\n',
                  "SE: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SE,' ',
                  "-log10 (P-value): ", gwas_surface_copy_use[i, j],' ',
                  "SNPID: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SNP,"\n",
                  "Gene =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$gene,"\n",
                  gsub('(.{1,90})(\\s|$)', '\\1\n', paste("SNP's in LD: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$snpsinld)),"\n",
                  "Linked SNP's ratio: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$snpsinldup,"\n",
                  gsub('(.{1,90})(\\s|$)', '\\1\n', paste("With other pheno: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$common)))
        }else{
          gwas_surface_copy_use[i, j] <-
            paste("Phenotype = ", rownames(gwas_surface_copy_use)[i], '\n',"kbp position range= ", colnames(gwas_surface_use)[j],'\n',
                  "BP =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$BP,"\n",
                  "Effect Size =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_use)[j],]$BETA,' ',
                  "SE =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SE,'\n',
                  "-log10 (P-value) = ", gwas_surface_copy_use[i, j],'\n',
                  "SNPID =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SNP,'\n',
                  "Gene =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$gene)
        }
    }
    gwas_surface_use <- gwas_surface
    upperlimit = max(gwas_surface_use, na.rm = TRUE)
    gwas_surface_prime_use <- gwas_surface_use
    gwas_surface_prime_use[,]=  upperlimit + 2
    colorheatmap = NULL

    if (geneview == TRUE) {
      print("GENE View is active")
      for (i in 1:nrow(gwas_surface_copy_gene_use)) {
        for (j in 1:ncol(gwas_surface_copy_gene_use)) gwas_surface_copy_gene_use[i,
                                                                                 j] <- paste(gwasmulti[gwasmulti$PHENO == rownames(gwas_surface_copy_gene_use)[i] &
                                                                                                         gwasmulti$lab == colnames(gwas_surface_copy_gene_use)[j],
                                                                                                       ]$gene)
      }
      xxx <- apply(gwas_surface_copy_gene_use, 2, function(x) unique(unlist(strsplit(x[!is.na(x)],
                                                                                     ";"))[duplicated(unlist(strsplit(x[!is.na(x)],
                                                                                                                      ";")))]))
      for (i in 1:nrow(gwas_surface_copy_gene_use)) {
        for (j in 1:ncol(gwas_surface_copy_gene_use)) {
          x <- gwas_surface_copy_gene_use[i, j]
          xx <- unlist(strsplit(x[!is.na(x)], ";"))
          if (length(intersect(xx, unlist(xxx[j]))) >
              0 & gwasmulti[gwasmulti$PHENO == rownames(gwas_surface_copy_gene_use)[i] &
                            gwasmulti$lab == colnames(gwas_surface_copy_gene_use)[j],
                            ]$logp > sliceval & !is.na(gwasmulti[gwasmulti$PHENO ==
                                                                 rownames(gwas_surface_copy_gene_use)[i] &
                                                                 gwasmulti$lab == colnames(gwas_surface_copy_gene_use)[j],
                                                                 ]$gene)) {
            gwas_surface_copy_gene_use[i, j] <- max(gwas_surface_use)
          }
          else {
            gwas_surface_copy_gene_use[i, j] <- sliceval
          }
        }
      }
      surfacecolourheat = gwas_surface_copy_gene_use
      colorheatmap = "Viridis"
    }
    else {
      surfacecolourheat = gwas_surface_use
    }

    p <- plot_ly() %>%
      layout(title = 'PheGWAS',scene = list(camera = list(eye = list(x=2, y=2, z=2)),yaxis = list(title = "Phenotypes",tickvals = rownames(gwas_surface_use))
                                            ,xaxis =list(title= "100 kbp divisions" ,
                                                         type = 'category',tickmode = "array",tickvals = colnames(gwas_surface_use), ticks = "outside"),
                                            zaxis =list(title= "-log 10(P)",range=c(sliceval,upperlimit + 2)),aspectmode = "manual",
                                            aspectratio = list( x = 2, y = 1, z = .8)))


    p <- p %>% add_surface(x=colnames(gwas_surface_use),y= rownames(gwas_surface_use),z = gwas_surface_use,opacity=.98,
                           surfacecolor = gwas_surface_use, cmin = sliceval ,cmax = upperlimit,
                           showscale = FALSE,text = gwas_surface_copy_use,hoverinfo = "text") %>%
      add_trace(x=colnames(gwas_surface_use),y= rownames(gwas_surface_use),z = gwas_surface_prime_use,
                type = "surface",surfacecolor = surfacecolourheat,showscale = FALSE,text = gwas_surface_copy_use,hoverinfo = "text",colorscale = colorheatmap)
    p
  }
}

######################################### ######################################### #########################################
############## FAST PROCESSPHEGWAS - THIS IS USED FOR PROCESISNG THE PHENOTYPE BEFORE THE LANDSCAPE FUNCTION
######################################### ######################################### #########################################
#' Prepare the dataframe to pass to landscape function
#'
#' @import tidyverse
#' @import tidyr
#' @import reshape2
#' @import factoextra
#' @import cluster
#' @import seriation
#' @importFrom data.table fread rbindlist
#' @param phenos Vector of names of dataframes that need to do iPheGWAS on.
#' @param dentogram to show structural differences
#' @return A processed dataframe to pass to PheGWAS landscape function
#' @details Make sure there are no duplicate rsid's in any of the dataframe, If there aremake sure to resolve it before passing it to this function.
#' @author George Gittu
#' @examples
#' \dontrun{
#' phenos <- c("HDL", "LDL", "TRIGS", "TOTALCHOLESTROL")
#' iphegwas(phenos)
#' iphegwas(phenos,dentogram = TRUE)
#' }
#' @export
iphegwas <- function(phenos,dentogram = FALSE){
list.df_pre = mget(phenos,envir = .GlobalEnv)
for (df in 1:length(list.df_pre)){
  if(length(grep("Z", colnames( list.df_pre[[df]]))) == 0){
    action3 = . %>% subset(rsid %in% skeltonsnps$rsid) %>% merge(skeltonsnps, by= "rsid", how='inner') %>%
      mutate(ZZ = beta/se)  %>%
      mutate(Z = ifelse(toupper(A1.x) == toupper(A1.y), ZZ * -1,ZZ)) %>% rename(FEATURE =rsid) %>%
      select("FEATURE","Z") %>%
      distinct(FEATURE, .keep_all= TRUE)
  list.df_pre[[df]] <-  list.df_pre[[df]] %>% action3()
  }else{
    action3 = . %>% drop_na(rsid) %>% subset(rsid %in% skeltonsnps$rsid) %>% merge(skeltonsnps, by= "rsid", how='inner') %>%
      mutate(Z = ifelse(toupper(A1.x) == toupper(A1.y), Z * -1,Z)) %>%rename(FEATURE =rsid) %>%
      select("FEATURE","Z") %>%
      distinct(FEATURE, .keep_all= TRUE)
    list.df_pre[[df]] <-  list.df_pre[[df]] %>% action3()
  }
  colnames(list.df_pre[[df]]) <- toupper(colnames(list.df_pre[[df]]))
}
xfull <- rbindlist(list.df_pre,idcol = "PHENOS")
gwas_surface_pre <- xfull[,c("PHENOS","FEATURE","Z")]
gwas_surface <- acast(gwas_surface_pre, PHENOS~FEATURE, value.var="Z") #>>>>>>>>>>>>>Z ang logcal
gwas_surface <- gwas_surface[ , colSums(is.na(gwas_surface)) == 0]
res.dist <- get_dist(gwas_surface, method = "pearson") #>>>>>>>>>>>>>>>>>>
correlation.coefficent <- as.matrix(1 - res.dist)[phenos,phenos]
A <- abs(correlation.coefficent)
AA <- as.dist(1- A)
res.hc1 <- hclust(d = AA,method = "mcquitty") #>>
res.hc1 <- reorder(res.hc1, AA, method = "OLO")
if(dentogram){
 fviz_dend(res.hc1, # Cut in four groups
                           cex = 2, # label size
                           k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
                           color_labels_by_k = TRUE, # color labels by groups
                           rect = TRUE,
                           main = "title"# Add rectangle around groups,
  )
}else{
  phenos[res.hc1$order]
}
}

######################################### ######################################### #########################################
############## FAST PROCESSPHEGWAS - THIS IS USED FOR PROCESISNG THE PHENOTYPE BEFORE THE LANDSCAPE FUNCTION
######################################### ######################################### #########################################
#' Prepare the dataframe to pass to landscape function
#' @import GenomicSEM
#' @param phenos Vector of names of dataframes that need to do iPheGWAS on.
#' @param dentogram to show structural differences
#' @return A processed dataframe to pass to PheGWAS landscape function
#' @details Make sure there are no duplicate rsid's in any of the dataframe, If there aremake sure to resolve it before passing it to this function.
#' @author George Gittu
#' @examples
#' \dontrun{
#' phenos <- c("HDL", "LDL", "TRIGS", "TOTALCHOLESTROL")
#' filennames <- c("HDL", "LDL", "TRIGS", "TOTALCHOLESTROL")
#' hm3 = "/Users/ggeorge/Desktop/ldscR/w_hm3.snplist"
#' ld = "/Users/ggeorge/Desktop/ldscR/MIX/eur_w_ld_chr"
#' N = c(53293,46350,51710,143677)
#' correlation.coeff <- ldscmod(phenos,filenames,N,hm3,ld)
#' }
#' @export
ldscmod <- function(phenos,filenames,N,hm3,ld){
  setwd(tempdir(check = TRUE))
  print(getwd())
  f <- list.files(getwd(), include.dirs = F, full.names = T, recursive = T)
  file.remove(f)
  for (i in 1:length(phenos)){
    munge(files =filenames[i],hm3 = hm3,N =N[i],
          trait.names =phenos[i])
  }
  filenames <- list.files(getwd(), pattern ="*.sumstats.gz", full.names=TRUE)
  f <- path_file(filenames)
  p <- sub('\\.sumstats.gz$', '', f)
  mm <- capture.output(ldsc(traits =f,
                            sample.prev = rep(NA, length(filenames)),
                            population.prev = rep(NA, length(filenames)) ,ld = ld,wld = ld,trait.names =p))
  preline<- grep("Genetic Correlation Results", mm, fixed = TRUE)
  postline<- grep("LDSC finished running", mm, fixed = TRUE)
  filedata <- mm[(preline+1):(postline-1)]
  filedata <- gsub('^.*between\\s*|\\s*\\(.*$', '', filedata)
  filedata <- sub( ":", "", filedata, fixed = TRUE)
  filedata <- sub( " and", "", filedata, fixed = TRUE)
  filedata_list<- strsplit(filedata, split = " ")
  new <- Reduce(rbind, filedata_list)
  names <- unique(c(new[,1], new[,2]))
  corrmat <- matrix(nrow =length(names),  ncol = (length(names)), dimnames = list(names, names))
  for (i in 1:length(filedata_list)){
    corrmat[filedata_list[[i]][1],filedata_list[[i]][2]] <- as.numeric(filedata_list[[i]][3])
    corrmat[filedata_list[[i]][2],filedata_list[[i]][1]] <- as.numeric(filedata_list[[i]][3])
  }
  diag(corrmat) <- 1
  return(corrmat)
}
