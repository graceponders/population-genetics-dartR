#Crayfish population genetic structure with dartR
    #Grace, October 2025 

#Analyzing SNPs on 80 individual crayfish from 4 locations.
#Using DartR (now dartRverse) following the iccb github tutorials.
#This code has been streamlined. For initial runs, troubleshooting and experimenting with the data, see "Mega_DartR_Script"

# 1: Setup and install -------------------------------------------------------
  
        # core packages:
    pacman::p_load (dartRverse,
                    dartR.base)
        
        #check available packages (all should be available): 
        dartRverse::dartRverse_install()
        
        #Extra packages
        pacman::p_load(ggplot2, 
                       tidyverse,
                       adegenet,
                       HardyWeinberg, 
                       pegas,
                       ape,
                       reshape2, 
                       ggtree, 
                       LEA, 
                       poppr,
                       hierfstat)
                    
# 2: Importing and handling genomic data -------------------------------------
        
      #following along with the tutorial https://green-striped-gecko.github.io/iccb2025/session01.html
      
      gl <- gl.read.dart("SNP.csv")
      
      #inspect.
      head(gl@other$ind.metrics)
      indNames(gl)[1:20]
      
      #need to assign populations 
      
      # Assign populations based on the prefix in sample names
      pop(gl) <- factor(gsub("", "", substr(indNames(gl), 1, 2)))
      
      # Check the new population assignments
      table(pop(gl))
      
      
  ###Filtering the gl ###
      #filter by callrate 
      gl_f <- gl.filter.callrate(gl, method = "loc", threshold = 0.75, verbose=0)
      
      #filter by individual callrate - this is a separate object as it removes the CD group.
      gl_fr <- gl.filter.callrate(gl_f, method = "ind", threshold = 0.75, verbose=0)
      
      #filter for monomorphs
      gl_f <- gl.filter.monomorphs(gl_f)
      gl_fr <- gl.filter.monomorphs(gl_fr)
      
      #filter for reproducibility
      gl_f <- gl.filter.reproducibility(gl_f, threshold = 0.98)
      gl_fr <- gl.filter.reproducibility(gl_fr, threshold = 0.98)
      
      #Filter for minor alleles, MAF:
      gl_f <- gl.filter.maf(gl_f, threshold = 0.05)
      gl_fr <- gl.filter.maf(gl_fr, threshold = 0.05)
      
      #inspect
      table(pop(gl_f))
      table(pop(gl_fr))
      gl.report.callrate(gl_fr)
      gl.compliance.check(gl_fr)
      
      #To recap:we now have 2 gl objects.
        #gl_f includes all 80 individuals - this is for neighbor joining trees, etc 
          #where I want to see the Cd outgroup for scale
        #gl_fr has also been filtered by individual call rate, removing the Cd outgroup. 
          #this is for pairwise analysis
            
            #heterozygosity reports; slightly lower heterozygosity; Inbreeding (FIS) is positive. 
            gl.report.heterozygosity(gl_f)
            gl.report.heterozygosity(gl_fr)
          
            
            #  After filtering, gl_fr contains 70 individuals and 7641 loci. 
                  #           gl_f contains 80 individuals and 13869 loci.
            
            #Reader note on filtering parameters: I ran this with many different filtering combinations with little to no impact on the output. 
              #Likely because each population is so distinct. 

# 3: Visualize structure with basic PCA ---------------------------------------------------------

      # Run PCA
      pca_f <- glPca(gl_f, nf=10)
      pca_fr  <- glPca(gl_fr, nf = 10)
      
      # Plot PCA
      s.class(pca_f$scores, pop(gl_f), col = rainbow(length(levels(pop(gl_f)))))
      s.class(pca_fr$scores, pop(gl_fr), col = rainbow(length(levels(pop(gl_fr))))) 
      
      #populations are very distinct, some RB adults have been misidentified.  
      #juveniles are split between RL and DC, where they were collected. 
      #Reader note: two new objects were made by manually reclassifying the juveniles by their collection site: gl2_f for all individuals and gl2_fr for robustus only. 
     
# 4: Sex linked markers ---------------------------------------------------

      #sex data
      sexlist <- read.csv("SEXLIST.csv")
      
      # add column named "id" to the gl
      gl_fr@other$ind.metrics$id <- indNames(gl_fr)
      gl_f@other$ind.metrics$id <- indNames(gl_f)
      names(gl_fr@other$ind.metrics)
      names(gl_f@other$ind.metrics)
      
      #combine sex info with a left join 
            gl_f@other$ind.metrics <- 
        left_join(gl_f@other$ind.metrics, sexlist, by = c("id" = "id"))
            
            gl_fr@other$ind.metrics <- 
              left_join(gl_fr@other$ind.metrics, sexlist, by = c("id" = "id"))

    #Run sex markers        
            sexlinked_report <- gl.report.sexlinked(gl_f, system = "zw", verbose = 3)      
            head(sexlinked_report$summary)
            
            sexlinked_report_fr <- gl.report.sexlinked(gl_fr, system = "zw", verbose = 3)      
            head(sexlinked_report_fr$summary)
            
        #No sexlinked markers detected (might need known sex >9 individuals)

# 5: Hardy-Weinberg equilibrium ----------------------------
         # usually you can use *gl.report.hwe* to do this easily, but ggtern is not updated, so I am using the pegas package instead.
            
            # Convert to GI
            gi <- gl2gi(gl_fr)
            
            # Run HWE for all populations
            hwe_results <- sapply(levels(pop(gi)), function(p) {
              hw.test(gi[pop(gi) == p, ], B = 0)
            }, simplify = FALSE)
            
            # Convert to data frame
            hwe_df <- do.call(rbind, lapply(names(hwe_results), function(p) {
              data.frame(
                locus = rownames(hwe_results[[p]]),
                pop = p,
                chi_squared = hwe_results[[p]][, "chi^2"],
                pval = hwe_results[[p]][, "Pr(chi^2 >)"]  )  }))
            
            # Remove NAs and apply FDR correction per population
            hwe_df <- hwe_df[!is.na(hwe_df$pval), ]
            hwe_df$pval_fdr <- ave(hwe_df$pval, hwe_df$pop, 
                                   FUN = function(x) p.adjust(x, method = "fdr"))
            
            # Check results 
            sum(hwe_df$pval_fdr < 0.05)  # total significant loci
            table(hwe_df[hwe_df$pval_fdr < 0.05, "pop"])  # failures per population
            
            #most failures are in the juvenile group
            #EXPECTED due to Wahlund effect (they were collected from 2 sites) 
            #remove juvs to create an "adult robustus" dataset (for pairwise analyses)
            
            gl_fra <- gl.keep.pop(gl_fr, pop.list = c("RB", "RL", "DC", "CO"))
            
            # Rerun HWE analysis
            hwe_df_adults <- hwe_df[hwe_df$pop != "JV", ]
            hwe_df_adults$pval_fdr <- ave(hwe_df_adults$pval, hwe_df_adults$pop, 
                                          FUN = function(x) p.adjust(x, method = "fdr"))
            
            sig_loci_adults <- hwe_df_adults[hwe_df_adults$pval_fdr < 0.05, ]
            loci_counts_adults <- table(sig_loci_adults$locus)
            print(table(loci_counts_adults))
            #No loci failures >2 sites. Clean.  

# 6: Population diversity metrics -----------------------------------------
            
            #Basic diversity report
            gl.report.diversity(gl2_fr)
            
            #Heterozygosity (observed vs expected)
            gl.report.heterozygosity(gl2_fr, method = "pop")
            
    #He vs Ho extracted into results table.
          
        #allelic richness
            #convert format genlight -> genind -> hierfstat 
            gi <- gl2gi(gl_fr)                          
            hf <- genind2hierfstat(gi, pop = pop(gi)) 
            
            #compute min.n (number of alleles to rarefy to).
                  #(default in allelic.richness is min genotyped * 2)
            n <- table(pop(gi))             # number of individuals per pop
            min_n <- 2 * min(n)             # diploid -> gene copies
            
            #allelic richness (rareified)
            AR <- allelic.richness(hf, min.n = min_n, diploid = TRUE)
            
            # AR$Ar is a matrix: rows = loci, cols = populations
            # per-population mean allelic richness across loci:
            AR_pop <- colMeans(AR$Ar, na.rm = TRUE)
            
            # overall mean / sd per-population 
            mean_AR <- mean(AR_pop)
            sd_AR   <- sd(AR_pop)
            
            # print tidy results
            print(AR_pop)
            print(mean_AR)
            print(sd_AR)

      #extract to results table       
            
            #pairwise private alleles 
            pa <- private_alleles(gi, report = "table", level = "population", count.alleles = TRUE)
            private_counts <- rowSums(pa > 0)
            private_counts
          
            
            #take all of this out into results table 
            
            
# 7: Pairwise Fst -------------------------------------------------------------

          # Calculate pairwise Fst
          fst <- gl.fst.pop(gl2_fr, nboots = 1000, nclusters = 1)
           
          # View the Fst matrix
          fst$Fsts
            
          # View the p-values
          fst$Pvalues
            
          
            #CO is the most isolated, but there is very high FST values everywhere. 
           #This indicates very limited gene flow between populations 
          
    # convert genetic distance to time since divergence 
    
          # Calculate Nei's genetic distance
          dist_nei <- gl.dist.pop(gl2_fr, method = "nei")
          
          # View distance matrix
          as.matrix(dist_nei)
          
          # Estimate divergence time
          # Nei's D = 2 × mutation rate × time
          # time = D / (2 × mutation rate)
          
          mutation_rate <- (1.3e-8)  #run with 7e−9, 1e-8, and 1.3e-8 
          loci <- 7641
          
          # Calculate divergence time
          divergence_time <-  as.matrix(dist_nei) / (2 * mutation_rate * loci)
    
          # View distance matrix
          as.matrix(divergence_time)
          
          
          
          #Checking if these numbers make sense by playing with effective population size:
          #Formula: T = -2 Ne ln(1-Fst)
          
          FSTs <- as.matrix(fst$Fsts)
          Ne <- 2000
          fst_div_times <- -Ne*log(1-FSTs)
          
          #Played around with this for a bit but Nes look feasible. Will calculate these properly later in the code. 
         
# 8: Ancestry -------------------------------------------------------------
           
   ###using sNMF, similar to STRUCTURE but faster
       
           # Convert genlight to geno format
           gl2geno(gl2_fr, outfile = "structure", outpath = getwd())
          
          #K should = 4 but we will confirm 
          testK <- snmf("structure.geno", K = 1:6, repetitions = 10, entropy = TRUE, project = "new") #this is computationally intense
          
          # Find best K (lowest cross-entropy)
          ce <- numeric(6)
          for(k in 1:6) {
            ce[k] <- min(cross.entropy(testK, K = k))}
          
          print(ce)
          #yes, best k is 4
          
          #run snmf with k=4
           str <- snmf("structure.geno", K = 4, entropy = TRUE, repetitions = 10)
           
           # Select best run
           best_run <- which.min(cross.entropy(str, K = 4))
           
           # Extract Q matrix
           Q <- Q(str, K = 4, run = best_run)
           
           # Combine Q with metadata
           df <- cbind(ind = indNames(gl2_fr),
                       pop = pop(gl2_fr),
                       as.data.frame(Q))
           
           # Long format for ggplot
           df_long <- tidyr::pivot_longer(df, cols = starts_with("V"),
                                          names_to = "Cluster",
                                          values_to = "Ancestry")
           
           # Order individuals by population
           df_long$ind <- factor(df_long$ind, levels = df$ind[order(df$pop)])
           
           cluster_colors <- c("V1"="#C09940","V2"="#7B997B",
                    "V3"="#5392A3","V4"="#706073")
           
           # Plot
           ggplot(df_long, aes(ind, Ancestry, fill = Cluster)) +
             geom_col(width = 0.8) +
             scale_fill_manual(values = cluster_colors) +
             facet_grid(~ pop, scales = "free_x", space = "free") +
             theme_void() +
             ylab("Ancestry") +
             xlab("Individuals")
           
          #bar chart of admixture
           


# 9: Neighbor joining tree -----------------------------------------------
    
           # Calculate pairwise distances between ALL INDIVIDUALS
           dist <- gl.dist.ind(gl_f, method = "euclidean")
           
           # Build NJ tree from individual distances
           tree <- nj(dist)
           
           # Check - should now have ~80+ tips (all individuals)
           tree
           tree$tip.label  # should show individual names like CD1, CD2, CO1, etc.
           
       
           pop_colors <- c("JV" = "#E77575", "DC" = "#D09940", "RB" = "#5392A3",
                           "RL" = "#7B997B", "CO" = "#706073", "CD" = "#98453A")
           
           ind_colors <- pop_colors[as.character(pop(gl_f))]
           
           # Plot tree WITHOUT tip labels, but WITH colored points
           plot(tree,
                type = "unrooted",
                show.tip.label = FALSE,  # hide labels
                edge.width = 1.5,
                no.margin = TRUE)   
           
           # Add colored points at tips
           tiplabels(pch = 19,          # filled circle
                     col = ind_colors,   # population colors
                     cex = 2.2,          # size of dots
                     frame = "none")     # no box around points
           # Add legend
           legend("topleft",
                  legend = names(pop_colors),
                  col = pop_colors,
                  pch = 19,
                  bty = "n",
                  cex = 1.3,
                  title = "Population")
           
           
           #done!
           