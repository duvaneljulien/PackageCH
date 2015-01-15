
####**********************************************************************
####  Written and Developed by: 
####**********************************************************************
####
####    Julien Duvanel, Copyright 2014
####    email: duvanel@stanford.edu
####
####**********************************************************************
####**********************************************************************

####**********************************************************************
####
####  Generate traits
####
####**********************************************************************

# TraitGenerator(plink.file = "/home/duvanel/BP/Data/GenoData1M/MendelAnalysis3",
#                export.path = "/home/duvanel/git/PackageCH/data/",
#                N = 5,
#                sbatch.file = "/home/duvanel/git/ClusterScript/ExtractSNPS.sbatch")

#' Extract randomly a few SNPs from genome and then generate fake traits
#' 
#' @title Generate fake traits
#' @param plink.file where are stored plink files
#' @param export.path where to export RData fake traits
#' @param N number of SNPs to extract from genome
#' @param sbatch.file where the sbatch file will be saved
#' @return nothing
#' @author Julien Duvanel 
#' @export
TraitGenerator <- function(plink.file, export.path, N = 20, sbatch.file = "ExtractSNPs.sbatch") {

    # datetime stamp (to save files)
    datetime.stamp <- format(Sys.time(), 
                             "%d%m%Y_%H%M%S")

    # Load plink bim data (contains SNPs list)
    bim.data <- read.table(file = paste0(plink.file, ".bim"))
    # Get the number of SNPs
    nbr.snps <- nrow(bim.data)
    
    # Get randomly N SNPs (their row id and name)
    rows.snps <- sample(x=seq(1, nbr.snps), size=N)
    name.snps <- bim.data$V2[rows.snps]  
    character.snps <- paste0(unlist(name.snps),collapse=",")
    
    # Write SBATCH and then execute it
    WriteSBATCH(plink.file, sbatch.file, character.snps, export.path, datetime.stamp)
    ExecuteSBATCH(sbatch.file)
    
    # Get results (stored into SNPS_datetime.stamp.raw)
    Phenotypes <- read.table(file = paste0(export.path, "SNPs_", datetime.stamp, ".raw"), header = TRUE)
    # Remove rows with NA's 
    Phenotypes <- Phenotypes[complete.cases(Phenotypes), ]
    # We must have "famid" and "id" as the two first columns
    colnames(Phenotypes) <- c("famid", "id")
    
    # Value of SNPs are stored in the last N columns
    snps <- Phenotypes[, (ncol(Phenotypes)-N+1):ncol(Phenotypes)]
    # Simulate alpha_i. The goal is to have
    # \sum_i alpha_i snps_i for i = 1,\dots,n
    alpha <- runif(N, -1, 1)
    
    # We create an alpha for the dominant effect
    alpha_d <- matrix(0, nrow = 3, ncol = N)
    alpha_d[1, ] <- runif(N, -1, 1)
    alpha_d[2, ] <- runif(N,  3, 4)
    alpha_d[3, ] <- runif(N,  15,20)
    
    nonlinear.matprod <- function(mat, alpha) {
        # Usually, as.matrix(mat) %*% alpha does the job
        # but here, the multiplication differs if mat[i,j] has 
        # a specific value.
        
        # For every row, we do the "special product"
        res <- apply(mat, 1, function(x) {
            ret <- 0
            for(i in 1:length(x)) {
                ret <- ret + x[i]*alpha[x[i]+1, i]
            }
            ret
        })
        
        as.matrix(res)
    }
    
    fake.add.trait <- as.matrix(snps) %*% alpha
    fake.dom.trait <- nonlinear.matprod(snps, alpha_d)
    
    # We create two traits. One is fully heritable and the other one not at all.
    P.raw <- cbind(Phenotypes[, 1:2], 
                   fake.add.trait, 
                   fake.add.trait + rnorm(n = nrow(Phenotypes), mean = 0, sd = 1),
                   fake.dom.trait,
                   fake.add.trait + fake.dom.trait,
                   rnorm(n = nrow(Phenotypes), mean = 0, sd = 1))
    
    # Save as P.raw (name is important)
    save(list = c("P.raw", "alpha"), file = paste0("data/Phenotype_Simulated_", datetime.stamp ,".RData"))
    
}

#' Write SBATCH
#'
#' @title Write SBATCH
#' @param plink.file where are stored plink files
#' @param sbatch.file where the sbatch file will be saved
#' @param character.snps list of SNPS (character separated with coma)
#' @param export.path where to export RData fake traits
#' @param datetime.stamp time stamp
#' @return nothing
#' @author Julien Duvanel
#' @export
WriteSBATCH <- function(plink.file, sbatch.file, character.snps, export.path, datetime.stamp) {
    
    # Write into a sbatch file
    sink(sbatch.file)
    cat("#!/bin/bash
        
        #set a job name
        #SBATCH --job-name=plinkExtractSNP
        
        #a file for job output, you can check job progress
        #SBATCH --output=plinkExtractSNP.out
        
        # a file for errors from the job
        #SBATCH --error=plinkExtractSNP.err
        
        #time you think you need; default is one hour
        #in minutes in this case
        #SBATCH --time=120:00
        
        #quality of service; think of it as job priority
        #SBATCH --qos=normal
        
        #number of nodes you are requesting
        #SBATCH --nodes=1
        
        #memory per node; default is 4000 MB per CPU
        #SBATCH --mem=4000
        
        #get emailed about job BEGIN, END, and FAIL 
        #SBATCH --mail-type=ALL 
        
        #who to send email to; please change to your email
        #SBATCH  --mail-user=duvanel@stanford.edu
        
        #task to run per node; each node has 16 cores
        #SBATCH --ntasks-per-node=16
        
        module load plink
        plink --bfile", plink.file," --snps ", character.snps, " --recodeA --out ", paste0(export.path, 
                                                                                           "SNPs_", 
                                                                                           datetime.stamp)
        ,"\n")
    sink()
}

#' Execute SBATCH
#'
#' @title Execute SBATCH
#' @param sbatch.file where the sbatch file is stored
#' @return nothing
#' @author Julien Duvanel
#' @export
ExecuteSBATCH <- function (sbatch.file) {
    
    # Give chmod +x
    system(paste0("chmod +x ", sbatch.file))
    
    # Execute using cluster
    system(paste0("sbatch ", sbatch.file))
    
    # We wait until it's finished (we hope so it is, in 10 sec)
    Sys.sleep(30)
    
}
