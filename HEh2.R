{
    args <- commandArgs(trailingOnly=T)

    GENE="APOL1"; #default gene to analyze
    INDIR="."; #default input directory is current working dir
    pref012=""; #prefix for the 012 file
    prefPOS=""; #prefix for the pos file
    suff012=""; #suffix for the 012 file
    suffPOS=""; #suffix for the pos file
    OUTDIR="."; #default output directory is current working dir
    OUTFILESUFF=""; #desired string to add to outfile
    APPEND="F";  #boolean whether to append to output file
    SNPLIST="";  #An optional file with particular SNPs to put into K bins (columns)
    SNPSET = array(,dim=0); #empty set to be filled if SNP sets are specified
    MINSNPS = 100; #minimum number of SNPs in a bin, otherwise merge adjacent
    WIDTH = "1000000"; #width of window around each gene
    QN = 1; #to perform quantile-normalization on phenotype or not
    Ktmp = 1; #number of MAF bins to include in analysis
    K.RANGE = c(1,2,5,10,20,50); #range of MAF bins to include in analysis
    maxMAC = -1; #maximum allele count to include in h2 inference; -1 = no max
    minMAC = 1; #minimum allele count to include in h2 inference
    minMAC.RANGE = c(1,2,5,36); #range of minimum allele count to include in h2 inference
    COVARCOLS = -1; #number of covariates to correct for, or list of columns
    COVARFILE = 0; #File with covariates in each column (individuals on rows)
    PHENOFILE = ""; #specify a phenotype file directly
    ngPC = 0; #number of genotype PCs to correct for 
    npPC = ngPC; #number of phenotype PCs to correct for
    DSAMP = 1.0;  #proportion of reads to downsample to for testing
    PERM = 0;  #permutation type (0,1,2,3,4,5)
    nPERM = 1;  #number of permutations to run
    nSIM = 1;  #number of simulations to run
    SIM = 0;  #simuation type (0=not a sim)
    Nc = 0;  #number of causal variants
    fracRare = 1;  #fraction of causal variants that are rare
    freqThresh = 0; #frequency threshold for causal variants
    causalmodel = 0;  #models for how causal variants are chosen
    betamodel = 0;  #models for how effect sizes are drawn
    betaALL = 0.05; #beta for all variants when betamodel=2
    h2 = 0.05;  #true value of h2 in simulation
    Kgen = 1;  #number of bins for generating causal variants
    Krem = 01;  #number of bins for remove from causal variants
    errorModel = 0;  #model for incorporating errors
    delta = 1.0;  #fraction of causal variants that have beta>0, only for betamodel=1,2,3.
    rho = 1.0;  #used in causalmodel=8 from Uricchio et al.
    tau = 1.0;  #used in causalmodel=8 from Uricchio et al.
    demog = "ten"; #can be either 'ten' for Tennessen or 'gravel' for Gravel
    dfe = "Boyko"; #can be either 'Boyko' or "Torg"
    errAdd = 0;  #rate of adding erroneous SNPs
    errRem = 0;  #rate of removing causal variants (erroneous)
    runHE = 1;  #indicator to run HE or not
    runLMM = 1; #indicator to run LMM or not
    RERUN = 0; #overwrite existing file
    LIMITSIMS = 0; #flag to limit the total number of simulations in a director
    SIMLIMIT = 500;#skip if there are more than this many output files
    NORESID = 0; #indicator to include covariates in HE regression instead of residualizing.
    ERESID = 0; #indicator to read pre-computed residual expression files.
    READRESID = -1; #indicator to read pre-computed residual genotype files.
    READGRM = 0; #indicator to read pre-computed GRM file.
    OUTPUTSTATS = 0; #indicator to print out statistics for regression
    SM = "";  #use strict mask 012 files

    iter=0; #
    
    RUNHELP=1;
    if(length(args) > 0){
        RUNHELP=0;
        for(i in 1:length(args)){
            foo=strsplit(args[i], "=");
            if(length(foo[[1]]) > 2){
                for(i in 3:length(foo[[1]])){
                    foo[[1]][2] = paste(foo[[1]][2],"=",foo[[1]][i],sep="")
                }
            }
            if(foo[[1]][1] == "GENE"){
                GENE=foo[[1]][2];
            }
            else if(foo[[1]][1] == "INDIR"){
                INDIR=foo[[1]][2];
            }
            else if(foo[[1]][1] == "pref012"){
                pref012=foo[[1]][2];
            }
            else if(foo[[1]][1] == "suff012"){
                suff012=foo[[1]][2];
            }
            else if(foo[[1]][1] == "prefPOS"){
                prefPOS=foo[[1]][2];
            }
            else if(foo[[1]][1] == "suffPOS"){
                suffPOS=foo[[1]][2];
            }
            else if(foo[[1]][1] == "OUTDIR"){
                OUTDIR=foo[[1]][2];
            }
            else if(foo[[1]][1] == "OUTFILESUFF"){
                OUTFILESUFF=foo[[1]][2];
                if(grep("^_", OUTFILESUFF, perl=T, invert=T)){ #add underscore
                    OUTFILESUFF = paste("_",OUTFILESUFF,sep="");
                }
            }
            else if(foo[[1]][1] == "SNPLIST"){
                SNPLIST=foo[[1]][2];
            }
            else if(foo[[1]][1] == "NORESID"){
                NORESID=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "ERESID"){
                ERESID=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "READRESID"){
                READRESID=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "READGRM"){
                READGRM=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "OUTPUTSTATS"){
                OUTPUTSTATS=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "MINSNPS"){
                MINSNPS=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "WIDTH"){
                WIDTH=foo[[1]][2];
            }
            else if(foo[[1]][1] == "QN"){
                QN=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "SM"){
                SM=".SM";
            }
            else if(foo[[1]][1] == "K.RANGE"){
                tmp=strsplit(foo[[1]][2], ",");
                K.RANGE=as.numeric(tmp[[1]][1])
                for(j in 2:length(tmp[[1]])){
                    K.RANGE[j] = as.numeric(tmp[[1]][j]);
                }
            }
            else if(foo[[1]][1] == "K"){
                Ktmp=as.numeric(foo[[1]][2]);
                K.RANGE=Ktmp;
            }
            else if(foo[[1]][1] == "maxMAC"){
                maxMAC=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "minMAC"){
                minMAC=as.numeric(foo[[1]][2]);
                minMAC.RANGE = minMAC;
            }
            else if(foo[[1]][1] == "minMAC.RANGE"){
                tmp=strsplit(foo[[1]][2], ",");
                minMAC.RANGE=as.numeric(tmp[[1]][1])
                for(j in 2:length(tmp[[1]])){
                    minMAC.RANGE[j] = as.numeric(tmp[[1]][j]);
                }
            }
            else if(foo[[1]][1] == "COVARCOLS"){
                COVARCOLS=foo[[1]][2];
                COVARCOLS = strsplit(COVARCOLS,",");
                COVARCOLS = as.numeric(unlist(COVARCOLS));
            }
            else if(foo[[1]][1] == "COVARFILE"){
                COVARFILE=foo[[1]][2];
                if(!file.exists(COVARFILE)){
                    cat("cannot find ",COVARFILE,"...\n",sep="");
                    quit('no');
                }
            }
            else if(foo[[1]][1] == "PHENOFILE"){
                PHENOFILE=foo[[1]][2];
                if(!file.exists(PHENOFILE)){
                    cat("cannot find ",PHENOFILE,"...\n",sep="");
                    quit('no');
                }
            }
            else if(foo[[1]][1] == "ngPC"){
                ngPC=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "npPC"){
                npPC=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "DSAMP"){
                DSAMP=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "PERM"){
                PERM=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "nPERM"){
                nPERM=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "nSIM"){
                nSIM=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "SIM"){
                SIM=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "Nc"){
                Nc=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "fracRare"){
                fracRare=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "freqThresh"){
                freqThresh=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "causalmodel"){
                causalmodel=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "betamodel"){
                betamodel=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "betaALL"){
                betaALL=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "h2"){
                h2=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "Kgen"){
                Kgen=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "Krem"){
                Krem=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "errorModel"){
                errorModel=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "errAdd"){
                errAdd=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "errRem"){
                errRem=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "delta"){
                delta=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "rho"){
                rho=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "tau"){
                tau=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "demog"){
                demog=foo[[1]][2];
                if(demog != "ten" & demog != "gravel"){
                    RUNHELP=1;
                    cat("demog must be either \"ten\" or \"gravel\".\n");
                }
            }
            else if(foo[[1]][1] == "dfe"){
                dfe=foo[[1]][2];
                if(dfe != "Boyko" & dfe != "Torg"){
                    RUNHELP=1;
                    cat("dfe must be either \"Boyko\" or \"Torg\".\n");
                }
            }
            else if(foo[[1]][1] == "runHE"){
                runHE=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "runLMM"){
                runLMM=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "RERUN"){
                RERUN=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "LIMITSIMS"){
                LIMITSIMS=as.numeric(foo[[1]][2]);
            }
            else if(foo[[1]][1] == "SIMLIMIT"){
                LIMITSIMS=1;
                SIMLIMIT=as.numeric(foo[[1]][2]);
            }
            else{
                if(foo[[1]][1] != "help" & foo[[1]][1] != "h"){
                    cat("I misunderstood option ",i,": ",args[i],"\n",sep="")
                }
                RUNHELP=1;
            }
        }
    }
    if(RUNHELP==1){
        cat("usage:  Rscript HE-reg_gene_PCadj_freqBin_width_rare_NZLMM.R <options>\n");
        cat("Options can be included in any order:\n");
        cat("\thelp [h]\tPrints out this help menu. For options below, default values are in parentheses.\n")
        cat("\tGENE=<gene> Gene to pull genotypes from  (",GENE,")\n",sep="");
        cat("\tINDIR=<INDIR> Directory that contains genotype data  (",INDIR,")\n",sep="");
        cat("\tpref012=<PREF> prefix genotype data file  (",pref012,")\n",sep="");
        cat("\tsuff012=<SUFF> suffix genotype data file  (",suff012,"); input file will be <pref012><GENE><suff012> (e.g. /some/dir/APOL1.012.txt.gz\n",sep="");
        cat("\tprefPOS=<PREF> prefix genotype data file  (",prefPOS,")\n",sep="");
        cat("\tsuffPOS=<SUFF> suffix positions data file  (",suffPOS,"); input file will be <prefPOS><GENE><suffPOS> (e.g. /some/dir/APOL1.POS.txt.gz\n",sep="");
        cat("\tOUTDIR=<OUTDIR> Root directory to put the output, files will be put in parameter-specific subdirectories  (",OUTDIR,")\n",sep="");
        cat("\tOUTPUTSTATS=<#> print statistics for regression. 0=no; -1=all bins; K=print only bin K  (",OUTDIR,")\n",sep="");
            cat("\tSNPLIST=<SNPfile> File with SNP positions in columns to put in each bin. Must specify K   (",SNPLIST,")\n",sep="");
        cat("\tNORESID=<0/1> Include covariats in HE regression instead of residualizing phenotypes and genotypes  (",NORESID,")\n",sep="");
        cat("\tERESID=<0/1> Read precomputed expression residuals, 0=no, 1=yes  (",ERESID,")\n",sep="");
        cat("\tREADRESID=<0/1> Read precomputed genotype residuals, 0=no, 1=yes, -1=try  (",READRESID,")\n",sep="");
        cat("\tMINSNPS=<MINSNPS> Minimum number of SNPs per bin for h2 analysis  (",MINSNPS,")\n",sep="");
        cat("\tWIDTH=<WIDTH> Width of the window around GENE in bp or CHR/GENOME  (",WIDTH,")\n",sep="");
        cat("\tQN=<0/1> Run quantile normalization [0/1]  (",QN,")\n",sep="");
        cat("\tK=<K> Number of allele frequency bins to infer  (",Ktmp,")\n",sep="");
        cat("\tmaxMAC=<minMAC> Maximum allele count to be included in h2 analysis; -1 = no max  (",minMAC,")\n",sep="");
        cat("\tminMAC=<minMAC> Minimum allele count to be included in h2 analysis  (",minMAC,")\n",sep="");
        cat("\tCOVARCOLS=<COVARCOLS> Number (or list) of covariates to control for  (",COVARCOLS,")\n",sep="");
        cat("\tCOVARFILE=<file> File with covariates in columns  (",COVARFILE,")\n",sep="");
        cat("\tngPC=<ngPC> Number of genotype PCs to control for  (",ngPC,")\n",sep="");
        cat("\tnpPC=<ngPC> Number of phenotype PCs to control for  (",npPC,")\n",sep="");
        cat("\tDSAMP=<proportion> The fraction of reads to downsample to for evaluating the effect of overall expression levels on inference  (",DSAMP,")\n",sep="");

        cat("\nPermutation parameters:\n\n");
        cat("\tPERM=<PERM> Run <nPERM> permutations of type: 0=none; 1=random; 2=permute within pop; 3=sample expression from random gene; 4=use file for permuted phenotypes; 5=permute expression values like 1 but maintain PC order   (",PERM,")\n",sep="");
        cat("\tnPERM=<nPERM> Run <nPERM> permutations of type ",PERM,"\n",sep="");
        
        cat("\nSimulation parameters:\n\n");
        cat("\tSIM=<SIM> Indicator turning on simulation machinery  (",SIM,")\n",sep="");
        cat("\tnSIM=<nSIM> Number of simulations to run  (",nSIM,")\n",sep="");

        cat("\tNc=<Nc> Number of causal variants   (",Nc,")\n",sep="");

        cat("\tcausalmodel=<causalmodel> Define how causals are distributed: 1-8  (",causalmodel,")\n\t1: <Nc> randomly draw causal variants, semi-favors rare; \n\t2: <Nc> causals sampled in proportion to frequency, favors common; \n\t3: <Nc> causals sampled unifromly from <Kgen> bins, exluding first <Krem> bins of rare variants; \n\t4: <Nc> causals sampled with inverse proportion to freq, linear decrease in sample prob;  \n\t5: <Nc> causals where <fracRare> have frequency <= <freqThresh>;  \n\t6: all variants are causal;  \n\t7: all variants with MAF <= <freqThresh> are causal;  \n\t8: <Nc> causals drawn from Uricchio et al 2016 distribution with given <demog>, <dfe>, <rho>, and <tau> files\n",sep="");
        
        cat("\tbetamodel=<betamodel> Define how effect sizes are distributed: 1-4  (",betamodel,")\n\t1: beta=0.4*log10(MAF);  \nn\t2: beta=<betaALL> for all causals;  \n\t3: beta=inverse variance, 1/(2*MAF*(1-MAF));  \n\t4: beta follows Uricchio et al 2016 using causalmodel=8",sep="");
        
        cat("\tbetaALL=<betaALL> beta for all variants when betamodel=2  (",betaALL,")\n",sep="");

        cat("\tfracRare=<fracRare> fraction of causal variants that are rare  (",fracRare,")\n",sep="");

        cat("\tfreqThresh=<freqThresh> Define frequency of rare variant for fracRare  (",freqThresh,")\n",sep="");
        cat("\trho=<rho> Correlation between fitness and effect size for betamodel=8  (",rho,")\n",sep="");
        cat("\ttau=<tau> Transformation of effect sizes for betamodel=8  (",tau,")\n",sep="");
        cat("\tdemog=<demog> Demographic model, can be either \"ten\" or \"gravel\" for betamodel=8  (",demog,")\n",sep="");
        cat("\tdfe=<dfe> Distribution of fitness effects, can be either \"Boyko\" or \"Torg\" for betamodel=8  (",dfe,")\n",sep="");
        cat("\th2=<h2> Specify desired h2 for simulation  (",h2,")\n",sep="");
        cat("\terrorModel=<errorModel> Define model for adding errors: 0-3  (",errorModel,")\n\t\t0:  nothing;  1: remove rpois(<errRem>*<Nc>) causal SNPs so that they contribute to phenotype but not GRM, only possible when simulating phenotype;  2: Add rpois(<errAdd>*#SNPs) non-causal SNPs to GRM, but not phenotype, with freq drawn from distribution of observed phenotypes, can be run on simulations or real data;  3: Seq error at causal SNPs, add rpois(<errAdd>*#causalAlleles) causal alleles at hom-ref/het sites and remove rpois(<errRem>*#causalAlleles) causal alleles at hom-alt/het sites, can only be run on simulations.\n",sep="");
            cat("\terrRem=<errRem> Fraction of erroneous genotypes removed  (",errRem,")\n",sep="");
        cat("\terrAdd=<errAdd> Fraction of erroneous genotypes added  (",errAdd,")\n",sep="");
        cat("\trunHE=<runHE> Indicator to run HE regression  (",runHE,")\n",sep="");
        cat("\trunLMM=<runLMM> Indicator to run LMM: OFF BY DEFAULT! (",runLMM,")\n",sep="");
        cat("\tRERUN=<RERUN> Indicator to rerun parameters and overwrite file  (",RERUN,")\n",sep="");
        quit("no");
    }

    ##this block of code will just see if the output files already exist and exit right away if nothing needs to be done
    if(RERUN == 0){
        DONE = 0;
        TODO = 0;
        
        for(minMAC in minMAC.RANGE){
            if(DONE < TODO){
                break; #already behind, run the analysis
            }
            for(Ktmp in K.RANGE){
                if(DONE < TODO){
                    break; #already behind, run the analysis
                }
                PRINTDIRCK = paste(OUTDIR,"/",WIDTH,sep="");
                if(!dir.exists(PRINTDIRCK)){
                    TODO = 1;
                    break; 
                }
                COVARFLAGS = paste("_ngPC=",ngPC,"_npPC=",npPC,sep="")
                
                if(COVARCOLS[1] >= 0){
                    foo = paste(COVARCOLS,collapse="-");
                    COVARFLAGS = paste(COVARFLAGS,"_CV=",foo,sep="")
                }
                
                if(SIM==0 & PERM==0){
                    if(maxMAC < 0){
                        PRINTDIRCK = paste(PRINTDIRCK,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,sep="");
                    }
                    else{
                        PRINTDIRCK = paste(PRINTDIRCK,"/K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,sep="");
                    }
                }
                else if(PERM>0){
                    PRINTDIRCK = paste(PRINTDIRCK,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,sep="");
                }
                else{
                    if(causalmodel == 8){
                        if(maxMAC < 0){
                            PRINTDIRCK = paste(PRINTDIRCK,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,sep="");
                        }
                        else{
                            PRINTDIRCK = paste(PRINTDIRCK,"/K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,sep="");
                        }
                    }
                    else{
                        if(maxMAC < 0){
                            PRINTDIRCK = paste(PRINTDIRCK,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,sep="");
                        }
                        else{
                            PRINTDIRCK = paste(PRINTDIRCK,"/K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,sep="");
                        }
                    }
                }
                if(!dir.exists(PRINTDIRCK)){
                    TODO = TODO+1;
                    break; #need to run the analysis
                }
                
                if(PERM==0 & SIM==0 & errorModel==0){
                    if(maxMAC < 0){
                        HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_reml.h2",sep="")
                    }
                    else{
                        HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_reml.h2",sep="")
                    }
                }
                else if(PERM>0){
                    if(nPERM==0){
                        cat("must specify number of permutations to run using nPERM=<nPERM>\n");
                        quit("no");
                    }
                    if(maxMAC < 0){
                        HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_reml.h2",sep="")
                    }
                    else{
                        HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_reml.h2",sep="")
                    }
                }
                else if(SIM>0){
                    if(nSIM==0){
                        cat("must specify number of simulations to run using nSIM=<nSIM>\n");
                        quit("no");
                    }
                    if(causalmodel == 8){
                        if(maxMAC < 0){
                            HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_reml.h2",sep="");
                        }
                        else{
                            HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_reml.h2",sep="");
                        }
                    }
                    else{
                        if(maxMAC < 0){
                            HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_reml.h2",sep="");
                        }
                        else{
                            HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_reml.h2",sep="");
                        }
                    }
                }
                else if(errorModel > 0){
                    if(maxMAC < 0){
                        HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_he.h2",sep="");
                        REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_reml.h2",sep="");
                    }
                    else{
                        HEOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_he.h2",sep="");
                        REMLOUTFILECK=paste(PRINTDIRCK,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_reml.h2",sep="");
                    }
                }
                if(PERM == 0){
                    HEBKCK = HEOUTFILECK;
                    HEBKCK1 = HEOUTFILECK;
                    REMLBKCK = REMLOUTFILECK;
                }
                else{
                    HEBKCK = HEOUTFILECK;
                    HEBKCK1 = HEOUTFILECK;
                    REMLBKCK = REMLOUTFILECK;
                }
                runHEtmp = runHE;
                runLMMtmp = runLMM;
                
                if(iter < nSIM){
                    TODO = 1;
                    break;
                }
                if(file.exists(paste(HEOUTFILECK,".gz",sep="")) |
                   file.exists(paste(HEBKCK,".gz",sep="")) |
                   file.exists(paste(HEBKCK,".gz",sep=""))){
                    ##ok, already done
                }
                else{
                    TODO = 1;
                    cat("1: ",(paste(HEOUTFILECK,".gz",sep="")),"=",file.exists(paste(HEOUTFILECK,".gz",sep="")),"\n",(paste(HEBKCK,".gz",sep="")),"=",file.exists(paste(HEBKCK,".gz",sep="")),"\n",(paste(HEBKCK,".gz",sep="")),"=",file.exists(paste(HEBKCK,".gz",sep="")),"\n");
                    break; #need to run analysis
                }
                if(runLMM==1){
                    if(file.exists(paste(REMLOUTFILECK,".gz",sep="")) |
                       file.exists(paste(REMLBKCK,".gz",sep=""))){
                        ##ok, already done
                    }
                    else{
                        TODO = 1;
                        cat("2: ",paste(REMLOUTFILECK,".gz",sep=""),"=",
                            file.exists(paste(REMLOUTFILECK,".gz",sep="")),"\n",
                            paste(REMLBKCK,".gz",sep=""),"=",file.exists(paste(REMLBKCK,".gz",sep="")),"\n");
                        break;
                    }
                }
            }
        }
        if(TODO == DONE){
            cat("All analyses done, exiting\n");
            quit("no");
        }
        else{
            cat("found a file not done: DONE=",DONE,"; TODO=",TODO,"\n");
        }
    }


    if(runLMM == 1){
        if(file.exists("reml_func.R")){
            source("reml_func.R");
        }
        else{
            cat("reml_func.R not in current working directory. Cannot run LMM\n");
            quit("no");
        }
    }

    if((errorModel > 0 | SIM > 0) & READRESID != 0){
        cat("Error, can only use precomputed residuals for data analysis, not simulations\n")
        quit("no");
    }

    if(SIM > 0 & causalmodel == 8){ #read data
        if(rho<0 | tau<0){
            cat("Must include rho and tau if causalmodel==8\n");
            quit("no");
        }
        betaMat = matrix(scan(paste(INDIR,"/OutBeta/",demog,dfe,"Sel_rho=",rho,"_tau=",tau,".txt",sep=""),quiet=T),byrow=T,nc=3)
        betaMat[,1] = apply(cbind(betaMat[,1],720-betaMat[,1]),1,min)
        betaMat = betaMat[order(betaMat[,1]),]
    }
    
    if(WIDTH == "CHR"){
        GENES=matrix(scan(paste(INDIR,"/matrix.eqtl/expressions.genes.meta.matrix.eqtl.txt",sep=""),sep="\t",what=""),nc=4,byrow=T)
        CHR = GENES[which(GENES[,1] == GENE)[1], 2];
        cat("found gene ",GENE," on ",CHR,"\n");
    }
    
    if(ngPC>0 | npPC>0){
        PCnames=matrix(scan(paste(INDIR,"/matrix.eqtl/covariates.eur.matrix.eqtl.txt",sep=""),nlines=1,what="character"),byrow=T,nr=1)
    
        phenonames=matrix(scan(paste(INDIR,"/GENOW/indiv.txt",sep=""),what="character"),byrow=T,nc=1);
    }
    
    if(PERM==2){
        POPS=read.table(paste(INDIR,"/GENOW/20140610_samples_to_analyse_v10.txt",sep=""))
        INDPOP = as.matrix(POPS[match(phenonames, POPS[,1]),2])            
    }
    else if(PERM==4){ #read permuted ID file
        PERMID = read.table(paste(INDIR,"/permPheno_withinPop.txt",sep=""));
    }
    
    if(ngPC>0 | npPC>0){
        pcs = read.table(paste(INDIR,"/matrix.eqtl/covariates.eur.matrix.eqtl.txt",sep=""),skip=1)
        pcs = pcs[,2:ncol(pcs)]
        pcs = pcs[,which(PCnames %in% phenonames)]
        pcs = as.matrix(pcs)
    }

    if(SNPLIST != ""){
        if(length(K.RANGE) != 1){
            cat("Must specify K when using SNPLIST\n");
            RUNHELP=1;
        }
        else{ #read SNP list
            cat("SNPLIST=<",SNPLIST,">\n",sep="");
            if(grep(".gz$",SNPLIST,perl=T)){
                GZIN = gzfile(SNPLIST,"r");
                maxCol=0;
                for(i in 1:Ktmp){
                    foo = as.numeric(scan(GZIN,nlines=1,quiet=T))
                    if(length(foo) > maxCol){
                        maxCol = length(foo);
                    }
                }
                close(GZIN);
                GZIN = gzfile(SNPLIST,"r");
                SNPSET = matrix(,nr=Ktmp,nc=maxCol);
                for(i in 1:Ktmp){
                    foo = as.numeric(scan(GZIN,nlines=1))
                    SNPSET[i,1:length(foo)] = foo;
                }
                close(GZIN);
            }
            else{
                GZIN = file(SNPLIST,"r");
                maxCol=0;
                for(i in 1:Ktmp){
                    foo = as.numeric(scan(GZIN,nlines=1))
                    if(length(foo) > maxCol){
                        maxCol = length(foo);
                    }
                }
                close(GZIN);
                GZIN = file(SNPLIST,"r");
                SNPSET = matrix(,nr=Ktmp,nc=maxCol);
                for(i in 1:Ktmp){
                    foo = as.numeric(scan(GZIN,nlines=1))
                    SNPSET[i,1:length(foo)] = foo;
                }
                close(GZIN);
            }
        }
        cat("read ",nrow(SNPSET)," rows of SNPs in SNPLIST\n");
    }
    
    ##get genotypes
    if(READGRM==0){
        if(READRESID != 1){
            if(WIDTH == "CHR"){
                if(READGRM==0){
                    if(pref012!="" | suff012!=""){
                        infile=paste(pref012,CHR,suff012,sep="");
                    }
                    else{
                        infile=paste(INDIR,"/GENOW/",WIDTH,"/",CHR,SM,".txt.012.gz",sep="")
                    }
                    GZIN=gzfile(infile,"r");
                    d=scan(GZIN,nlines=1)
                    close(GZIN);
                    GZIN=gzfile(infile,"r");
                    genosBK = matrix(scan(GZIN),byrow=T,nc=length(d))
                    close(GZIN)
                    cat("dim(genosBK) =",dim(genosBK),"\n")
                    notmono = array(,dim=0)
                    for(i in 1:ncol(genosBK)){
                        if(length(table(genosBK[,i])) > 1){
                            notmono=c(notmono,i)
                        }
                    }
                    genosBK = genosBK[,notmono]
                    cat("dim(genosBK) ->",dim(genosBK),"\n")
                }
            }
            else{
                if(pref012!="" | suff012!=""){
                    infile=paste(pref012,GENE,suff012,sep="");
                }
                else{
                    infile=paste(INDIR,"/GENOW/",WIDTH,"/",GENE,".",WIDTH,SM,".txt.012.gz",sep="")
                }
                cat("reading",infile,"\n");
                GZIN=gzfile(infile,"r");
                d=scan(GZIN,nlines=1)
                close(GZIN);
                GZIN=gzfile(infile,"r");
                genosBK = matrix(scan(GZIN),byrow=T,nc=length(d))
                close(GZIN)
                cat("dim(genosBK) =",dim(genosBK),"\n")
                notmono = array(,dim=0)
                for(i in 1:ncol(genosBK)){
                    if(length(table(genosBK[,i])) > 1){
                        notmono=c(notmono,i)
                    }
                }
                genosBK = genosBK[,notmono]
                cat("dim(genosBK) ->",dim(genosBK),"\n")
            }
            if(SNPLIST != ""){ #subsampling SNPs based on file            
                if(prefPOS!="" | suffPOS!=""){
                    infile=paste(prefPOS,GENE,suffPOS,sep="");
                }
                else{
                    infile=paste(INDIR,"/GENOW/",WIDTH,"/",GENE,".",WIDTH,".txt.pos.gz",sep="");
                }
                GZIN=gzfile(infile,"r");
                SNPPOSBK=scan(GZIN);
                close(GZIN);
                SNPPOSBK = SNPPOSBK[notmono];
            }

            efreqsBK  = apply(genosBK,2,function(x){mean(x,na.rm=T)})/2
            emafsBK =  apply(cbind(efreqsBK,1-efreqsBK),1,min)
            sefreqsBK = sort(emafsBK,index.return=T)
            MBK = ncol(genosBK)
            NBK = nrow(genosBK)
            
            if(MBK < MINSNPS){
                cat("not enough SNPs read for this gene, read",MBK,", want at least ",MINSNPS,"; exiting\n");
                quit("no");
            }
            
            ## double check to see if the genotype residual files are present
            if(READRESID != 0){
                READRESID = 0; 
                if(WIDTH=="CHR"){
                    RESIDFILE=paste(INDIR,"/GENOW/",WIDTH,"/",CHR,"_ngPC=",ngPC,"_npPC=",npPC,".txt.resid.gz",sep="")
                }
                else{
                    RESIDFILE=paste(INDIR,"/GENOW/",WIDTH,"/",GENE,"_",WIDTH,"_ngPC=",ngPC,"_npPC=",npPC,".txt.resid.gz",sep="")
                }
                ##check if file exists
                if(file.exists(RESIDFILE)){
                    ##if it exists, read it to make sure it is ok
                    GZIN=gzfile(RESIDFILE,"r");
                    foo = try(matrix(scan(GZIN),byrow=T,nc=MBK));
                    close(GZIN);
                    if(class(foo) == "try-error"){ #don't do anything, just procede as normal
                    }
                    else{ # genotype residual files exist, double check they are the right size
                        cat(RESIDFILE," exists!\n");
                        if(nrow(foo) == NBK+1 & ncol(foo) == MBK){
                            cat("genotype residual file ok, let's use it!!\n");
                            efreqsBK=foo[1,];
                            genosBK = foo[2:nrow(foo),];
                            MBK = ncol(genosBK)
                            NBK = nrow(genosBK)
                            efreqsBK = efreqsBK/NBK/2
                            emafsBK =  apply(cbind(efreqsBK,1-efreqsBK),1,min)
                            sefreqsBK = sort(emafsBK,index.return=T)
                            READRESID = 1;
                        }
                    }
                }
            }
        }
        else{
            if(READGRM!=0){
                cat("cannot READRESID and READGRM...\n");
                quit("no");
            }
            if(WIDTH=="CHR"){
                RESIDFILE=paste(INDIR,"/GENOW/",WIDTH,"/",CHR,"_ngPC=",ngPC,"_npPC=",npPC,".txt.resid.gz",sep="")
            }
            else{
                RESIDFILE=paste(INDIR,"/GENOW/",WIDTH,"/",GENE,"_",WIDTH,"_ngPC=",ngPC,"_npPC=",npPC,".txt.resid.gz",sep="")
            }
            if(!file.exists(RESIDFILE)){
                cat("error, residfile DNE, expecting ",RESIDFILE,"\n");
                quit("no");
            }
            GZIN=gzfile(RESIDFILE,"r");
            efreqsBK=scan(file=GZIN,nlines=1)
            close(GZIN);
            GZIN=gzfile(RESIDFILE,"r");
            genosBK = matrix(scan(file=GZIN, skip=1),byrow=T,nc=length(efreqsBK))
            close(GZIN);
        
            MBK = ncol(genosBK)
            NBK = nrow(genosBK)
            efreqsBK = efreqsBK/NBK/2
            emafsBK =  apply(cbind(efreqsBK,1-efreqsBK),1,min)
            sefreqsBK = sort(emafsBK,index.return=T)
        }
    }
    
    MISSINGDATA = 0;
    if(is.na(min(genosBK))){
        MISSINGDATA = 1;
    }
    iter=1;
    HEpassed=0;
    REMLpassed=0;
    while(iter <= nSIM){ #loop
        cat("iter=",iter,"/",nSIM,"\n",sep="")

        if(READGRM==0){
            genos = genosBK
            efreqs = efreqsBK
            emafs = emafsBK
            sefreqs = sort(emafs,index.return=T)
        }
        
        iter = iter+1
        if(SIM == 0){ ##real data analysis
            ##get phenotypes
            if(ERESID == 1){ #read expression residuals
                phenos = scan(paste(INDIR,"/matrix.eqtl/log2/",GENE,".eresid.pheno",sep=""))
            }
            else{
                if(PHENOFILE == ""){
                    phenos = scan(paste(INDIR,"/matrix.eqtl/log2/",GENE,".log2.pheno",sep=""))
                }
                else{
                    if(grep(".gz$",PHENOFILE,perl=T)){
                        GZIN = gzfile(PHENOFILE,"r");
                        phenos = scan(GZIN)
                        close(GZIN);
                    }
                    else{
                        phenos = scan(PHENOFILE)
                    }
                }
            }
            if(PERM == 3){ #choose a random gene
                GENES=matrix(scan(paste(INDIR,"/matrix.eqtl/expressions.genes.meta.matrix.eqtl.txt",sep=""),sep="\t",what=""),nc=4,byrow=T)
                
                ##sample a gene, make sure phenos exist, and uncorrelated (trans) to original
                FOUND=0;
                tries = 0;
                while(FOUND == 0){
                    tries = tries+1;
                    G = sample(GENES[,1],1);
                    if(ERESID == 1){
                        PHfile = paste(INDIR,"/matrix.eqtl/log2/",G,".eresid.pheno",sep="")
                    }
                    else{
                        PHfile = paste(INDIR,"/matrix.eqtl/log2/",G,".log2.pheno",sep="")
                    }
                    if(file.exists(PHfile)){
                        tmpphenos = scan(PHfile,quiet=T)
                        cor=cor.test(phenos,tmpphenos,method="kendall"); #rank correlation
                        if(cor$p.value > 0.05){
                            cat("sampled expression GENE=",G,": cor.test=",cor$estimate,"; p=",cor$p.value,"; tries=",tries,"\n",sep="");
                            phenos = tmpphenos;
                            FOUND = 1;
                        }
                    }
                }
            }
            else if(PERM == 5){ #shuffle phenotypes, but keep everything else the same
                phenos = sample(phenos);
            }
            
            ##Down sample reads if desired
            if(DSAMP < 1.0){
                phenos = log2(rgamma(length(phenos), DSAMP*(2^phenos), 1));
            }
        }
        else{  #simulate phenotypes!
            if(READGRM != 0){
                cat("cannot run simulations with READGRM=1\n");
                quit("no");
            }
            ##select causal variants
            if (causalmodel == 1) { #sample Nc from all variants
                causals = sort(sample(1:MBK,Nc))
            }
            else if (causalmodel == 2) { #sample Nc with prob given by MAF
                psamp = emafsBK/sum(emafsBK)
                causals = sample(1:MBK,Nc,prob=psamp)
            }
            else if (causalmodel == 3) {
                st = 1+(Krem)*floor(MBK/Kgen) #divide SNPs into Kgen bins, exclude first Krem. This is the index of the first kept SNP
                bsnps = sefreqsBK$ix[st:MBK] #get the index from the freq-sorted array
                causals = bsnps[sample(1:length(bsnps),Nc)]
            }
            else if (causalmodel == 4) {
                psamp = (0.5-emafsBK)/sum(0.5-emafsBK)
                causals = sample(1:MBK,Nc,prob=psamp)
            }
            else if(causalmodel == 5) {
                set1 = which(efreqsBK > 0 & efreqsBK < 1 & efreqsBK <= freqThresh) #set of "rare" SNPs
                set2 = which(efreqsBK > 0 & efreqsBK < 1 & efreqsBK > freqThresh) #set of "common" SNPs
                if(length(set1) <= round(fracRare*Nc) |
                   length(set2) <= Nc-round(fracRare*Nc)){
                    cat("not enough causals...\n");
                    SKIP = 1;
                    break;
                }
                ##cat("fracRare=",fracRare,"; Nc=",Nc,"; freqThresh=",freqThresh,"\n");
                ##cat("l(set1)=",length(set1)," E=",round(fracRare*Nc),"; l(set2)=",length(set2)," E=",Nc-round(fracRare*Nc),"\n");
                causals = c(sample(set1, round(fracRare*Nc)));
                ##cat("causals1:",efreqsBK[causals],sep="\n");
                causals = c(causals, sample(set2, Nc-round(fracRare*Nc)));
                ##cat("causals2:",efreqsBK[causals],sep="\n");
                ##quit("no")
            }
            else if(causalmodel == 6){
                causals = 1:MBK; #all SNPs are causal
            }
            else if(causalmodel == 7){
                causals = which(emafsBK <= freqThresh); #all "rare" SNPs are causal
            }
            else if(causalmodel == 8){ #follows Uricchio et al 2016 distribution
                ##idea:  assign selection coefficient based on sampling allele at same frequency from Uricchio model, then sample effect size based on that fitness
                causals = array(,dim=0);
                if(Nc >= MBK){ #all SNPs are causal
                    betas=0;
                    for(bi in 1:MBK){
                        causals[bi] = bi;
                        betas[bi] = betaMat[sample(which(betaMat[,1]==round(emafsBK[bi]*2*NBK)), 1),3]; #sample a random SNP with same frequency, get simulated fitness
                    }
                }
                else{ #sample Nc causals according to causal frequency distribution from Uricchio model.
                    causals = array(,dim=0);
                    betas = array(,dim=0);
                    ccol = sort(sample(1:nrow(betaMat),Nc)); #sample Nc causal freqs
                    b = matrix(betaMat[ccol,c(1,3)], ncol=2); #get (freq, effect size) pairs
                    b[,1] = apply(cbind(b[,1],2*NBK-b[,1]),1, min); #convert to MAC
                    ufreq = unique(b[,1]) #unique causal frequencies
                    for(bi in 1:length(ufreq)){ #now sample SNPs at each causal frequency
                        sampset=which(round(2*NBK*emafsBK)==ufreq[bi]);
                        foocnt=0;
                        ## there is a chance no SNPs at right frequency. get something close
                        while(length(sampset) <= sum(b[,1]==ufreq[bi])){
                            foocnt = foocnt+1;
                            sampset=which(round(2*NBK*emafsBK)>=ufreq[bi]-foocnt &
                                              round(2*NBK*emafsBK)<=ufreq[bi]+foocnt);
                        }
                        causals = c(causals, sample(sampset,sum(b[,1]==ufreq[bi])));
                        betas = c(betas, sample(b[which(b[,1]==ufreq[bi]),2],sum(b[,1]==ufreq[bi])));
                    }
                }
            }
            else{
                cat("Must enter causalmodel = 1,2,...,8.  Read ",causalmodel,"\n");
                quit("no");
            }
            
            ##choose betas for the causal SNPs
            if(causalmodel != 8){
                if (betamodel==1) {
                    betas =  sign(delta-runif(length(causals)))*(-0.4*log(emafs[causals]))
                }
                else if(betamodel==2){
                    betas = sign(delta-runif(length(causals)))*rep(betaALL, Nc);
                }
                else if(betamodel==3){
                    betas = sign(delta-runif(length(causals)))*1/(2*emafs[causals]*(1-emafs[causals]))
                }
                else if (betamodel==4){
                    ##specified above in causalmodel == 8
                }
                else{
                    cat("Must enter betamodel = 1,2,3, or 4.  Read ",betamodel,"\n");
                    quit("no");
                }
            }
            
            ##generate genetic component of phenotype
            ##handling R stupidity (or mine) of dealing with 1 versus multiple causals
            if (Nc == 1) {
                phenosg = genosBK[,causals]*betas	
            }
            else{
                phenosg = as.matrix(genosBK[,causals])%*%betas	
            }
            
            ##generate phenotypes with specified heritability
            cat("var(phenosg)=",var(phenosg),"; NBK=",NBK,"; h2=",h2,"\n",sep="");
            phenos = phenosg + rnorm(NBK,0,sqrt(var(phenosg)*(1/h2-1)))
        }   

        if(errorModel!=0 & READGRM!=0){
            cat("cannot run error model simulations with READGRM=1\n");
            quit("no");
        }
        if(errorModel == 1){
            if(errRem > 0 | SIM!=0){ # remove some causal sites from genotype data
                remCaus = sample(causals, rpois(1,errRem*length(causals)));
                good = which(!(1:MBK %in% remCaus))
                genos = genosBK[,good];
                M = ncol(genos)
                N = nrow(genos)
                efreqs  = apply(genos,2,function(x){mean(x,na.rm=T)})/2
                emafs =  apply(cbind(efreqs,1-efreqs),1,min)
                sefreqs = sort(emafs,index.return=T)
            }
            else{
                cat("errorModel=1, but SIM=",SIM," and errRem=",errRem,"??\n",sep="");
                quit("no");
            }
        }
        else if(errorModel == 2){
            if(errAdd > 0){ # add some non-causal sites to genotype data
                added = rpois(1,M*errAdd);
                cat("adding ",added," SNPs to GRM\n",sep="");
                addFreqs = sample(efreqs, added, replace=T);
                for(i in 1:added){
                    g = array(0,dim=N);
                    j=1;
                    while(j <= round(addFreqs[i]*N*2)){
                        indiv=sample(1:N,1)
                        if(g[indiv] < 2){
                            g[indiv] = g[indiv]+1;
                            j=j+1;
                        }
                    }
                    genos = cbind(genos,g);
                }
                M = ncol(genos)
                N = nrow(genos)
                efreqs  = apply(genos,2,function(x){mean(x,na.rm=T)})/2
                emafs =  apply(cbind(efreqs,1-efreqs),1,min)
                sefreqs = sort(emafs,index.return=T)
                cat("there are now ",M," SNPs\n",sep="");
            }
        }
        else if(errorModel == 3){
            if(SIM==1 & (errAdd > 0 | errRem > 0)){ # add noise to genotypes at causal sites
                ##idea:  find total causal alleles, and increase this by <errAdd>% and decrease original number of causal alleles by <errRem>%.  If <errRem>=<errAdd> should result in no change in causal alleles, but will be shuffled.
                
                causHet  = which(genos[,causals] == 1);
                causHom  = which(genos[,causals] == 2);
                nonCaus = which(genos[,causals] == 0);
                numCaus = length(causHet) + 2*length(causHom);
                remCaus = round(numCaus*errRem);
                addCaus = round(numCaus*errAdd);
                
                while(remCaus > 0){
                    indiv = sample(c(causHet, causHom, causHom), 1);
                    if(genos[,causals][indiv] > 0){
                        genos[,causals][indiv] = genos[,causals][indiv]-1;
                        remCaus = remCaus-1;
                    }
                    
                }
                while(addCaus > 0){
                    indiv = sample(c(causHet, nonCaus, nonCaus), 1);
                    if(genos[,causals][indiv] < 2){
                        genos[,causals][indiv] = genos[,causals][indiv]+1;
                        addCaus = addCaus-1;
                    }
                }
            }
        }

        if(READGRM==0){
            if(READRESID != 1){
                genosTMP=genos;
                M = ncol(genosTMP)
                N = nrow(genosTMP)
                efreqsTMP  = apply(genosTMP,2,function(x){mean(x,na.rm=T)})/2
                emafsTMP =  apply(cbind(efreqsTMP,1-efreqsTMP),1,min)
                sefreqsTMP = sort(emafsTMP,index.return=T)
            }
            else{
                genosTMP=genosBK;
                M = ncol(genosTMP)
                N = nrow(genosTMP)
                efreqsTMP  = efreqsBK
                emafsTMP =  emafsBK
                sefreqsTMP = sefreqsBK
            }
        }

        if(NORESID == 0){ #residualize genotypes and phenotypes
            ##adjust phenotypes for PCs
            if(ngPC > 0){
                gpcs=read.table(paste(INDIR,"/matrix.eqtl/GEUVADIS_EUR.mac_2.pca.loadings.txt",sep=""),skip=1);
                ##take subset to desired #PCs
                gpcs=as.matrix(gpcs[match(gpcs[,1], phenonames), 2:(ngPC+1)])
                if(SIM == 0){ 
                    phenos = summary(lm(phenos~gpcs))$resid
                }
            }            
            if(npPC > 0 & SIM == 0){
                if(npPC == 1){
                    phenos = summary(lm(phenos~pcs[1,]))$residuals
                }
                else{
                    phenos = summary(lm(phenos~t(pcs[1:npPC,])))$residuals
                }
            }
        }
        
        ##normalize phenotypes to mean 0 and variance 1
        nphenos = scale(phenos)
        if(QN == 1){
            sphenos = sort(nphenos,index.return=T)
            phenost = qnorm((1:length(nphenos))/(length(nphenos)+1))
            nphenos[sphenos$ix] = phenost
        }

        if(PERM == 4){ #use precomputed sample ID
            tmpphenos = nphenos;
            nphenos = tmpphenos[PERMID[,1]]
        }

        
        ##product of corresponding pairs of phenotypes
        pp = nphenos%*%t(nphenos)
        pp = pp[upper.tri(pp)]

        ##filter out SNPs with MAF < minMAC
        for(minMAC in minMAC.RANGE){
            REMLFAILED = 0;
            if(READGRM==0){
                if(SNPLIST != "" & (minMAC>1 | maxMAC>0)){
                    cat("can only use minMAC=1, maxMAC=-1 (defaults) when defining SNPLIST\n");
                    quit("no");
                }
                cat("minMAC=",minMAC,"; maxMAC=",maxMAC,"\n",sep="");
                if(maxMAC < 0){
                    good = which(emafsTMP>((minMAC-0.5)/(2*NBK)))
                }
                else{
                    good = which(emafsTMP>((minMAC-0.5)/(2*NBK)) & emafsTMP<((maxMAC+0.5)/(2*NBK)))
                }
                genos = genosTMP[,good];
                if(SNPLIST != ""){ #subsampling SNPs based on file
                    SNPPOS = SNPPOSBK[good];
                    if(length(SNPPOSBK) != length(emafsTMP)){
                        cat("oops... length(SNPPOSBK)=",length(SNPPOSBK)," != length(emafsTMP)=",length(emafsTMP),sep="");
                        quit("no");
                    }
                }
                M = ncol(genos)
                N = nrow(genos)
                efreqs = efreqsTMP[good]
                emafs = emafsTMP[good]
                sefreqs = sort(emafs,index.return=T)
                
                if(M < MINSNPS){
                    cat("not enough SNPs read for this gene, read",MBK,", want at least",MINSNPS," (",M < MINSNPS,"); exiting\n",sep="");
                    next;
                }

                if(READRESID != 1){
                    if(NORESID == 0){
                        ##adjust genotypes for PCs
                        if(ngPC > 0){
                            gpcs=read.table(paste(INDIR,"/matrix.eqtl/GEUVADIS_EUR.mac_2.pca.loadings.txt",sep=""),skip=1);
                            ##take subset to desired #PCs
                            gpcs=as.matrix(gpcs[match(gpcs[,1], phenonames), 2:(ngPC+1)])
                            
                            for (i in 1:M){
                                genos[,i] = summary(lm(genos[,i]~gpcs[,1]))$resid 
                            }
                        }
                        if(npPC > 0 & SIM == 0){
                            for (i in 1:M){
                                if(npPC == 1){
                                    genos[,i] = summary(lm(genos[,i]~pcs[1,]))$resid 
                                }
                                else{
                                    genos[,i] = summary(lm(genos[,i]~t(pcs[1:npPC,])))$resid 
                                }
                            }
                        }
                    }
                    ##normalize genotypes to mean 0 and variance 1
                    zgenos = scale(genos)
                }
                else{
                    zgenos = genos;
                }
            }
            for(Ktmp in K.RANGE){
                PRINTDIR = paste(OUTDIR,"/",WIDTH,sep="");
                ifelse(!dir.exists(PRINTDIR), dir.create(PRINTDIR), FALSE)
                COVARFLAGS = paste("_ngPC=",ngPC,"_npPC=",npPC,sep="")

                if(COVARCOLS[1] >= 0){
                    foo = paste(COVARCOLS,collapse="-");
                    COVARFLAGS = paste(COVARFLAGS,"_CV=",foo,sep="")
                }
                
                if(SIM==0 & PERM==0){
                    if(maxMAC < 0){
                        PRINTDIR = paste(PRINTDIR,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,sep="");
                    }
                    else{
                        PRINTDIR = paste(PRINTDIR,"/K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,sep="");
                    }
                }
                else if(PERM>0){
                    PRINTDIR = paste(PRINTDIR,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,sep="");
                }
                else{
                    if(causalmodel == 8){
                        if(maxMAC < 0){
                            PRINTDIR = paste(PRINTDIR,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,sep="");
                        }
                        else{
                            PRINTDIR = paste(PRINTDIR,"/K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,sep="");
                        }
                    }
                    else{
                        if(maxMAC < 0){
                            PRINTDIR = paste(PRINTDIR,"/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,sep="");
                        }
                        else{
                            PRINTDIR = paste(PRINTDIR,"/K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,sep="");
                        }
                    }
                }
                ifelse(!dir.exists(PRINTDIR), dir.create(PRINTDIR), FALSE)
                cat("there are ",length(list.files(path=PRINTDIR, pattern="he.h2.gz")),
                    " simulations already done.\n",sep="");
                if(LIMITSIMS==1 &
                   length(list.files(path=PRINTDIR, pattern="he.h2.gz")) >= SIMLIMIT){
                    cat("Skipping\n");
                    next; #skip to next parameter
                }
                
                if(PERM==0 & SIM==0 & errorModel==0){
                    nSIM=1; #all analyses in a loop, this makes the loop run once.
                    if(maxMAC < 0){
                        HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_reml.h2",sep="")
                    }
                    else{
                        HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,OUTFILESUFF,"_reml.h2",sep="")
                    }
                }
                else if(PERM>0){
                    if(nPERM==0){
                        cat("must specify number of permutations to run using nPERM=<nPERM>\n");
                        quit("no");
                    }
                    if(maxMAC < 0){
                        HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_reml.h2",sep="")
                    }
                    else{
                        HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_he.h2",sep="")
                        REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_PERM=",PERM,OUTFILESUFF,"_reml.h2",sep="")
                    }
                }
                else if(SIM>0){
                    if(nSIM==0){
                        cat("must specify number of simulations to run using nSIM=<nSIM>\n");
                        quit("no");
                    }
                    if(causalmodel == 8){
                        if(maxMAC < 0){
                            HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_reml.h2",sep="");
                        }
                        else{
                            HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_delta=",delta,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,"_",demog,dfe,OUTFILESUFF,"_reml.h2",sep="");
                        }
                    }
                    else{
                        if(maxMAC < 0){
                            HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_reml.h2",sep="");
                        }
                        else{
                            HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_he.h2",sep="");
                            REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_CM=",causalmodel,"_BM=",betamodel,"_h2=",h2,"_Nc=",Nc,"_FR=",fracRare,"_FT=",freqThresh,"_Kgen=",Kgen,"_Krem=",Krem,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,"_rho=",rho,"_tau=",tau,OUTFILESUFF,"_reml.h2",sep="");
                        }
                    }
                }
                else if(errorModel > 0){
                    if(maxMAC < 0){
                        HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_he.h2",sep="");
                        REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_reml.h2",sep="");
                    }
                    else{
                        HEOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_he.h2",sep="");
                        REMLOUTFILE=paste(PRINTDIR,"/",GENE,"_",WIDTH,"_PCadj_freqBin_K=",Ktmp,"_MAC=",minMAC,"-",maxMAC,COVARFLAGS,"_QN=",QN,"_DSAMP=",DSAMP,"_SIM=",SIM,"_errMod=",errorModel,"_errAdd=",errAdd,"_errRem=",errRem,OUTFILESUFF,"_reml.h2",sep="");
                    }
                }
                if(PERM == 0){
                    HEBK = HEOUTFILE;
                    HEBK1 = HEOUTFILE;
                    REMLBK = REMLOUTFILE;
                }
                else{
                    HEBK = HEOUTFILE;
                    HEBK1 = HEOUTFILE;
                    REMLBK = REMLOUTFILE;
                }
                runHEtmp = runHE;
                runLMMtmp = runLMM;
                
                if(RERUN == 0 & (file.exists(paste(HEOUTFILE,".gz",sep="")) |
                                                  file.exists(paste(HEBK,".gz",sep="")) |
                                                      file.exists(paste(HEBK1,".gz",sep="")))){
                    cat(HEOUTFILE,".gz exists\n",sep="");
                    runHEtmp = 0;  #no need to rerun this
                    runMLLtmp = 0; #no need to rerun this either...
                }
                if(RERUN == 0 & (file.exists(paste(REMLOUTFILE,".gz",sep="")) |
                                     file.exists(paste(REMLBK,".gz",sep="")))){
                    cat(REMLOUTFILE,".gz exists\n",sep="");
                    runLMMtmp = 0; #no need to rerun this
                }
                if(runHEtmp == 0 && (runLMMtmp == 0 | REMLFAILED == 1)){
                    next;
                }

                cat("K=",Ktmp,"\n",sep="");
                ##check to see if GRM already exists...
                if(WIDTH=="CHR"){
                    GRMFILE=paste(INDIR,"GENOW/CHR/",CHR,"_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,".txt.grm.gz",sep="");
                    if(file.exists(GRMFILE)){
                        READGRM=1;
                        cat("Found GRM: ",GRMFILE,"  Let's use it.\n");
                    }
                }
                else if(WIDTH=="GENOME"){
                    GRMFILE=paste(INDIR,"GENOW/GENOME/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,".txt.grm.gz",sep="");
                    if(file.exists(GRMFILE)){
                        READGRM=1;
                        cat("Found GRM: ",GRMFILE,"  Let's use it.\n");
                    }
                }
                
                if(READGRM == 1){
                    if(WIDTH=="CHR"){
                        GRMFILE=paste(INDIR,"GENOW/CHR/",CHR,"_K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,".txt.grm.gz",sep="");
                        if(file.exists(GRMFILE)){
                            GZIN = gzfile(GRMFILE,"r");
                            IBSg = matrix(scan(GZIN),byrow=T,nc=Ktmp);
                            cat(dim(IBSg),"\n");
                            close(GZIN);
                        }
                        else{
                            cat("GRM file ",GRMFILE," does not exist\n");
                            quit("no");
                        }
                    }
                    else if(WIDTH=="GENOME"){
                        GRMFILE=paste(INDIR,"GENOW/GENOME/K=",Ktmp,"_MAC=",minMAC,COVARFLAGS,".txt.grm.gz",sep="");
                        if(file.exists(GRMFILE)){
                            GZIN = gzfile(GRMFILE,"r");
                            IBSg = matrix(scan(GZIN),byrow=T,nc=Ktmp);
                            cat(dim(IBSg),"\n");
                            close(GZIN);
                        }
                        else{
                            cat("GRM file ",GRMFILE," does not exist\n");
                            quit("no");
                        }
                    }
                    else{
                        cat("sorry, only CHR and GENOME GRM files are supported\n");
                        quit("no");
                    }
                }
                cat("OUT =",HEOUTFILE,"\n");
                Ktest=Ktmp;

                if(READGRM==0){
                    ##generate kinship matrices that we are using to test with stratified by frequency in Ktest bins
                    Ks = list()
                    Kt = list()
                    FBIN = matrix(,nc=Ktest,nr=3) #by row: 1) upper frequency, 2) number of SNPs, 3) expected h2 if SIM==1
                    st = 1

                    if(SNPLIST == ""){ #no SNP partitions provided, do it by frequency
                        for (bin in 1:Ktest) {  #break SNPs up into bins
                            cat("working on pulling out bin ",bin,"\n")
                            if(st <= M){
                                en = st+floor((M-st+1)/(Ktest-bin+1))-1  #pull out first quantile
                                en = max(which(sefreqs$x == sefreqs$x[en])); #keep SNPs with same MAF together
                                FREQ = sefreqs$x[en]
                                
                                while(en-st+1 < MINSNPS & en < M){
                                    FREQ = FREQ+1/2/N;
                                    en = max(which(sefreqs$x <= FREQ));
                                }
                                if(M-en < MINSNPS){
                                    en = M;
                                }
                                
                                Mb = en-st+1
                                if(MISSINGDATA){
                                    Ks[[bin]]=apply(matrix(1:N,nr=N),1,function(x){apply(matrix(1:N,nr=N),1,function(y){sum(zgenos[x,sefreqs$ix[st:en]]*zgenos[y,sefreqs$ix[st:en]],na.rm=T)/sum(!is.na(zgenos[x,sefreqs$ix[st:en]]*zgenos[y,sefreqs$ix[st:en]]))})})
                                    cat("missing data: N=",N,"; (N,2)=",choose(N,2),"; dim(Ks[",bin,"]=",dim(Ks[[bin]]),"\n",sep="");
                                }
                                else{
                                    Ks[[bin]] = zgenos[,sefreqs$ix[st:en]]%*%t(zgenos[,sefreqs$ix[st:en]])/(Mb)
                                }
                                
                                FBIN[1,bin] = sefreqs$x[en]
                                FBIN[2,bin] = Mb
                                
                                if(SIM == 1){ #simulation, calculate expected h2 in each bin
                                    FBIN[3,bin] = 0;
                                    bcs2 = subset(causals,efreqsBK[causals]>=sefreqs$x[st] &
                                                      efreqsBK[causals]<=sefreqs$x[en])
                                    bidx2 = which(efreqsBK[causals]>=sefreqs$x[st] &
                                                      efreqsBK[causals]<=sefreqs$x[en])
                                    
                                    if (length(bcs2)==0) {
                                        FBIN[3,bin] = 0
                                    }
                                    if (length(bcs2)==1) {
                                        FBIN[3,bin]=var(genosBK[,bcs2]*betas[bidx2])
                                    }
                                    if (length(bcs2)>1) {
                                        FBIN[3,bin]=var(as.matrix(genosBK[,bcs2])%*%betas[bidx2])
                                    }                                    
                                }
                                st = en+1;
                            }
                            else{
                                Ktest = bin-1;
                                break;
                            }
                        }
                    }
                    else{
                        for(bin in 1:Ktest) {  #partition SNPs up into bins
                            cat("working on pulling out bin ",bin,"\n")
                            if(length(SNPPOS) != length(sefreqs$ix)){
                                cat("oops... length(SNPPOS)=",length(SNPPOS)," != length(sefreqs$ix)=",length(sefreqs$ix),sep="");
                                quit("no")
                            }
                            set = which(SNPPOS %in% SNPSET[bin,]) #index of SNPSET SNPs
                            seix = which(sefreqs$ix %in% set)
                            Mb = length(set); #number of SNPs in SNPSET
                            
                            Ks[[bin]] = zgenos[,sefreqs$ix[seix]]%*%t(zgenos[,sefreqs$ix[seix]])/(Mb)
                            FBIN[1,bin] = max(sefreqs$x[seix])
                            FBIN[2,bin] = Mb
                            
                            cat("max MAF=",FBIN[1,bin],"  #SNPs =",FBIN[2,bin],"\n");
                        }
                    }
                    if(SIM == 1){
                        FBIN[3,] = h2*FBIN[3,]/sum(FBIN[3,])
                    }
                    
                    cat("got SNP bins, now K matrix\n");
                    ##pull elements of kinship matrix for HE regression.
                    IBSg = NULL
                    for (bin in 1:Ktest){
                        Kt[[bin]] = Ks[[bin]]
                        foo = Ks[[bin]]
                        IBSg = cbind(IBSg,foo[upper.tri(foo)])
                    }
                    cat("dim(IBSg)=",dim(IBSg),"\n");
                }
                
                perm=0;
                while(perm < nPERM){ #will run at least once
                    if(PERM == 3 & nPERM > 1){
                        cat("you broke this code... PERM=3 only works if nPERM==1\n");
                        quit('n');
                    }
                    if(PERM == 1){ #randomly permute all indivs
                        SAMPORD = sample(1:N);
                        nphenos = nphenos[SAMPORD];
                        
                        ##product of corresponding pairs of phenotypes
                        pp = nphenos%*%t(nphenos)
                        pp = pp[upper.tri(pp)]
                        cat("permuted phenotypes: ",perm,"/",nPERM,"\n",sep="");
                        ##permpcs = pcs[sample(. . . )]
                    }
                    else if(PERM == 2){ #permute within each pop group only maintain structure
                        SAMPORD = 1:N;
                        for(foo in unique(INDPOP)){
                            SAMPORD[which(INDPOP==foo)] = sample(which(INDPOP==foo));
                        }
                        nphenos = nphenos[SAMPORD];
                        
                        ##product of corresponding pairs of phenotypes
                        pp = nphenos%*%t(nphenos)
                        pp = pp[upper.tri(pp)]
                        cat("permuted phenotypes: ",perm,"/",nPERM,"\n",sep="");
                    }
                    if(runHE & (RERUN==1 | !file.exists(paste(HEOUTFILE,".gz",sep=""))
                                | perm>0)){
                        cat("running HE\n");
                        ##H-E regression woohoo!  put a try catch here in case of error
                        if(NORESID == 0 | (COVARCOLS[1]==-1 & ngPC==0 & npPC==0)){
                            HEmodel = try(lm(pp~IBSg));
                            cat("no covars\n");
                        }
                        else{
                            nPCS = npPC + ngPC + length(COVARCOLS);
                            PCS = matrix(0,nrow=length(pp),ncol=nPCS)
                            
                            PCScol = 0;
                            
                            if(COVARCOLS[1] >= 0){
                                cat("controling for ",length(COVARCOLS)," covariates: ",COVARCOLS,"\n");
                                COVAR = read.table(COVARFILE,skip=1,sep="\t");
                                cat("dim(COVAR)=",dim(COVAR),"; pulling out ",COVARCOLS,"\n");
                                for(p in COVARCOLS){
                                    if(is.numeric(COVAR[,p])){
                                        ##can only include numeric covariates
                                        PCScol = PCScol+1;
                                        pc = COVAR[,p];
                                        pcp = scale(pc)%*%t(scale(pc))
                                        pcp = pcp[upper.tri(pcp)]
                                        PCS[,PCScol] = pcp
                                    }
                                    else{
                                        cat("Can only include numeric covriates...\n");
                                        cat("column ",p," is no not numeric: ",COVAR[,p],"\n");
                                        quit("no");
                                    }
                                }
                            }
                            if(npPC > 0){
                                pcs = read.table(paste(INDIR,"/matrix.eqtl/covariates.eur.matrix.eqtl.txt",sep=""),skip=1)
                                pcs = pcs[,2:ncol(pcs)]
                                pcs = pcs[,which(PCnames %in% phenonames)]
                                pcs = as.matrix(pcs)
                                for(p in 1:npPC) {
                                    PCScol = PCScol+1;
                                    pc = pcs[p,]
                                    pcp = scale(pc)%*%t(scale(pc))
                                    pcp = pcp[upper.tri(pcp)]
                                    PCS[,PCScol] = pcp
                                }
                            }

                            if(ngPC > 0){
                                gpcs=read.table(paste(INDIR,"/matrix.eqtl/GEUVADIS_EUR.mac_2.pca.loadings.txt",sep=""),skip=1);
                                ##take subset to desired #PCs
                                        #gpcs=as.matrix(gpcs[match(gpcs[,1], phenonames), 2:(ngPC+1)])
                                for(p in 1:ngPC) {
                                    PCScol = PCScol+1;
                                    pc = pcs[match(gpcs[,1], phenonames),p]
                                    pcp = scale(pc)%*%t(scale(pc))
                                    pcp = pcp[upper.tri(pcp)]
                                    PCS[,PCScol] = pcp
                                }
                            }
                            cat("dim(PCS)=",dim(PCS),"\n")
                            HEmodel = try(lm(pp~IBSg + PCS));
                        }

                        if(OUTPUTSTATS != 0){
                            cat("pp=c(",pp[1],sep="");
                            for(t in 2:length(pp)){
                                cat(", ",pp[t],sep="");
                            }
                            cat(")\n");
                            for(tk in 1:Ktmp){
                                if(OUTPUTSTATS==tk){
                                    cat("IBSg=c(",IBSg[1,tk],sep="");
                                }
                                else if(OUTPUTSTATS==-1){
                                    cat("IBSg[",tk,"]=c(",IBSg[1,tk],sep="");
                                }
                                if(OUTPUTSTATS==tk | OUTPUTSTATS==-1){
                                    for(t in 2:nrow(IBSg)){
                                        cat(", ",IBSg[t,tk],sep="");
                                    }
                                    cat(")\n");
                                }
                            }
                        }
                        if(class(HEmodel) == "try-error"){
                            cat("HE failed for outfile ",HEOUTFILE,"\n");
                        }
                        else{
                            HEpassed = 1;
                            coef = summary(HEmodel)$coef
                            ##cat("dim(coef) = ",dim(coef),"\nROWNAMES=",rownames(coef),"\nCOLNAMES=",colnames(coef),"\n",apply(coef,1,function(x){cat(x,"\n")}),"\n");
                            sres = t(coef[2:(Ktest+1),1]);
                            srse = t(coef[2:(Ktest+1),2]);
                            srp = t(coef[2:(Ktest+1),4]);
                            if(Ktest < Ktmp){
                                sres = c(sres[1,],rep(NA,Ktmp-Ktest))
                                srse = c(srse[1,],rep(NA,Ktmp-Ktest))
                                srp = c(srp[1,],rep(NA,Ktmp-Ktest))
                            }
                            if(nrow(coef) > Ktest+1){
                                sres = c(sres,coef[(Ktest+2):nrow(coef),1]);
                                srse = c(srse,coef[(Ktest+2):nrow(coef),2]);
                                srp = c(srp,coef[(Ktest+2):nrow(coef),4]);
                            }
                            
                            if(nPERM > 1 & perm > 0){
                                APPEND="T";
                                system(paste("gunzip -q ",paste(HEOUTFILE,".gz",sep="")))

                            }
                            if(nSIM == 1){
                                if(SIM==0){ #not a simulation
                                    if(READGRM==1){
                                        write(t(c(sum(sres,na.rm=T),sres,srse,srp)),file=HEOUTFILE, ncolumns=3*Ktmp+1+ifelse(exists("nPCS"),3*nPCS,0), append=APPEND);
                                        system(paste("gzip -f ",HEOUTFILE))
                                    }
                                    else{
                                        write(t(c(sum(sres,na.rm=T),sres,FBIN[1,],FBIN[2,],srse,srp)),file=HEOUTFILE, ncolumns=5*Ktmp+1+ifelse(exists("nPCS"),3*nPCS,0), append=APPEND);
                                        system(paste("gzip -f ",HEOUTFILE))
                                    }
                                }
                                else{ #include expected h2 in each bin
                                    write(t(c(sum(sres,na.rm=T),sres,FBIN[1,],FBIN[2,],h2,FBIN[3,],srse,srp)),file=HEOUTFILE, ncolumns=6*Ktmp+2+ifelse(exists("nPCS"),3*nPCS,0));
                                    system(paste("gzip -f ",HEOUTFILE))
                                }
                            }
                            else{
                                if(SIM == 0){ #not a simulation
                                    if(READGRM==0){
                                        if(file.exists(HEOUTFILE)){
                                            system(paste("gunzip ",HEOUTFILE))
                                        }
                                        write(t(c(sum(sres,na.rm=T),sres,srse,srp)),file=HEOUTFILE, ncolumns=3*Ktmp+1+ifelse(exists("nPCS"),3*nPCS,0),append=T);
                                        system(paste("gzip -f ",HEOUTFILE))
                                    }
                                    else{
                                        if(file.exists(HEOUTFILE)){
                                            system(paste("gunzip ",HEOUTFILE))
                                        }
                                        write(t(c(sum(sres,na.rm=T),sres,FBIN[1,],FBIN[2,],srse,srp)),file=HEOUTFILE, ncolumns=5*Ktmp+1+ifelse(exists("nPCS"),3*nPCS,0),append=T);
                                        system(paste("gzip -f ",HEOUTFILE))
                                    }
                                }
                                else{
                                    if(file.exists(paste(HEOUTFILE,".gz",sep=""))){
                                        system(paste("gunzip ",HEOUTFILE,".gz",sep=""))
                                    }
                                    write(t(c(sum(sres,na.rm=T),sres,FBIN[1,],FBIN[2,],h2,FBIN[3,],srse,srp)),file=HEOUTFILE, ncolumns=6*Ktmp+2+ifelse(exists("nPCS"),3*nPCS,0),append=T);
                                    system(paste("gzip -f ",HEOUTFILE))
                                }
                            }
                        }
                    }

                    if(runLMM & REMLFAILED==0 &
                       (RERUN==1 | !file.exists(paste(REMLOUTFILE,".gz",sep="")))){
                        if(READGRM != 0){
                            cat("Cannot run REML with READGRM yet...\n");
                            quit("no");
                        }
                        cat("running LMM\n");
                        
                        ##LMM ESTIMATION: put a try catch here as well, much more likely to get convergence problems especially with many bins
                        Kt2 = Ktest+1
                        lmmRes <- try(aiREML(Kt,nphenos,rep(1,N),rep(1/Kt2,Kt2),verbose=F))
                        if(class(lmmRes) == "try-error"){
                            cat("aiREML failed for outfile ",REMLOUTFILE,"\n");
                            REMLFAILED = 1;
                        }
                        else{
                            REMLpassed=1;
                            remlOUT = c(sum(lmmRes$vars[1:Ktest])/sum(lmmRes$vars),
                                lmmRes$vars[1:Ktest]/sum(lmmRes$vars));
                            if(Ktest < Ktmp){
                                remlOUT = c(remlOUT, rep(NA,Ktmp-Ktest))
                            }
                            
                            if(nPERM > 1){
                                APPEND="T";
                                system(paste("gunzip -f ",paste(REMLOUTFILE,".gz",sep="")))
                            }
                            
                            if(nSIM == 1 & nPERM==0){
                                if(SIM==0){ #not a simulation
                                    write(t(c(remlOUT,FBIN[1,],FBIN[2,])),file=REMLOUTFILE, ncolumns=3*Ktmp+1);
                                    system(paste("gzip -f ",REMLOUTFILE))
                                    
                                }
                                else{ #simulation, include expected h2 in output
                                    write(t(c(remlOUT,FBIN[1,],FBIN[2,],h2,FBIN[3,])),file=REMLOUTFILE, ncolumns=4*Ktmp+2);
                                    system(paste("gzip -f ",REMLOUTFILE))
                                }
                            }
                            else{
                                if(SIM == 0){ #not a simulation
                                    if(file.exists(REMLOUTFILE)){
                                        system(paste("gunzip ",REMLOUTFILE))
                                    }
                                    write(t(c(remlOUT,FBIN[1,],FBIN[2,])),file=REMLOUTFILE, ncolumns=3*Ktmp+1,append=T);
                                    system(paste("gzip -f ",REMLOUTFILE))
                                }
                                else{
                                    if(file.exists(REMLOUTFILE)){
                                        system(paste("gunzip ",REMLOUTFILE))
                                    }
                                    write(t(c(remlOUT,FBIN[1,],FBIN[2,],h2,FBIN[3,])),file=REMLOUTFILE, ncolumns=4*Ktmp+2,append=T);
                                    system(paste("gzip -f ",REMLOUTFILE))
                                }
                            }
                        }
                    }
                    perm=perm+1;
                }
            }
        }
    }
}


