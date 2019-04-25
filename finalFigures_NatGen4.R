setwd("~/Dropbox/Research/Data/GEUVADIS/INFW_HEvREML5/OUTSUM");
FIGDIR="~/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures6/"
require(RColorBrewer)
expGenes = read.table("~/Dropbox/Research/Data/GEUVADIS/GENESETS/goodExp.txt")
pcol=brewer.pal(9, "Set1")
pcol=c(pcol,"black")
tmp=pcol[1];
pcol[1] = pcol[7];
pcol[7] = pcol[5];
pcol[5] = tmp;

library(viridis)

## Notes for Nat Gen figure requirements:
##
##  1-column figure: w=3.46457 in h=5.11811-8.66142 in
w1col=3.46457
##
##  2-column figure: w=7.08661 in h=7.28346-8.85827 in
w2col=7.08661
##
##
##
##
##

{
                                        #MAIN FIGURES
    {  #FIGURE 1!!  Main simulation results
        NSIMTYPE=440#1760; 
        PLOTIT=1
        PLOTSUP=0;
        yl = c(0,0); #pseudo global variables.  Need to run block below twice.
        xl = 0;
        { # 1A-G: Violin plots of relative bias per SNP bin
            h2=c(0.2); #c(0.02,0.05,0.1,0.5)#c(0.2)
            K=c(1,2,5,10,20)
            minMAC=1;#c(1,2,5,36);
            FT=c(0,0.01,0.05,0.1)
            FR=c(0.01,0.05,0.1,0.5,1)
            Nc=c(1,10,100,1000)
            BM=c(1,2)
            rho=c(0,0.5,0.8,0.9,0.95,1.0)
            tau=c(0.5,0.8,1.0,1.5)
            CM=c(1,5,8)

            SIMOUT = matrix(,nr=1000,nc=18);
            nSIMOUT = 0;
            Eh2 = array(,dim=c(NSIMTYPE,20));
            MAF = 0;
            CLvec = 0;
            PCHvec = 0;
            
            if(PLOTIT){
                pdf(paste(FIGDIR,"/Main/1B-F.pdf",sep=""),height=2,width=w2col,pointsize=10)
                par(mfrow=c(1,length(K)),mar=c(1,1.5,.5,.5),oma=c(1.7,2.5,0,0))
            }
            panel = 0;
            totCombos=0;
            totSims=0;
            for(h in h2){ #rows
                toth2 = array(,dim=c(0,NSIMTYPE));
                STAT.he = array(,dim=c(length(K),NSIMTYPE));
                Kcnt=0;
                for(k in K){ #columns
                    Kcnt = Kcnt+1;
                    panel=panel+1;
                    ##get data for each panel
                    DATA.he = matrix(,nc=4*max(K)+3,nr=0);
                    hedatapoints=0;
                    dp=0;
                    for(mm in minMAC){
                        for(cm in CM){
                            for(bm in BM){
                                for(nc in Nc){
                                    for(fr in FR){
                                        for(ft in FT){
                                            for(r in rho){
                                                for(t in tau){
                                                    for(type in c("","_tenBoyko")){
                                                        INFILEhe=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_he.h2.gz",sep="");
                                                        
                                                        if(file.exists(INFILEhe)){
                                                            con=gzfile(INFILEhe);
                                                            d=matrix(scan(file=con,quiet=T),byrow=T,nc=4*k+2);
                                                            close(con);
                                                            if(nrow(d)>10){
                                                                dp=dp+1;
                                                                PCHvec[dp] = log10(nc)+1;
                                                                CLvec[dp]=2; #default color
                                                                if(cm==5 & nc==1000 & fr==1 & ft<=0.1){
                                                                    CLvec[dp]=3; #overestimates
                                                                }
                                                                else if(cm==8){#(cm==5 & nc==1 & fr==1 & ft==0.01) | (cm==8 & ((nc==1 & t==1.5) | (nc==10 & t==0.8)))){
                                                                    CLvec[dp]=1; #underestimates
                                                                }

                                                                DATA.he = rbind(DATA.he,rep(NA,4*max(K)+3));
                                                                DATA.he[dp,1:(4*k+3)] = c(apply(d,2,function(x){mean(x,na.rm=T)}),nrow(d));
                                                                if(k == 20){
                                                                    Eh2[dp,] = apply(d[,62+1:20],2,function(x){mean(x,na.rm=T)});
                                                                    if(length(MAF) == 1){
                                                                        MAF = apply(d[,21+1:20],2,function(x){mean(x,na.rm=T)});
                                                                    }
                                                                }
                                                                x = d[,1];
                                                                if(mean(x>h)>0.90){#mean(x)-2*sd(x)/sqrt(length(x))>h){#
                                                                    STAT.he[Kcnt,dp] = 1;
                                                                    cat(totCombos,": K=",k,"; x>",h,"=",mean(x>h),": mm=",mm,"; cm=",cm,"; bm=",bm,"; nc=",nc,"; fr=",fr,"; ft=",ft,"; r=",r,"; t=",t,"\n",sep="");
                                                                }
                                                                else if(mean(x<h)>0.90){#mean(x)+2*sd(x)/sqrt(length(x))<h){#
                                                                    STAT.he[Kcnt,dp] = -1;
                                                                    cat(totCombos,": K=",k,"; x<",h,"=",mean(x>h),": mm=",mm,"; cm=",cm,"; bm=",bm,"; nc=",nc,"; fr=",fr,"; ft=",ft,"; r=",r,"; t=",t,"\n",sep="");
                                                                }
                                                                else{
                                                                    STAT.he[Kcnt,dp] = 0;
                                                                }
                                                                
                                                                se=sd(x,na.rm=T)/sqrt(sum(!is.na(x)));
                                                                ##nSIMOUT = nSIMOUT+1;
                                                                ##SIMOUT[nSIMOUT,1,k,
                                                                
                                                                ##if(DATA.he[dp,1]-1.96*se<h & DATA.he[dp,1]+se>h){
                                                                ##  STAT.he[Kcnt,dp] = 0;
                                                                ##}
                                                                ##else if(DATA.he[dp,1]-1.96*se>h){
                                                                ##    STAT.he[Kcnt,dp] = 1;
                                                                ##}
                                                                ##else{
                                                                ##  STAT.he[Kcnt,dp] = -1;
                                                                ##}
                                                                hedatapoints = hedatapoints+1;
                                                                totCombos=totCombos+1;
                                                                totSims=totSims + nrow(d);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    cat(totCombos,": K=",k,"; h2=",h,": HE=",hedatapoints,"\n",sep="");
                    if(nrow(DATA.he) > 0){
                        for(i in 2:(k+1)){
                            rb.he=(DATA.he[,i]-DATA.he[,3*k+1+i])#/DATA.he[,3*k+1+i]
                            d = density(rb.he)
                            xl = max(xl, d$y);
                            yl = range(yl, d$x);
                        }
                        rb.he=(DATA.he[,1]-DATA.he[,3*k+2])#/DATA.he[,3*k+2]
                        toth2 = rbind(toth2,rb.he);
                        d = density(rb.he)
                        if(k==1){
                            d1=d;
                        }
                        else if(k==2){
                            d2=d;
                        }
                        else if(k==5){
                            d5=d;
                        }
                        else if(k==10){
                            d10=d;
                        }
                        else if(k==20){
                            d20=d;
                        }

                        if(PLOTIT){
                            plot(1,1,pch=16,xlab="",ylab="",xaxt='n',yaxt='n',ylim=yl,xlim=c(0.5,k+0.5),type='n')
                            if(k<10){
                                axis(side=1,at=1:k, labels=c(1:k),padj=-1, cex.axis=0.8)
                            }
                            else{
                                axis(side=1,at=1:k, labels=NA, padj=-1, cex.axis=0.8)
                                axis(side=1,at=c(1,5,10,15,20), padj=-1, cex.axis=0.8)
                            }
                            
                            axis(side=2,padj=1, cex.axis=0.8)
                            axis(side=2,padj=1,at=seq(-1,2,by=0.1),labels=NA, cex.axis=0.8)
                            abline(h=0,col='grey')
                            mtext(side=1,"SNP bin",line=1.5)
                            if(k==1){
                                mtext(side=2,expression(paste("Bias of ",{widehat(h^2)}[bin],sep="")),line=1.2)
                            }
                            mtext(side=3,paste("K=",k,sep=""),line=-1.5);
                            for(i in 2:(k+1)){
                                rb.he=(DATA.he[,i]-DATA.he[,3*k+i+1])#/DATA.he[,3*k+1+i]
                                d = density(rb.he)
                                polygon(x=i-1+c(-d$y,rev(d$y))/max(d$y)/3, y=c(d$x,rev(d$x)))
                            }
                        }
                    }
                }
                if(PLOTIT){
                    dev.off();
                }

                if(PLOTSUP){
                {
                    pdf(paste(FIGDIR,"SUPP/1A_h2=",h,"_minMAC=",mm,".pdf",sep=""),height=1,width=w2col/4,pointsize=7)
                    par(mar=c(2.6,3,0.75,0.5))
                    plot(1,1,xlim=c(0.5,5.5),ylim=range(toth2,na.rm=T), type='n',xaxt='n',yaxt='n',xlab="",ylab="",main="")
                    axis(side=1,at=1:5, labels=c(1,2,5,10,20),padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"Number of SNP bins",line=1.5)
                    ##mtext(side=2,expression({hat(h)}^2),line=1.2)
                    mtext(side=2,expression(paste("Bias of ",{widehat(h^2)}[tot],sep="")),line=1.2)
                    abline(h=0,col='grey',lwd=2)
                    boxplot(t(toth2),outline=F,add=T,xaxt='n',yaxt='n',xlab="",ylab="",lwd=2)
                    if(1){
                        for(i in 1:ncol(toth2)){
                            CL=ifelse(STAT.he[1,i]==0,rgb(0,0,1,0.2),ifelse(STAT.he[1,i]<0,rgb(1,0,0,0.2),rgb(1,0,0,0.2)))
                            points(1+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d1$y[which.min(abs(toth2[1,i]-d1$x))]/max(d1$y,na.rm=T)/3), toth2[1,i],pch=16,col=CL,lwd=0.5)
                            CL=ifelse(STAT.he[2,i]==0,rgb(0,0,1,0.2),ifelse(STAT.he[2,i]<0,rgb(1,0,0,0.2),rgb(1,0,0,0.2)))
                            points(2+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d2$y[which.min(abs(toth2[2,i]-d2$x))]/max(d2$y,na.rm=T)/3), toth2[2,i],pch=16,col=CL,lwd=0.5)
                            CL=ifelse(STAT.he[3,i]==0,rgb(0,0,1,0.2),ifelse(STAT.he[3,i]<0,rgb(1,0,0,0.2),rgb(1,0,0,0.2)))
                            points(3+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d5$y[which.min(abs(toth2[3,i]-d5$x))]/max(d5$y,na.rm=T)/3), toth2[3,i],pch=16,col=CL,lwd=0.5)
                            CL=ifelse(STAT.he[4,i]==0,rgb(0,0,1,0.2),ifelse(STAT.he[4,i]<0,rgb(1,0,0,0.2),rgb(1,0,0,0.2)))
                            points(4+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d10$y[which.min(abs(toth2[4,i]-d10$x))]/max(d10$y,na.rm=T)/3), toth2[4,i],pch=16,col=CL,lwd=0.5)
                            CL=ifelse(STAT.he[5,i]==0,rgb(0,0,1,0.2),ifelse(STAT.he[5,i]<0,rgb(1,0,0,0.2),rgb(1,0,0,0.2)))
                            points(5+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d20$y[which.min(abs(toth2[5,i]-d20$x))]/max(d20$y,na.rm=T)/3), toth2[5,i],pch=16,col=CL,lwd=0.5)
                        }
                    }
                    if(0){
                        polygon(x=1+c(-d1$y,rev(d1$y))/max(d1$y)/3, y=c(d1$x,rev(d1$x)))
                        polygon(x=2+c(-d2$y,rev(d2$y))/max(d2$y)/3, y=c(d2$x,rev(d2$x)))
                        polygon(x=3+c(-d5$y,rev(d5$y))/max(d5$y)/3, y=c(d5$x,rev(d5$x)))
                        polygon(x=4+c(-d10$y,rev(d10$y))/max(d10$y)/3, y=c(d10$x,rev(d10$x)))
                        polygon(x=5+c(-d20$y,rev(d20$y))/max(d20$y)/3, y=c(d20$x,rev(d20$x)))
                    }
                    dev.off()
                }
                if(0){
                    pdf(paste(FIGDIR,"SUPP/simTraj_overUnder.pdf",sep=""),height=2,width=w2col,pointsize=10)
                    par(mfrow=c(1,length(K)),mar=c(1,1.5,.5,.5),oma=c(1.7,2.5,0,0))
                    for(i in 1:5){
                        plot(1,1,type='n',xlim=range(MAF),ylim=c(0,1),log='x')
                        for(j in which(STAT.he[i,]==0)){
                            lines(MAF,cumsum(Eh2[j,])/sum(Eh2[j,]),col=rgb(0,0,0,1))
                        }
                        for(j in which(STAT.he[i,]==-1)){
                            lines(MAF,cumsum(Eh2[j,])/sum(Eh2[j,]),col=rgb(0,0,1,1))
                        }
                        for(j in which(STAT.he[i,]==1)){
                            lines(MAF,cumsum(Eh2[j,])/sum(Eh2[j,]),col=rgb(1,0,0,1))
                        }
                    }
                    dev.off()
                }
                }                
            }
            cat("total param combos = ",totCombos,"\n",sep="");
            cat("total number sims  = ",totSims,"\n",sep="");
            
            if(PLOTIT){
                pdf(paste(FIGDIR,"/Main/1A.pdf",sep=""),height=2,width=w2col,pointsize=10)
                par(mar=c(2.6,3,0.75,0.5))
                plot(1,1,xlim=c(0.5,5.5),ylim=range(toth2), type='n',xaxt='n',yaxt='n',xlab="",ylab="",main="")
                axis(side=1,at=1:5, labels=c(1,2,5,10,20),padj=-1, cex.axis=0.8)
                axis(side=2,padj=1, cex.axis=0.8)
                mtext(side=1,"Number of SNP bins (K)",line=1.5)
                ##mtext(side=2,expression({hat(h)}^2),line=1.2)
                mtext(side=2,expression(paste("Bias of ",{widehat(h^2)}[tot],sep="")),line=1.2)
                abline(h=0,col='grey',lwd=2)
                boxplot(t(toth2),outline=F,add=T,xaxt='n',yaxt='n',xlab="",ylab="",lwd=1)
                if(1){
                    for(foo in 1:2){
                        for(i in 1:ncol(toth2)){
                            ##cl = col2rgb(pcol[8])/255;
                            cl = col2rgb(pcol[c(8,6,2)][CLvec[i]])/255;

                            ##cl = col2rgb(viridis(length(unique(CLvec)),option="C")[CLvec[i]])/255;
                            if(0){
                                if(CLvec[i] == 1){
                                    cl = col2rgb(pcol[8])/255;
                                }
                                else{
                                    cl = col2rgb(pcol[2])/255;
                                }
                            }
                            CL=rgb(cl[1],cl[2],cl[3],0.2)
                            PCH=20
                            ##PCH=c(20,17,18,15)[PCHvec[i]]
                            
                            if(foo==1){
                                if((toth2[1,i]>0 & CLvec[i]==2) | (toth2[1,i]<0 & CLvec[i]==1)){
                                    points(1+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d1$y[which.min(abs(toth2[1,i]-d1$x))]/max(d1$y)/3), toth2[1,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }
                            else{
                                if(CLvec[i]==3 | (toth2[1,i]<0 & CLvec[i]==2) | (toth2[1,i]>0 & CLvec[i]==1)){
                                    points(1+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d1$y[which.min(abs(toth2[1,i]-d1$x))]/max(d1$y)/3), toth2[1,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }

                            if(foo==1){
                                if((toth2[2,i]>0 & CLvec[i]==2) | (toth2[2,i]<0 & CLvec[i]==1)){
                                    points(2+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d2$y[which.min(abs(toth2[2,i]-d2$x))]/max(d2$y)/3), toth2[2,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }
                            else{
                                if(CLvec[i]==3 | (toth2[2,i]<0 & CLvec[i]==2) | (toth2[2,i]>0 & CLvec[i]==1)){
                                    points(2+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d2$y[which.min(abs(toth2[2,i]-d2$x))]/max(d2$y)/3), toth2[2,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }

                            if(foo==1){
                                if((toth2[3,i]>0 & CLvec[i]==2) | (toth2[3,i]<0 & CLvec[i]==1)){
                                    points(3+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d5$y[which.min(abs(toth2[3,i]-d5$x))]/max(d5$y)/3), toth2[3,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }
                            else{
                                if(CLvec[i]==3 | (toth2[3,i]<0 & CLvec[i]==2) | (toth2[3,i]>0 & CLvec[i]==1)){
                                    points(3+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d5$y[which.min(abs(toth2[3,i]-d5$x))]/max(d5$y)/3), toth2[3,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }

                            if(foo==1){
                                if((toth2[4,i]>0 & CLvec[i]==2) | (toth2[4,i]<0 & CLvec[i]==1)){
                                    points(4+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d10$y[which.min(abs(toth2[4,i]-d10$x))]/max(d10$y)/3), toth2[4,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }
                            else{
                                if(CLvec[i]==3 | (toth2[4,i]<0 & CLvec[i]==2) | (toth2[4,i]>0 & CLvec[i]==1)){
                                    points(4+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d10$y[which.min(abs(toth2[4,i]-d10$x))]/max(d10$y)/3), toth2[4,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }
                            
                            if(foo==1){
                                if((toth2[5,i]>0 & CLvec[i]==2) | (toth2[5,i]<0 & CLvec[i]==1)){
                                    points(5+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d20$y[which.min(abs(toth2[5,i]-d20$x))]/max(d20$y)/3), toth2[5,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }
                            else{
                                if(CLvec[i]==3 | (toth2[5,i]<0 & CLvec[i]==2) | (toth2[5,i]>0 & CLvec[i]==1)){
                                    points(5+ifelse(runif(1)<0.5,-1,1)*runif(1,max=d20$y[which.min(abs(toth2[5,i]-d20$x))]/max(d20$y)/3), toth2[5,i],pch=PCH,col=CL,lwd=0.5)
                                }
                            }

                        }
                    }
                }
                if(0){
                    polygon(x=1+c(-d1$y,rev(d1$y))/max(d1$y)/3, y=c(d1$x,rev(d1$x)))
                    polygon(x=2+c(-d2$y,rev(d2$y))/max(d2$y)/3, y=c(d2$x,rev(d2$x)))
                    polygon(x=3+c(-d5$y,rev(d5$y))/max(d5$y)/3, y=c(d5$x,rev(d5$x)))
                    polygon(x=4+c(-d10$y,rev(d10$y))/max(d10$y)/3, y=c(d10$x,rev(d10$x)))
                    polygon(x=5+c(-d20$y,rev(d20$y))/max(d20$y)/3, y=c(d20$x,rev(d20$x)))
                }

                cl = col2rgb(pcol[c(8,6,2)])/255;
                CL1=rgb(cl[1,1],cl[2,1],cl[3,1],0.6)
                CL2=rgb(cl[1,2],cl[2,2],cl[3,2],0.6)
                CL3=rgb(cl[1,3],cl[2,3],cl[3,3],0.6)
                legend("topright",c("Many rare causals with large, positive effects","Evolutionary models","Other models"), col=c(CL3,CL1,CL2), bg=rgb(1,1,1,0.8),pch=20,box.col=rgb(1,1,1,0.8), cex=0.8)
                box()
                dev.off()
            }           
        }
    }
    


    { # Figure 2.
        ## This is based on reading *.meanQuant files
        ##cum h2:  dat.0[,1:K]
        ##mean maxMAF per bin: dat.0[,K+(gPC+pPC)+1:K]
        ##Number of SNPs:  dat.0[,2*K+(gPC+pPC)+1:K]
        
        W="1000000"
        K=20
        gPC=10
        pPC=10
        QN=1
        
        {##get data
            MAF = array(,dim=c(5,20))
            SNPs= array(,dim=c(5,20))
            QMAT = array(,dim=c(2,5,20))
            H2MATS = array(,dim=c(5,20))
            MM = c(36,5,2,1);
            for(i in 1:length(MM)){
                ##first read in data h2 (cumsum)
                filename=paste("merge_1000000_K=20_MAC=",MM[i],"_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="");
                foo=read.table(filename);
                MAF[i,1:K] = unlist(foo[1,K+(gPC+pPC)+1:K])
                SNPs[i,1:K] = unlist(foo[1,2*K+(gPC+pPC)+1:K])
                h2 = foo[1,1:K];
                ql = foo[2,1:K];
                qu = foo[3,1:K];
                ##transform to h2 per bin
                for(j in K:2){
                    h2[j] = h2[j] - h2[j-1];
                    ql[j] = ql[j] - ql[j-1];
                    qu[j] = qu[j] - qu[j-1];
                }
                ##now get permuted values
                filename=paste("merge_1000000_K=20_MAC=",MM[i],"_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="");
                foo=read.table(filename);
                h2p = foo[1,1:K];
                qlp = foo[2,1:K];
                qup = foo[3,1:K];
                ##transform to h2 per bin
                for(j in K:2){
                    h2p[j] = h2p[j] - h2p[j-1];
                    qlp[j] = qlp[j] - qlp[j-1];
                    qup[j] = qup[j] - qup[j-1];
                }
                
                H2MATS[i,1:K] = unlist(h2-h2p)
                QMAT[1,i,1:K] = unlist(ql-qlp) #lower
                QMAT[2,i,1:K] = unlist(qu-qup) #upper
            }

            ##now get gnomAD partition data
            filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_ALLgnomAD.he.meanQuant";
            foo=read.table(paste("ALLgnomAD_SFS_K=",K,".txt",sep=""));
            MAF[5,1:K] = foo[,1];
            MAF[5,1] = 1e-5
            foo=read.table(filename);
            h2 = foo[1,1:K];
            ql = foo[2,1:K];
            qu = foo[3,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2[j] = h2[j] - h2[j-1];
                ql[j] = ql[j] - ql[j-1];
                qu[j] = qu[j] - qu[j-1];
            }
            ##now get permuted values
            filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_ALLgnomAD.he.meanQuant";
            foo=read.table(filename);
            h2p = foo[1,1:K];
            qlp = foo[2,1:K];
            qup = foo[3,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2p[j] = h2p[j] - h2p[j-1];
                qlp[j] = qlp[j] - qlp[j-1];
                qup[j] = qup[j] - qup[j-1];
            }
            
            H2MATS[5,1:K] = unlist(h2-h2p)
            QMAT[1,5,1:K] = unlist(ql-qlp) #lower
            QMAT[2,5,1:K] = unlist(qu-qup) #upper

            GETCI=0;
            NBOOT=1000;
            if(GETCI == 1){
            {##  generate bootstraps for CIs for Fig 2A: cis-trans-corr
                cat("Getting bootstraps for Fig 2A\n");
                    Q2A = array(,dim=c(2,3,20));
                    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.h2.gz",sep="");
                    F=gzfile(filename);
                    h2=read.table(F);
                    ##close(F);
                    dim(h2)
                    h2 = subset(h2,h2[,1] %in% expGenes[,1])
                    dim(h2)
                    
                    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz",sep="");
                    F=gzfile(filename);
                    h2p=read.table(F);
                    ##close(F);
                    dim(h2p)
                    h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
                    dim(h2p);
                    
                    B=array(,dim=c(NBOOT,K));
                    BP=array(,dim=c(NBOOT,K));
                    BD=array(,dim=c(NBOOT,K)); #difference
                    
                    for(i in 1:NBOOT){
                        h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:K];
                        B[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:K];
                        BP[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        BD[i,] = B[i,]-BP[i,];
                    }
                    for(i in 1:K){
                        Q2A[1,1,i] = quantile(B[,i],0.975) #upper Q raw
                        Q2A[2,1,i] = quantile(B[,i],0.025) #lower Q raw
                        Q2A[1,2,i] = quantile(BP[,i],0.975) #upper Q trans
                        Q2A[2,2,i] = quantile(BP[,i],0.025) #lower Q trans
                        Q2A[1,3,i] = quantile(BD[,i],0.975) #upper Q corr
                        Q2A[2,3,i] = quantile(BD[,i],0.025) #lower Q corr
                    }
                }

                
                {## generate bootstraps for CIs for Fig 2B
                    Q2B = array(,dim=c(2,20));
                    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.h2.gz",sep="");
                    F=gzfile(filename);
                    h2=read.table(F);
                    ##close(F);
                    dim(h2)
                    h2 = subset(h2,h2[,1] %in% expGenes[,1])
                    dim(h2)
                    
                    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz",sep="");
                    F=gzfile(filename);
                    h2p=read.table(F);
                    ##close(F);
                    dim(h2p)
                    h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
                    dim(h2p);
                    
                    B=array(,dim=c(NBOOT,K));
                    BP=array(,dim=c(NBOOT,K));
                    BD=array(,dim=c(NBOOT,K)); #difference
                    BR=array(,dim=c(NBOOT,K)); #ratio for proportion
                    
                    for(i in 1:NBOOT){
                        h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:K];
                        B[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:K];
                        BP[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        BD[i,] = B[i,]-BP[i,];
                        BR[i,] = BD[i,]/sum(BD[i,]);
                    }
                    for(i in 1:K){
                        Q2B[1,i] = quantile(BR[,i],0.975) #upper Q
                        Q2B[2,i] = quantile(BR[,i],0.025) #lower Q
                    }
                }

                {## generate bootstraps for CIs for Fig 2Binset
                    Q2Bi = array(,dim=c(2,20));
                    filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_SINGgnomAD.he.h2.gz";
                    F=gzfile(filename);
                    h2=read.table(F);
                    ##close(F);
                    dim(h2)
                    h2 = subset(h2,h2[,1] %in% expGenes[,1])
                    dim(h2)
                    
                    filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_SINGgnomAD.he.h2.gz";
                    F=gzfile(filename);
                    h2p=read.table(F);
                    ##close(F);
                    dim(h2p)
                    h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
                    dim(h2p);
                    
                    B=array(,dim=c(NBOOT,K));
                    BP=array(,dim=c(NBOOT,K));
                    BD=array(,dim=c(NBOOT,K)); #difference
                    BR=array(,dim=c(NBOOT,K)); #ratio for proportion
                    
                    for(i in 1:NBOOT){
                        h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:K];
                        B[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:K];
                        BP[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        BD[i,] = B[i,]-BP[i,];
                        BR[i,] = BD[i,]/sum(BD[i,]);
                    }
                    for(i in 1:K){
                        Q2Bi[1,i] = quantile(BR[,i],0.975) #upper Q
                        Q2Bi[2,i] = quantile(BR[,i],0.025) #lower Q
                    }
                }

                
                {## CIs for 2C: cumulative h2 for each minMAC
                    K=20;
                    nSB=6;
                    Kt=K+nSB-1
                    Q2C = array(,dim=c(2,5,K+nSB-1));
                    for(i in 1:4){
                        filename=paste("merge_1000000_K=20_MAC=",MM[i],"_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.h2.gz",sep="");
                        F=gzfile(filename);
                        h2=read.table(F);
                        ##close(F);
                        dim(h2)
                        h2 = subset(h2,h2[,1] %in% expGenes[,1])
                        dim(h2)
                        
                        filename=paste("merge_1000000_K=20_MAC=",MM[i],"_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz",sep="");
                        F=gzfile(filename);
                        h2p=read.table(F);
                        ##close(F);
                        dim(h2p)
                        h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
                        dim(h2p);
                        
                        B=array(,dim=c(NBOOT,K));
                        BP=array(,dim=c(NBOOT,K));
                        BD=array(,dim=c(NBOOT,K)); #difference
                        
                        for(j in 1:NBOOT){
                            h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:K];
                            B[j,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                            h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:K];
                            BP[j,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                            BD[j,] = cumsum(B[j,]-BP[j,]);
                        }
                        
                        for(j in 1:K){
                            Q2C[1,i,j] = quantile(BD[,j],0.975,na.rm=T) #upper Q
                            Q2C[2,i,j] = quantile(BD[,j],0.025,na.rm=T) #lower Q
                        }
                    }
                    
                    filename=paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_gnomAD.SINGpart.he.h2.gz",sep="")
                    F=gzfile(filename);
                    h2=read.table(F);
                    ##close(F);
                    dim(h2)
                    h2 = subset(h2,h2[,1] %in% expGenes[,1])
                    dim(h2)
                    
                    filename=paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_gnomAD.SINGpart.he.h2.gz",sep="");
                    F=gzfile(filename);
                    h2p=read.table(F);
                    ##close(F);
                    dim(h2p)
                    h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
                    dim(h2p);
                    
                    B=array(,dim=c(NBOOT,Kt));
                    BP=array(,dim=c(NBOOT,Kt));
                    BD=array(,dim=c(NBOOT,Kt)); #difference
                    
                    for(j in 1:NBOOT){
                        h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:Kt];
                        B[j,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:Kt];
                        BP[j,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        BD[j,] = cumsum(B[j,]-BP[j,]);
                    }
                    
                    for(j in 1:Kt){
                        Q2C[1,5,j] = quantile(BD[,j],0.975,na.rm=T) #upper Q
                        Q2C[2,5,j] = quantile(BD[,j],0.025,na.rm=T) #lower Q
                    }
                }
                
                {##bootstraps for CIs for Fig 2D: various variant types
                    Q2D = array(,dim=c(3,4));
                    filename="merge_1000000_K=4_MAC=1_ngPC=0_npPC=0_QN=1_DSAMP=1_singSVT.he.h2.gz";
                    F=gzfile(filename);
                    h2=read.table(F);
                    dim(h2)
                    h2 = subset(h2,h2[,1] %in% expGenes[,1])
                    dim(h2)
                    SNPs2D = apply(h2[,11:14],2,mean)
                    
                    filename="merge_1000000_K=4_MAC=1_ngPC=0_npPC=0_QN=1_DSAMP=1_PERM=3_singSVT.he.h2.gz";
                    F=gzfile(filename);
                    h2p=read.table(F);
                    ##close(F);
                    dim(h2p)
                    h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
                    dim(h2p);
                    
                    B=array(,dim=c(NBOOT,4));
                    BP=array(,dim=c(NBOOT,4));
                    BD=array(,dim=c(NBOOT,4)); #difference
                    
                    for(i in 1:NBOOT){
                        h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:4];
                        B[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:4];
                        BP[i,] = apply(h2t,2,function(x){mean(x,na.rm=T)});
                        BD[i,] = B[i,]-BP[i,];
                    }
                    for(i in 1:4){
                        Q2D[1,i] = quantile(BD[,i],0.975) #upper Q
                        Q2D[2,i] = quantile(BD[,i],0.025) #lower Q
                    }
                    Q2D[3,] = apply(h2[,2+1:4],2,function(x){mean(x,na.rm=T)})-apply(h2p[,2+1:4],2,function(x){mean(x,na.rm=T)});
                }
            }
        }

        
        ## PLOT DATA!! 
        {## Figure 2A: cis vs trans - h2 per bin for observed vs permuted
            K=20
            ##first read in data h2 (cumsum)
            filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="");
            foo=read.table(filename);
            h2 = foo[1,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2[j] = h2[j] - h2[j-1];
            }
            ##now get permuted values
            filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="");
            foo=read.table(filename);
            h2p = foo[1,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2p[j] = h2p[j] - h2p[j-1];
            }

            {
                pdf(paste(FIGDIR,"/Main/2A_raw_perm_cor_h2.pdf",sep=""),height=1.5,width=w2col/2,pointsize=10)
                par(mar=c(2,2.8,0.3,0.3))
                plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.022),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                axis(side=2,padj=1.3,cex.axis=0.7,at=seq(0,0.02,by=0.005),labels=NA)
                axis(side=2,padj=1.3,cex.axis=0.7,at=c(0,0.01,0.02))
                mtext(side=1,"Minor Allele Frequency",line=1)
                mtext(side=2,expression(paste({widehat(h^2)}[bin],sep="")),line=1.15)

                abline(h=seq(0,0.02,by=0.005),col=rgb(0,0,0,0.1),lwd=0.75)
                abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                
                legend("topright",c("Raw PCA-corrected","trans-permutation","Permutation-corrected"), col=c("blue",pcol[8],pcol[4]), lty=1,lwd=c(.5,.5,2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                box()

                polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(Q2A[1,1,],rev(Q2A[2,1,])), col=rgb(0,0,1,0.1),border=NA)
                pci = pcol[8]
                rgbi = col2rgb(pci)/255;
                polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(Q2A[1,2,],rev(Q2A[2,2,])), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.1),border=NA)
                pci = pcol[4]
                rgbi = col2rgb(pci)/255;
                polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(Q2A[1,3,],rev(Q2A[2,3,])), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.1),border=NA)
                points(MAF[4,],h2,pch=20,col=rgb(0,0,1,1))
                lines(MAF[4,],h2,col=rgb(0,0,1,1),lwd=0.5)
                points(MAF[4,],h2p,pch=20,col=pcol[8])
                lines(MAF[4,],h2p,col=pcol[8],lwd=0.5)

                points(MAF[4,],h2-h2p,pch=16,col=pcol[4])
                lines(MAF[4,],h2-h2p,lwd=2,col=pcol[4])
                dev.off()

                STAGES=0;
                if(STAGES){
                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_0.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.022),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7,at=seq(0,0.02,by=0.005),labels=NA)
                    axis(side=2,padj=1.3,cex.axis=0.7,at=c(0,0.01,0.02))
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste({widehat(h^2)}[bin],sep="")),line=1.15)

                    abline(h=seq(0,0.02,by=0.005),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    dev.off()

                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_0cum.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7)
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste("Cum ",{widehat(h^2)}[bin],sep="")),line=1.15)

                    abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    dev.off()



                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_1.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.022),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7,at=seq(0,0.02,by=0.005),labels=NA)
                    axis(side=2,padj=1.3,cex.axis=0.7,at=c(0,0.01,0.02))
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste({widehat(h^2)}[bin],sep="")),line=1.15)

                    abline(h=seq(0,0.02,by=0.005),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    
                    legend("topright",c("Raw PCA-corrected"), col=c("blue",pcol[8],pcol[4]), lty=1,lwd=c(.5,.5,2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                    box()

                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(Q2A[1,1,],rev(Q2A[2,1,])), col=rgb(0,0,1,0.1),border=NA)
                    pci = pcol[8]
                    rgbi = col2rgb(pci)/255;
                    points(MAF[4,],h2,pch=20,col=rgb(0,0,1,1))
                    lines(MAF[4,],h2,col=rgb(0,0,1,1),lwd=0.5)
                    dev.off()

                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_1cum.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7)
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste("Cum ",{widehat(h^2)}[bin],sep="")),line=1.15)

                    abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    
                    legend("topleft",c("Raw PCA-corrected"), col=c("blue",pcol[8],pcol[4]), lty=1,lwd=c(.5,.5,2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                    box()

                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(cumsum(Q2A[1,1,]),rev(cumsum(Q2A[2,1,]))), col=rgb(0,0,1,0.1),border=NA)
                    points(MAF[4,],cumsum(t(h2)),pch=20,col=rgb(0,0,1,1))
                    lines(MAF[4,],cumsum(t(h2)),col=rgb(0,0,1,1),lwd=0.5)
                    dev.off()

                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_1cumNorm.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(0,1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7)
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste("Cum. Prop. ",{widehat(h^2)},sep="")),line=1.15)

                    abline(h=seq(0,1,by=0.1),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    
                    legend("topleft",c("Perm-corrected"), col=c("blue",pcol[8],pcol[4]), lty=1,lwd=c(.5,.5,2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                    box()

                    points(MAF[4,],cumsum(t(h2-h2p))/sum(h2-h2p),pch=20,col=rgb(0,0,1,1))
                    lines(MAF[4,],cumsum(t(h2-h2p))/sum(h2-h2p),col=rgb(0,0,1,1),lwd=0.5)
                    dev.off()

                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_2.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.022),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7,at=seq(0,0.02,by=0.005),labels=NA)
                    axis(side=2,padj=1.3,cex.axis=0.7,at=c(0,0.01,0.02))
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste({widehat(h^2)}[bin],sep="")),line=1.15)

                    abline(h=seq(0,0.02,by=0.005),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    
                    legend("topright",c("Raw PCA-corrected","trans-permutation"), col=c("blue",pcol[8],pcol[4]), lty=1,lwd=c(.5,.5,2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                    box()

                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(Q2A[1,1,],rev(Q2A[2,1,])), col=rgb(0,0,1,0.1),border=NA)
                    pci = pcol[8]
                    rgbi = col2rgb(pci)/255;
                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(Q2A[1,2,],rev(Q2A[2,2,])), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.1),border=NA)
                    points(MAF[4,],h2,pch=20,col=rgb(0,0,1,1))
                    lines(MAF[4,],h2,col=rgb(0,0,1,1),lwd=0.5)
                    points(MAF[4,],h2p,pch=20,col=pcol[8])
                    lines(MAF[4,],h2p,col=pcol[8],lwd=0.5)
                    dev.off()

                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_2cum.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7)
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste("Cum ",{widehat(h^2)}[bin],sep="")),line=1.15)

                    abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    
                    legend("topleft",c("Raw PCA-corrected","trans-permutation"), col=c("blue",pcol[8],pcol[4]), lty=1,lwd=c(.5,.5,2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                    box()

                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(cumsum(Q2A[1,1,]),rev(cumsum(Q2A[2,1,]))), col=rgb(0,0,1,0.1),border=NA)
                    pci = pcol[8]
                    rgbi = col2rgb(pci)/255;
                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(cumsum(Q2A[1,2,]),rev(cumsum(Q2A[2,2,]))), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.1),border=NA)
                    points(MAF[4,],cumsum(t(h2)),pch=20,col=rgb(0,0,1,1))
                    lines(MAF[4,],cumsum(t(h2)),col=rgb(0,0,1,1),lwd=0.5)
                    points(MAF[4,],cumsum(t(h2p)),pch=20,col=pcol[8])
                    lines(MAF[4,],cumsum(t(h2p)),col=pcol[8],lwd=0.5)
                    dev.off()
                    
                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_3cum.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,0.1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7)
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste("Cum ",{widehat(h^2)}[bin],sep="")),line=1.15)

                    abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    
                    legend("topleft",c("Raw PCA-corrected","trans-permutation","Permutation-corrected"), col=c("blue",pcol[8],pcol[4]), lty=1,lwd=c(.5,.5,2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                    box()

                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(cumsum(Q2A[1,1,]),rev(cumsum(Q2A[2,1,]))), col=rgb(0,0,1,0.1),border=NA)
                    pci = pcol[8]
                    rgbi = col2rgb(pci)/255;
                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(cumsum(Q2A[1,2,]),rev(cumsum(Q2A[2,2,]))), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.1),border=NA)
                    pci = pcol[4]
                    rgbi = col2rgb(pci)/255;
                    polygon(x=c(MAF[4,],rev(MAF[4,])), y=c(cumsum(Q2A[1,3,]),rev(cumsum(Q2A[2,3,]))), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.1),border=NA)
                    points(MAF[4,],cumsum(t(h2)),pch=20,col=rgb(0,0,1,1))
                    lines(MAF[4,],cumsum(t(h2)),col=rgb(0,0,1,1),lwd=0.5)
                    points(MAF[4,],cumsum(t(h2p)),pch=20,col=pcol[8])
                    lines(MAF[4,],cumsum(t(h2p)),col=pcol[8],lwd=0.5)

                    points(MAF[4,],cumsum(t(h2-h2p)),pch=16,col=pcol[4])
                    lines(MAF[4,],cumsum(t(h2-h2p)),lwd=2,col=pcol[4])
                    dev.off()

                    pdf(paste(FIGDIR,"/NotUsed/2A_raw_perm_cor_h2_cdf.pdf",sep=""),height=2,width=2.75,pointsize=10)
                    par(mar=c(2,2.8,0.3,0.3))
                    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-0.001,1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                    axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                    axis(side=2,padj=1.3,cex.axis=0.7)
                    mtext(side=1,"Minor Allele Frequency",line=1)
                    mtext(side=2,expression(paste("CDF ",{widehat(h^"2'")}[bin],sep="")),line=1.15)

                    abline(h=seq(0,1,by=0.1),col=rgb(0,0,0,0.1),lwd=0.75)
                    abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                    
                    legend("topleft",c("Permutation-corrected"), col=c(pcol[4]), lty=1,lwd=c(2), bg=rgb(1,1,1,0.8),box.col="white",pch=c(20,20,16),cex=0.8)
                    box()

                    pci = pcol[4]
                    rgbi = col2rgb(pci)/255;
                    points(MAF[4,],cumsum(t(h2-h2p))/sum(h2-h2p),pch=16,col=pcol[4])
                    lines(MAF[4,],cumsum(t(h2-h2p))/sum(h2-h2p),lwd=2,col=pcol[4])
                    dev.off()
                }
            }
        }

        
        { ##plot figure 2B: prop h2 per bin for minMAC=1
            pdf(paste(FIGDIR,"/Main/2B_prop_h2_vs_MAF_minMAC=1.pdf",sep=""),height=1.5,width=w2col/2,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(-.03,0.28),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,at=c(1e-3,1e-2,1e-1,0.5),labels=c("1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7,at=seq(0,0.25,by=0.05),labels=NA)
            axis(side=2,padj=1.3,cex.axis=0.7,at=c(0,.1,.2),labels=c("0","0.1","0.2"))
            axis(side=2,padj=1.3,cex.axis=0.7,at=c(0.05,.15,.25),labels=c("0.05","0.15","0.25"))
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Prop. of  ",{widehat(h^"2'")}[total],sep="")),line=1.15)
            abline(h=seq(0,1,by=0.05),col=rgb(0,0,0,0.1),lwd=0.75)
            abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
            box()
            for(i in 4){
                pci = pcol[i]
                rgbi = col2rgb(pci)/255;
                
                lines(MAF[i,1:K],H2MATS[i,1:K]/sum(H2MATS[i,1:K]),col=pcol[i],lwd=0.5)
                for(j in 1:K){
                    segments(x0=MAF[i,j],y0=Q2B[1,j],y1=Q2B[2,j],col=pcol[i])
                }
                points(MAF[i,1:K],H2MATS[i,1:K]/sum(H2MATS[i,1:K]),col=pcol[i],pch=20)
            }
            
            dev.off()
        }
        
        { ##plot figure 2B_inset: prop h2_singleton per bin for gnomAD
            K=20
            SING=array(,dim=K);
            ##now get gnomAD partition data
            filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_SINGgnomAD.he.meanQuant"
            foo=read.table(filename);
            h2 = foo[1,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2[j] = h2[j] - h2[j-1];
            }
            ##now get permuted values
            filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_SINGgnomAD.he.meanQuant";
            foo=read.table(filename);
            h2p = foo[1,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2p[j] = h2p[j] - h2p[j-1];
            }
            SING[1:K] = unlist(h2-h2p)
            
            pdf(paste(FIGDIR,"/Main/2Binset_prop_h2_vs_MAF_minMAC=1.pdf",sep=""),height=0.8,width=1.7,pointsize=7)
            par(mar=c(1,1.1,0.1,0.1),oma=c(.7,1,0.05,0.05))
            plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=c(-0.05,0.9),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=NA,padj=-2, cex.axis=0.72)
            axis(side=1,at=c(1e-5,1e-4,1e-2,0.5),labels=c(paste("*"),"1e-4","0.01","0.5"),padj=-2, cex.axis=0.72)
            axis(side=2,padj=1.3,cex.axis=0.72)
            mtext(side=2,expression(paste("Prop ",{widehat(h^"2'")}[singleton],sep="")),line=.8,cex=0.7)
            mtext(side=1,"gnomAD MAF",line=.7,cex=0.7)
            abline(h=seq(0,1,by=0.2),col=rgb(0,0,0,0.1),lwd=0.75)
            abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
            box()
            for(i in 5){
                pci = pcol[i]
                rgbi = col2rgb(pci)/255;
                
                lines(MAF[i,1:K],SING/sum(SING),col=pcol[i],lwd=0.5)
                for(j in 1:K){
                    segments(x0=MAF[i,j],y0=Q2Bi[1,j],y1=Q2Bi[2,j],col=pcol[i])
                }
                ##points(MAF[i,1:K],SING/sum(SING),col=pcol[i],pch=20,lwd=0.01)
                symbols(x=MAF[i,1:K], y=SING/sum(SING), circles=rep(1,20), bg=pcol[i], fg=NULL,add=T,inches=.015)
            }

            dev.off()
        }

        
        { ##plot figure 2C: cum h2 vs MAF across minMAC
            K=20
            nSB=6;
            Kt=K+nSB-1;
            pdf(paste(FIGDIR,"/Main/2C_cum_h2_vs_MAF_by_minMAC.pdf",sep=""),height=1.5,width=w2col/2,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=c(0,0.09),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Cumulative ",{widehat(h^"2'")},sep="")),line=1.15)
            abline(h=seq(0,0.1,by=0.02),col=rgb(0,0,0,0.1),lwd=0.75)
            abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
            
            legend("topleft",c("gnomAD",expression("MAC" >= 1),expression("MAC" >= 2),expression("MAC" >= 5),expression("MAF" >= "5%")), col=pcol[5:1], lty=1,lwd=2,bg=rgb(1,1,1,0.8),ncol=2,box.col="white",pch=20,cex=0.8)
            box()
            for(i in 1:4){
                pci = pcol[i]
                rgbi = col2rgb(pci)/255;
                
                polygon(x=c(MAF[i,1:K], rev(MAF[i,1:K])), y=c(Q2C[1,i,1:K],rev(Q2C[2,i,1:K])), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.2),border=NA)
                lines(MAF[i,1:K],cumsum(H2MATS[i,1:K]),col=pcol[i],lwd=0.5)
                points(MAF[i,1:K],cumsum(H2MATS[i,1:K]),col=pcol[i],pch=20)
            }
            
            d1=read.table(paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_gnomAD.SINGpart.meanQuant",sep=""))
            d2=read.table(paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_gnomAD.SINGpart.meanQuant",sep=""))
            SINGMAF = c(1e-5, 6.46e-5, 0.0002, 0.0004, 1e-3, 1/720);
            gMAF = c(SINGMAF[1:nSB], MAF[4,2:K])
            t1=unlist(d1[1,1:Kt])
            tb1 = t1;
            t2=unlist(d2[1,1:Kt])
            tb2 = t2;
            for(i in Kt:2){
                tb1[i] = t1[i]-t1[i-1];
                tb2[i] = t2[i]-t2[i-1];
            }
            pci = pcol[5]
            rgbi = col2rgb(pci)/255;
            polygon(x=c(gMAF, rev(gMAF)), y=c(Q2C[1,5,1:Kt],rev(Q2C[2,5,1:Kt])), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.2),border=NA)
            lines(gMAF, cumsum(tb1-tb2),col=pcol[5])
            points(gMAF, cumsum(tb1-tb2),col=pcol[5],pch=20)
            dev.off()

            tb3=cumsum(tb1-tb2);
        }

        { ## plot h2/SNP for gnomAD singleton partition
            d1=read.table(paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_gnomAD.SINGpart.meanQuant",sep=""))
            d2=read.table(paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_gnomAD.SINGpart.meanQuant",sep=""))
            SINGMAF = c(1e-5, 6.46e-5, 0.0002, 0.0004, 1e-3, 1/720);
            gMAF = c(SINGMAF[1:nSB], MAF[4,2:K])
            t1=unlist(d1[1,1:Kt])
            tb1 = t1;
            t2=unlist(d2[1,1:Kt])
            tb2 = t2;
            for(i in Kt:2){
                tb1[i] = t1[i]-t1[i-1];
                tb2[i] = t2[i]-t2[i-1];
            }
            nSNPs=d1[1,70+1:25]
            h2pS=(tb1-tb2)/nSNPs;
            h2pSn=((tb1-tb2)/sum((tb1-tb2)))/(nSNPs/sum(nSNPs));
            pdf(paste(FIGDIR,"/SUPP/S23_h2perSNP_gnomAD.pdf",sep=""),height=1.5,width=2.75,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=range(h2pS),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste({widehat(h^"2'")}," per SNP",sep="")),line=1.15)
            abline(h=seq(0,1.5e-5,length=4),col=rgb(0,0,0,0.1),lwd=0.75)
            abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
            
            lines(gMAF, h2pS,col=pcol[5])
            points(gMAF, h2pS,col=pcol[5],pch=20)
            dev.off()

            pdf(paste(FIGDIR,"/SUPP/S23_h2perSNP_gnomAD_rel.pdf",sep=""),height=1.5,width=2.75,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=range(h2pSn),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste({widehat(h^"2'")}," enrichment",sep="")),line=1.15)
            abline(h=1,col=rgb(0,0,0,0.1),lwd=0.75)
            ##abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)

            lines(gMAF, h2pSn,col=pcol[5])
            points(gMAF, h2pSn,col=pcol[5],pch=20)
            dev.off()

        }


        {##Figure 2D: h2_singleton for multiple classes
            LEG = c("Neanderthal","Globally rare","Other Singletons","INDELs/SVs");
            pdf(paste(FIGDIR,"/Main/2D_h2_vs_numSing_SVcat.pdf",sep=""),height=1.5,width=w2col/2,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=range(SNPs2D,na.rm=T), ylim=c(-0.001,0.045),xaxt='n',yaxt='n',xlab="",ylab="",log='x')
            abline(h=seq(-.04,.06,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
            axis(side=1,padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Mean singletons per gene (log scale)",line=1)
            mtext(side=2,expression(paste({widehat(h^"2'")}[singleton],sep="")),line=1.15)
            
            PLOTORD = c(1,4,3,2);
            colord=c(9,7,4,5)
            for(j in 1:4){
                segments(x0=SNPs2D[PLOTORD[j]],y0=Q2D[1,PLOTORD[j]],y1=Q2D[2,PLOTORD[j]],col=pcol[colord[j]])
            }
            points(SNPs2D[PLOTORD],Q2D[3,PLOTORD],col=pcol[colord],pch="-")
            legend("topleft",rev(LEG[PLOTORD]), col=rev(pcol[colord]), bg=rgb(1,1,1,0.8),box.col="white",pch="+",cex=0.8)
            box()
            dev.off()           

        }
    }



    { #Figure 3!
        library(viridis)
        ##library(image.plot)
        GETDATA=1
        if(GETDATA){
            read.table(paste(FIGDIR,"/LawrenceData/lsq_all.txt",sep=""))->t
            read.table(paste(FIGDIR,"LawrenceData/real_infer.5000.txt",sep=""))->r

            ##tau: true=t[,2]; inf=t[,8]
            ##rho: true=t[,1]; inf=t[,7]
            ##p: true=t[,3]; inf=t[,9]

            BREAKS = 50;
            Bx = seq(0, 1.001, length=BREAKS+1)
            By = seq(0, 1.001, length=BREAKS+1)

            TAU=matrix(,nr=BREAKS,nc=BREAKS);
            RHO=matrix(,nr=BREAKS,nc=BREAKS);
            P=matrix(,nr=BREAKS,nc=BREAKS);
            JOINT=matrix(,nr=BREAKS,nc=BREAKS);

            for(i in 1:BREAKS){
                cat("pulling out i=",i,"\n");
                for(j in 1:BREAKS){
                    RHO[i,j] = sum(t[,1]>=Bx[i] & t[,1]<Bx[i+1] & t[,7]>=By[j] & t[,7]<By[j+1])
                    TAU[i,j] = sum(t[,2]>=Bx[i] & t[,2]<Bx[i+1] & t[,8]>=By[j] & t[,8]<By[j+1])
                    P[i,j] = sum(t[,3]>=Bx[i] & t[,3]<Bx[i+1] & t[,9]>=By[j] & t[,9]<By[j+1])
                    JOINT[i,j] = sum(r[,1]>=Bx[i] & r[,1]<Bx[i+1] & r[,2]>=By[j] & r[,2]<By[j+1])
                }
            }
            t = read.table("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures5/LawrenceData/real_infer.5000.txt");
            al<-function(m,s) {
                return( (m^2 - m^3 -m*(s^2))/s^2)
            }
            
            be <-function(m,s) {
                return( (m-1)*(m^2-m+s^2)/(s^2))
            }
            
            ## sample marginal tau
            marg_tau<-c()
            numE = 0;
            for (i in seq(1,5000)) {
                numE = numE+1;
                ##cat(numE,": ",t$V1[i], t$V2[i], t$V4[i], t$V5[i],"\n");
                marg_tau<-c(marg_tau,rbeta(100,al(t$V2[i],t$V5[i]),be(t$V2[i],t$V5[i])))
            }
            
            ## samples marginal rho
            marg_rho<-c()
            numS = 0;
            for (i in seq(1,5000)) {
                numS = numS+1;
                marg_rho<-c(marg_rho,rbeta(100,al(t$V1[i],t$V4[i]),be(t$V1[i],t$V4[i])))
            }
            
            ## samples marginal phi
            marg_phi<-c()
            numS = 0;
            for (i in seq(1,5000)) {
                numS = numS+1;
                marg_phi<-c(marg_phi,rbeta(100,al(t$V3[i],t$V6[i]),be(t$V3[i],t$V6[i])))
            }
            
            
            read.table(paste(FIGDIR,"LawrenceData/real_data.posterior.txt",sep=""))->real_infer
            drho = density(real_infer$V1,adjust=.5);
            dtau = density(real_infer$V2,adjust=.5);
            dphi = density(real_infer$V3,adjust=.5);
            BR = seq(0,1,by=0.02)
            hrho = hist(marg_rho,breaks=BR,plot=FALSE)
            htau = hist(marg_tau,breaks=BR,plot=FALSE)
            hphi = hist(marg_phi,breaks=BR,plot=FALSE)
            pci = pcol[1]
            ri = col2rgb(pci)/255;


        }

        if(1){ #3A: marginal posteriors, 3B: joint rho/tau posterior, 3C out of sample draws
            {
                pdf(paste(FIGDIR,"/Main/3A_ridgeline_Marg_posterior.pdf",sep=""),height=1.5,width=w2col/3,pointsize=10)
                par(mar=c(2.15,1.5,0,0.4))
                plot(1,1,type='n',xlim=c(0,1),ylim=c(1,4),xaxt='n',yaxt='n',xlab="",ylab="",bty='n')
                mtext(side=2,"Posterior",line=0.5)
                axis(side=1,padj=-1.5,cex=0.8)
                
                
                pci = pcol[1]
                ri = col2rgb(pci)/255;
                polygon(x=c(dtau$x,min(dtau$x)),y=1+c((dtau$y/max(dtau$y)),0)*0.95,col=rgb(ri[1],ri[2],ri[3],0.1),border="white")
                pci = pcol[2]
                ri = col2rgb(pci)/255;
                polygon(x=c(drho$x,min(drho$x)),y=2+c((drho$y/max(drho$y)),0)*0.95,col=rgb(ri[1],ri[2],ri[3],0.1),border="white")
                pci = pcol[4]
                ri = col2rgb(pci)/255;
                polygon(x=c(dphi$x,min(dphi$x)),y=3+c((dphi$y/max(dphi$y)),0)*0.95,col=rgb(ri[1],ri[2],ri[3],0.1),border="white")
                lines(dtau$x,1+(dtau$y/max(dtau$y))*0.95,col=rgb(0,0,0,0.5))
                lines(drho$x,2+(drho$y/max(drho$y))*0.95,col=rgb(0,0,0,0.5))
                lines(dphi$x,3+(dphi$y/max(dphi$y))*0.95,col=rgb(0,0,0,0.5))


                for(i in 2:length(BR)){
                    pci = pcol[1]
                    r1 = col2rgb(pci)/255;
                    polygon(x=BR[c(i-1,i-1,i,i)], y=1+c(0,htau$counts[i-1],htau$counts[i-1],0)/max(htau$counts)*0.95,col=rgb(1,1,1,.5),border="white")
                    polygon(x=BR[c(i-1,i-1,i,i)], y=1+c(0,htau$counts[i-1],htau$counts[i-1],0)/max(htau$counts)*0.95,col=rgb(r1[1],r1[2],r1[3],0.5),border="white")

                    pci = pcol[2]
                    r2 = col2rgb(pci)/255;
                    polygon(x=BR[c(i-1,i-1,i,i)], y=2+c(0,hrho$counts[i-1],hrho$counts[i-1],0)/max(hrho$counts)*0.95,col=rgb(1,1,1,.5),border="white")
                    polygon(x=BR[c(i-1,i-1,i,i)], y=2+c(0,hrho$counts[i-1],hrho$counts[i-1],0)/max(hrho$counts)*0.95,col=rgb(r2[1],r2[2],r2[3],0.5),border="white")

                    pci = pcol[4]
                    r3 = col2rgb(pci)/255;
                    polygon(x=BR[c(i-1,i-1,i,i)], y=3+c(0,hphi$counts[i-1],hphi$counts[i-1],0)/max(hphi$counts)*0.95,col=rgb(1,1,1,.5),border="white")
                    polygon(x=BR[c(i-1,i-1,i,i)], y=3+c(0,hphi$counts[i-1],hphi$counts[i-1],0)/max(hphi$counts)*0.95,col=rgb(r3[1],r3[2],r3[3],0.5),border="white")
                }

                
                abline(h=1:3,col=rgb(0,0,0,0.2))
                yb = c(htau$counts[1]/max(htau$counts)*.95,hrho$counts[1]/max(hrho$counts)*.95,hphi$counts[1]/max(hphi$counts)*.95)
                text(x=-0.02,y=1:3+yb,las=1,labels=c(expression(tau),expression(rho),expression(phi)), adj=0.5,offset=0.5,col=c(rgb(r1[1],r1[2],r1[3],0.7),rgb(r2[1],r2[2],r2[3],0.7),rgb(r3[1],r3[2],r3[3],0.7)))

                xf=c(0.25,0.3,0.35)
                xb = c(which.min(abs(dtau$x-xf[1])), which.min(abs(drho$x-xf[2])),which.min(abs(dphi$x-xf[3])))
                yb = c(dtau$y[xb[1]]/max(dtau$y)*0.9, drho$y[xb[2]]/max(drho$y)*0.9, dphi$y[xb[3]]/max(dphi$y)*0.8)
                text(x=xf,y=1:3+yb,las=1,labels=c(expression(bar(tau)),expression(bar(rho)),expression(bar(phi))), pos=3,offset=0.5,col=rgb(0,0,0,0.5))
                dev.off()
            }

            {
                pdf(paste(FIGDIR,"/Main/3B_joint_rho_tau.pdf",sep=""),height=1.5,width=w2col/3,pointsize=10)
                par(mar=c(2.15,2,.5,.5))
                image(Bx[1:BREAKS], By[1:BREAKS], log(JOINT),xlim=c(0,1),ylim=c(0,1),
                      col = viridis(BREAKS, option = "C"),
                      xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
                abline(0,1,col=rgb(0,0,0,0.1))

                axis(side=1,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
                axis(side=1,padj=-2,at=c(0,0.5,1.0),cex.axis=0.7)
                axis(side=2,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
                axis(side=2,padj=1.3,at=c(0,0.5,1.0), cex.axis=0.7)
                mtext(side=1,expression(bar(rho)),line=1.1)
                mtext(side=2,expression(bar(tau)),line=1)
                box()
                dev.off()

                pdf(paste(FIGDIR,"/Main/3B_joint_leg.pdf",sep=""),height=1,width=0.4,pointsize=10)
                par(mar=c(1,.6,1,.6),oma=c(0,0,0,0))
                tcol = rev(viridis(5,option="C"))
                image(t(matrix(seq(1,max(JOINT),length=BREAKS),nr=BREAKS,nc=1)),
                      col = viridis(BREAKS, option = "C"),
                      xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
                ticks = c(1,2,3,5,10)
                text(x=0,y=0,labels=ticks[1],cex=0.7,col="black",adj=c(0.5,0))
                text(x=0,y=seq(0,1,length=5)[2:4],labels=ticks[2:4],cex=0.7,col="black",adj=c(0.5,0.5))
                text(x=0,y=1,labels=ticks[5],cex=0.7,col="black",adj=c(0.5,1))
                mtext(side=3,"Count",cex=0.7)
                dev.off()
            }
            {
                ##plotting the out of sample simulations

                library(scales)
                par(mfrow=c(1,1))

                read.table(paste(FIGDIR,"LawrenceData/resamp_sims.txt",sep=""))->out_samp
                read.table(paste(FIGDIR,"LawrenceData/transformed.real.data",sep=""))->real
                LMAF=c(1,2,3,5,10,20,60,120,180,240,360)/720;

                pdf(paste(FIGDIR,"/Main/3C_outofSampleDraws.pdf",sep=""),height=1.5,width=w2col/3,pointsize=10)
                par(mar=c(2.15,2.5,.5,.5))
                plot(LMAF,t(out_samp[1,7:17]),col=alpha('darkgray',0.02),log='x',type='l',xlab='',ylab='',xlim=c(0.001,0.5),ylim=c(0,1),xaxt='n', yaxt='n')
                axis(side=1,padj=-1.75,at=c(0.001,0.005,0.02,0.1,0.5),labels=c("0.001","0.005","0.02","0.1","0.5"),cex.axis=0.7)
                axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.1,0.9,by=0.1),labels=NA,tck=-0.04)
                axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,1,by=0.5),labels=c("0","0.5","1"))
                abline(h=seq(0,1,by=0.1),col=rgb(0,0,0,0.1),lwd=0.75)
                mtext(side=1,"Minor Allele Frequency",line=1)
                mtext(side=2,expression(paste("Cum. Prop. ",{widehat(h^2)},sep="")),line=.9)
                legend("topleft",c("Data","Posterior Draws","Neutral Model"),col=c("black","darkgray",pcol[8]),lty=c(1,1,2),pch=c(20,-1,-1),bg=rgb(1,1,1,0.8),box.col=rgb(1,1,1,0.8), cex=0.8)
                box()
                lines(LMAF,LMAF*2, col=pcol[8],lty=2)
                for (i in seq(2,length(out_samp$V1))) {
                    lines(LMAF,cumsum(t(out_samp[i,7:17])),col=alpha('darkgray',0.1))
                }
                lines(LMAF,cumsum(t(real)),lty=1,pch=20,type='o')
                dev.off()

            }
        }
        
        ## histogram of inferred values
        read.table(paste(FIGDIR,"LawrenceData/real_data.posterior.txt",sep=""))->real_infer
        drho = density(real_infer$V1,adjust=.5);
        dtau = density(real_infer$V2,adjust=.5);
        dphi = density(real_infer$V3,adjust=.5);
        BR = seq(0,1,by=0.02)
        hrho = hist(real_infer$V1,breaks=BR,plot=FALSE)
        htau = hist(real_infer$V2,breaks=BR,plot=FALSE)
        hphi = hist(real_infer$V3,breaks=BR,plot=FALSE)
        pci = pcol[1]
        ri = col2rgb(pci)/255;
        {
            pdf(paste(FIGDIR,"/NotUsed/3D_ridgeline_posterior.pdf",sep=""),height=2,width=2.375,pointsize=10)
            par(mar=c(1,1.5,0.5,0.4))
            plot(1,1,type='n',xlim=c(0,1),ylim=c(1,4),xaxt='n',yaxt='n',xlab="",ylab="",bty='n')
            mtext(side=2,"Posterior Distribution",line=0.5)
            axis(side=1,padj=-1.5,cex=0.8)
            pci = pcol[1]
            ri = col2rgb(pci)/255;
            polygon(x=c(dtau$x,min(dtau$x)),y=1+c((dtau$y/max(dtau$y)),0)*0.95,col=rgb(ri[1],ri[2],ri[3],0.1),border="white")
            pci = pcol[2]
            ri = col2rgb(pci)/255;
            polygon(x=c(drho$x,min(drho$x)),y=2+c((drho$y/max(drho$y)),0)*0.95,col=rgb(ri[1],ri[2],ri[3],0.1),border="white")
            pci = pcol[4]
            ri = col2rgb(pci)/255;
            polygon(x=c(dphi$x,min(dphi$x)),y=3+c((dphi$y/max(dphi$y)),0)*0.95,col=rgb(ri[1],ri[2],ri[3],0.1),border="white")
            for(i in 2:length(BR)){
                pci = pcol[1]
                ri = col2rgb(pci)/255;
                polygon(x=BR[c(i-1,i-1,i,i)], y=1+c(0,htau$counts[i-1],htau$counts[i-1],0)/max(htau$counts)*0.95,col=rgb(ri[1],ri[2],ri[3],0.2),border="white")

                pci = pcol[2]
                ri = col2rgb(pci)/255;
                polygon(x=BR[c(i-1,i-1,i,i)], y=2+c(0,hrho$counts[i-1],hrho$counts[i-1],0)/max(hrho$counts)*0.95,col=rgb(ri[1],ri[2],ri[3],0.2),border="white")

                pci = pcol[4]
                ri = col2rgb(pci)/255;
                polygon(x=BR[c(i-1,i-1,i,i)], y=3+c(0,hphi$counts[i-1],hphi$counts[i-1],0)/max(hphi$counts)*0.95,col=rgb(ri[1],ri[2],ri[3],0.2),border="white")
            }
            lines(dtau$x,1+(dtau$y/max(dtau$y))*0.95,col=rgb(0,0,0,0.5))
            lines(drho$x,2+(drho$y/max(drho$y))*0.95,col=rgb(0,0,0,0.5))
            lines(dphi$x,3+(dphi$y/max(dphi$y))*0.95,col=rgb(0,0,0,0.5))
            abline(h=1:3,col=rgb(0,0,0,0.2))
            text(x=-0.02,y=1:3+.3,las=1,labels=c(expression(tau),expression(rho),expression(phi)), pos=3,offset=0.5)
            dev.off()
            
        }
        if(0){
            pdf(paste(FIGDIR,"/NotUsed/3D-F.pdf",sep=""),height=1.5,width=4.75,pointsize=10)
            par(mfrow=c(1,3),mar=c(2.5,1,.5,.5),oma=c(0,2,0,0))
            hist(real_infer$V1,breaks=BR,col='gray',xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(0,1),main="",ylim=c(0,300))
            axis(side=1,padj=-2,cex.axis=0.7)
            axis(side=2,padj=1.3, cex.axis=0.7)
            mtext(side=1,expression(rho),line=1.5)
            mtext(side=2,"Posterior Draws",line=1.5)
            
            hist(real_infer$V2,breaks=BR,col='gray',xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(0,1),main="",ylim=c(0,450))
            axis(side=1,padj=-2,cex.axis=0.7)
            axis(side=2,padj=1.3, cex.axis=0.7)
            mtext(side=1,expression(tau),line=1.5)

            hist(real_infer$V3,breaks=BR,col='gray',xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(0,1),main="",ylim=c(0,150))
            axis(side=1,padj=-2,cex.axis=0.7)
            axis(side=2,padj=1.3, cex.axis=0.7)
            mtext(side=1,expression(phi),line=1.5)
            dev.off()
        }
        ##plotting the out of sample simulations

        library(scales)
        par(mfrow=c(1,1))

        read.table(paste(FIGDIR,"LawrenceData/resamp_sims.txt",sep=""))->out_samp
        read.table(paste(FIGDIR,"LawrenceData/transformed.real.data",sep=""))->real
        LMAF=c(1,2,3,5,10,20,60,120,180,240,360)/720;

        pdf(paste(FIGDIR,"/NotUsed/3E_outofSampleDraws.pdf",sep=""),height=2,width=2.375,pointsize=10)
        par(mar=c(2,2.5,.5,.5))
        plot(LMAF,t(out_samp[1,7:17]),col=alpha('darkgray',0.02),log='x',type='l',xlab='',ylab='',xlim=c(0.001,0.5),ylim=c(0,1),xaxt='n', yaxt='n')
        axis(side=1,padj=-1.75,at=c(0.001,0.005,0.02,0.1,0.5),labels=c("0.001","0.005","0.02","0.1","0.5"),cex.axis=0.7)
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.1,0.9,by=0.1),labels=NA,tck=-0.04)
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,1,by=0.5),labels=c("0","0.5","1"))
        abline(h=seq(0,1,by=0.1),col=rgb(0,0,0,0.1),lwd=0.75)
        mtext(side=1,"Minor Allele Frequency",line=1)
        mtext(side=2,expression(paste("Cum. Prop. ",{widehat(h^2)},sep="")),line=.9)
        legend("topleft",c("Data","Model"),col=c("black","darkgray"),lty=1,pch=c(20,-1),bg=rgb(1,1,1,0.8),box.col=rgb(1,1,1,0.8), cex=0.8)
        box()
        for (i in seq(2,length(out_samp$V1))) {
            lines(LMAF,cumsum(t(out_samp[i,7:17])),col=alpha('darkgray',0.02))
        }
        lines(LMAF,cumsum(t(real)),lty=1,pch=20,type='o')
        dev.off()
    }
}


{
#### SUPPLEMENTAL FIGURES
    { ##Fig S2 and S3
        {
            ##complicated sim plot:  panel columns have increasing K, rows have increasing h2. Each panel has relative bias for total h2 as function of Nc. Adjacent points for each minMAF will be HE (filled) and REML (open) points.  Color indicates FR/FT (green/blu scale, square) or rho/tau (red/blue scale circle) combination. Different beta models will be offset from each other.
            
            K=c(1,2,5,10,20)
            h2=c(0.02,0.05,0.1,0.2,0.5)
            FT=c(0,0.01,0.05,0.1)
            FR=c(0.01,0.05,0.1,0.5,1)
            Nc=c(1,10,100,1000)
            BM=c(1,2)
            rho=c(0,0.5,0.8,0.9,0.95,1.0)
            tau=c(0.5,0.8,1.0)
            CM=c(1,5,8)
            minMAC=c(1,2,5,36);
            
            pdf(paste(FIGDIR,"/SUPP/S2_allSims_relbias_h2tot.pdf",sep=""),height=4.48,width=w1col,pointsize=10)
            par(mfrow=c(length(h2),length(K)),mar=c(.5,.75,.5,0),oma=c(1,2,1,2))
            panel = 0;
            totCombos=0;
            totSims=0;
            for(h in h2){ #rows
                for(k in K){ #columns

                    ##get data for each panel
                    DATA.he = matrix(,nc=4*k+3,nr=0);
                    DATA.reml = matrix(,nc=4*k+3,nr=0);
                    hedatapoints=0;
                    remldatapoints=0;
                    dp=0;
                    for(mm in minMAC){
                        for(cm in CM){
                            for(ft in FT){
                                for(fr in FR){
                                    for(nc in Nc){
                                        for(bm in BM){
                                            for(r in rho){
                                                for(t in tau){
                                                    INFILEhe=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_he.h2.gz",sep="");
                                                    INFILEreml=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_reml.h2.gz",sep="")
                                                    
                                                    if(file.exists(INFILEhe) | file.exists(INFILEreml)){
                                                        ##cat(INFILEhe,"=",file.exists(INFILEhe),"; ",INFILEreml,"=",file.exists(INFILEreml),"\n");
                                                        dp=dp+1;
                                                        DATA.he = rbind(DATA.he,rep(NA,4*k+3));
                                                        DATA.reml = rbind(DATA.reml,rep(NA,4*k+3));

                                                        if(file.exists(INFILEhe)){
                                                            con=gzfile(INFILEhe);
                                                            d=matrix(scan(file=con,quiet=T),byrow=T,nc=4*k+2);
                                                            if(nrow(d)>10){
                                                                DATA.he[dp,] = c(apply(d,2,function(x){mean(x,na.rm=T)}),nrow(d));
                                                                hedatapoints = hedatapoints+1;
                                                                totCombos=totCombos+1;
                                                                totSims=totSims + nrow(d);
                                                            }
                                                            close(con);
                                                        }
                                                        if(file.exists(INFILEreml)){
                                                            con=gzfile(INFILEreml);
                                                            d=matrix(scan(file=con,quiet=T),byrow=T,nc=4*k+2);
                                                            if(nrow(d)>10){
                                                                DATA.reml[dp,] = c(apply(d,2,function(x){mean(x,na.rm=T)}),nrow(d));
                                                                remldatapoints = remldatapoints+1;
                                                                totCombos=totCombos+1;
                                                                totSims=totSims + nrow(d);
                                                            }
                                                            close(con);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    cat("K=",k,"; h2=",h,": HE=",hedatapoints,"; REML=",remldatapoints,"\n",sep="");
                    panel=panel+1;
                    if(nrow(DATA.he) > 0 & nrow(DATA.reml) > 0){
                        rb.he=(DATA.he[,1]-h)/h
                        rb.reml=(DATA.reml[,1]-h)/h;
                        plot(rb.he,pch=16,xlab="",ylab="",xaxt='n',yaxt='n',ylim=c(-5,5))
                        points(rb.reml,col='blue')
                        abline(h=0,col='red')
                        if(panel <= length(K)){
                            mtext(side=3,paste("K=",k,sep=""));
                        }
                        if(panel %% length(K) == 0){
                            mtext(side=4,paste("h2=",h,sep=""),line=0.5);
                        }
                        if(panel %% length(K) == 1){
                            axis(side=2,padj=1)
                            mtext(side=2,"Rel Bias",line=1.2)
                        }
                        else{
                            axis(side=2,padj=1,labels=NA)

                        }
                        if(panel > (length(h2)-1)*length(K)){
                            mtext(side=1,"Sim",line=.5)
                        }
                    }
                    else{
                        plot(1,1,type='n',xlab="",ylab="",xaxt='n',yaxt='n')
                    }
                }
            }
            dev.off();
            cat("total param combos = ",totCombos,"\n",sep="");
            cat("total number sims  = ",totSims,"\n",sep="");
        }

        {
            ##complicated sim plot:  panel columns have increasing K, rows have increasing h2. Each panel has relative bias for total h2 as function of Nc. Adjacent points for each minMAF will be HE (filled) and REML (open) points.  Color indicates FR/FT (green/blu scale, square) or rho/tau (red/blue scale circle) combination. Different beta models will be offset from each other.
            
            K=c(1,2,5,10,20)
            h2=c(0.02,0.05,0.1,0.2,0.5)
            FT=c(0,0.01,0.05,0.1)
            FR=c(0.01,0.05,0.1,0.5,1)
            Nc=c(1,10,100,1000)
            BM=c(1,2)
            rho=c(0,0.5,0.8,0.9,0.95,1.0)
            tau=1.0#c(0.5,0.8,1.0)
            CM=c(1,5,8)
            minMAC=c(1,2,5,36);
            
            pdf(paste(FIGDIR,"/SUPP/S3_allSims_relbias_h2tot_nsims.pdf",sep=""),height=4.48,width=w1col,pointsize=10)
            par(mfrow=c(length(h2),length(K)),mar=c(.5,.75,.5,0),oma=c(1,2,1,2))
            panel = 0;
            totCombos=0;
            totSims=0;
            for(h in h2){ #rows
                for(k in K){ #columns
                    
                    ##get data for each panel
                    DATA.he = 0;
                    DATA.reml = 0
                    hedatapoints=0;
                    remldatapoints=0;
                    dp=0;
                    for(mm in minMAC){
                        for(cm in CM){
                            for(ft in FT){
                                for(fr in FR){
                                    for(nc in Nc){
                                        for(bm in BM){
                                            for(r in rho){
                                                for(t in tau){
                                                    INFILEhe=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_he.h2.gz",sep="");
                                                    INFILEreml=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_reml.h2.gz",sep="")
                                                    
                                                    if(file.exists(INFILEhe) | file.exists(INFILEreml)){
                                                        ##cat(INFILEhe,"=",file.exists(INFILEhe),"; ",INFILEreml,"=",file.exists(INFILEreml),"\n");
                                                        dp=dp+1;
                                                        DATA.he[dp] = 0
                                                        DATA.reml[dp] = 0

                                                        if(file.exists(INFILEhe)){
                                                            con=gzfile(INFILEhe);
                                                            d=matrix(scan(file=con,quiet=T),byrow=T,nc=4*k+2);
                                                            DATA.he[dp] = nrow(d);
                                                            hedatapoints = hedatapoints+1;
                                                            totCombos=totCombos+1;
                                                            totSims=totSims + nrow(d);
                                                            close(con);
                                                        }
                                                        if(file.exists(INFILEreml)){
                                                            con=gzfile(INFILEreml);
                                                            d=matrix(scan(file=con,quiet=T),byrow=T,nc=4*k+2);
                                                            DATA.reml[dp] = nrow(d);
                                                            remldatapoints = remldatapoints+1;
                                                            totCombos=totCombos+1;
                                                            totSims=totSims + nrow(d);
                                                            close(con);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    cat("K=",k,"; h2=",h,": HE=",hedatapoints,"; REML=",remldatapoints,"\n",sep="");
                    panel=panel+1;
                    if(length(DATA.he) > 0 & length(DATA.reml) > 0){
                        rb.he=DATA.he
                        rb.reml=DATA.reml
                        plot(rb.he,pch=16,xlab="",ylab="",xaxt='n',yaxt='n',ylim=c(0,500))
                        points(rb.reml,col='blue')
                        abline(h=0,col='red')
                        if(panel <= length(K)){
                            mtext(side=3,paste("K=",k,sep=""));
                        }
                        if(panel %% length(K) == 0){
                            mtext(side=4,paste("h2=",h,sep=""),line=0.5);
                        }
                        if(panel %% length(K) == 1){
                            axis(side=2,padj=1)
                            mtext(side=2,"Num Conv",line=1.2)
                        }
                        else{
                            axis(side=2,padj=1,labels=NA)
                        }
                        if(panel > (length(h2)-1)*length(K)){
                            mtext(side=1,"Sim",line=.5)
                        }
                    }
                    else{
                        plot(1,1,type='n',xlab="",ylab="",xaxt='n',yaxt='n')
                    }
                }
            }
            dev.off();
            cat("total param combos = ",totCombos,"\n",sep="");
            cat("total number sims  = ",totSims,"\n",sep="");
        }

        for(MAFRARE in c(0.01,0.025,0.05)){ #frequency threshold for rare
        {  ##REPEATED FOR h2_rare only!!
            ##complicated sim plot:  panel columns have increasing K, rows have increasing h2. Each panel has relative bias for total h2 as function of Nc. Adjacent points for each minMAF will be HE (filled) and REML (open) points.  Color indicates FR/FT (green/blu scale, square) or rho/tau (red/blue scale circle) combination. Different beta models will be offset from each other.
            K=c(1,2,5,10,20)
            h2=c(0.02,0.05,0.1,0.2,0.5)
            FT=c(0,0.01,0.05,0.1)
            FR=c(0.01,0.05,0.1,0.5,1)
            Nc=c(1,10,100,1000)
            BM=c(1,2)
            rho=c(0,0.5,0.8,0.9,0.95,1.0)
            tau=c(0.5,0.8,1.0)
            CM=c(1,5,8)
            minMAC=c(1,2,5,36);
            
            pdf(paste(FIGDIR,"/SUPP/allSims_relbias_h2rare_",MAFRARE,".pdf",sep=""),height=11,width=8.5,pointsize=10)
            par(mfrow=c(length(h2),length(K)),mar=c(1,1.5,.5,.5),oma=c(1,1,1,1.5))
            panel = 0;
            totCombos=0;
            totSims=0;
            for(h in h2){ #rows
                for(k in K){ #columns
                    
                    ##get data for each panel
                    DATA.he = matrix(,nc=3,nr=0);
                    DATA.reml = matrix(,nc=3,nr=0);
                    hedatapoints=0;
                    remldatapoints=0;
                    dp=0;
                    for(mm in minMAC){
                        for(cm in CM){
                            for(ft in FT){
                                for(fr in FR){
                                    for(nc in Nc){
                                        for(bm in BM){
                                            for(r in rho){
                                                for(t in tau){
                                                    INFILEhe=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_he.h2.gz",sep="");
                                                    INFILEreml=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_reml.h2.gz",sep="")
                                                    
                                                    if(file.exists(INFILEhe) | file.exists(INFILEreml)){
                                                        ##cat(INFILEhe,"=",file.exists(INFILEhe),"; ",INFILEreml,"=",file.exists(INFILEreml),"\n");
                                                        dp=dp+1;
                                                        DATA.he = rbind(DATA.he,rep(NA,3));
                                                        DATA.reml = rbind(DATA.reml,rep(NA,3));

                                                        if(file.exists(INFILEhe)){
                                                            con=gzfile(INFILEhe);
                                                            d=matrix(as.numeric(gsub("NaN",0,scan(file=con,quiet=T))),byrow=T,nc=4*k+2);
                                                            if(nrow(d)>10){
                                                                if(sum(d[1,k+1+1:k]<=MAFRARE, na.rm=T) > 0){
                                                                    DATA.he[dp,] = c(mean(apply(d[,2:ncol(d)],1,function(x){sum(x[which(x[k+1:k]<=MAFRARE)],na.rm=T)}),na.rm=T),mean(apply(d[,2:ncol(d)],1,function(x){sum(x[3*k+1+which(x[k+1:k]<=MAFRARE)],na.rm=T)}),na.rm=T),nrow(d));
                                                                    hedatapoints = hedatapoints+1;
                                                                    totCombos=totCombos+1;
                                                                    totSims=totSims + nrow(d);
                                                                }
                                                                else{
                                                                    DATA.he[dp,] = c(NA,NA,0);
                                                                }
                                                            }
                                                            close(con);
                                                        }
                                                        if(file.exists(INFILEreml)){
                                                            con=gzfile(INFILEreml);
                                                            d=matrix(as.numeric(gsub("NaN",0,scan(file=con,quiet=T))),byrow=T,nc=4*k+2);
                                                            if(nrow(d)>10){
                                                                if(sum(d[1,k+1+1:k]<=MAFRARE, na.rm=T) > 0){
                                                                    DATA.reml[dp,] = c(mean(apply(d[,2:ncol(d)],1,function(x){sum(x[which(x[k+1:k]<=MAFRARE)],na.rm=T)}),na.rm=T),mean(apply(d[,2:ncol(d)],1,function(x){sum(x[3*k+1+which(x[k+1:k]<=MAFRARE)],na.rm=T)}),na.rm=T),nrow(d));
                                                                    remldatapoints = remldatapoints+1;
                                                                    totCombos=totCombos+1;
                                                                    totSims=totSims + nrow(d);
                                                                }
                                                                else{
                                                                    DATA.reml[dp,] = c(NA,NA,0);
                                                                }
                                                            }
                                                            close(con);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    cat("K=",k,"; h2=",h,": HE=",hedatapoints,"; REML=",remldatapoints,"\n",sep="");
                    panel=panel+1;
                    if(nrow(DATA.he) > 0 & nrow(DATA.reml) > 0){
                        hTRUE = DATA.he[,2];
                        maxBias = min(1, 5*h,na.rm=T);
                        rb.he=(DATA.he[,1]-hTRUE)
                        for(i in which(rb.he > maxBias)){
                            rb.he[i] = maxBias;
                        }
                        for(i in which(rb.he < -maxBias)){
                            rb.he[i] = -maxBias;
                        }
                        rb.reml=(DATA.reml[,1]-hTRUE);
                        for(i in which(rb.reml > maxBias)){
                            rb.reml[i] = maxBias;
                        }
                        for(i in which(rb.reml < -maxBias)){
                            rb.reml[i] = -maxBias;
                        }
                        plot(rb.he,pch=16,xlab="",ylab="",xaxt='n',yaxt='n',ylim=c(-maxBias,maxBias))
                        points(rb.reml,col='blue')
                        axis(side=2,padj=1)
                        abline(h=0,col='red')
                        if(panel <= length(K)){
                            mtext(side=3,paste("K=",k,sep=""));
                        }
                        if(panel %% length(K) == 0){
                            mtext(side=4,paste("h2=",h,sep=""),line=0.5);
                        }
                        if(panel %% length(K) == 1){
                            mtext(side=2,"Bias (Inf-True)",line=1)
                        }
                        if(panel > (length(h2)-1)*length(K)){
                            mtext(side=1,"Sim",line=1)
                        }
                    }
                    else{
                        plot(1,1,type='n',xlab="",ylab="",xaxt='n',yaxt='n')
                    }
                }
            }
            dev.off();
            cat("total param combos = ",totCombos,"\n",sep="");
            cat("total number sims  = ",totSims,"\n",sep="");
        }

        {  ##REPEATED FOR h2_common only!!
            ##complicated sim plot:  panel columns have increasing K, rows have increasing h2. Each panel has relative bias for total h2 as function of Nc. Adjacent points for each minMAF will be HE (filled) and REML (open) points.  Color indicates FR/FT (green/blu scale, square) or rho/tau (red/blue scale circle) combination. Different beta models will be offset from each other.
            K=c(1,2,5,10,20)
            h2=c(0.02,0.05,0.1,0.2,0.5)
            FT=c(0,0.01,0.05,0.1)
            FR=c(0.01,0.05,0.1,0.5,1)
            Nc=c(1,10,100,1000)
            BM=c(1,2)
            rho=c(0,0.5,0.8,0.9,0.95,1.0)
            tau=c(0.5,0.8,1.0)
            CM=c(1,5,8)
            minMAC=c(1,2,5,36);
            
            pdf(paste(FIGDIR,"/SUPP/allSims_relbias_h2common_",MAFRARE,".pdf",sep=""),height=11,width=8.5,pointsize=10)
            par(mfrow=c(length(h2),length(K)),mar=c(1,1.5,.5,.5),oma=c(1,1,1,1.5))
            panel = 0;
            totCombos=0;
            totSims=0;
            for(h in h2){ #rows
                for(k in K){ #columns
                    
                    ##get data for each panel
                    DATA.he = matrix(,nc=3,nr=0);
                    DATA.reml = matrix(,nc=3,nr=0);
                    hedatapoints=0;
                    remldatapoints=0;
                    dp=0;
                    for(mm in minMAC){
                        for(cm in CM){
                            for(ft in FT){
                                for(fr in FR){
                                    for(nc in Nc){
                                        for(bm in BM){
                                            for(r in rho){
                                                for(t in tau){
                                                    INFILEhe=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_he.h2.gz",sep="");
                                                    INFILEreml=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/OUTSIM_HEvREML_SUM/merge_K=",k,"_MAC=",mm,"_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=",cm,"_BM=",bm,"_h2=",h,"_Nc=",nc,"_FR=",fr,"_FT=",ft,"_Kgen=1_Krem=1_errMod=1_errAdd=0_errRem=0_rho=",r,"_tau=",t,"_reml.h2.gz",sep="")
                                                    
                                                    if(file.exists(INFILEhe) | file.exists(INFILEreml)){
                                                        ##cat(INFILEhe,"=",file.exists(INFILEhe),"; ",INFILEreml,"=",file.exists(INFILEreml),"\n");
                                                        dp=dp+1;
                                                        DATA.he = rbind(DATA.he,rep(NA,3));
                                                        DATA.reml = rbind(DATA.reml,rep(NA,3));

                                                        if(file.exists(INFILEhe)){
                                                            con=gzfile(INFILEhe);
                                                            d=matrix(as.numeric(gsub("NaN",0,scan(file=con,quiet=T))),byrow=T,nc=4*k+2);
                                                            if(nrow(d)>10){
                                                                if(sum(d[1,k+1+1:k]>MAFRARE, na.rm=T) > 0){
                                                                    DATA.he[dp,] = c(mean(apply(d[,2:ncol(d)],1,function(x){sum(x[which(x[k+1:k]>MAFRARE)],na.rm=T)}),na.rm=T),mean(apply(d[,2:ncol(d)],1,function(x){sum(x[3*k+1+which(x[k+1:k]>MAFRARE)],na.rm=T)}),na.rm=T),nrow(d));
                                                                    hedatapoints = hedatapoints+1;
                                                                    totCombos=totCombos+1;
                                                                    totSims=totSims + nrow(d);
                                                                }
                                                                else{
                                                                    DATA.he[dp,] = c(NA,NA,0);
                                                                }
                                                            }
                                                            close(con);
                                                        }
                                                        if(file.exists(INFILEreml)){
                                                            con=gzfile(INFILEreml);
                                                            d=matrix(as.numeric(gsub("NaN",0,scan(file=con,quiet=T))),byrow=T,nc=4*k+2);
                                                            if(nrow(d)>10){
                                                                if(sum(d[1,k+1+1:k]>MAFRARE, na.rm=T) > 0){
                                                                    DATA.reml[dp,] = c(mean(apply(d[,2:ncol(d)],1,function(x){sum(x[which(x[k+1:k]>MAFRARE)],na.rm=T)}),na.rm=T),mean(apply(d[,2:ncol(d)],1,function(x){sum(x[3*k+1+which(x[k+1:k]>MAFRARE)],na.rm=T)}),na.rm=T),nrow(d));
                                                                    remldatapoints = remldatapoints+1;
                                                                    totCombos=totCombos+1;
                                                                    totSims=totSims + nrow(d);
                                                                }
                                                                else{
                                                                    DATA.reml[dp,] = c(NA,NA,0);
                                                                }
                                                            }
                                                            close(con);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    cat("K=",k,"; h2=",h,": HE=",hedatapoints,"; REML=",remldatapoints,"\n",sep="");
                    panel=panel+1;
                    if(nrow(DATA.he) > 0 & nrow(DATA.reml) > 0){
                        hTRUE = DATA.he[,2];
                        maxBias = min(1,5*h)
                        rb.he=(DATA.he[,1]-hTRUE)
                        rb.reml=(DATA.reml[,1]-hTRUE);
                        for(i in which(rb.he > maxBias)){
                            rb.he[i] = maxBias;
                        }
                        for(i in which(rb.he < -maxBias)){
                            rb.he[i] = -maxBias;
                        }
                        for(i in which(rb.reml > maxBias)){
                            rb.reml[i] = maxBias;
                        }
                        for(i in which(rb.reml < -maxBias)){
                            rb.reml[i] = -maxBias;
                        }
                        plot(rb.he,pch=16,xlab="",ylab="",xaxt='n',yaxt='n',ylim=c(-maxBias,maxBias))
                        points(rb.reml,col='blue')
                        axis(side=2,padj=1)
                        abline(h=0,col='red')
                        if(panel <= length(K)){
                            mtext(side=3,paste("K=",k,sep=""));
                        }
                        if(panel %% length(K) == 0){
                            mtext(side=4,paste("h2=",h,sep=""),line=0.5);
                        }
                        if(panel %% length(K) == 1){
                            mtext(side=2,"Bias (Inf-True)",line=1)
                        }
                        if(panel > (length(h2)-1)*length(K)){
                            mtext(side=1,"Sim",line=1)
                        }
                    }
                    else{
                        plot(1,1,type='n',xlab="",ylab="",xaxt='n',yaxt='n')
                    }
                }
            }
            dev.off();
            cat("total param combos = ",totCombos,"\n",sep="");
            cat("total number sims  = ",totSims,"\n",sep="");
        }
        }
    }

    { # Fig S4: plot noah sims
        d=read.table("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures5/SUPP/NoahSimData.txt");
        pdf(paste(FIGDIR,"/SUPP/S4_Noah_fracCausalSims.pdf",sep=""),height=1.5,width=2.375,pointsize=10)
        par(mar=c(2,2.8,0.3,0.3))
        plot(1,1,type='n',xlab="",ylab="",xaxt='n',yaxt='n',main="",xlim=c(1,4),ylim=c(0,0.8));
        h2r = c(0.1,0.25,0.5,0.75);
        fcr = c(0.001,0.01,0.05,0.5);
        axis(side=1,at=1:4,labels=fcr,padj=-2, cex.axis=0.7)
        axis(side=2,at=h2r,padj=1.3, cex.axis=0.7)
        mtext(side=1,"Frac. Causal",line=1)
        mtext(side=2,expression({widehat(h^2)}[total]),line=1.15)
        abline(h=h2r,col=rgb(0,0,0,0.1),lwd=0.75)
        cnt = 0;
        for(h2 in 1:length(h2r)){
            for(fc in 1:length(fcr)){
                cnt=cnt+1;
                points(fc, d[cnt,3],pch="-")
                segments(x0=fc,y0=d[cnt,3]-2*d[cnt,4], y1=d[cnt,3]+2*d[cnt,4])
            }
        }
        dev.off()

    }

    { #Figure S5:  simulation figure with adding genotyping errors

        K=20;
        errAdd = c(0,0.01,0.03,0.1,0.3,0.9);
        errRem = c(0,0.01,0.03,0.1,0.3,0.9);

        QN=1;
        h2tot1 = matrix(,nr=length(errAdd),nc=length(errRem));
        for(i1 in 1:length(errAdd)){
            for(i2 in 1:length(errRem)){
                ea = errAdd[i1];
                er = errRem[i2];
                filename=paste("merge_1000000_K=20_MAC=1_ngPC=0_npPC=0_QN=1_DSAMP=1_SIM=1_CM=8_BM=0_h2=0.2_Nc=100_FR=1_FT=0_delta=1_Kgen=1_Krem=1_errMod=3_errAdd=",ea,"_errRem=",er,"_rho=0.9_tau=0.5_tenBoyko.he.h2.gz",sep="");
                GZIN=gzfile(filename,"r");
                foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                close(GZIN);
                GZIN=gzfile(filename,"r");
                he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                close(GZIN);
                h2tot1[i1,i2] = sum(apply(he[,2+1:K],2,function(x){mean(as.numeric(x),na.rm=T)}));
            }
        }

        QN=0;
        h2tot0 = matrix(,nr=length(errAdd),nc=length(errRem));
        h2se0 = matrix(,nr=length(errAdd),nc=length(errRem));
        for(i1 in 1:length(errAdd)){
            for(i2 in 1:length(errRem)){
                ea = errAdd[i1];
                er = errRem[i2];
                filename=paste("merge_1000000_K=20_MAC=1_ngPC=0_npPC=0_QN=",QN,"_DSAMP=1_SIM=1_CM=8_BM=0_h2=0.2_Nc=100_FR=1_FT=0_delta=1_Kgen=1_Krem=1_errMod=3_errAdd=",ea,"_errRem=",er,"_rho=0.9_tau=0.5_tenBoyko.he.h2.gz",sep="");
                GZIN=gzfile(filename,"r");
                foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                close(GZIN);
                GZIN=gzfile(filename,"r");
                he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                close(GZIN);
                h2tot0[i1,i2] = sum(apply(he[,2+1:K],2,function(x){mean(as.numeric(x),na.rm=T)}));
                h2se0[i1,i2] = sum(apply(he[,2+1:K],2,function(x){sd(as.numeric(x),na.rm=T)/sqrt(length(x))}));
            }
        }

        {
            pdf(paste(FIGDIR,"/SUPP/S5_errSims.pdf",sep=""),height=1.5,width=2.375,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xaxt='n',yaxt='n',xlab="",ylab="",main="",xlim=c(0.5,6.5),ylim=c(0.10,0.23))
            abline(h=0.2,col='grey')
            ##abline(h=seq(0,0.24,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
            ##abline(v=errAdd,col=rgb(0,0,0,0.1),lwd=0.75)
            axis(side=1,at=1:6,padj=-2, cex.axis=0.7,labels=errRem)
            axis(side=2,padj=1.3, cex.axis=0.7)
            mtext(side=1,"Error rate (missing)",line=1)
            mtext(side=2,expression({widehat(h^2)}[total]),line=1.15)
            
            tcol=c("red","orange","yellow","green","blue","purple")
            legend("bottom",c("errAdd=0","errAdd=0.01","errAdd=0.03","errAdd=0.1","errAdd=0.3","errAdd=0.9"),col=tcol,pch=20, cex=0.7, bg=rgb(1,1,1,0.8),box.col=rgb(1,1,1,0.2),ncol=2)
            box()
            
            for(i in 1:length(errAdd)){
                segments(x0=1:6+(i-3.5)/10,y0=h2tot0[i,]-2*h2se0[i,],y1=h2tot0[i,]+2*h2se0[i,],col=tcol[i],pch=20);
                points(1:6+(i-3.5)/10,h2tot0[i,],col=tcol[i],pch=20);
            }
            dev.off()
        }
    }

    { # Figure S6 & S7
        MAC="1";
        K=20
        for(BM in c("CONST", "NORM")){
            for(QN in c(0,1)){
                FIGNUM="";
                if(BM=="NORM" & QN==0){
                    FIGNUM="S6_";
                }
                if(BM=="NORM" & QN==1){
                    FIGNUM="S7_";
                }
                pdf(paste(FIGDIR,"/SUPP/",FIGNUM,"h2_simSingCausal_by_h2_fracCaus_beta=",BM,"_QN=",QN,"_MAC=",MAC,".pdf",sep=""),height=4,width=6,pointsize=10);
                par(mfrow=c(4,4),mar=c(0.5,1.5,1,0),oma=c(1.5,2,1.5,2))
                CNT=0;
                for(h2 in c(0.1,0.25,0.5,0.75)){
                    for(fc in c(0.001,0.01,0.05,0.5)){
                        CNT=CNT+1;
                        ANAL1="ALLtgp"
                        ANAL2="ALLGEU"
                        SRC1="TGP"
                        SRC2="GEU"
                        heh2=matrix(,nr=4,nc=20);
                        if(1){
                            filename=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/INFW_SIM/OUTSUM/merge_1000000_K=",K,"_MAC=",MAC,"_ngPC=0_npPC=0_QN=",QN,"_DSAMP=1_",ANAL1,"_fC=",fc,"_h2=",h2,"_src=",SRC1,"_BM=",BM,".he.h2.gz",sep="");
                            GZIN=gzfile(filename,"r");
                            foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                            close(GZIN);
                            GZIN=gzfile(filename,"r");
                            he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                            close(GZIN);
                            d1tot=density(as.numeric(he[,2]));
                            d1sing=density(as.numeric(he[,3]));
                            d1tot$mean=mean(as.numeric(he[,2]));
                            d1tot$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                            d1sing$mean=mean(as.numeric(he[,2]));
                            d1sing$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                                        #b1tot= boxplot(as.numeric(he[,2]))
                                        #b1sing= boxplot(as.numeric(he[,3]))
                                        #heh2[1,] = apply(he[,1:20+2],2,function(x){mean(as.numeric(x),na.rm=T)})
                            
                            filename=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/INFW_SIM/OUTSUM/merge_1000000_K=",K,"_MAC=",MAC,"_ngPC=0_npPC=0_QN=",QN,"_DSAMP=1_",ANAL2,"_fC=",fc,"_h2=",h2,"_src=",SRC1,"_BM=",BM,".he.h2.gz",sep="");
                            GZIN=gzfile(filename,"r");
                            foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                            close(GZIN);
                            GZIN=gzfile(filename,"r");
                            he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                            close(GZIN);
                            d2tot=density(as.numeric(he[,2]));
                            d2sing=density(as.numeric(he[,3]));
                            d2tot$mean=mean(as.numeric(he[,2]));
                            d2tot$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                            d2sing$mean=mean(as.numeric(he[,2]));
                            d2sing$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                                        #b2tot= boxplot(as.numeric(he[,2]))
                                        #b2sing= boxplot(as.numeric(he[,3]))
                                        #heh2[2,] = apply(he[,1:20+2],2,function(x){mean(as.numeric(x),na.rm=T)})
                            
                            filename=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/INFW_SIM/OUTSUM/merge_1000000_K=",K,"_MAC=",MAC,"_ngPC=0_npPC=0_QN=",QN,"_DSAMP=1_",ANAL1,"_fC=",fc,"_h2=",h2,"_src=",SRC2,"_BM=",BM,".he.h2.gz",sep="");
                            GZIN=gzfile(filename,"r");
                            foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                            close(GZIN);
                            GZIN=gzfile(filename,"r");
                            he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                            close(GZIN);
                            d3tot=density(as.numeric(he[,2]));
                            d3sing=density(as.numeric(he[,3]));
                            d3tot$mean=mean(as.numeric(he[,2]));
                            d3tot$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                            d3sing$mean=mean(as.numeric(he[,2]));
                            d3sing$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                                        #b3tot= boxplot(as.numeric(he[,2]))
                                        #b3sing= boxplot(as.numeric(he[,3]))
                                        #heh2[3,] = apply(he[,1:20+2],2,function(x){mean(as.numeric(x),na.rm=T)})
                            
                            filename=paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/INFW_SIM/OUTSUM/merge_1000000_K=",K,"_MAC=",MAC,"_ngPC=0_npPC=0_QN=",QN,"_DSAMP=1_",ANAL2,"_fC=",fc,"_h2=",h2,"_src=",SRC2,"_BM=",BM,".he.h2.gz",sep="");
                            GZIN=gzfile(filename,"r");
                            foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                            close(GZIN);
                            GZIN=gzfile(filename,"r");
                            he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                            close(GZIN);
                            d4tot=density(as.numeric(he[,2]));
                            d4sing=density(as.numeric(he[,3]));
                            d4tot$mean=mean(as.numeric(he[,2]),na.rm=T);
                            d4tot$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                            d4sing$mean=mean(as.numeric(he[,2]));
                            d4sing$se=sd(as.numeric(he[,2]))/sqrt(nrow(he));
                                        #b4tot= boxplot(as.numeric(he[,2]),notch=T)
                                        #b4sing= boxplot(as.numeric(he[,3]),notch=T)
                                        #heh2[4,] = apply(he[,1:20+2],2,function(x){mean(as.numeric(x),na.rm=T)})
                        }
                        {
                            YLIM=range(0,d1tot$mean,d2tot$mean,d3tot$mean,d4tot$mean,h2);
                            if(CNT<=4){
                                YLIM=range(YLIM,0.14)
                            }
                            else if(CNT<=8){
                                YLIM=range(YLIM,0.4)
                            }
                            else if(CNT<=12){
                                YLIM=range(YLIM,0.8)
                            }
                            else{
                                YLIM=range(YLIM,1)
                            }
                            plot(1,1,type='n',xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(0.5,4.5),ylim=YLIM);
                            if(CNT<=4){
                                mtext(side=3,paste("fracCausal=",fc,sep=""),line=0.5)
                            }
                            if(CNT%%4 == 0){
                                mtext(side=4,paste("h2=",h2,sep=""),line=0.5)
                            }
                            if(CNT%%4 == 1){
                                mtext(side=2,expression({widehat(h^2)}[tot]),line=1.2)
                            }
                            if(CNT>=13){
                                mtext(side=1,"TGP",adj=0.2,cex=0.8,line=0.5)
                                mtext(side=1,"GEUVADIS",adj=.95,cex=0.8,line=0.5)
                            }
                            
                            axis(side=2,padj=1)
                            polygon(x=c(0,2.5,2.5,0),y=c(-1,-1,10,10),col=rgb(1,1,0,.1))
                            polygon(x=c(2.5,5,5,2.5),y=c(-1,-1,10,10),col=rgb(0,0,0,.1))
                            abline(h=0,col='grey')
                            abline(h=h2,col='black',lty=2)
                            points(c(1:4),c(d1tot$mean,d2tot$mean,d3tot$mean,d4tot$mean),col=pcol[1:2],pch=16)
                            segments(x0=1,y0=d1tot$mean-2*d1tot$se,y1=d1tot$mean+2*d1tot$se,col=pcol[1]);
                            segments(x0=2,y0=d2tot$mean-2*d2tot$se,y1=d2tot$mean+2*d2tot$se,col=pcol[2]);
                            segments(x0=3,y0=d3tot$mean-2*d3tot$se,y1=d3tot$mean+2*d3tot$se,col=pcol[1]);
                            segments(x0=4,y0=d4tot$mean-2*d4tot$se,y1=d4tot$mean+2*d4tot$se,col=pcol[2]);
                            ##legend("bottomright",c("TGP MAF partition","GEUVADIS MAF partition"),col=pcol[1:2],pch=20,bty='n')
                        }
                    }
                }
                dev.off()
            }
        }
    }

    { #investigate distribution of expression: Fig S8 and S9
        {
            ex = scan("../../matrix.eqtl/pooled.log2.pheno")
            x=seq(from=min(ex),to=max(ex),length=1000)
            cdf = 0;
            for(i in 1:length(x)){
                cdf[i] = mean(ex >= x[i]);
            }
            pdf(paste(FIGDIR,"/SUPP/S8_dist_expn_all_ind_all_genes.pdf",sep=""),height=3.5,width=4.75,pointsize=8)
            par(mfrow=c(3,1), mar=c(1,1,1,1), oma=c(2,2,0,0))
            foo=hist(ex,1000,plot=FALSE)
            plot(foo$mids, foo$counts, type='h', ,xlab="",ylab="",xaxt='n', yaxt='n',main="")
            axis(side=1,padj=-1,cex.axis=.8,at=0:15)
            axis(side=2,padj=1,cex.axis=.8)
            mtext(side=2,"Counts", line=1.5)
            plot(foo$mids, foo$counts, type='h', ,xlab="",ylab="",xaxt='n', yaxt='n',main="",log='y')
            axis(side=1,padj=-1,cex.axis=.8,at=0:15)
            axis(side=2,padj=1,cex.axis=.8)
            mtext(side=2,"Counts (log-scale)",line=1.5)

            plot(x,cdf,type='l',lwd=2,xaxt='n',yaxt='n')
            mtext(side=1,expression(log[2](FPKM)),line=2)
            mtext(side=2,"Fraction >= x", line=1.5)
            axis(side=1,padj=-1,cex.axis=.8,at=0:15)
            axis(side=2,padj=1,cex.axis=.8)

            for(val in seq(0,0.6,by=.1)){
                foo=which.min(abs(cdf-val))
                segments(-10,val,x[foo],val)
                segments(x[foo],val,x[foo],-1)
            }

            dev.off()
        }

        {
            F = scan("../../matrix.eqtl/files.log2.pheno", what="")
            emat = matrix(,nr=length(F),nc=360);
            i=0;
            for(f in F){
                i=i+1;
                emat[i,]=scan(paste("../../matrix.eqtl/log2/",f,sep=""),quiet=T)
            }
            Y=c(1,1.5,2,log2(5),3,4);
            X=seq(0,1,by=.1)
            P=matrix(,nr=length(X), nc=length(Y));
            for(i in 1:length(X)){
                for(j in 1:length(Y)){
                    P[i,j] = sum(apply(emat,1,function(y){mean(y>Y[j])}) > X[i]);
                }
            }
            ##P=
            ##        1   1.5     2  2.32     3    4
            ##0   15499 14182 13108 12463 11122 8783
            ##0.1 12625 11711 10867 10338  9067 6757
            ##0.2 12168 11295 10480  9923  8684 6349
            ##0.3 11865 11037 10201  9677  8413 6057
            ##0.4 11636 10775  9981  9466  8181 5802
            ##0.5 11386 10553  9764  9253  7967 5546
            ##0.6 11164 10366  9555  9045  7731 5328
            ##0.7 10922 10130  9342  8803  7454 5089
            ##0.8 10661  9885  9097  8528  7133 4779
            ##0.9 10295  9497  8717  8126  6709 4441

            ##Choose X=0.5 Y=1.5

            x=0.5;
            y=1.5;
            expGenes2 = array(,dim=0);
            foo=apply(emat,1,function(y){mean(y>1.5)});
            expGenes2=F[which(foo>=0.5)];
            expGenes2 = cbind(expGenes2, apply(emat[which(foo>=0.5),],1,mean));
            if(0){
                write(t(expGenes2), file="../GENESETS/goodExp.txt",ncolumns=2);
                for(x in X){
                    for(y in Y){
                        expGenes2 = array(,dim=0);
                        foo=apply(emat,1,function(z){mean(z>y)});
                        expGenes2=F[which(foo>=x)];
                        expGenes2 = cbind(expGenes2, apply(emat[which(foo>=x),],1,mean));
                        write(t(expGenes2), file=paste("../GENESETS/goodExp_X=",x,"_Y=",y,".txt",sep=""),ncolumns=2);
                    }
                }
            }
        }

        {
            W="1000000"
            K=20
            PC=10
            QN=1
            X=seq(0,1,by=.1)
            Y=c(0,1,1.5,2,log2(5),3,4);
            h2p = array(,dim=c(length(X), length(Y), K));
            Ns = array(,dim=c(length(X),length(Y)));
            
            for(PERMSTRING in c("","_PERM=3")){
                FILE=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1",PERMSTRING,"_NORESID.he.h2.gz",sep="");
                d = as.matrix(read.table(FILE)); #total dataset
                F = scan("../../matrix.eqtl/files.log2.pheno", what="")
                emat = matrix(,nr=length(F),nc=360);
                cnt=0;
                for(f in F){
                    cnt=cnt+1;
                    emat[cnt,]=scan(paste("../../matrix.eqtl/log2/",f,sep=""),quiet=T)
                    F[cnt]=strsplit(f, ".log2.pheno")[1];
                }
                
                My = 1;
                my=0;
                { #plot CDF(h2) for range of X% indivs with log2(expn)>Y
                    pdf(paste(FIGDIR,"/SUPP/S9_CDFh2vMAF_by_XY",PERMSTRING,".pdf",sep=""),height=5,width=4,pointsize=8)
                    par(mfrow=c(length(X),length(Y)), mar=c(.1,.1,.1,.1), oma=c(3,3.5,1.5,1.5))
                    
                    for(i in 1:length(X)){
                        for(j in 1:length(Y)){
                            if(1){
                                expGenes2 = array(,dim=0);
                                foo=apply(emat,1,function(y){mean(y>=Y[j])});
                                expGenes2=F[which(foo>=X[i])];
                                expGenes2 = cbind(expGenes2, apply(emat[which(foo>=X[i]),],1,mean));
                                
                                d2 = subset(d,d[,1] %in% expGenes2[,1]) #subset of genes that have "good" expression
                                d2m=apply(d2[,2:ncol(d2)],2,function(x){mean(as.numeric(x),na.rm=T)}); #mean of columns
                                d2c=cumsum(d2m[1+1:K]);
                                
                                if(max(d2c/d2c[K]) > My){
                                    My = max(d2c/d2c[K]);
                                }
                                if(min(d2c/d2c[K]) < my){
                                    my = min(d2c/d2c[K]);
                                }
                            }
                            plot(d2m[1+K+(2*PC)+1:K],d2c/d2c[K],ylim=c(my,My),xlim=c(0.001, 0.5),type='b',log='x',xaxt='n',yaxt='n',xlab="",ylab="",main="")
                            if(PERMSTRING == ""){
                                h2p[i,j,1:K] = d2c;
                                Ns[i,j] = nrow(d2);
                            }
                            else{
                                h2p[i,j,1:K] = h2p[i,j,1:K] - d2c;
                                h2p[i,j,1:K] = h2p[i,j,1:K]; #norm
                            }
                            if(i == 6 && j == 3){
                                polygon(x=c(0.00001,1,1,0.00001), y=c(-1,-1,2,2),col=rgb(1,1,0,0.2))
                            }
                            text(0.001,0.85, substitute(paste("N=",s1,sep=""), list(s1 = nrow(d2))), pos=4, offset=0.2)
                            
                            if(d2c[1]/d2c[K]<0.22){
                                text(0.001,0.65, substitute(paste(h^2,"=",s3,sep=""), list(s3 = format(d2c[K],digits=2))), pos=4, offset=0.2)
                            }
                            else{
                                text(0.5,0.25, substitute(paste(h^2,"=",s3,sep=""), list(s3 = format(d2c[K],digits=2))), pos=2, offset=0.2)
                            }
                            text(0.5,0.05, substitute(paste(psi,"=",s2,sep=""), list(s2 = format(d2c[1]/d2c[K],digits=2))), pos=2, offset=0.2)
                            
                            if(j==1){
                                axis(side=2,at=seq(0,1,by=0.25), padj=1, cex.axis=0.9,labels=NA)
                                if((i+1) %% 2 == 0){
                                    axis(side=2,at=c(0,0.5,1.0),padj=1, cex.axis=0.9)
                                }
                                if(i==6){
                                    mtext(side=2,expression(paste("CDF ",widehat(h^2),sep="")),line=1.5)
                                }
                            }
                            else{
                                axis(side=2,at=seq(0,1,by=.25), padj=1,labels=NA, cex.axis=0.8)
                                if(j==length(Y)){
                                    mtext(side=4,paste("X=",X[i],sep=""),line=0.22)
                                }
                            }
                            if(i == 1){
                                if(j == 5){
                                    mtext(side=3,expression(paste("Y=",log[2](5),sep="")),line=0);
                                }
                                else{
                                    mtext(side=3,paste("Y=",Y[j],sep=""),line=0);
                                }
                            }
                            if(i == length(X)){
                                axis(side=1,labels=NA, padj=-1, cex.axis=0.9)
                                if((j+1) %% 2 == 0){
                                    axis(side=1,at=c(0.001,0.5),labels=c("0.001","0.5"), padj=-1, cex.axis=0.9)
                                    axis(side=1,at=c(0.02),padj=-1, cex.axis=0.9)
                                }
                                mtext(side=1,"MAF",line=1.5)
                            }
                            else{
                                axis(side=1,padj=-1,labels=NA, cex.axis=0.8)
                            }
                        }
                        ## break;
                    }
                    dev.off()
                }

                { #plot perm-corrected CDF(h2)
                    pdf(paste(FIGDIR,"/SUPP/S9_permCorr_CDFh2vMAF_by_XY",PERMSTRING,".pdf",sep=""),height=5,width=4,pointsize=8)
                    par(mfrow=c(length(X),length(Y)), mar=c(.1,.1,.1,.1), oma=c(3,3.5,1.5,1.5))
                    
                    for(i in 1:length(X)){
                        for(j in 1:length(Y)){
                            plot(d2m[1+K+(2*PC)+1:K],h2p[i,j,]/h2p[i,j,K],ylim=range(0,1,h2p[i,j,]/h2p[i,j,K]),xlim=c(0.001, 0.5),type='b',log='x',xaxt='n',yaxt='n',xlab="",ylab="",main="")
                            if(i == 6 && j == 3){
                                polygon(x=c(0.00001,1,1,0.00001), y=c(-1,-1,2,2),col=rgb(1,1,0,0.2))
                            }
                            text(0.001,0.85, substitute(paste("N=",s1,sep=""), list(s1 = Ns[i,j])), pos=4, offset=0.2)
                            if(d2c[1]/d2c[K]<0.22){
                                text(0.001,0.65, substitute(paste(h^"2'","=",s3,sep=""), list(s3 = format(h2p[i,j,K],digits=2))), pos=4, offset=0.2)
                            }
                            else{
                                text(0.5,0.25, substitute(paste(h^"2'","=",s3,sep=""), list(s3 = format(h2p[i,j,K],digits=2))), pos=2, offset=0.2)
                            }
                            text(0.5,0.05, substitute(paste(psi,"=",s2,sep=""), list(s2 = format(h2p[i,j,1]/h2p[i,j,K],digits=2))), pos=2, offset=0.2)
                            
                            if(j==1){
                                axis(side=2,at=seq(0,1,by=0.25), padj=1, cex.axis=0.9,labels=NA)
                                if((i+1) %% 2 == 0){
                                    axis(side=2,at=c(0,0.5,1.0),padj=1, cex.axis=0.9)
                                }
                                if(i==6){
                                    mtext(side=2,expression(paste("CDF ",widehat(h^"2'"),sep="")),line=1.5)
                                }
                            }
                            else{
                                axis(side=2,at=seq(0,1,by=.25), padj=1,labels=NA, cex.axis=0.8)
                                if(j==length(Y)){
                                    mtext(side=4,paste("X=",X[i],sep=""),line=0.22)
                                }
                            }
                            if(i == 1){
                                if(j == 5){
                                    mtext(side=3,expression(paste("Y=",log[2](5),sep="")),line=0);
                                }
                                else{
                                    mtext(side=3,paste("Y=",Y[j],sep=""),line=0);
                                }
                            }
                            if(i == length(X)){
                                axis(side=1,labels=NA, padj=-1, cex.axis=0.9)
                                if((j+1) %% 2 == 0){
                                    axis(side=1,at=c(0.001,0.5),labels=c("0.001","0.5"), padj=-1, cex.axis=0.9)
                                    axis(side=1,at=c(0.02),padj=-1, cex.axis=0.9)
                                }
                                mtext(side=1,"MAF",line=1.5)
                            }
                            else{
                                axis(side=1,padj=-1,labels=NA, cex.axis=0.8)
                            }
                        }
                        ## break;
                    }
                    dev.off()
                }
            }
        }
    }

    { #Fig S10: h2 as a function of mean expression quantile (tot=rare+common)
        W="1000000"
        K=20;
        PC=10;
        QN=1;
        minMAC=1

        for(PERMSTRING in c("","_PERM=3")){
            ge = read.table("../../GENO/log2_exp_summary_gene.txt");
            GZIN=gzfile(paste("merge_",W,"_K=",K,"_MAC=",minMAC,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1",PERMSTRING,"_NORESID.he.h2.gz",sep=""),"r");
            he=scan(file=GZIN,nlines=1,what="",quiet=T)
            close(GZIN);
            GZIN=gzfile(paste("merge_",W,"_K=",K,"_MAC=",minMAC,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1",PERMSTRING,"_NORESID.he.h2.gz",sep=""),"r");
            he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(he))
            close(GZIN);
            G=which(he[,1] %in% expGenes[,1])
            he2 = he[G,]

            QT = c(1.0,0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
            GAP=0.1;
            smry=matrix(,nr=length(QT),nc=14);
            smry2=matrix(,nr=length(QT),nc=14);
            for(i in 1:length(QT)){
                q=QT[i];
                s = match(ge[which(ge[,2]<=quantile(ge[,2],q,na.rm=T) & ge[,2]>=quantile(ge[,2],q-GAP,na.rm=T)),1],he[,1])
                t = match(ge[which(ge[,2]<=quantile(ge[,2],q,na.rm=T)),1],he[,1])
                h2R = 0;
                h2C = 0;
                h2Rt = 0;
                h2Ct = 0;
                for(j in 1:length(s)){
                    if(!is.na(s[j])){
                        h2R[j] = sum(as.numeric(he[s[j],which(as.numeric(he[s[j],2*K+2+1:K]) <= 0.05)+2]))
                        h2C[j] = sum(as.numeric(he[s[j],which(as.numeric(he[s[j],2*K+2+1:K]) > 0.05)+2]))
                    }
                }
                
                smry[i,1]=QT[i]
                smry[i,2]=mean(as.numeric(he[s,2]),na.rm=T)
                smry[i,3]=mean(as.numeric(he[t,2]),na.rm=T)
                smry[i,4]=mean(h2R,na.rm=T)
                smry[i,5]=mean(h2C,na.rm=T)
                ##smry[i,6]=mean(h2Rt,na.rm=T)
                ##smry[i,7]=mean(h2Ct,na.rm=T)
                smry[i,8]=1.96*sd(as.numeric(he[s,52]),na.rm=T)/sqrt(length(he[s,52]))
                ##smry[i,9]=1.96*sd(as.numeric(he[t,52]),na.rm=T)/sqrt(length(he[t,52]))
                smry[i,10]=1.96*sd(h2R,na.rm=T)/sqrt(length(h2R))
                smry[i,11]=1.96*sd(h2C,na.rm=T)/sqrt(length(h2C))
                ##smry[i,12]=1.96*sd(h2Rt,na.rm=T)/sqrt(length(h2Rt))
                ##smry[i,13]=1.96*sd(h2Ct,na.rm=T)/sqrt(length(h2Ct))
                smry[i,14]=1.96*sd(h2R+h2C, na.rm=T)/sqrt(length(h2R+h2C))
                cat(smry[i,],"\n")
            }
            for(i in 1:length(QT)){
                
                q=QT[i];
                s = match(ge[which(ge[,2]<=quantile(ge[,2],q,na.rm=T) & ge[,2]>=quantile(ge[,2],q-GAP,na.rm=T)),1],he2[,1])
                t = match(ge[which(ge[,2]<=quantile(ge[,2],q,na.rm=T)),1],he2[,1])
                h2R = 0;
                h2C = 0;
                h2Rt = 0;
                h2Ct = 0;
                for(j in 1:length(s)){
                    if(!is.na(s[j])){
                        h2R[j] = sum(as.numeric(he2[s[j],which(as.numeric(he2[s[j],2*K+2+1:K]) <= 0.05)+2]),na.rm=T)
                        h2C[j] = sum(as.numeric(he2[s[j],which(as.numeric(he2[s[j],2*K+2+1:K]) > 0.05)+2]),na.rm=T)
                        ##h2Rt[j] = sum(as.numeric(he2[t[j],which(as.numeric(he2[s[j],K+2+1:K]) <= 0.01)+2]))
                        ##h2Ct[j] = sum(as.numeric(he2[t[j],which(as.numeric(he2[s[j],K+2+1:K]) > 0.01)+2]))
                    }
                }
                
                smry2[i,1]=QT[i]
                smry2[i,2]=mean(as.numeric(he2[s,2]),na.rm=T)
                smry2[i,3]=mean(as.numeric(he2[t,2]),na.rm=T)
                if(!is.na(smry2[i,2])){
                    smry2[i,4]=mean(h2R,na.rm=T)
                    smry2[i,5]=mean(h2C,na.rm=T)
                }
                ##smry2[i,6]=mean(h2Rt,na.rm=T)
                ##smry2[i,7]=mean(h2Ct,na.rm=T)
                smry2[i,8]=1.96*sd(as.numeric(he2[s,52]),na.rm=T)/sqrt(length(he2[s,52]))
                ##smry2[i,9]=1.96*sd(as.numeric(he2[t,52]),na.rm=T)/sqrt(length(he2[t,52]))
                smry2[i,10]=1.96*sd(h2R,na.rm=T)/sqrt(length(h2R))
                smry2[i,11]=1.96*sd(h2C,na.rm=T)/sqrt(length(h2C))
                ##smry2[i,12]=1.96*sd(h2Rt,na.rm=T)/sqrt(length(h2Rt))
                ##smry2[i,13]=1.96*sd(h2Ct,na.rm=T)/sqrt(length(h2Ct))
                smry2[i,14]=1.96*sd(h2R+h2C, na.rm=T)/sqrt(length(h2R+h2C))
                cat(smry2[i,],"\n")
            }
            ##smry=smry[1:9,]
            
            {
                pdf(paste(FIGDIR,"/SUPP/S10_h2_tot_fn_MeanExpnQuant_K=",K,"_PC=",PC,"_QN=",QN,PERMSTRING,".pdf",sep=""),height=2.5,width=7.25/3,pointsize=10)
                par(mar=c(2.6,2.7,0.5,0.5))

                c1=col2rgb(pcol[1])/255;
                c2=col2rgb(pcol[2])/255;
                
                plot(smry[,1],(smry[,4]+smry[,5]), xlim=range(smry[1:(nrow(smry)),1]),ylim=c(-.02,.13),xaxt='n',yaxt='n',xlab="",ylab="",main="",pch=16,type='n')
                abline(h=0,col='grey')
                axis(side=1,padj=-1,at=seq(0.1,0.9,by=0.2),labels=NA, cex.axis=0.8)
                axis(side=1,padj=-1, cex.axis=0.8)
                axis(side=2,padj=1, cex.axis=0.8)
                mtext(side=1,"Mean Expression Quantile",line=1.5)
                if(PERMSTRING==""){
                    mtext(side=2,expression(widehat(h^2)),line=1.2)
                }
                else{
                    mtext(side=2,expression(widehat(h^2)[perm]),line=1.2)
                }
                
                lines(smry[,1],(smry[,4]+smry[,5]),col=rgb(0,0,0,0.2))
                points(smry[,1],(smry[,4]+smry[,5]),pch=16,col=rgb(0,0,0,0.2))
                lines(smry[,1],(smry[,4]),pch=16,col=rgb(c1[1],c1[2],c1[3],0.2))
                lines(smry[,1],(smry[,5]), pch=16,col=rgb(c2[1],c2[2],c2[3],0.2))
                points(smry[,1],(smry[,4]),pch=16,col=rgb(c1[1],c1[2],c1[3],0.2))
                points(smry[,1],(smry[,5]), pch=16,col=rgb(c2[1],c2[2],c2[3],0.2))
                for(i in 1:nrow(smry)){
                    segments(x0=smry[i,1],x1=smry[i,1],y0=smry[i,4]+smry[i,5]-smry[i,14],y1=smry[i,4]+smry[i,5]+smry[i,14],col=rgb(0,0,0,0.2))
                    segments(x0=smry[i,1]-.02,x1=smry[i,1]+.02,y0=smry[i,4]+smry[i,5]-smry[i,14],y1=smry[i,4]+smry[i,5]-smry[i,14],col=rgb(0,0,0,0.2))
                    segments(x0=smry[i,1]-.02,x1=smry[i,1]+.02,y0=smry[i,4]+smry[i,5]+smry[i,14],y1=smry[i,4]+smry[i,5]+smry[i,14],col=rgb(0,0,0,0.2))
                    
                    segments(x0=smry[i,1],x1=smry[i,1],y0=smry[i,4]-smry[i,10],y1=smry[i,4]+smry[i,10], col=rgb(c1[1],c1[2],c1[3],0.2))
                    segments(x0=smry[i,1]-.02,x1=smry[i,1]+.02,y0=smry[i,4]-smry[i,10],y1=smry[i,4]-smry[i,10], col=rgb(c1[1],c1[2],c1[3],0.2))
                    segments(x0=smry[i,1]-.02,x1=smry[i,1]+.02,y0=smry[i,4]+smry[i,10],y1=smry[i,4]+smry[i,10], col=rgb(c1[1],c1[2],c1[3],0.2))

                    segments(x0=smry[i,1],x1=smry[i,1],y0=smry[i,5]-smry[i,11],y1=smry[i,5]+smry[i,11], col=rgb(c2[1],c2[2],c2[3],0.2))
                    segments(x0=smry[i,1]-.02,x1=smry[i,1]+.02,y0=smry[i,5]-smry[i,11],y1=smry[i,5]-smry[i,11], col=rgb(c2[1],c2[2],c2[3],0.2))
                    segments(x0=smry[i,1]-.02,x1=smry[i,1]+.02,y0=smry[i,5]+smry[i,11],y1=smry[i,5]+smry[i,11], col=rgb(c2[1],c2[2],c2[3],0.2))
                }
                
                lines(smry2[,1],(smry2[,4]+smry2[,5]))
                points(smry2[,1],(smry2[,4]+smry2[,5]),pch=16)
                lines(smry2[,1],(smry2[,4]),pch=16,col=pcol[1])
                lines(smry2[,1],(smry2[,5]), pch=16,col=pcol[2])
                points(smry2[,1],(smry2[,4]),pch=16,col=pcol[1])
                points(smry2[,1],(smry2[,5]), pch=16,col=pcol[2])
                for(i in 1:nrow(smry2)){
                    segments(x0=smry2[i,1],x1=smry2[i,1],y0=smry2[i,4]+smry2[i,5]-smry2[i,14],y1=smry2[i,4]+smry2[i,5]+smry2[i,14])
                    segments(x0=smry2[i,1]-.02,x1=smry2[i,1]+.02,y0=smry2[i,4]+smry2[i,5]-smry2[i,14],y1=smry2[i,4]+smry2[i,5]-smry2[i,14])
                    segments(x0=smry2[i,1]-.02,x1=smry2[i,1]+.02,y0=smry2[i,4]+smry2[i,5]+smry2[i,14],y1=smry2[i,4]+smry2[i,5]+smry2[i,14])
                    
                    segments(x0=smry2[i,1],x1=smry2[i,1],y0=smry2[i,4]-smry2[i,10],y1=smry2[i,4]+smry2[i,10], col=pcol[1])
                    segments(x0=smry2[i,1]-.02,x1=smry2[i,1]+.02,y0=smry2[i,4]-smry2[i,10],y1=smry2[i,4]-smry2[i,10], col=pcol[1])
                    segments(x0=smry2[i,1]-.02,x1=smry2[i,1]+.02,y0=smry2[i,4]+smry2[i,10],y1=smry2[i,4]+smry2[i,10], col=pcol[1])

                    segments(x0=smry2[i,1],x1=smry2[i,1],y0=smry2[i,5]-smry2[i,11],y1=smry2[i,5]+smry2[i,11], col=pcol[2])
                    segments(x0=smry2[i,1]-.02,x1=smry2[i,1]+.02,y0=smry2[i,5]-smry2[i,11],y1=smry2[i,5]-smry2[i,11], col=pcol[2])
                    segments(x0=smry2[i,1]-.02,x1=smry2[i,1]+.02,y0=smry2[i,5]+smry2[i,11],y1=smry2[i,5]+smry2[i,11], col=pcol[2])
                }
                foo=legend("topleft",c(expression({h^2}["total"]),expression({h^2}["rare"]),expression({h^2}["common"])),pch=16,lty=1,col=c("black",pcol[1],pcol[2]),bty='n')
                text(foo$text$x,foo$rect$top-foo$rect$h-0.01,"Final Gene Set",col=rgb(0,0,0,1),pos=3,cex=0.7)
                text(foo$text$x,foo$rect$top-foo$rect$h-0.02,"All Genes",col=rgb(0,0,0,0.2),pos=3,cex=0.7)
                dev.off()
            }
        }
    }

    { #Figure S11: multi-panel plot showing cdf(h2) vs MAF acros K and minMAC
        PC=10
        gPC = PC
        pPC = PC
        Krange=c(1,2,5,10,20,50);
        minMACrange = c(36,5,2,1)
        
        for(QN in c(0,1)){
            pdf(paste(FIGDIR,"SUPP/S11_he_h2_byKminMAC_PC=",PC,"_QN=",QN,".pdf",sep=""),height=2.5,width=7.5,pointsize=10)
            par(mfrow=c(1,length(minMACrange)), mar=c(2.7,1,1.5,0.5), oma=c(0,3.3,0,0))
            PANEL=0;
            for(mm in minMACrange){
                PANEL=PANEL+1;
                
                if(1){
                    heh2 = matrix(,nr=length(Krange),nc=max(Krange));
                    heh2cdf = matrix(,nr=length(Krange),nc=max(Krange));
                    heh2p = matrix(,nr=length(Krange),nc=max(Krange));
                    heh2cdfp = matrix(,nr=length(Krange),nc=max(Krange));
                    heMAC = matrix(,nr=length(Krange),nc=max(Krange));
                    
                    CNT=0;
                    for(K in Krange){
                        CNT=CNT+1;

                        filename=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1_NORESID.he.h2.gz",sep="");
                        if(file.exists(filename)){
                            GZIN=gzfile(filename,"r");
                            foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                            close(GZIN);
                            GZIN=gzfile(filename,"r");
                            he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                            close(GZIN);
                            
                            G=which(he[,1] %in% expGenes[,1])
                            he = he[G,]
                        }
                        else{
                            cat("File DNE: ",filename,"\n");
                        }

                        filename=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1_PERM=3_NORESID.he.h2.gz",sep="");
                        if(file.exists(filename)){
                            GZIN=gzfile(filename,"r");
                            foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                            close(GZIN);
                            GZIN=gzfile(filename,"r");
                            hep=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                            close(GZIN);
                            
                            G=which(hep[,1] %in% expGenes[,1])
                            hep = hep[G,]
                        }
                        else{
                            cat("File DNE: ",filename,"\n");
                        }
                        
                        if(K>1){
                            heh2[CNT,1:K] = apply(he[,1:K+2],2,function(x){mean(as.numeric(x),na.rm=T)});
                            heh2cdf[CNT,1:K] = cumsum(heh2[CNT,1:K])
                            heh2p[CNT,1:K] = apply(hep[,1:K+2],2,function(x){mean(as.numeric(x),na.rm=T)});
                            heh2cdfp[CNT,1:K] = cumsum(heh2p[CNT,1:K])
                            heMAC[CNT,1:K] = apply(he[,1:K+K+2+gPC+pPC],2,function(x){mean(as.numeric(x),na.rm=T)});
                        }
                        else{
                            heh2[CNT,1] = mean(as.numeric(he[,3]), na.rm=T);
                            heh2cdf[CNT,1] = heh2[CNT,1];
                            heh2p[CNT,1] = mean(as.numeric(hep[,3]), na.rm=T);
                            heh2cdfp[CNT,1] = heh2p[CNT,1];
                            heMAC[CNT,1] = 0.5;
                        }
                    }
                }
                
                plot(1,1,type='n',xlab="",ylab="",xaxt="n",yaxt='n',log='x',ylim=c(0,0.12),xlim=c(0.001,0.5))
                tcol=pcol[c(1,2,3,5,4,9)];
                axis(side=1,padj=-1,labels=NA, cex.axis=0.8)
                axis(side=1,padj=-1,at=c(0.001,0.005,0.02,0.1,0.5),labels=c("0.001","0.005","0.02","0.1","0.5"), cex.axis=0.8)
                axis(side=2,padj=0.9, cex.axis=0.8)
                abline(h=seq(0,0.12,by=0.02),col=rgb(0,0,0,0.1))
                mtext(side=1,"Minor Allele Frequency",line=1.5)
                if(PANEL==1){
                    mtext(side=2,expression(paste("Cumulative ",{widehat(h^2)},sep="")),line=1.25)
                    mtext(side=3,expression("MAF" >= "5%"),line=0);
                }
                if(PANEL==2){
                    mtext(side=3,expression("MAC" >= 5),line=0);
                }
                if(PANEL==3){
                    mtext(side=3,expression("MAC" >= 2),line=0);
                }
                if(PANEL==4){
                    mtext(side=3,"All SNPs",line=0);
                }
                for(i in length(Krange):1){
                    lines(heMAC[i,], heh2cdf[i,],col=tcol[i])
                    points(heMAC[i,], heh2cdf[i,],col=tcol[i], pch=16)

                    ci=col2rgb(tcol[i])/255;
                    lines(heMAC[i,], heh2cdfp[i,],col=rgb(ci[1],ci[2],ci[3],0.25))
                    points(heMAC[i,], heh2cdfp[i,],col=rgb(ci[1],ci[2],ci[3],0.25), pch=16)
                }
                if(mm==36){
                    foo=legend("top",c("K=1","K=10","K=2","K=20","K=5","K=50"),lty=1,pch=16,col=tcol,ncol=3,bty="n")

                    text(10^foo$text$x[1],foo$rect$top-foo$rect$h-0.01,"Raw PCA-corrected",col=rgb(0,0,0,1),pos=4,cex=0.7)
                    text(10^foo$text$x[1],foo$rect$top-foo$rect$h-0.02,"trans-permutation",col=rgb(0,0,0,0.2),pos=4,cex=0.7)

                }
            }
            dev.off()
        }
    }
    
    { #Figure S12: effect of window length
        Kr=c(5,10,20);
        Wr=c("20000","100000","1000000");
        SUFF=c("","","_NORESID");
        he=array(,dim=c(3,20));
        QS12=array(,dim=c(2,3,20));
        QS12c=array(,dim=c(2,3,20));
        hep=array(,dim=c(3,20))
        QS12p=array(,dim=c(2,3,20));
        MAF=array(,dim=c(3,20));
        for(CNT in 1:3){
            K=Kr[CNT];
            W=Wr[CNT];
            suff=SUFF[CNT];
            he[CNT,1:K] = scan(paste("merge_",W,"_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1",suff,".he.meanQuant",sep=""),nlines=1)[1:K];
            hep[CNT,1:K] = scan(paste("merge_",W,"_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3",suff,".he.meanQuant",sep=""),nlines=1)[1:K];
            MAF[CNT,1:K] = scan(paste("merge_",W,"_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1",suff,".he.meanQuant",sep=""),nlines=1)[1:K+K+gPC+pPC];
            filename=paste("merge_",W,"_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1",suff,".he.h2.gz",sep="");
            F=gzfile(filename);
            h2=read.table(F);
            h2 = subset(h2,h2[,1] %in% expGenes[,1])
            filename=paste("merge_",W,"_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3",suff,".he.h2.gz",sep="");
            F=gzfile(filename);
            h2p=read.table(F);
            h2p = subset(h2p,h2p[,1] %in% expGenes[,1])

            NBOOT=1000
            B=array(,dim=c(NBOOT,K));
            Bp=array(,dim=c(NBOOT,K));
            Bc=array(,dim=c(NBOOT,K));
            for(i in 1:NBOOT){
                h2t=h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:K];
                h2t=t(apply(h2t,1,function(x){for(j in which(is.na(x))){x[j]=0}; x}))
                B[i,] = apply(t(apply(h2t,1,cumsum)),2,mean);
                h2t=h2p[sample(1:nrow(h2p),nrow(h2p),replace=T),2+1:K];
                h2t=t(apply(h2t,1,function(x){for(j in which(is.na(x))){x[j]=0}; x}))
                Bp[i,] = apply(t(apply(h2t,1,cumsum)),2,mean);
                Bc[i,] = B[i,]-Bp[i,];
            }
            QS12[1,CNT,1:K] = apply(B,2,function(x){quantile(x,0.975,na.rm=T)})
            QS12[2,CNT,1:K] = apply(B,2,function(x){quantile(x,0.025,na.rm=T)})
            QS12p[1,CNT,1:K] = apply(Bp,2,function(x){quantile(x,0.975,na.rm=T)})
            QS12p[2,CNT,1:K] = apply(Bp,2,function(x){quantile(x,0.025,na.rm=T)})
            QS12c[1,CNT,1:K] = apply(Bc,2,function(x){quantile(x,0.975,na.rm=T)})
            QS12c[2,CNT,1:K] = apply(Bc,2,function(x){quantile(x,0.025,na.rm=T)})
        }

        {
            pdf(paste(FIGDIR,"SUPP/S12_window.pdf",sep=""),height=2,width=2.25,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(0,0.1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            axis(side=1,at=c(1e-3,1e-2,1e-1,0.5),labels=c("1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Cumulative ",{widehat(h^2)},sep="")),line=1.15)
            
            abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
            
            cols = c("blue","green",pcol[4]);
            foo=legend("topleft",c("W=1Mb (K=20)","W=100kb (K=10)","W=20kb (K=5)"), col=rev(cols), lty=1, bg=rgb(1,1,1,0.8),box.col="white",cex=0.8)
            legend(10^foo$rect$left, foo$text$y[3],c("Raw PCA-corrected","trans-permutation","perm-corrected"), lty=c(2,3,1), bg=rgb(1,1,1,0.8),box.col=NA,cex=0.8)
            legend("topleft",c("W=1Mb (K=20)","W=100kb (K=10)","W=20kb (K=5)"), col=rev(cols), lty=1, box.col=NA,cex=0.8)
            box()
            for(i in 1:3){
                K=Kr[i];
                pci=col2rgb(cols[i])/255;
                polygon(x=c(MAF[i,1:K],rev(MAF[i,1:K])), y=c(QS12p[1,i,1:K],rev(QS12p[2,i,1:K])),col=rgb(pci[1],pci[2],pci[3],0.1),border=NA)
                lines(MAF[i,1:K],hep[i,1:K],col=cols[i],lty=3)
            }
            for(i in 1:3){
                K=Kr[i];
                pci=col2rgb(cols[i])/255;
                polygon(x=c(MAF[i,1:K],rev(MAF[i,1:K])), y=c(QS12[1,i,1:K],rev(QS12[2,i,1:K])),col=rgb(pci[1],pci[2],pci[3],0.1),border=NA)
                lines(MAF[i,1:K],he[i,1:K],col=cols[i],lty=2)
            }
            for(i in 1:3){
                K=Kr[i];
                pci=col2rgb(cols[i])/255;
                polygon(x=c(MAF[i,1:K],rev(MAF[i,1:K])), y=c(QS12c[1,i,1:K],rev(QS12c[2,i,1:K])),col=rgb(pci[1],pci[2],pci[3],0.2),border=NA)
                lines(MAF[i,1:K],he[i,1:K]-hep[i,1:K],col=cols[i],lty=1)
                ##points(MAF[i,1:K],he[i,1:K]-hep[i,1:K],col=cols[i],pch=20)
            }
            dev.off()
        }
    }

    { #Fig S13: plot all vs eQTL genes
        { #get data
            RUNBOOT=1;
            K=20;
            PERMstring=""
            FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",gPC,"_npPC=",pPC,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.h2.gz",sep="");
            if(!file.exists(FILE)){
                cat("error. file ",FILE," DNE\n");
            }
            gz = gzfile(FILE,"r");
            foo = scan(file=gz,nlines=1,what="");
            close(gz);
            gz = gzfile(FILE,"r");
            dat.0 = matrix(scan(file=gz,what=""),byrow=T,nc=length(foo));
            close(gz)
            dat.0 = dat.0[which(dat.0[,1] %in% expGenes[,1]),]
            MAFS13 = apply(dat.0[,42+1:K],2,function(x){mean(as.numeric(x),na.rm=T)})
            dat.0 = t(apply(dat.0[,c(1,2+1:K)],1,function(x){for(j in which(is.na(x))){x[j] = 0;}; c(x[1],cumsum(x[2:length(x)]))}));
            
            h2.0=apply(dat.0[,2:ncol(dat.0)], 2, function(x){mean(as.numeric(x))})

            PERMstring="_PERM=3"
            FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",gPC,"_npPC=",pPC,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.h2.gz",sep="");
            if(!file.exists(FILE)){
                cat("error. file ",FILE," DNE\n");
            }
            gz = gzfile(FILE,"r");
            foo = scan(file=gz,nlines=1,what="");
            close(gz);
            gz = gzfile(FILE,"r");
            dat.0p = matrix(scan(file=gz,what=""),byrow=T,nc=length(foo));
            close(gz)
            dat.0p = dat.0p[which(dat.0p[,1] %in% expGenes[,1]),]
            dat.0p = t(apply(dat.0p[,c(1,2+1:K)],1,function(x){for(j in which(is.na(x))){x[j] = 0;}; c(x[1],cumsum(x[2:length(x)]))}));
            
            h2.0p=apply(dat.0p[,2:ncol(dat.0p)], 2, function(x){mean(as.numeric(x))})

            ##eQTL genes
            eQTL=read.table("../../GENESETS/Geuvadis.MatrixEqtl.merged.TopSnp.RhoMaf.EUR.N373.7Columns.tsv",skip=1);
            eQTL = eQTL[which(eQTL[,5]<1e-5),2];
            eQTL=as.matrix(eQTL);
            for(i in 1:nrow(eQTL)){
                eQTL[i,1]=unlist(strsplit(eQTL[i,1],"\\."))[1];
            }
            geneMap=matrix(scan("../../GENESETS/geneMap.txt",what="character",skip=1,sep="\t"),byrow=T,nc=2);
            gm = match(eQTL,geneMap[,2])
            gm = gm[which(!is.na(gm))];
            eQTL = geneMap[gm,1]
            x=which(dat.0[,1] %in% eQTL)
            dat.e = dat.0[x,]
            h2.e = apply(dat.e[,2:ncol(dat.e)], 2, function(x){x=as.numeric(x);mean(x,na.rm=T)});
            x=which(dat.0p[,1] %in% eQTL)
            dat.ep = dat.0p[x,]
            h2.ep = apply(dat.ep[,2:ncol(dat.ep)], 2, function(x){x=as.numeric(x);mean(x,na.rm=T)});

            ##non-eQTL genes
            x=which(!(dat.0[,1] %in% eQTL))
            dat.ne = dat.0[x,]
            h2.ne = apply(dat.ne[,2:ncol(dat.ne)], 2, function(x){x=as.numeric(x);mean(x,na.rm=T)})
            x=which(!(dat.0p[,1] %in% eQTL))
            dat.nep = dat.0p[x,]
            h2.nep = apply(dat.nep[,2:ncol(dat.nep)], 2, function(x){x=as.numeric(x);mean(x,na.rm=T)})

            if(RUNBOOT){
                QS13 = array(,dim=c(2,3,K));
                QS13p = array(,dim=c(2,3,K));
                QS13c = array(,dim=c(2,3,K));
                NBOOT = 1000;
                B = array(,dim=c(NBOOT,3,K))
                Bp = array(,dim=c(NBOOT,3,K))
                Bc = array(,dim=c(NBOOT,3,K))
                for(i in 1:NBOOT){
                    B[i,1,] = apply(dat.0[sample(1:nrow(dat.0),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #all genes
                    B[i,2,] = apply(dat.e[sample(1:nrow(dat.e),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #eQTL genes
                    B[i,3,] = apply(dat.ne[sample(1:nrow(dat.ne),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #non-eQTL genes

                    Bp[i,1,] = apply(dat.0p[sample(1:nrow(dat.0p),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #all genes
                    Bp[i,2,] = apply(dat.ep[sample(1:nrow(dat.ep),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #eQTL genes
                    Bp[i,3,] = apply(dat.nep[sample(1:nrow(dat.nep),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #non-eQTL genes

                    for(j in 1:3){
                        Bc[i,j,] = B[i,j,]-Bp[i,j,];
                    }
                }
                for(i in 1:3){
                    QS13[1,i,1:K] = apply(B[,i,1:K],2,function(x){quantile(x,0.975)})
                    QS13[2,i,1:K] = apply(B[,i,1:K],2,function(x){quantile(x,0.025)})

                    QS13p[1,i,1:K] = apply(Bp[,i,1:K],2,function(x){quantile(x,0.975)})
                    QS13p[2,i,1:K] = apply(Bp[,i,1:K],2,function(x){quantile(x,0.025)})

                    QS13c[1,i,1:K] = apply(Bc[,i,1:K],2,function(x){quantile(x,0.975)})
                    QS13c[2,i,1:K] = apply(Bc[,i,1:K],2,function(x){quantile(x,0.025)})
                }
            }

        }
        { #plot CDFs with envelope
            cols = c(pcol[4],"green","blue");
            CVstring="";
            pdf(paste(FIGDIR,"SUPP/S13_CDF_env_trueVperm+eQTL_NORESID_gPC=",gPC,"_pPC=",pPC,"_QN=",QN,CVstring,".pdf",sep=""),height=2,width=2.25,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1, type='n', log='x', xlim=c(0.001,0.5), ylim=c(0,0.14), xaxt='n', yaxt='n', pch=16, col="black",xlab="",ylab="")
            axis(side=1,at=c(1e-3,1e-2,1e-1,0.5),labels=c("1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Cumulative ",{widehat(h^"2")},sep="")),line=1.15)
            abline(h=seq(0,0.14,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
            foo=legend("topleft",c("eQTL genes (n=4322)","all genes (n=10078)","non-eQTL genes (n=5756)"), col=c(cols[2],cols[1],cols[3]), lty=1, bg=rgb(1,1,1,0.8),box.col="white",cex=0.7)
            legend(10^foo$rect$left, foo$text$y[3],c("Raw PCA-corrected","trans-permutation","perm-corrected"), lty=c(2,3,1), bg=rgb(1,1,1,0.8),box.col=NA,cex=0.7)
            legend("topleft",c("eQTL genes (n=4322)","all genes (n=10078)","non-eQTL genes (n=5756)"), col=c(cols[2],cols[1],cols[3]), lty=1, box.col=NA,cex=0.7)
            box()

            for(i in 1:3){
                pci = col2rgb(cols[i])/255;
                polygon(x=c(MAFS13,rev(MAFS13)), y=c(QS13p[1,i,1:K], rev(QS13p[2,i,1:K])), col=rgb(pci[1], pci[2], pci[3], 0.1), border=NA)
            }
            lines(MAFS13,h2.0p,col=cols[1],lty=3)
            lines(MAFS13,h2.ep,col=cols[2],lty=3)
            lines(MAFS13,h2.nep,col=cols[3],lty=3)

            for(i in 1:3){
                pci = col2rgb(cols[i])/255;
                polygon(x=c(MAFS13,rev(MAFS13)), y=c(QS13[1,i,1:K], rev(QS13[2,i,1:K])), col=rgb(pci[1], pci[2], pci[3], 0.1), border=NA)
            }
            lines(MAFS13,h2.0,col=cols[1],lty=2)
            lines(MAFS13,h2.e,col=cols[2],lty=2)
            lines(MAFS13,h2.ne,col=cols[3],lty=2)

            for(i in 1:3){
                pci = col2rgb(cols[i])/255;
                polygon(x=c(MAFS13,rev(MAFS13)), y=c(QS13c[1,i,1:K], rev(QS13c[2,i,1:K])), col=rgb(pci[1], pci[2], pci[3], 0.2), border=NA)
            }
            lines(MAFS13,h2.0-h2.0p,col=cols[1],lty=1,lwd=2)
            lines(MAFS13,h2.e-h2.ep,col=cols[2],lty=1,lwd=2)
            lines(MAFS13,h2.ne-h2.nep,col=cols[3],lty=1,lwd=2)

            dev.off()
        }
    }


    { #Fig S14: plot cdf h2 for different PCs
        { #get data
            NBOOT = 1000;
            W="1000000";
            K=20;
            PERMstring=""
            PCr = c(0,1,5,10,20);

            QS14 = array(,dim=c(3,length(PCr),K));
            QS14p = array(,dim=c(3,length(PCr),K));
            QS14c = array(,dim=c(3,length(PCr),K));

            for(i in 1:length(PCr)){
                PC=PCr[i];
                PERMstring=""
                FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.h2.gz",sep="");
                if(file.exists(FILE)){
                    gz = gzfile(FILE,"r");
                    foo = scan(file=gz,nlines=1,what="");
                    close(gz);
                    gz = gzfile(FILE,"r");
                    dat.0 = matrix(scan(file=gz,what=""),byrow=T,nc=length(foo));
                    close(gz)
                    dat.0 = dat.0[which(dat.0[,1] %in% expGenes[,1]),]
                    MAF = apply(dat.0[,42+1:K],2,function(x){mean(as.numeric(x),na.rm=T)})
                    dat.0 = t(apply(dat.0[,c(1,2+1:K)],1,function(x){for(j in which(is.na(x))){x[j] = 0;}; c(x[1],cumsum(x[2:length(x)]))}));
                    
                    QS14[3,i,1:K]=apply(dat.0[,2:ncol(dat.0)], 2, function(x){mean(as.numeric(x))})
                    
                    PERMstring="_PERM=3"
                    FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",gPC,"_npPC=",pPC,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.h2.gz",sep="");
                    if(!file.exists(FILE)){
                        cat("error. file ",FILE," DNE\n");
                    }
                    gz = gzfile(FILE,"r");
                    foo = scan(file=gz,nlines=1,what="");
                    close(gz);
                    gz = gzfile(FILE,"r");
                    dat.0p = matrix(scan(file=gz,what=""),byrow=T,nc=length(foo));
                    close(gz)
                    dat.0p = dat.0p[which(dat.0p[,1] %in% expGenes[,1]),]
                    dat.0p = t(apply(dat.0p[,c(1,2+1:K)],1,function(x){for(j in which(is.na(x))){x[j] = 0;}; c(x[1],cumsum(x[2:length(x)]))}));
                    
                    QS14p[3,i,1:K]=apply(dat.0p[,2:ncol(dat.0p)], 2, function(x){mean(as.numeric(x))})
                    
                    QS14c[3,i,1:K]=QS14[3,i,1:K]-QS14p[3,i,1:K];

                    if(PC == 10){
                        B = array(,dim=c(NBOOT,K))
                        Bp = array(,dim=c(NBOOT,K))
                        Bc = array(,dim=c(NBOOT,K))
                        for(j in 1:NBOOT){
                            B[j,] = apply(dat.0[sample(1:nrow(dat.0),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #all genes
                            Bp[j,] = apply(dat.0p[sample(1:nrow(dat.0p),replace=T),1+1:K], 2, function(x){mean(as.numeric(x),na.rm=T)}) #all genes
                            Bc[j,] = B[j,]-Bp[j,];
                        }
                        
                        QS14[2,i,1:K] = apply(B, 2, function(x){quantile(x,0.025)})
                        QS14p[2,i,1:K] = apply(Bp, 2, function(x){quantile(x,0.025)})
                        QS14c[2,i,1:K] = apply(Bc, 2, function(x){quantile(x,0.025)})
                        
                        QS14[1,i,1:K] = apply(B, 2, function(x){quantile(x,0.975)})
                        QS14p[1,i,1:K] = apply(Bp, 2, function(x){quantile(x,0.975)})
                        QS14c[1,i,1:K] = apply(Bc, 2, function(x){quantile(x,0.975)})
                    }
                }
                else{
                    PERMstring=""
                    FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.meanQuant",sep="");
                    if(!file.exists(FILE)){
                        cat("error. file ",FILE," DNE\n");
                    }
                    QS14[3,i,1:K]=matrix(scan(file=FILE,nlines=1)[1:K],byrow=T,nc=K);
                    PERMstring="_PERM=3"
                    FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.meanQuant",sep="");
                    if(!file.exists(FILE)){
                        cat("error. file ",FILE," DNE\n");
                    }
                    QS14p[3,i,1:K]=matrix(scan(file=FILE,nlines=1)[1:K],byrow=T,nc=K);
                    QS14c[3,i,1:K] = QS14[3,i,1:K]-QS14p[3,i,1:K];
                }
            }
        }

        { #plot CDFs with envelope
            cols = c(pcol[4],"green","blue");
            CVstring="";
            pdf(paste(FIGDIR,"SUPP/S14_CDF_env_h2vsPC.pdf",sep=""),height=2.25,width=2.25,pointsize=10)
            par(mar=c(2,2.7,0.3,0.3))
            plot(1,1, type='n', log='x', xlim=c(0.001,0.5), ylim=c(0.002,0.1), xaxt='n', yaxt='n', pch=16, col="black",xlab="",ylab="")
            axis(side=1,at=c(1e-3,1e-2,1e-1,0.5),labels=c("1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Cumulative ",{widehat(h^"2")},sep="")),line=1.15)
            abline(h=seq(0,0.14,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
            foo=legend("topleft",c("PC=0","PC=1","PC=5","PC=10","PC=20"), col=pcol[1:length(PCr)], lty=1, bg=rgb(1,1,1,0.8),box.col="white",cex=0.7)
            legend(10^foo$rect$left, foo$text$y[length(PCr)],c("Raw PCA-corrected","trans-permutation","perm-corrected"), lty=c(2,3,1), bg=rgb(1,1,1,0.8),box.col=NA,cex=0.7)
            legend("topleft",c("PC=0","PC=1","PC=5","PC=10","PC=20"), col=pcol[1:length(PCr)], lty=1,box.col=NA,cex=0.7)
            box()
            cols=pcol;

            for(i in 1:length(PCr)){
                pci = col2rgb(cols[i])/255;
                lines(MAF,QS14p[3,i,1:K],col=cols[i],lty=3)
                lines(MAF,QS14[3,i,1:K],col=cols[i],lty=2)
                lines(MAF,QS14c[3,i,1:K],col=cols[i],lty=1)
            }
            i=4;
            pci = col2rgb(cols[i])/255;
            polygon(x=c(MAF,rev(MAF)), y=c(QS14p[1,i,1:K], rev(QS14p[2,i,1:K])), col=rgb(pci[1], pci[2], pci[3], 0.5), border=NA)
            polygon(x=c(MAF,rev(MAF)), y=c(QS14[1,i,1:K], rev(QS14[2,i,1:K])), col=rgb(pci[1], pci[2], pci[3], 0.5), border=NA)
            polygon(x=c(MAF,rev(MAF)), y=c(QS14c[1,i,1:K], rev(QS14c[2,i,1:K])), col=rgb(pci[1], pci[2], pci[3], 0.5), border=NA)


            dev.off()
        }
    }
    
    {##Figure S15: h2_bin as a function of genotype likelihood
        K=4
        GL=array(,dim=K);
        GLp=array(,dim=K);
        GLc=array(,dim=K);
        ##now get gnomAD partition data
        filename="merge_1000000_K=4_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_GLpart.he.meanQuant"
        foo=read.table(filename);
        h2 = foo[1,1:K];
        ##transform to h2 per bin
        for(j in K:2){
            h2[j] = h2[j] - h2[j-1];
        }
        ##now get permuted values
        filename="merge_1000000_K=4_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_GLpart.he.meanQuant";
        foo=read.table(filename);
        h2p = foo[1,1:K];
        ##transform to h2 per bin
        for(j in K:2){
            h2p[j] = h2p[j] - h2p[j-1];
        }
        GL[1:K] = unlist(h2)
        GLp[1:K] = unlist(h2p)
        GLc[1:K] = unlist(h2-h2p)
        
        pdf(paste(FIGDIR,"SUPP/S15_h2_vs_GL_SING.pdf",sep=""),height=2,width=2.25,pointsize=10)
        par(mar=c(2,2.8,0.3,0.3))
        plot(1,1,type='n',xlim=c(1,4), ylim=c(-0.04,0.06),xaxt='n',yaxt='n',xlab="",ylab="")
        abline(h=seq(-.04,.06,by=0.02),col=rgb(0,0,0,0.1),lwd=0.75)
        axis(side=1,padj=-2,at=1:4,cex.axis=0.7,tck=-0.03)
        axis(side=2,padj=1.3, cex.axis=0.7)
        axis(side=2,padj=1.3, cex.axis=0.7,at=-0.02)
        mtext(side=1,"Quartile",line=1.025)
        mtext(side=1, adj=0, "Low GL", cex=0.7, line=.7)
        mtext(side=1, adj=1, "High GL", cex=0.7, line=.7)
        mtext(side=2,expression(paste({widehat(h^2)}[singleton],sep="")),line=1.15)
        legend("topleft",c("Raw PCA-corrected","trans-permutation","perm-corrected"), bg=rgb(1,1,1,0.8),box.col=NA,cex=0.7,pch=c(6,1,20),lty=c(2,3,1),lwd=c(1,1,2))
        box()
        points(1:4,GL,col=pcol[10],pch=6)
        points(1:4,GLp,col=pcol[10],pch=1)
        points(1:4,GLc,col=pcol[10],pch=20)
        lines(1:4,GL,col=pcol[10],lty=2)
        lines(1:4,GLp,col=pcol[10],lty=3)
        lines(1:4,GLc,col=pcol[10],lty=1,lwd=2)
        dev.off()           
    }
    
    { #Figure S16: strict mask analysis with SNPs and exons (eg for RNAseq read mapping Q)
        W="1000000"
        CVstring="";
        K=20
        gPC=10
        pPC=10
        QN=1
        mm=1
        GORD = read.table("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/matrix.eqtl/gencode.v12.exonStrictMaskOverlap.txt")
        
        { ##get data
            QS16 = array(,dim=c(7,3,3,K)); #[data subset, raw/perm/corr, quantiles, 1:K]
            MAF=array(,dim=c(3,K));
            PS=c("","_PERM=3");
            for(p in 1:length(PS)){
                PERMSTRING=PS[p];
                FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",gPC,"_npPC=",pPC,"_QN=",QN,"_DSAMP=1",PERMSTRING,"_NORESID.he.h2.gz",sep="");
                if(file.exists(FILE)){
                    gz = gzfile(FILE,"r");
                    foo = scan(file=gz,nlines=1,what="");
                    close(gz);
                    gz = gzfile(FILE,"r");
                    dat.1 = matrix(scan(file=gz,what=""),byrow=T,nc=length(foo));
                    dat.tot = dat.1;
                    close(gz)
                    dat.1 = dat.1[which(dat.1[,1] %in% expGenes[,1]),]
                    QS16[1,p,3,1:K]=apply(t(apply(dat.1[,2+1:K],1,cumsum)), 2, function(x){x=as.numeric(x); mean(x,na.rm=T)})
                    if(p==1){
                        MAF[p,1:K] = apply(dat.1[,1:K+2+K+gPC+pPC], 2, function(x){x=as.numeric(x); mean(x,na.rm=T)});
                    }
                }
                else{
                    cat("error. file ",FILE," DNE\n");
                }
                
                SM="_SM";
                FILE=paste("merge_",W,"_K=",K,"_MAC=",mm,"_ngPC=",gPC,"_npPC=",pPC,"_QN=",QN,"_DSAMP=1",PERMSTRING,"_NORESID",SM,".he.meanQuant",sep="");
                if(file.exists(FILE)){
                    dat.2 = read.table(file=FILE);
                    QS16[2,p,3,1:K]=scan(file=FILE,nlines=1)[1:K]
                }
                else{
                    cat("error. file ",FILE," DNE\n");
                }
                
                GSUB = which(GORD[,1] %in% dat.1[,1]) 
                GSUB2 = which(!(GORD[,1] %in% dat.1[,1]))
                GSUB3 = which(!(dat.tot[,1] %in% dat.1[,1]))
                QS16[7,p,3,1:K]=apply(t(apply(dat.tot[GSUB3,2+1:K],1,cumsum)), 2, function(x){x=as.numeric(x); mean(x,na.rm=T)})
                
                QNT = c(0,0.25,0.5,0.75,1);
                for(i in 2:5){
                    GQ = array(GORD[round((QNT[i-1]*nrow(GORD)):(QNT[i]*nrow(GORD))),1]);
                    QS16[1+i,p,3,1:K]=apply(t(apply(dat.1[which(dat.1[,1] %in% GQ),2+1:K],1,cumsum)), 2, function(x){x=as.numeric(x); mean(x,na.rm=T)})
                }
            }
            for(i in 1:7){
                QS16[i,3,3,1:K] = QS16[i,1,3,1:K]-QS16[i,2,3,1:K];
            }
        }

        { #plot CDFs
            pdf(paste(FIGDIR,"SUPP/S16_CDF_strictMask_EXONS_NORESID_gPC=",gPC,"_pPC=",pPC,"_QN=",QN,CVstring,".pdf",sep=""),height=2.25,width=2.25,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1, type='n', log='x', xlim=c(0.001,0.5), ylim=c(-.04,0.1), xaxt='n', yaxt='n', pch=16, col="black",xlab="",ylab="")
            axis(side=1,at=c(1e-3,1e-2,1e-1,0.5),labels=c("1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Cumulative ",{widehat(h^2)},sep="")),line=1.15)
            abline(h=seq(-.04,.06,by=0.02),col=rgb(0,0,0,0.1),lwd=0.75)

            colord=c(4,1,2,3,5,6,7);
            legend("topleft",c("All SNPs","Strict Mask SNPs","SME Q1","SME Q2","SME Q3","SME Q4","Excluded"),col=pcol[colord],lty=1, bty="n", cex=0.7)
            for(i in 1:7){
                pci=col2rgb(pcol[colord[i]])/255;
                lines(MAF[1,1:K],QS16[i,2,3,1:K],col=rgb(pci[1],pci[2],pci[3],0.5),lty=3)
                lines(MAF[1,1:K],QS16[i,1,3,1:K],col=rgb(pci[1],pci[2],pci[3],0.5),lty=2)
                lines(MAF[1,1:K],QS16[i,3,3,1:K],col=rgb(pci[1],pci[2],pci[3],1),lty=1,lwd=2)
            }
            dev.off()
        }
    }
    
    { #Fig S17: compare HE vs REML: plot number of genes converged as fn(K).
        { ## get data
            W="1000000";
            Krange=c(1,2,5,10,20,50);
            PC=0;
            minMAC = 1;
            
            numG.1=matrix(0,nr=4,nc=length(Krange));
            
            QN=1;
            CNT=0;
            for(K in Krange){
                CNT=CNT+1;
                GZIN=gzfile(paste("merge_",W,"_K=",K,"_MAC=",minMAC,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1_NORESID.he.h2.gz",sep=""),"r");
                he=scan(file=GZIN,nlines=1,what="",quiet=T)
                close(GZIN);
                GZIN=gzfile(paste("merge_",W,"_K=",K,"_MAC=",minMAC,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1_NORESID.he.h2.gz",sep=""),"r");
                he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(he))
                close(GZIN);
                
                numG.1[1,CNT] = nrow(he);
                
                G=which(he[,1] %in% expGenes[,1])
                he = he[G,]
                numG.1[3,CNT] = nrow(he);
                
                GZIN=gzfile(paste("merge_",W,"_K=",K,"_MAC=",minMAC,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1_NORESID.reml.h2.gz",sep=""),"r");
                reml=scan(file=GZIN,nlines=1,what="",quiet=T)
                close(GZIN);
                GZIN=gzfile(paste("merge_",W,"_K=",K,"_MAC=",minMAC,"_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1_NORESID.reml.h2.gz",sep=""),"r");
                reml=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(reml))
                close(GZIN);
                
                numG.1[2,CNT] = nrow(reml);
                
                G=which(reml[,1] %in% expGenes[,1])
                
                if(length(G) > 1){
                    reml = reml[G,]
                    numG.1[4,CNT] = nrow(reml);
                }
            }
        }
        
        {## number of genes analyzed
            pdf(paste(FIGDIR,"SUPP/S17_heVSreml_numG_byK_minMAC=",minMAC,"_PC=",PC,"_QN=",QN,".pdf",sep=""),height=1.75,width=2.25,pointsize=10)
            par(mar=c(2.3,2.3,0.8,0.5))
            plot(1,1,type='n',xlab="",ylab="",xaxt="n",yaxt='n',log='x',ylim=range(numG.1),xlim=range(Krange),xlog=T)
            axis(side=1,padj=-1, cex.axis=0.8)
            axis(side=2,padj=1, cex.axis=0.8)
            mtext(side=1,"Number of SNP bins (K)",line=1.2)
            mtext(side=2,"# Genes",line=1.2)
            c1=col2rgb(pcol[2])/255;
            c2=col2rgb(pcol[1])/255;
            lines(Krange, numG.1[1,], col=rgb(c1[1],c1[2],c1[3],0.3), type='b', pch=16)
            lines(Krange, numG.1[2,], col=rgb(c2[1],c2[2],c2[3],0.3), type='b', pch=16)
            lines(Krange, numG.1[3,], col=pcol[2], type='b', pch=16)
            lines(Krange, numG.1[4,], col=pcol[1], type='b', pch=16)
            
            text(30,21000,"All Genes",pos=3,cex=0.7,col=rgb(0,0,0,0.3))
            text(30,6000,"Final Subset",pos=3,cex=0.7)
            
            legend("bottomleft",c("HE","REML"),col=c(pcol[2:1]),bty='n',lty=1,pch=16,cex=0.7,ncol=1)
            dev.off()
        }
    }


    { #Fig S18: shadow effects, rare variants tagging common and vice versa
        if(1){ #get data
            K=20;
            H2MAT = matrix(, nr=9, nc=K);
            MAFS18 = matrix(, nr=9, nc=K+1);
            PC=10;
            gPC=PC
            pPC=PC
            CVstring=""
            W="1000000"
            QN=1
            PERMstring=""
            
            CNT = 1;
            filename = "";
            filename[1]="merge_1000000_K=1_MAC=1-1_ngPC=10_npPC=10_QN=1_DSAMP=1";
            K[1] = 1;
            filename[2]="merge_1000000_K=5_MAC=1-8_ngPC=10_npPC=10_QN=1_DSAMP=1";
            K[2] = 5;
            filename[3]="merge_1000000_K=8_MAC=1-23_ngPC=10_npPC=10_QN=1_DSAMP=1";
            K[3] = 8;
            filename[4]="merge_1000000_K=11_MAC=1-72_ngPC=10_npPC=10_QN=1_DSAMP=1";
            K[4] = 11;
            filename[5]="merge_1000000_K=15_MAC=1-161_ngPC=10_npPC=10_QN=1_DSAMP=1";
            K[5] = 15;
            filename[6]="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1";
            K[6] = 20;
            
            for(i in 1:6){
                if(i<6){
                    FILE=paste(filename[i],".he.h2.gz",sep="");
                }
                else{
                    FILE=paste(filename[i],"_NORESID.he.h2.gz",sep="");
                }
                GZIN=gzfile(FILE,"r");
                foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                close(GZIN);
                GZIN=gzfile(FILE,"r");
                he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                close(GZIN);
                dim(he)
                G=which(he[,1] %in% expGenes[,1])
                he = he[G,]
                dim(he)

                if(i<6){
                    FILE=paste(filename[i],"_PERM=3.he.h2.gz",sep="");
                }
                else{
                    FILE=paste(filename[i],"_PERM=3_NORESID.he.h2.gz",sep="");
                }
                GZIN=gzfile(FILE,"r");
                foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                close(GZIN);
                GZIN=gzfile(FILE,"r");
                he2=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                close(GZIN);
                dim(he2)
                G=which(he2[,1] %in% expGenes[,1])
                he2 = he2[G,]
                dim(he2)

                if(K[i] > 1){
                    H2MAT[i,1:K[i]] = apply(t(apply(he[,2+1:K[i]],1,cumsum)),2,function(x){mean(as.numeric(x),na.rm=T)})-apply(t(apply(he2[,2+1:K[i]],1,cumsum)),2,function(x){mean(as.numeric(x),na.rm=T)})
                    MAFS18[i,1:K[i]] = apply(he[,2+K[i]+(2*PC)+1:K[i]],2,function(x){mean(as.numeric(x),na.rm=T)})
                }
                else{
                    H2MAT[i,1] = mean(as.numeric(he[,3]),na.rm=T)-mean(as.numeric(he2[,3]),na.rm=T)
                    MAFS18[i,1] = mean(as.numeric(he[,2+K[i]+(2*PC)+1:K[i]]),na.rm=T)
                }

                if(i==6){
                    ##GET CI and subset CI
                    BKS=c(1,2,4,9)
                    QS18 = array(,dim=c(2,4,20));
                    NBOOT=1000;
                    for(j in 1:length(BKS)){
                        BOOT = matrix(,nr=NBOOT,nc=20);
                        for(k in 1:NBOOT){
                            t=sample(1:nrow(he),nrow(he),replace=T);
                            t2=sample(1:nrow(he2),nrow(he2),replace=T);
                            BOOT[k,BKS[j]:20] = apply(t(apply(he[t,2+BKS[j]:20],1,cumsum)),2,function(x){mean(as.numeric(x),na.rm=T)})-apply(t(apply(he2[t2,2+BKS[j]:20],1,cumsum)),2,function(x){mean(as.numeric(x),na.rm=T)});
                        }
                        QS18[1,j,] = apply(BOOT,2,function(x){quantile(x,0.025,na.rm=T)})
                        QS18[2,j,] = apply(BOOT,2,function(x){quantile(x,0.975,na.rm=T)})
                    }
                }
            }


            CNT=6;
            for(mmt in c(2,5,36)){
                CNT = CNT+1;
                K[CNT] = 20;
                PERMstring="";
                FILE=paste("merge_",W,"_K=",K[CNT],"_MAC=",mmt,"_ngPC=",gPC,"_npPC=",pPC,CVstring,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.h2.gz",sep="");
                if(!file.exists(FILE)){
                    cat("error. file ",FILE," DNE\n");
                }
                gz = gzfile(FILE,"r");
                foo = scan(file=gz,nlines=1,what="");
                close(gz);
                gz = gzfile(FILE,"r");
                he = matrix(scan(file=gz,what=""),byrow=T,nc=length(foo));
                dim(he)
                close(gz)
                he = he[which(he[,1] %in% expGenes[,1]),]
                dim(he)

                PERMstring="_PERM=3";
                FILE=paste("merge_",W,"_K=",K[CNT],"_MAC=",mmt,"_ngPC=",gPC,"_npPC=",pPC,CVstring,"_QN=",QN,"_DSAMP=1",PERMstring,"_NORESID.he.h2.gz",sep="");
                if(!file.exists(FILE)){
                    cat("error. file ",FILE," DNE\n");
                }
                gz = gzfile(FILE,"r");
                foo = scan(file=gz,nlines=1,what="");
                close(gz);
                gz = gzfile(FILE,"r");
                he2 = matrix(scan(file=gz,what=""),byrow=T,nc=length(foo));
                dim(he2)
                close(gz)
                he2 = he2[which(he2[,1] %in% expGenes[,1]),]
                dim(he2)

                
                H2MAT[CNT,1:K[CNT]] = apply(t(apply(he[,2+1:K[CNT]],1,cumsum)),2,function(x){mean(as.numeric(x),na.rm=T)})-apply(t(apply(he2[,2+1:K[CNT]],1,cumsum)),2,function(x){mean(as.numeric(x),na.rm=T)})
                MAFS18[CNT,1:K[CNT]] = apply(he[,2+K[CNT]+(2*PC)+1:K[CNT]],2,function(x){mean(as.numeric(x),na.rm=T)})
            }
        }

        tmpcol = pcol[c(1:3,5:6,4,7:9)]
        

        

        {# plot it on single panel
            pdf(paste(FIGDIR,"SUPP/S18_taggingEffects_1panel+leg_2.pdf",sep=""),height=2,width=4,pointsize=10)
            par(mar=c(2,2.8,0.3,5.5))
            plot(1,1,type='n',ylim=range(0,0.065),xlim=range(0.001,MAFS18,na.rm=T),xaxt='n',yaxt='n',xlab="",ylab="",log='x')
            abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
            
            tcol = pcol
            tcol[6]=pcol[3]
            tcol[3]=pcol[6]
            rcol = col2rgb(tcol[6])/255;

            axis(side=1,padj=-2,labels=NA, cex.axis=0.7)
            axis(side=1,padj=-2,at=c(0.001,0.005,0.02,0.1,0.5),labels=c("0.001","0.005","0.02","0.1","0.5"),cex.axis=0.7)
            axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.01,0.09,by=0.02),labels=NA,tck=-0.04)
            axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,0.1,by=0.02),labels=c("0","0.02","0.04","0.06","0.08","0.1"))

            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Cumulative ",{widehat(h^"2'")},sep="")),line=1.15)
            pci = col2rgb(tmpcol[6])/255; #use same shading for all envelopes
            for(i in 1:9){
                if(i>=6){
                    tmaf=MAFS18[6,BKS[i-5]:K[i]];
                    polygon(x=c(tmaf,rev(tmaf)), y=c(QS18[1,i-5,BKS[i-5]:K[i]],rev(QS18[2,i-5,BKS[i-5]:K[i]])),col=rgb(pci[1],pci[2],pci[3],0.25),border=NA)
                    
                }
                lines(MAFS18[i,1:K[i]],H2MAT[i,1:K[i]], col=tmpcol[i])
                points(MAFS18[i,1:K[i]],H2MAT[i,1:K[i]], col=tmpcol[i], pch=20)
            }

            par(xpd=T)
            legend(0.75,0.065,c("Singletons",expression("MAF" <= "1%"),expression("MAF" <= "3%"),expression("MAF" <= "10%"),expression("MAF" <= "5%"),"ALL SNPs",expression("MAF" >= "0.2%"),expression("MAF" >= "1.8%"),expression("MAF" >= "5%")), col=tmpcol, lty=c(0,rep(1,8)),lwd=2,ncol=1,bty='n',pch=20,cex=0.7)
            dev.off()

        }

    }

    {## Fig S19: reverse cumuqlative h2 for each minMAC
        { ##get data
            K=20;
            nSB=6;
            Kt=K+nSB-1
            QS19 = array(,dim=c(3,2,Kt));

            MAF = array(,dim=c(1,K))
            filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="");
            foo=read.table(filename);
            MAF[1:K] = unlist(foo[1,K+(gPC+pPC)+1:K])
            SINGMAF = c(1e-5, 6.46e-5, 0.0002, 0.0004, 1e-3, 1/720);
            gMAF = c(SINGMAF[1:nSB], MAF[2:K])



            ##first get data with GEUVADIS MAF partition
            filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.h2.gz",sep="");
            F=gzfile(filename);
            h2=read.table(F);
            dim(h2)
            h2 = subset(h2,h2[,1] %in% expGenes[,1])
            dim(h2)
            
            filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz",sep="");
            F=gzfile(filename);
            h2p=read.table(F);
            dim(h2p)
            h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
            dim(h2p);

            NBOOT=1000
            B=array(,dim=c(NBOOT,K));
            BP=array(,dim=c(NBOOT,K));
            BD=array(,dim=c(NBOOT,K)); #difference
            
            for(j in 1:NBOOT){
                h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:K];
                B[j,] = apply(t(apply(h2t,1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});
                h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:K];
                BP[j,] = apply(t(apply(h2t,1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});
                BD[j,] = (B[j,]-BP[j,]);
            }
            
            QS19[3,1,1:K] = apply(t(apply(h2[,2+1:K],1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});
            QS19[3,1,1:K] = QS19[3,1,1:K]-apply(t(apply(h2p[,2+1:K],1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});

            for(j in 1:K){
                QS19[1,1,j] = quantile(BD[,j],0.975,na.rm=T) #upper Q
                QS19[2,1,j] = quantile(BD[,j],0.025,na.rm=T) #lower Q
            }
            
            ## now get data for gnomad partition
            Kt=25;
            filename=paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_gnomAD.SINGpart.he.h2.gz",sep="")
            F=gzfile(filename);
            h2=read.table(F);
            ##close(F);
            dim(h2)
            h2 = subset(h2,h2[,1] %in% expGenes[,1])
            dim(h2)
            
            filename=paste("merge_1000000_K=",Kt,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_gnomAD.SINGpart.he.h2.gz",sep="");
            F=gzfile(filename);
            h2p=read.table(F);
            ##close(F);
            dim(h2p)
            h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
            dim(h2p);
            
            B2=array(,dim=c(NBOOT,Kt));
            BP2=array(,dim=c(NBOOT,Kt));
            BD2=array(,dim=c(NBOOT,Kt)); #difference
            
            for(j in 1:NBOOT){
                h2t = h2[sample(1:nrow(h2),nrow(h2),replace=T),2+1:Kt];
                B2[j,] = apply(t(apply(h2t,1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});
                h2t = h2p[sample(1:nrow(h2),nrow(h2p),replace=T),2+1:Kt];
                BP2[j,] = apply(t(apply(h2t,1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});
                BD2[j,] = (B2[j,]-BP2[j,]);
            }
            
            QS19[3,2,1:Kt] = apply(t(apply(h2[,2+1:Kt],1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});
            QS19[3,2,1:Kt] = QS19[3,2,1:Kt]-apply(t(apply(h2p[,2+1:Kt],1,function(x){rev(cumsum(rev(x)))})),2,function(x){mean(x,na.rm=T)});

            for(j in 1:Kt){
                QS19[1,2,j] = quantile(BD2[,j],0.975,na.rm=T) #upper Q
                QS19[2,2,j] = quantile(BD2[,j],0.025,na.rm=T) #lower Q
            }
        }

        
        { ## now plot it
            
            pdf(paste(FIGDIR,"SUPP/S19_revcumh2.pdf",sep=""),height=1.75,width=4,pointsize=10)
            par(mar=c(2.,2.6,0.3,0.3))
            plot(1,1,type='n',ylim=range(0,0.09),xlim=range(gMAF,MAF,na.rm=T),xaxt='n',yaxt='n',xlab="",ylab="",log='x')
            abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
            axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.01,0.09,by=0.02),labels=NA,tck=-0.04)
            axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,0.1,by=0.02),labels=c("0","0.02","0.04","0.06","0.08","0.1"))

            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Rev. Cum. ",{widehat(h^"2'")},sep="")),line=1.15)

            rcol = col2rgb(pcol[5])/255;
            polygon(x=c(gMAF,rev(gMAF)),y=c(QS19[1,2,1:Kt],rev(QS19[2,2,1:Kt])), col=rgb(rcol[1],rcol[2],rcol[3],0.3),border=NA)

            rcol = col2rgb(pcol[4])/255;
            polygon(x=c(MAF,rev(MAF)),y=c(QS19[1,1,1:K],rev(QS19[2,1,1:K])), col=rgb(rcol[1],rcol[2],rcol[3],0.3),border=NA)
            
            lines(gMAF,QS19[3,2,1:Kt],col=pcol[5])
            lines(MAF,QS19[3,1,1:K],col=pcol[4])
            
            points(gMAF,QS19[3,2,1:Kt],col=pcol[5],pch=20)
            points(MAF,QS19[3,1,1:K],col=pcol[4],pch=20)

            legend("topright",c("Pooled Singletons","gnomAD MAF partition"), col=pcol[4:5], lty=1,lwd=2,ncol=1,bty='n',pch=20,cex=0.7)

            dev.off();
        }
    }

    {## Fig S20: h2 as a function of GERP score
        { ##get data
            K=20;
            QS20 = array(,dim=c(3,2,K));
            
            GERP = array(,dim=c(2,K))
            GERP[1,] = read.table("SINGgerp_SFS_K=20.txt")[,1]
            GERP[2,] = read.table("ALLgerp_SFS_K=20.txt")[,1]
            GERP[1,1] = -GERP[1,K];
            GERP[2,1] = -GERP[2,K];
            filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_SINGgerp.he.meanQuant",sep="");
            foo=read.table(filename);
            filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_ALLgerp.he.meanQuant",sep="");
            foo2=read.table(filename);
            QS20[3,1,1] = foo[1,1]
            QS20[3,2,1] = foo2[1,1]
            QS20[1,1,1] = foo[2,1]
            QS20[1,2,1] = foo2[2,1]
            QS20[2,1,1] = foo[3,1]
            QS20[2,2,1] = foo2[3,1]            
            for(i in 2:K){
                QS20[3,1,i] = foo[1,i]-foo[1,i-1];
                QS20[3,2,i] = foo2[1,i]-foo2[1,i-1];
                QS20[1,1,i] = foo[2,i]-foo[2,i-1];
                QS20[1,2,i] = foo2[2,i]-foo2[2,i-1];
                QS20[2,1,i] = foo[3,i]-foo[3,i-1];
                QS20[2,2,i] = foo2[3,i]-foo2[3,i-1];
            }
            
        }

        
        { ## now plot it
            
            pdf(paste(FIGDIR,"SUPP/S20_h2_vs_GERP.pdf",sep=""),height=2,width=4,pointsize=10)
            par(mar=c(2,2.8,0.4,0.3))
            plot(1,1,type='n',ylim=range(-.005,0.025),xlim=range(GERP,-8,8,na.rm=T),xaxt='n',yaxt='n',xlab="",ylab="")
            abline(h=seq(0,0.1,by=0.005),col=rgb(0,0,0,0.1),lwd=0.75)
            axis(side=1,padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3, cex.axis=0.7)
            
            mtext(side=1,"GERP score",line=1)
            mtext(side=2,expression(paste({widehat(h^2)[bin]},sep="")),line=1.15)

            lines(GERP[1,],QS20[3,1,])
            lines(GERP[2,],QS20[3,2,],col='blue')
            points(GERP[1,],QS20[3,1,],pch=20)
            points(GERP[2,],QS20[3,2,],col='blue',pch=20)

            for(i in 1:20){
                segments(x0=GERP[1,i],y0=QS20[1,1,i],y1=QS20[2,1,i],col='black')
                segments(x0=GERP[2,i],y0=QS20[1,2,i],y1=QS20[2,2,i],col='blue')
            }
            
            legend("topright",c("Singletons Only","All Variants"), col=c("black","blue"), lty=1,lwd=2,ncol=1,bty='n',pch=20,cex=0.7)

            dev.off();
        }
    }
    
    
    { ## Fig S21: permutation figure with just singletons

        files = c("merge_1000000_K=20_MAC=1_ngPC=0_npPC=0_QN=1_DSAMP=1_PERM=1_NORESID.he.h2.gz", "merge_1000000_K=20_MAC=1_ngPC=0_npPC=0_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz", "merge_1000000_K=20_MAC=1_ngPC=1_npPC=1_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz", "merge_1000000_K=20_MAC=1_ngPC=5_npPC=5_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz", "merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz", "merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.h2.gz");
        QS21 = array(,dim=c(3,length(files)));
        for(i in 1:length(files)){
            foo=read.table(files[i])
            foo=foo[which(foo[,1] %in% expGenes[,1]),3]
            QS21[3,i] = mean(foo,na.rm=T);
            NBOOT=1000;
            B = array(,dim=NBOOT);
            for(j in 1:NBOOT){
                B[j] = mean(sample(foo,replace=T),na.rm=T)
            }
            QS21[1,i] = quantile(B,0.025)
            QS21[2,i] = quantile(B,0.975)
        }
        
        
        pdf(paste(FIGDIR,"SUPP/S21_h2_vs_permType.pdf",sep=""),height=2,width=2.25,pointsize=10)
        par(mar=c(6.8,2.8,0.4,0.3))
        plot(1,1,type='n',ylim=range(-.002,0.022),xlim=range(1:6,na.rm=T),xaxt='n',yaxt='n',xlab="",ylab="")
        abline(h=seq(0,0.1,by=0.005),col=rgb(0,0,0,0.1),lwd=0.75)
        axis(side=1,at=1:6,labels=c("Random Perm","trans PC=0","trans PC=1","trans PC=5","trans PC=10","Data"),las=3,cex=0.7)
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,0.02,by=0.005),labels=NA)
        axis(side=2,padj=1.3, cex.axis=0.7,at=c(0,0.01,0.02),labels=c(0,0.01,0.02))
        
        mtext(side=2,expression(paste({widehat(h^2)[singletons]},sep="")),line=1.15)

        points(1:length(files),QS21[3,],pch=20);
        for(i in 1:length(files)){
            segments(x0=i,y0=QS21[1,i],y1=QS21[2,i])
        }
        dev.off()
    }
    
    { ## fig S22: h2_singleton for various CGI subsets
        if(1){ ##get data
            K=4;
            QS22 = array(,dim=c(3,K)); #{lowQ,highQ,mean}{bins}
            file="merge_1000000_K=4_MAC=1_ngPC=0_npPC=0_QN=1_DSAMP=1_CGIsub.he.h2.gz"
            zfile = gzfile(file,"r");
            foo=scan(zfile,nlines=1,what="")
            close(zfile);
            zfile = gzfile(file,"r");
            d=matrix(scan(zfile,what=""),ncol=length(foo),byrow=T)
            close(zfile);
            d2=matrix(as.numeric(d[,1:K+2]),ncol=K);
            QS22[3,] = apply(d2,2,function(x){mean(x,na.rm=T)});

            ##get avg number of SNPs per bin
            M = apply(d[,10+1:K],2,function(x){round(mean(as.numeric(x),na.rm=T))});
            
            ##now do boot straps
            NBOOT=1000;
            B = array(,dim=c(NBOOT,K));
            for(i in 1:NBOOT){
                b = d2[sample(1:nrow(d2),replace=T),]
                B[i,] = apply(b,2,function(x){mean(x,na.rm=T)});
            }
            for(j in 1:K){
                QS22[1,j] = quantile(B[,j],0.025);
                QS22[2,j] = quantile(B[,j],0.975);
            }
        }
        
        plotOrder=c(2,1,4,3);
        
        pdf(paste(FIGDIR,"SUPP/S22_h2_vs_CGIsub_singletons.pdf",sep=""),height=2.5,width=3,pointsize=10)
        par(mar=c(5,3,.5,.5))
        plot(1,1,type='n',xlim=c(0.75,4.25),ylim=c(-0.01,0.0325),xlab="",ylab="",xaxt='n', yaxt='n')
        abline(h=0,col=rgb(0,0,0,0.1))
        LABELS=c(paste("+CGI\n+gnomAD\n(n=",M[1],")",sep=""),paste("+CGI\n-gnomAD\n(n=",M[2],")",sep=""),paste("-CGI\n+gnomAD\n(n=",M[3],")",sep=""),paste("-CGI\n-gnomAD\n(n=",M[4],")",sep=""))
        axis(side=1,at=1:4,labels=LABELS[plotOrder],las=3,cex=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=2,expression(paste({widehat(h^2)[singletons]},sep="")),line=1.15)
        abline(h=0,col=rgb(0,0,0,0.1))
        for(i in 1:4){
            segments(x0=i,y0=QS22[1,plotOrder[i]], y1=QS22[2,plotOrder[i]], col="black")
            points(i,QS22[3,plotOrder[i]],col="black", pch=16)
        }
        dev.off()
    }

    { #fig S25: performance of rejection sampler
        GETDATA=0;
        if(GETDATA){
            read.table(paste(FIGDIR,"/LawrenceData/lsq_all.txt",sep=""))->t
            read.table(paste(FIGDIR,"LawrenceData/real_infer.5000.txt",sep=""))->r

            ##tau: true=t[,2]; inf=t[,8]
            ##rho: true=t[,1]; inf=t[,7]
            ##p: true=t[,3]; inf=t[,9]

            BREAKS = 50;
            Bx = seq(0, 1.001, length=BREAKS+1)
            By = seq(0, 1.001, length=BREAKS+1)

            TAU=matrix(,nr=BREAKS,nc=BREAKS);
            RHO=matrix(,nr=BREAKS,nc=BREAKS);
            P=matrix(,nr=BREAKS,nc=BREAKS);
            JOINT=matrix(,nr=BREAKS,nc=BREAKS);

            for(i in 1:BREAKS){
                cat("pulling out i=",i,"\n");
                for(j in 1:BREAKS){
                    RHO[i,j] = sum(t[,1]>=Bx[i] & t[,1]<Bx[i+1] & t[,7]>=By[j] & t[,7]<By[j+1])
                    TAU[i,j] = sum(t[,2]>=Bx[i] & t[,2]<Bx[i+1] & t[,8]>=By[j] & t[,8]<By[j+1])
                    P[i,j] = sum(t[,3]>=Bx[i] & t[,3]<Bx[i+1] & t[,9]>=By[j] & t[,9]<By[j+1])
                    JOINT[i,j] = sum(r[,1]>=Bx[i] & r[,1]<Bx[i+1] & r[,2]>=By[j] & r[,2]<By[j+1])
                }
            }
        }
        
        pdf(paste(FIGDIR,"SUPP/S25A.pdf",sep=""),height=1.5,width=4.75/3,pointsize=10)
        par(mar=c(2,2,.5,2))
        image(Bx[1:BREAKS], By[1:BREAKS], log(RHO),xlim=c(0,1),ylim=c(0,1),
              col = viridis(BREAKS, option = "C"),
              xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
        abline(0,1,col=rgb(0,0,0,0.1))

        axis(side=1,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
        axis(side=1,padj=-2,at=c(0,0.5,1.0),cex.axis=0.7)
        axis(side=2,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
        axis(side=2,padj=1.3,at=c(0,0.5,1.0), cex.axis=0.7)
        mtext(side=1,expression("True " * rho),line=1)
        mtext(side=2,expression("Inferred " * rho),line=.9)
        dev.off()

        pdf(paste(FIGDIR,"SUPP/S25Aleg.pdf",sep=""),height=1.5,width=0.4,pointsize=10)
        par(mar=c(3,.6,1,.6),oma=c(0,0,0,0))
        tcol = rev(viridis(5,option="C"))
        image(t(matrix(seq(1,max(RHO),length=BREAKS),nr=BREAKS,nc=1)),
              col = viridis(BREAKS, option = "C"),
              xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
        ticks = c(1,5,25,125,425)
        text(x=0,y=0,labels=ticks[1],cex=0.7,col="black",adj=c(0.5,0))
        text(x=0,y=seq(0,1,length=5)[2:4],labels=ticks[2:4],cex=0.7,col="black",adj=c(0.5,0.5))
        text(x=0,y=1,labels=ticks[5],cex=0.7,col="black",adj=c(0.5,1))
        mtext(side=3,"Count",cex=0.7)
        dev.off()

        pdf(paste(FIGDIR,"SUPP/S25B.pdf",sep=""),height=1.5,width=4.75/3,pointsize=10)
        par(mar=c(2,2,.5,2))
        image(Bx[1:BREAKS], By[1:BREAKS], log(TAU),xlim=c(0,1),ylim=c(0,1),
              col = viridis(BREAKS, option = "C"),
              xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
        abline(0,1,col=rgb(0,0,0,0.1))

        axis(side=1,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
        axis(side=1,padj=-2,at=c(0,0.5,1.0),cex.axis=0.7)
        axis(side=2,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
        axis(side=2,padj=1.3,at=c(0,0.5,1.0), cex.axis=0.7)
        mtext(side=1,expression("True " * tau),line=1)
        mtext(side=2,expression(phantom(rho) * "Inferred " * tau * phantom(rho)),line=.9)
        dev.off()

        pdf(paste(FIGDIR,"SUPP/S25Bleg.pdf",sep=""),height=1.5,width=0.4,pointsize=10)
        par(mar=c(3,.6,1,.6),oma=c(0,0,0,0))
        tcol = rev(viridis(5,option="C"))
        image(t(matrix(seq(1,max(TAU),length=BREAKS),nr=BREAKS,nc=1)),
              col = viridis(BREAKS, option = "C"),
              xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
        ticks = c(1,5,25,125,625)
        text(x=0,y=0,labels=ticks[1],cex=0.7,col="black",adj=c(0.5,0))
        text(x=0,y=seq(0,1,length=5)[2:4],labels=ticks[2:4],cex=0.7,col="black",adj=c(0.5,0.5))
        text(x=0,y=1,labels=ticks[5],cex=0.7,col="black",adj=c(0.5,1))
        mtext(side=3,"Count",cex=0.7)
        dev.off()

        pdf(paste(FIGDIR,"SUPP/S25C.pdf",sep=""),height=1.5,width=4.75/3,pointsize=10)
        par(mar=c(2,2,.5,2))
        image(Bx[1:BREAKS], By[1:BREAKS], log(P),xlim=c(0,1),ylim=c(0,1),
              col = viridis(BREAKS, option = "C"),
              xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
        abline(0,1,col=rgb(0,0,0,0.1))

        axis(side=1,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
        axis(side=1,padj=-2,at=c(0,0.5,1.0),cex.axis=0.7)
        axis(side=2,padj=-2,at=seq(0,1,by=.1),cex.axis=0.7,labels=NA,tck=-0.04)
        axis(side=2,padj=1.3,at=c(0,0.5,1.0), cex.axis=0.7)
        mtext(side=1,expression("True " * phi),line=1)
        mtext(side=2,expression("Inferred " * phi),line=.9)
        dev.off()

        pdf(paste(FIGDIR,"SUPP/S25Cleg.pdf",sep=""),height=1.5,width=0.4,pointsize=10)
        par(mar=c(3,.6,1,.6),oma=c(0,0,0,0))
        tcol = rev(viridis(5,option="C"))
        image(t(matrix(seq(1,max(P),length=BREAKS),nr=BREAKS,nc=1)),
              col = viridis(BREAKS, option = "C"),
              xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
        ticks = c(1,4,16,64,256)
        text(x=0,y=0,labels=ticks[1],cex=0.7,col="black",adj=c(0.5,0))
        text(x=0,y=seq(0,1,length=5)[2:4],labels=ticks[2:4],cex=0.7,col="black",adj=c(0.5,0.5))
        text(x=0,y=1,labels=ticks[5],cex=0.7,col="black",adj=c(0.5,1))
        mtext(side=3,"Count",cex=0.7)
        dev.off()
    }

    { #Fig S27: "true singleton MAF
        
        SIMD=matrix(scan(file="/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures2/Lawrence/h2_by_tf.postEst.txt"),byrow=T,nc=2);
        SIMD = SIMD[which(SIMD[,1]<=0.5),]

        K=20;
        
        MAFS27=array(,dim=c(2,K))
        foo=read.table(paste("SINGgnomAD_SFS_K=",K,".txt",sep=""));
        MAFS27[1,1:K] = foo[,1];
        MAFS27[1,1] = 1e-5;

        foo=read.table(paste("SINGtgp_SFS_K=",K,".txt",sep=""));
        MAFS27[2,1:K] = foo[,1];

        heh2=array(,dim=c(2,K));
        ##now get gnomAD partition data
        filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_SINGgnomAD.he.meanQuant"
        foo=read.table(filename);
        h2 = foo[1,1:K];
        ##transform to h2 per bin
        for(j in K:2){
            h2[j] = h2[j] - h2[j-1];
        }
        ##now get permuted values
        filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_SINGgnomAD.he.meanQuant";
        foo=read.table(filename);
        h2p = foo[1,1:K];
        ##transform to h2 per bin
        for(j in K:2){
            h2p[j] = h2p[j] - h2p[j-1];
        }
        heh2[1,1:K] = unlist(h2-h2p)

        ##now get gnomAD partition data
        filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_SINGtgp.he.meanQuant"
        foo=read.table(filename);
        h2 = foo[1,1:K];
        ##transform to h2 per bin
        for(j in K:2){
            h2[j] = h2[j] - h2[j-1];
        }
        ##now get permuted values
        filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_SINGtgp.he.meanQuant";
        foo=read.table(filename);
        h2p = foo[1,1:K];
        ##transform to h2 per bin
        for(j in K:2){
            h2p[j] = h2p[j] - h2p[j-1];
        }
        heh2[2,1:K] = unlist(h2-h2p)

        
        pdf(paste(FIGDIR,"SUPP/S27_h2sing_MAF_TGP_gnomAD_Sim_K=",K,"_PC=",gPC,".pdf",sep=""),height=2,width=2,pointsize=10)
        par(mar=c(2,2.5,1,.5))
        plot(1,1,type='n',xlim=c(1e-6,0.5),ylim=range(1e-9,1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
        abline(h=seq(0,1,by=0.1),col=rgb(0,0,0,0.1),lwd=0.75)
        for(i in 1:2){
            lines(MAFS27[i,],heh2[i,1:K]/sum(heh2[i,1:K],na.rm=T),col=pcol[i])
            points(MAFS27[i,],heh2[i,1:K]/sum(heh2[i,1:K],na.rm=T),col=pcol[i],pch=20)
            ##segments(x0=MAFS[i,],y0=Ql[i,], y1=Qh[i,], col=pcol[i])
        }
        lines(SIMD[,1],(SIMD[,2]),col=pcol[9],lwd=2)

        axis(side=1,padj=-1.75,cex.axis=0.7,at=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c("1e-6","1e-5","1e-4","0.001","0.01","0.1","0.5"))
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.1,0.9,by=0.1),labels=NA,tck=-0.04)
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,1,by=0.5),labels=c("0","0.5","1"))
        mtext(side=1,"\"True\" Singleton MAF",line=1)
        mtext(side=2,expression(paste("Prop. of  ",{widehat(h^2)[singleton]},sep="")),line=1.1)
        legend("topright",c(expression("TGP"),expression("gnomAD"),"Model"), col=pcol[c(2:1,9)], bg=rgb(1,1,1,0.8),pch=20,box.col=rgb(1,1,1,0.8), cex=0.8)
        box();
        dev.off()
    }
}
    

    

    
    {
        
        pdf(paste(FIGDIR,"SUPP/2A_violin_h2tot_obsONLY_K=",K,"_PC=",gPC,".pdf",sep=""),height=1.5,width=2.375,pointsize=10)
        par(mar=c(2,2.8,0.3,0.3))
        plot(1,1,type='n',ylim=range(0,0.1),xlim=c(2.66,6.34),xaxt='n',yaxt='n',xlab="",ylab="")
        abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
        axis(side=1,at=3:6,labels=c("1e-3","3e-3","7e-3","0.05"),padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3, cex.axis=0.7,at=c(1:4/100,6:9/100),labels=NA,tck=-0.04)
        axis(side=2,padj=1.3, cex.axis=0.7,at=c(0,0.05,0.1),labels=c("0","0.05","0.1"))
        mtext(side=1,"Minimum MAF",line=1)
        mtext(side=2,expression({widehat(h^2)}[total]),line=1.15)

        legend("bottom",c(expression(paste(MAC >= "1")),expression(paste(MAC >= "2")),expression(paste(MAC >= "5")),expression(paste(MAF >= "5%"))),col=pcol[3:6],lty=1, cex=0.8,bg=rgb(1,1,1,0.8),lwd=2,box.col=rgb(1,1,1,0.2),ncol=2)
        box()
        for(i in 3:6){
            d = density(QMAT[i,1,])
            pci = pcol[i]
            ri = col2rgb(pci)/255;
            polygon(x=i+c(d$y,-rev(d$y))/max(d$y)/3, y=c(d$x,rev(d$x)), col="white",border=NA)
            polygon(x=i+c(d$y,-rev(d$y))/max(d$y)/3, y=c(d$x,rev(d$x)), col=rgb(ri[1],ri[2],ri[3],0.1),border=pci)
        }
        dev.off()
    }

    {
        pdf(paste(FIGDIR,"2A_points_h2tot_obsONLY_K=",K,"_PC=",gPC,".pdf",sep=""),height=1.5,width=2.375,pointsize=10)
        par(mar=c(1.5,2.8,0.3,0.3))
        plot(1,1,type='n',ylim=range(0.04,0.1),xlim=c(2.66,6.34),xaxt='n',yaxt='n',xlab="",ylab="")
        abline(h=seq(0,0.1,by=0.01),col=rgb(0,0,0,0.1),lwd=0.75)
        axis(side=1,at=3:6,labels=c(expression(paste(MAC >= "1")),expression(paste(MAC >= "2")),expression(paste(MAC >= "5")),expression(paste(MAF >= "5%"))),padj=-1.5, cex.axis=0.7)
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.05,0.09,by=0.02),labels=NA,tck=-0.04)
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.04,0.1,by=0.02),labels=c("0.04","0.06","0.08","0.1"))
        mtext(side=2,expression({widehat(h^2)}[total]),line=1.15)

        for(i in 3:6){
            d = density(QMAT[i,1,])
            pci = pcol[i]
            ri = col2rgb(pci)/255;
            #polygon(x=i+c(d$y,-rev(d$y))/max(d$y)/3, y=c(d$x,rev(d$x)), col="white",border=NA)
            #polygon(x=i+c(d$y,-rev(d$y))/max(d$y)/3, y=c(d$x,rev(d$x)), col=rgb(ri[1],ri[2],ri[3],0.1),border=pci)
                        H = quantile(QMAT[i,1,],0.995)
            L = quantile(QMAT[i,1,],0.005)
            pci = pcol[i]
            ri = col2rgb(pci)/255;
            points(i,sum(H2MATS[i,1+1:K]),col=pcol[i],pch=20)
            segments(x0=i,y0=L,y1=H,col=pcol[i])
        }
        dev.off()
    }
    
    { ##proportion of h2 per bin
        pdf(paste(FIGDIR,"2B_h2_propByMAC_NORESID_gPC=",gPC,"_pPC=",pPC,"_QN=",QN,CVstring,".pdf",sep=""),height=1.5,width=2.375,pointsize=10)
        par(mar=c(2,2.8,0.3,0.3))
        plot(1,1, log='x', xlim=c(0.001,0.5), ylim=range(0,hq,lq,0.3), xaxt='n', yaxt='n',xlab="",ylab="",main="",type='n')
        abline(h=seq(0,0.3,by=0.05),col=rgb(0,0,0,0.1),lwd=0.75)
        hq = 0;
        lq = 0;
        p = H2MATS[3,1+1:K]/sum(H2MATS[3,1+1:K])
        lines(MAF[3,1+1:K], p)
        for(i in 1:K){
            hq[i] = quantile(QMAT[3,1+i,]/QMAT[3,1,],0.995)
            lq[i] = quantile(QMAT[3,1+i,]/QMAT[3,1,],0.005)
            points(MAF[3,1+i], p[i],pch=20)
            segments(x0=MAF[3,1+i], y0=lq[i],y1=hq[i])
        }

        axis(side=1,padj=-2,labels=NA, cex.axis=0.7)
        axis(side=1,padj=-2,at=c(0.001,0.005,0.02,0.1,0.5),labels=c("0.001","0.005","0.02","0.1","0.5"),cex.axis=0.7)

        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.05,0.25,by=0.1),labels=NA,tck=-0.04)
        axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,0.3,by=0.1),labels=c("0","0.1","0.2","0.3"))

        mtext(side=1,"Minor Allele Frequency",line=1)
        mtext(side=2,expression(paste("Prop. of  ",{widehat(h^2)[total]},sep="")),line=1.15)
        dev.off()
    }

    { # 2C: proportion of h2 for singletons partitioned by gnomAD/TGP freq
        if(0){#get data
            PART=c("SINGgnomAD", "SINGtgp");
            K = 20
            W="1000000";
            PC=10;
            QN=1;
            MAFS=matrix(,nr=length(PART),nc=K);

            heh2 = matrix(,nr=length(PART),nc=2*K+1);
            heh2se = matrix(,nr=length(PART),nc=2*K+1);
            NBOOT = 1000;
            Qprop=array(,dim=c(2,K,NBOOT));
            Qh = matrix(,nr=2,nc=K)
            cumH = matrix(,nr=2,nc=K)
            cumL = matrix(,nr=2,nc=K);
            Ql = matrix(,nr=2,nc=K)
            Qx = matrix(,nr=length(PART), nc=K);
            Xlim = matrix(,nr=length(PART), nc=2);
            Ylim = matrix(,nr=length(PART), nc=2);
            CNT=0;
            for(part in PART){
                CNT=CNT+1;
                
                foo=read.table(paste(part,"_SFS_K=",K,".txt",sep=""));
                MAFS[CNT,1:K] = foo[,1];
                if(MAFS[CNT,1] < 1e-5){
                    MAFS[CNT,1] = 1e-5;
                }
                filename=paste("merge_",W,"_K=",K,"_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=",QN,"_DSAMP=1_",part,".he.h2.gz",sep="");
                GZIN=gzfile(filename,"r");
                foo=scan(file=GZIN,nlines=1,what="",quiet=T)
                close(GZIN);
                GZIN=gzfile(filename,"r");
                he=matrix(scan(file=GZIN,what="",quiet=T),byrow=T,nc=length(foo))
                close(GZIN);
                G=which(he[,1] %in% expGenes[,1])
                he = he[G,]
                heh2[CNT,1:(K+1)] = apply(he[,1:(K+1)+1],2,function(x){mean(as.numeric(x))});
                heh2[CNT,K+1+1:K] = apply(he[,2*K+2+1:K],2,function(x){mean(as.numeric(x))});
                for(i in 1:NBOOT){
                    d=sample(1:nrow(he),nrow(he),replace=T);
                    h=apply(he[d,2+1:K],2,function(x){mean(as.numeric(x),na.rm=T)});
                    for(j in 1:K){
                        Qprop[CNT,j,i] = h[j]/sum(h,na.rm=T);
                    }
                }
                for(j in 1:K){
                    Qh[CNT,j] = quantile(Qprop[CNT,j,],0.995);
                    Ql[CNT,j] = quantile(Qprop[CNT,j,],0.005);
                }
                heh2se[CNT,1:(K+1)] = apply(he[,1:(K+1)+1],2,function(x){x=as.numeric(x); sd(x)/sqrt(length(which(!is.na(x))))});
                heh2se[CNT,K+1+1:K] = apply(he[,2*K+2+1:K],2,function(x){x=as.numeric(x); sd(x)/sqrt(length(which(!is.na(x))))});
                
                
                Ylim[CNT,] = range(Ylim[CNT,], heh2[CNT,1:K+1]+2*heh2se[CNT,1:K+1],heh2[CNT,1:K+1]-2*heh2se[CNT,1:K+1],na.rm=T);
                
                foo=read.table(paste(part,"_SFS_K=",K,".txt",sep=""));
                Qx[CNT,1:K] = foo[,1];
                if(Qx[CNT,1] == 0){
                    Qx[CNT,1] = 1e-05;
                }
                else if(Qx[CNT,1] == -100){
                    Qx[CNT,1] = NA;
                }
                Xlim[CNT,]=range(Xlim[CNT,], Qx[CNT,1:K], na.rm=T);
            }
            
        }

        {
            pdf(paste(FIGDIR,"2C_prop_h2tot_GEUsing_K=",K,"_PC=",gPC,".pdf",sep=""),height=2,width=2,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=c(-.01,1),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            ##abline(h=0,col=rgb(0,0,0,0.1),lwd=0.75)
            abline(h=seq(0,1,by=0.1),col=rgb(0,0,0,0.1),lwd=0.75)
                                        #abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
            
            legend("topright",c(expression("TGP"),expression("gnomAD")), col=pcol[2:1], bg=rgb(1,1,1,0.8),pch=20,box.col=rgb(1,1,1,0.8), cex=0.8)
            box()
            
            for(i in 1:2){
                lines(MAFS[i,],heh2[i,1+1:K]/sum(heh2[i,1+1:K]),col=pcol[i])
                points(MAFS[i,],heh2[i,1+1:K]/sum(heh2[i,1+1:K]),col=pcol[i],pch=20)
                segments(x0=MAFS[i,],y0=Ql[i,], y1=Qh[i,], col=pcol[i])
            }

            axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=1,at=c(1e-3,0.5),labels=c("1e-3","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0.1,0.9,by=0.2),labels=NA,tck=-0.04)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Singleton Global MAF",line=1)
            mtext(side=2,expression(paste("Prop. of  ",{widehat(h^2)[singleton]},sep="")),line=1.1)
            dev.off()
        }
    }  

    {
        pdf(paste(FIGDIR,"2D_cum_h2tot_obsVtgpVgnomAD_K=",K,"_PC=",gPC,".pdf",sep=""),height=2,width=2.75,pointsize=10)
        par(mar=c(2,2.8,0.3,0.3))
        plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=c(0,0.22),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
        abline(h=seq(0,0.25,by=0.05),col=rgb(0,0,0,0.1),lwd=0.75)
        abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)

        legend("topleft",c("gnomAD","TGP",expression("MAC" >= 1),expression("MAC" >= 2),expression("MAC" >= 5),expression("MAF" >= "5%")), col=pcol, lty=1,lwd=2,bg=rgb(1,1,1,0.8),ncol=2,box.col="white",pch=20,cex=0.8)
        box()
        for(i in 1:6){
            pci = pcol[i]
            rgbi = col2rgb(pci)/255;
            hq = 0;
            lq = 0;
            CQ = apply(t(QMAT[i,1+1:K,]),1,function(x){cumsum(x)})

            for(j in 1:K){
                hq[j] = quantile(CQ[j,],0.995)
                lq[j] = quantile(CQ[j,],0.005)
            }
            polygon(x=c(MAF[i,1+1:K], rev(MAF[i,1+1:K])), y=c(hq,rev(lq)), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.2),border=NA)
            lines(MAF[i,1+1:K],cumsum(H2MATS[i,1+1:K]),col=pcol[i])
            points(MAF[i,1+1:K],cumsum(H2MATS[i,1+1:K]),col=pcol[i],pch=20)
        }

        axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=1,"Minor Allele Frequency",line=1)
        mtext(side=2,expression(paste("Cumulative ",{widehat(h^2)},sep="")),line=1.15)
        dev.off()
    }
    
    {
        pdf(paste(FIGDIR,"SUPP/2A_violin_h2tot_obsVtgpVgnomAD_K=",K,"_PC=",gPC,".pdf",sep=""),height=2,width=2.375,pointsize=10)
        par(mar=c(2.3,2.8,0.5,0.3))
        plot(1,1,type='n',ylim=range(0,QMAT[i,1,],0.2),xlim=c(0.66,6.34),xaxt='n',yaxt='n',xlab="",ylab="")
        abline(h=seq(0,0.2,by=0.05),col=rgb(0,0,0,0.1),lwd=0.75)
        abline(v=1:6,col=rgb(0,0,0,0.1),lwd=0.75)
        axis(side=1,at=1:6,labels=c(paste("\"0\""),"2e-4","1e-3","3e-3","7e-3","0.05"),padj=-2, cex.axis=0.7)
        ##axis(side=1,at=c(3,5),labels=c("1.4e-3","7e-3"),padj=-1.5, cex.axis=0.7,tick="NA")
        axis(side=2,padj=1.3, cex.axis=0.8)
        mtext(side=1,"Minimum MAF",line=1.15)
        mtext(side=2,expression({widehat(h^2)}[total]),line=1.15)

        legend("topright",c(expression(gnomAD["n=15k"]),expression(TGP["n=2505"]),expression(paste(MAC["n=360"] >= "1")),expression(paste(MAC["n=360"] >= "2")),expression(paste(MAC["n=360"] >= "5")),expression(paste(MAF["n=360"] >= "5%"))),col=pcol,lty=1, cex=0.7,bg=rgb(1,1,1,0.8),lwd=2,box.col="white")
        box()
        for(i in 1:6){
            d = density(QMAT[i,1,])
            pci = pcol[i]
            ri = col2rgb(pci)/255;
            polygon(x=i+c(d$y,-rev(d$y))/max(d$y)/3, y=c(d$x,rev(d$x)), col="white",border=NA)
            polygon(x=i+c(d$y,-rev(d$y))/max(d$y)/3, y=c(d$x,rev(d$x)), col=rgb(ri[1],ri[2],ri[3],0.1),border=pci)
        }
        dev.off()
    }

    
    
    
    { ##ridgeline plot
        pdf(paste(FIGDIR,"2A_ridgeline_h2tot_K=",K,"_PC=",gPC,".pdf",sep=""),height=2,width=2.375,pointsize=10)
        par(mar=c(2.5,0.4,0.5,0.4))
        plot(1,1,type='n',ylim=c(0,3),xlim=range(-.17,QMAT[i,1,],0.2),xaxt='n',yaxt='n',xlab="",ylab="",bty="n")
        axis(side=1,padj=-2.5, cex.axis=0.8, at=c(0,0.05,0.1,0.15,0.2),tick=F)
        mtext(side=1,expression({widehat(h^2)}[total]),line=1.5,adj=0.8)
        ##axis(side=2,at=0:5/2+.2,las=1,labels=c(expression(paste(MAF[n=360] >= "5%")),expression(paste(MAC[n=360] >= "5")),expression(paste(MAC[n=360] >= "2")),expression(paste(MAC[n=360] >= "1")),expression(TGP["n=2505"]),expression(gnomAD["n=15k"])),tick=F, hadj=0.2)
        text(x=0,y=0:5/2+.15,las=1,labels=c(expression(paste(MAF["n=360"] >= "5%")),expression(paste(MAC["n=360"] >= "5")),expression(paste(MAC["n=360"] >= "2")),expression(paste(MAC["n=360"] >= "1")),expression(TGP["n=2505"]),expression(gnomAD["n=15k"])), pos=2)
        abline(h=0:5/2,col=rgb(0,0,0,0.1))
        abline(v=seq(0,.2,by=0.05),col=rgb(0,0,0,0.1))
        for(i in 1:6){
            d = density(QMAT[i,1,],adjust=0.5)
            pci = pcol[i]
            ri = col2rgb(pci)/255;
            polygon(x=d$x,y=(3-i/2)+d$y/max(d$y)/2.25, col="white",border=NA)
            polygon(x=d$x,y=(3-i/2)+d$y/max(d$y)/2.25, col=rgb(ri[1],ri[2],ri[3],0.1),border="NA")
            lines(d$x,(3-i/2)+d$y/max(d$y)/2.25,col=pcol[i],lwd=0.75)
        }
        dev.off()
    }
}

{
    { #compare genome-wide MAC GRMs
        FREQ=c("1-1","2-2","3-7","8-360"); #SAMPLE
        FREQT=c("1-1","2-50","51-2504"); #TGP
        jpeg(paste(FIGDIR,"SUPP/compGRMs_MAC_SAMPLEvsTGP.jpg",sep=""),height=7,width=7,pointsize=10,units="in",res=300)
        par(mfrow=c(length(FREQ),length(FREQT)),mar=c(2,2,.5,.5))
        for(i in 1:length(FREQ)){
            for(j in 1:length(FREQT)){
                c1 = -22;
                foo=gzfile(paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/GRM/GRM_chr",c1,"_MAC=",FREQ[i],"_SAMPLE.txt.gz",sep=""),"r");
                d1 = read.table(foo)
                d1 = d1[upper.tri(d1)]
                close(foo);
                
                c2=c1;
                foo=gzfile(paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/GRM/GRM_chr",c2,"_MAC=",FREQT[j],"_TGP.txt.gz",sep=""),"r");
                d2 = read.table(foo)
                d2 = d2[upper.tri(d2)]
                close(foo);
                plot(1,1,type='n',xlim=range(d1),ylim=range(d2),xaxt='n',yaxt='n',xlab="",ylab="")
                abline(0,1)
                axis(side=1,padj=-2, cex.axis=0.7)
                axis(side=2,padj=1.3, cex.axis=0.7)
                mtext(side=1,paste("MAC=",FREQ[i],sep=""),line=1.15)
                mtext(side=2,paste("MAC=",FREQT[j],sep=""),line=1.15)
                points(unlist(d1),unlist(d2),pch=20,col=rgb(0,0,0,0.25))
            }
        }
        dev.off();   
    }

    { #compare genome-wide MAC GRMs
        FREQ=c("1-1","2-2","3-7","8-360"); #SAMPLE
                                        #FREQT=c("1-1","2-50","51-2504"); #TGP
        jpeg(paste(FIGDIR,"SUPP/compGRMs_MAC_SAMPLE.jpg",sep=""),height=7,width=7,pointsize=10,units="in",res=300)
        par(mfrow=c(length(FREQ),length(FREQ)),mar=c(2,2,.5,.5))
        for(i in 1:length(FREQ)){
            for(j in 1:length(FREQ)){
                c1 = -22;
                foo=gzfile(paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/GRM/GRM_chr",c1,"_MAC=",FREQ[i],"_SAMPLE.txt.gz",sep=""),"r");
                d1 = read.table(foo)
                d1 = d1[upper.tri(d1)]
                close(foo);
                
                c2=c1;
                if(i == j){
                    c2 = -1;
                }
                foo=gzfile(paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/GRM/GRM_chr",c2,"_MAC=",FREQ[j],"_SAMPLE.txt.gz",sep=""),"r");
                d2 = read.table(foo)
                d2 = d2[upper.tri(d2)]
                close(foo);
                plot(1,1,type='n',xlim=range(d1),ylim=range(d2),xaxt='n',yaxt='n',xlab="",ylab="")
                abline(0,1)
                axis(side=1,padj=-2, cex.axis=0.7)
                axis(side=2,padj=1.3, cex.axis=0.7)
                mtext(side=1,paste("MAC=",FREQ[i],sep=""),line=1.15)
                mtext(side=2,paste("MAC=",FREQ[j],sep=""),line=1.15)
                points(unlist(d1),unlist(d2),pch=20,col=rgb(0,0,0,0.25))
            }
        }
        dev.off();   
    }

    { #compare genome-wide MAC GRMs
                                        #FREQ=c("1-1","2-2","3-7","8-360"); #SAMPLE
        FREQ=c("1-1","2-50","51-2504"); #TGP
        jpeg(paste(FIGDIR,"SUPP/compGRMs_MAC_TGP.jpg",sep=""),height=7,width=7,pointsize=10,units="in",res=300)
        par(mfrow=c(length(FREQ),length(FREQ)),mar=c(2,2,.5,.5))
        for(i in 1:length(FREQ)){
            for(j in 1:length(FREQ)){
                c1 = -22;
                foo=gzfile(paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/GRM/GRM_chr",c1,"_MAC=",FREQ[i],"_TGP.txt.gz",sep=""),"r");
                d1 = read.table(foo)
                d1 = d1[upper.tri(d1)]
                close(foo);
                
                c2=c1;
                if(i == j){
                    c2 = -1;
                }
                foo=gzfile(paste("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/GRM/GRM_chr",c2,"_MAC=",FREQ[j],"_TGP.txt.gz",sep=""),"r");
                d2 = read.table(foo)
                d2 = d2[upper.tri(d2)]
                close(foo);
                plot(1,1,type='n',xlim=range(d1),ylim=range(d2),xaxt='n',yaxt='n',xlab="",ylab="")
                abline(0,1)
                axis(side=1,padj=-2, cex.axis=0.7)
                axis(side=2,padj=1.3, cex.axis=0.7)
                mtext(side=1,paste("MAC=",FREQ[i],sep=""),line=1.15)
                mtext(side=2,paste("MAC=",FREQ[j],sep=""),line=1.15)
                points(unlist(d1),unlist(d2),pch=20,col=rgb(0,0,0,0.25))
            }
        }
        dev.off();   
    }
}


{
    setwd("~/Dropbox/Research/Data/GEUVADIS/INFW_HEvREML5/OUTSUM");
    FIGDIR="~/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures5/"
    GRMs = c("1-1_SAMPLE", "2-2_SAMPLE", "3-7_SAMPLE", "8-360_SAMPLE", "1-1_TGP", "2-50_TGP", "51-2504_TGP");
    PART = c("GEUVADISMAF", "SINGtgp", "SINGgnomAD", "ALLtgp", "ALLgnomAD");
    pcol=brewer.pal(9, "Set1")
    pcol[6] = "black";
    
    YLIM = matrix(,nr=c(length(GRMs)),nc=2)
    for(i in 1:length(GRMs)){
        YLIM[i,]=c(-.02,0.1);
    }
    {
        {
            for(CORTYPE in c("SAMPLE","TGP")){
                if(CORTYPE == "SAMPLE"){ 
                    GRM = c(1:4);
                }
                else{
                    GRM = c(5:7);
                }
                for(PC in c(0,10)){
                    pdf(paste(FIGDIR,"HEh2_GRMcor_PC=",PC,"_",CORTYPE,".pdf",sep=""),height=2*length(GRM),width=2*length(GRM),pointsize=10)
                    par(mfrow=c(length(GRM),length(GRM)),mar=c(1,.5,1,.5),oma=c(1,2.5,0.5,1));
                    panel = 0;
                    for(g1 in 1:length(GRM)){
                        grm1 = GRM[g1];
                        for(g2 in 1:length(GRM)){
                            panel = panel+1;
                            grm2 = GRM[g2];
                            if(grm1 == grm2){ #diagonal, just single GRM
                                heP=matrix(,nr=length(PART)+1,nc=20);
                                heT=matrix(,nr=length(PART)+1,nc=20);
                                for(part in 1:length(PART)){
                                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_GRM_",GRMs[grm1],"_",PART[part],".he.meanQuant",sep="")
                                    if(file.exists(FILE)){
                                        heP[part,] = scan(FILE,nlines=1)[1:20]
                                    }
                                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_GRM_",GRMs[grm1],"_",PART[part],".he.meanQuant",sep="")
                                    if(file.exists(FILE)){
                                        heT[part,] = scan(FILE,nlines=1)[1:20]
                                    }
                                }
                                FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="")
                                if(file.exists(FILE)){
                                    heT[nrow(heT),] = scan(FILE,nlines=1)[1:20]
                                }
                                FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="")
                                if(file.exists(FILE)){
                                    heP[nrow(heP),] = scan(FILE,nlines=1)[1:20]
                                }
                                
                                YLIM[grm1,] = c(min(YLIM[grm1,1],heP,heT,na.rm=T),max(YLIM[grm1,2],heP,heT,na.rm=T))
                                cat(YLIMD,"\n");
                                
                                plot(1,1,type='n',ylim=YLIM[grm1,],xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                                abline(h=0)
                                for(i in 1:nrow(heT)){
                                    points(1:20,heT[i,],pch=1,col=pcol[i])
                                    lines(1:20,heT[i,],col=pcol[i],lty=1)
                                    points(1:20,heP[i,],pch=2,col=pcol[i])
                                    lines(1:20,heP[i,],col=pcol[i],lty=2)
                                }
                            }
                            else if(grm1 < grm2){ #upper tri, plot both true and permutation h2
                                heP=matrix(,nr=length(PART)+1,nc=20);
                                heT=matrix(,nr=length(PART)+1,nc=20);
                                for(part in 1:length(PART)){
                                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_GRM_",GRMs[grm1],"_",GRMs[grm2],"_",PART[part],".he.meanQuant",sep="")
                                    if(file.exists(FILE)){
                                        heP[part,] = scan(FILE,nlines=1)[1:20]
                                    }
                                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_GRM_",GRMs[grm1],"_",GRMs[grm2],"_",PART[part],".he.meanQuant",sep="")
                                    if(file.exists(FILE)){
                                        heT[part,] = scan(FILE,nlines=1)[1:20]
                                    }
                                }
                                FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="")
                                if(file.exists(FILE)){
                                    heT[nrow(heT),] = scan(FILE,nlines=1)[1:20]
                                }
                                FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="")
                                if(file.exists(FILE)){
                                    heP[nrow(heP),] = scan(FILE,nlines=1)[1:20]
                                }
                                
                                YLIM[grm1,] = c(min(YLIM[grm1,1],heP,heT,na.rm=T),max(YLIM[grm1,2],heP,heT,na.rm=T))
                                
                                plot(1,1,type='n',ylim=YLIM[grm1,],xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                                abline(h=0)
                                for(i in 1:nrow(heT)){
                                    points(1:20,heT[i,],pch=1,col=pcol[i])
                                    lines(1:20,heT[i,],col=pcol[i],lty=1)
                                    points(1:20,heP[i,],pch=2,col=pcol[i])
                                    lines(1:20,heP[i,],col=pcol[i],lty=2)
                                }
                                
                            }
                            else{ #lower tri, plot perm-corrected h2
                                heP=matrix(,nr=length(PART)+1,nc=20);
                                heT=matrix(,nr=length(PART)+1,nc=20);
                                for(part in 1:length(PART)){
                                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_GRM_",GRMs[grm2],"_",GRMs[grm1],"_",PART[part],".he.meanQuant",sep="")
                                    if(file.exists(FILE)){
                                        heP[part,] = scan(FILE,nlines=1)[1:20]
                                    }
                                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_GRM_",GRMs[grm2],"_",GRMs[grm1],"_",PART[part],".he.meanQuant",sep="")
                                    if(file.exists(FILE)){
                                        heT[part,] = scan(FILE,nlines=1)[1:20]
                                    }
                                }
                                FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="")
                                if(file.exists(FILE)){
                                    heT[nrow(heT),] = scan(FILE,nlines=1)[1:20]
                                }
                                FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="")
                                if(file.exists(FILE)){
                                    heP[nrow(heP),] = scan(FILE,nlines=1)[1:20]
                                }
                                
                                YLIM[grm1,] = c(min(YLIM[grm1,1],heP,heT,na.rm=T),max(YLIM[grm1,2],heP,heT,na.rm=T))
                                
                                plot(1,1,type='n',ylim=YLIM[grm1,],xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")

                                abline(h=0)
                                for(i in 1:nrow(heT)){
                                    heP2 = heP[i,];
                                    for(j in 20:2){
                                        heP2[j] = heP2[j]-heP[i,j-1];
                                    }
                                    points(1:20,heT[i,]-heP2,pch=1,col=pcol[i])
                                    lines(1:20,heT[i,]-heP2,col=pcol[i],lty=1)
                                }
                                
                            }
                            axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                            if(g2 == 1){
                                axis(side=2,padj=1, cex.axis=0.8)
                            }
                            else{
                                axis(side=2,labels=NA)
                            }
                            if(g1==length(GRM)){
                                mtext(side=1,"MAF bin",line=1)
                            }
                            if(g2==1){
                                mtext(side=2,expression({widehat(h^2)}),line=1)
                            }
                            if((CORTYPE == "SAMPLE" & panel == 8) | (CORTYPE == "TGP" & panel == 9)){
                                legend("bottom",ncol=2,c(PART,"original"),col=pcol,lty=1,cex=0.9,bg=rgb(1,1,1,0.8),box.col=rgb(0,0,0,0))
                                box()
                            }
                            if(g1==1){
                                mtext(side=3,GRMs[grm2])
                            }
                            if(g2==length(GRM)){
                                mtext(side=4,GRMs[grm1],line=0.25)
                            }
                        }
                    }
                    dev.off()
                }
            }
        }

        {
            ##OK! Figured it out! Want 1-1_SAMPLE and 8-360_SAMPLE GRMs to correct with 10 PCs!!
            ## now create 3-panel plot with:
            ## A: data/permuted h2_bin
            ## B: perm-corrected h2_bin
            ## C: perm-corrected cummulative h2
            
            CORTYPE = "SAMPLE";
            if(CORTYPE == "SAMPLE"){ 
                GRM = c(1:4);
            }
            else{
                GRM = c(5:7);
            }
            g1 = 1;
            grm1 = GRM[g1];
            g2 = length(GRM);
            grm2 = GRM[g2];

            for(PC in c(0,10)){
                heP=matrix(,nr=length(PART),nc=20);
                heP2=matrix(,nr=length(PART),nc=20);
                heT=matrix(,nr=length(PART),nc=20);
                heT2=matrix(,nr=length(PART),nc=20);
                heTC=matrix(,nr=length(PART),nc=20);
                heTN=matrix(,nr=length(PART),nc=20);
                for(part in 1:length(PART)){
                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_GRM_",GRMs[grm1],"_",GRMs[grm2],"_",PART[part],".he.meanQuant",sep="")
                    if(file.exists(FILE)){
                        heP[part,] = scan(FILE,nlines=1)[1:20]
                    }
                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_GRM_",GRMs[grm1],"_",GRMs[grm2],"_",PART[part],".he.meanQuant",sep="")
                    if(file.exists(FILE)){
                        heT[part,] = scan(FILE,nlines=1)[1:20]
                    }

                    heP2[part,] = heP[part,];
                    heT2[part,] = heT[part,];
                    for(i in 20:2){
                        heP2[part,i] = heP[part,i]-heP[part,i-1];
                        heT2[part,i] = heT[part,i]-heT[part,i-1];

                    }
                    heTC[part,] = heT[part,]-heP2[part,];
                    for(i in 20:2){
                        heTC[part,i] = heTC[part,i]-heTC[part,i-1];
                    }
                    heTN[part,] = heT[part,]/heT[part,20]
                }

                {
                    pdf(paste(FIGDIR,"HEh2_GRMcor_PC=",PC,"_",GRMs[grm1],"_",GRMs[grm2],".pdf",sep=""),height=2,width=6,pointsize=10)
                    par(mfrow=c(1,4),mar=c(2,3,0.5,.5),oma=c(0,0.1,0,0))
                    LEGORD = c(5,4,1,2,3)
                    plot(1,1,type='n',ylim=range(heT2,heP2,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression({widehat(h^2)}[bin]),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heT2[i,],pch=1,col=pcol[i])
                        lines(1:20,heT2[i,],col=pcol[i],lty=1)
                        points(1:20,heP2[i,],pch=2,col=pcol[i])
                        lines(1:20,heP2[i,],col=pcol[i],lty=2)
                    }

                    plot(1,1,type='n',ylim=range(heTC,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Corrected ",{widehat(h^2)}[bin])),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heTC[i,],pch=1,col=pcol[i])
                        lines(1:20,heTC[i,],col=pcol[i],lty=1)
                    }

                    plot(1,1,type='n',ylim=range(0,0.2,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Cor. Prop. ",{widehat(h^2)}[bin])),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heTC[i,]/sum(heTC[i,]),pch=1,col=pcol[i])
                        lines(1:20,heTC[i,]/sum(heTC[i,]),col=pcol[i],lty=1)
                    }

                    plot(1,1,type='n',ylim=range(heTC,0.1,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Cor. Cum. ",{widehat(h^2)})),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heT[i,]-heP2[i,],pch=1,col=pcol[i])
                        lines(1:20,heT[i,]-heP2[i,],col=pcol[i],lty=1)
                    }
                    legend("topleft",PART[LEGORD],col=pcol[LEGORD],pch=1,cex=0.9,bg=rgb(1,1,1,0.8),box.col=rgb(0,0,0,0))
                    box()
                    dev.off()
                }
            }
        }

        {
            ##OK! Figured it out! Want 1-1_SAMPLE and 8-360_SAMPLE GRMs to correct with 10 PCs!!
            ## now create 3-panel plot with:
            ## A: data/permuted cum h2
            ## B: perm-corrected cum h2
            ## C: perm-corrected prop cum h2
            
            CORTYPE = "SAMPLE";
            if(CORTYPE == "SAMPLE"){ 
                GRM = c(1:4);
            }
            else{
                GRM = c(5:7);
            }
            g1 = 1;
            grm1 = GRM[g1];
            g2 = length(GRM);
            grm2 = GRM[g2];

            for(PC in c(0)){
                heP=matrix(,nr=length(PART),nc=20);
                heP2=matrix(,nr=length(PART),nc=20);
                heT=matrix(,nr=length(PART),nc=20);
                heT2=matrix(,nr=length(PART),nc=20);
                heTC=matrix(,nr=length(PART),nc=20);
                heTN=matrix(,nr=length(PART),nc=20);
                heTCN=matrix(,nr=length(PART),nc=20);
                for(part in 1:length(PART)){
                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_GRM_",GRMs[grm1],"_",GRMs[grm2],"_",PART[part],".he.meanQuant",sep="")
                    if(file.exists(FILE)){
                        heP[part,] = scan(FILE,nlines=1)[1:20]
                    }
                    FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_GRM_",GRMs[grm1],"_",GRMs[grm2],"_",PART[part],".he.meanQuant",sep="")
                    if(file.exists(FILE)){
                        heT[part,] = scan(FILE,nlines=1)[1:20]
                    }

                    heP2[part,] = heP[part,];
                    heT2[part,] = heT[part,];
                    for(i in 20:2){
                        heP2[part,i] = heP[part,i]-heP[part,i-1];
                        heT2[part,i] = heT[part,i]-heT[part,i-1];

                    }
                    heTC[part,] = heT[part,]-heP2[part,];
                    heTCN[part,] = heT[part,];
                    for(i in 20:2){
                        heTC[part,i] = heTC[part,i]-heTC[part,i-1];
                    }
                    heTN[part,] = heT[part,]/heT[part,20]
                    for(i in 1:20){
                        heTCN[part,i] = heTCN[part,i]-heP2[part,i];
                    }
                    for(i in 1:nrow(heTCN)){
                        heTCN[part,] = heTCN[part,]/max(heTCN[part,])
                    }
                }

                {
                    pdf(paste(FIGDIR,"HEh2_GRMcor_PC=",PC,"_",GRMs[grm1],"_",GRMs[grm2],"_3panel.pdf",sep=""),height=2,width=6,pointsize=10)
                    par(mfrow=c(1,3),mar=c(2,3,0.5,.5),oma=c(0,0.1,0,0))
                    LEGORD = c(5,4,1,2,3)
                    plot(1,1,type='n',ylim=range(heT,heP,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression({widehat(h^2)}[cum]),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heT[i,],pch=1,col=pcol[i])
                        lines(1:20,heT[i,],col=pcol[i],lty=1)
                        points(1:20,heP[i,],pch=2,col=pcol[i])
                        lines(1:20,heP[i,],col=pcol[i],lty=2)
                    }
                    legend("topleft",PART[LEGORD],col=pcol[LEGORD],pch=1,cex=0.9,bg=rgb(1,1,1,0.8),box.col=rgb(0,0,0,0))
                    box()

                    plot(1,1,type='n',ylim=range(c(0,apply(heTC,1,sum)),na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Corrected ",{widehat(h^2)}[cum])),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,cumsum(heTC[i,]),pch=1,col=pcol[i])
                        lines(1:20,cumsum(heTC[i,]),col=pcol[i],lty=1)
                    }

                    plot(1,1,type='n',ylim=range(0.02,1,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="",log='y')
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Cor. Prop. ",{widehat(h^2)}[cum])),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heTCN[i,],pch=1,col=pcol[i])
                        lines(1:20,heTCN[i,],col=pcol[i],lty=1)
                    }
                    dev.off()
                }
            }
        }

        {
            ##original analysis
            
            for(PC in c(0,10)){
                heP=matrix(,nr=length(PART),nc=20);
                heP2=matrix(,nr=length(PART),nc=20);
                heT=matrix(,nr=length(PART),nc=20);
                heT2=matrix(,nr=length(PART),nc=20);
                heTC=matrix(,nr=length(PART),nc=20);
                for(part in 1:length(PART)){
                    if(part==5){
                        FILE=paste("../../INFW_HEvREML4/OUTSUM/merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_",PART[part],".he.meanQuant",sep="")
                    }
                    else{
                        FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_PERM=3_",PART[part],".he.meanQuant",sep="")
                    }
                    if(file.exists(FILE)){
                        heP[part,] = scan(FILE,nlines=1)[1:20]
                    }
                    if(part==5){
                        FILE=paste("../../INFW_HEvREML4/OUTSUM/merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_",PART[part],".he.meanQuant",sep="")
                    }
                    else{
                        FILE=paste("merge_1000000_K=20_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1_",PART[part],".he.meanQuant",sep="")
                    }
                    if(file.exists(FILE)){
                        heT[part,] = scan(FILE,nlines=1)[1:20]
                    }

                    heP2[part,] = heP[part,];
                    heT2[part,] = heT[part,];
                    for(i in 20:2){
                        heP2[part,i] = heP[part,i]-heP[part,i-1];
                        heT2[part,i] = heT[part,i]-heT[part,i-1];
                    }
                    heTC[part,] = heT2[part,]-heP2[part,];
                }

                {
                    pdf(paste(FIGDIR,"HEh2_noGRMcor_PC=",PC,".MIX.pdf",sep=""),height=2,width=6,pointsize=10)
                    par(mfrow=c(1,4),mar=c(2,3,0.5,.5),oma=c(0,0.1,0,0))
                    LEGORD = c(5,4,1,2,3)
                    plot(1,1,type='n',ylim=range(heT2,heP2,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression({widehat(h^2)}[bin]),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heT2[i,],pch=1,col=pcol[i])
                        lines(1:20,heT2[i,],col=pcol[i],lty=1)
                        points(1:20,heP2[i,],pch=2,col=pcol[i])
                        lines(1:20,heP2[i,],col=pcol[i],lty=2)
                    }

                    plot(1,1,type='n',ylim=range(heTC,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Corrected ",{widehat(h^2)}[bin])),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heTC[i,],pch=1,col=pcol[i])
                        lines(1:20,heTC[i,],col=pcol[i],lty=1)
                    }

                    plot(1,1,type='n',ylim=range(0,1,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Cor. Prop. ",{widehat(h^2)}[bin])),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,heTC[i,]/sum(heTC[i,]),pch=1,col=pcol[i])
                        lines(1:20,heTC[i,]/sum(heTC[i,]),col=pcol[i],lty=1)
                    }

                    plot(1,1,type='n',ylim=range(heTC,0.1,na.rm=T),xlim=c(1,20),xaxt='n',yaxt='n',ylab="",xlab="",main="")
                    axis(side=1,at=1:20, labels=NA,padj=-1, cex.axis=0.8)
                    axis(side=2,padj=1, cex.axis=0.8)
                    mtext(side=1,"MAF bin",line=1)
                    mtext(side=2,expression(paste("Cor. Cum. ",{widehat(h^2)})),line=1.1,cex=0.8)
                    abline(h=0)
                    for(i in rev(LEGORD)){
                        points(1:20,cumsum(heTC[i,]),pch=1,col=pcol[i])
                        lines(1:20,cumsum(heTC[i,]),col=pcol[i],lty=1)
                    }
                    legend("topleft",PART[LEGORD],col=pcol[LEGORD],pch=1,cex=0.9,bg=rgb(1,1,1,0.8),box.col=rgb(0,0,0,0))
                    box()
                    dev.off()
                }
            }
        }
    }
}

{
    ##find singleton partition
    minSING=c(1e-5, 6.46e-5,0.0002, 0.0004, 0.001, 0.5)
    
    part="SINGgnomAD";
    fooG=read.table(paste(part,"_SFS_K=",K,".txt",sep=""));
    part="SINGtgp";
    fooT=read.table(paste(part,"_SFS_K=",K,".txt",sep=""));

    N=matrix(,nr=length(minSING),nc=2)
    pmin = -1
    for(i in 1:length(minSING)){
        N[i,1] = sum(fooG[which(fooG[,1]<=minSING[i] & fooG[,1]>pmin),2])/sum(fooG[,2]);
        N[i,2] = sum(fooT[which(fooT[,1]<=minSING[i] & fooT[,1]>pmin),2])/sum(fooT[,2]);
        pmin = minSING[i];
    }
    
    bin=c(1,2,5,9,14,20); s=foo[1,2]; cat(1,"\t",1,"\t",foo[1,2],"\n");for(i in 2:length(bin)){s[i]=sum(foo[(bin[i-1]+1):bin[i],2]); cat(i," ",bin[i]," ",foo[bin[i],1]," ",s[i],"\n")}; plot(foo[bin,1],s/sum(s),log='x',ylim=range(0,s/sum(s)),type='b');abline(v=2/720)

    bin=c(1,2,5,9,14,20); s=foo[1,2]; cat(1,"\t",1,"\t",foo[1,2],"\n");for(i in 2:length(bin)){s[i]=sum(foo[(bin[i-1]+1):bin[i],2]); cat(i," ",bin[i]," ",foo[bin[i],1]," ",s[i],"\n")}; plot(foo[bin,1],s/sum(s),log='x',ylim=range(0,s/sum(s)),type='b');abline(v=2/720)
    
}

{ ## plot h2 as a function of MAF for different singleton partitions
    SINGgMAF = c(1e-5, 6.46e-5, 0.0002, 0.0004, 0.001, 0.0017);
    GDMAC = c(2, 3, 5, 8, 11, 17, 25, 36, 50, 68, 89, 114, 143, 174, 207, 242, 280, 319, 360)/720;
    xval = c(SINGgMAF,GDMAC);
    PC=10;
    
    h2=array(,dim=c(2,2,6,25)); # [PERM 0/3] [TGP/gnomAD] [singBins=1:6] [bin=1:25]
    cumh2=array(,dim=c(2,2,6,25)); # [PERM 0/3] [TGP/gnomAD] [singBins=1:6] [bin=1:25]
    rcumh2=array(,dim=c(2,2,6,25)); # [PERM 0/3] [TGP/gnomAD] [singBins=1:6] [bin=1:25]
    PERM=c("","_PERM=3");
    PART=c("TGP","gnomAD");
    for(p1 in 1:length(PERM)){
        perm = PERM[p1];
        for(p2 in 1:length(PART)){
            part = PART[p2];
            for(SB in 1:6){
                FILE=paste("merge_1000000_K=",19+SB,"_MAC=1_ngPC=",PC,"_npPC=",PC,"_QN=1_DSAMP=1",perm,"_",part,".SINGpart.meanQuant",sep="")
                if(file.exists(FILE)){
                    d=scan(FILE,nlines=1)[1:(19+SB)];
                    cumh2[p1,p2,SB,(25-(19+SB)+1):25] = d;
                    h2[p1,p2,SB,25-(19+SB)+1] = cumh2[p1,p2,SB,25-(19+SB)+1];
                    for(i in (25-(19+SB)+2):25){
                        h2[p1,p2,SB,i] = cumh2[p1,p2,SB,i]-cumh2[p1,p2,SB,i-1];
                    }
                    rcumh2[p1,p2,SB,25] = h2[p1,p2,SB,25];
                    for(i in 24:(25-(19+SB)+1)){
                        rcumh2[p1,p2,SB,i] = h2[p1,p2,SB,i]+rcumh2[p1,p2,SB,i+1];
                    }
                }
            }
        }
    }

    ##single plot with 9 panels: columns: data, perm, corr;  rows: h2, cumh2, rcumh2
    {
        pdf(paste(FIGDIR,"SUPP/HEh2_SINGpart_PC=",PC,".all.pdf",sep=""),height=6,width=6,pointsize=10)
        {
            par(mfrow=c(3,3),mar=c(1,1,.5,.5),oma=c(2,2,1.5,0.5))
            {
                for(p1 in 1:2){
                    plot(1,1,type='n',xlim=range(SINGgMAF,GDMAC),ylim=range(h2,na.rm=T),log='x',xaxt='n',yaxt='n')
                    abline(h=0,col='grey')
                    if(p1==1){
                        mtext(side=2,expression(widehat(h^2)[bin]), line=1.1,cex=0.8)
                    }
                    if(p1==1){
                        mtext(side=3, "Obs Data")
                    }
                    else{
                        mtext(side=3, "trans-Permutation")
                    }
                    axis(side=1, padj=-1.5, cex.axis=0.8)
                    axis(side=2, padj=1, cex.axis=0.8)
                    for(p2 in 1:2){
                        for(SB in 2:6){
                            lines(xval, h2[p1,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]),lty=p1)
                            points(xval, h2[p1,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                        }
                        SB=1;
                        lines(xval, h2[p1,p2,SB,], col="black",lty=p1)
                        points(xval, h2[p1,p2,SB,], col="black")

                    }
                }
                plot(1,1,type='n',xlim=range(SINGgMAF,GDMAC),ylim=range(h2,na.rm=T),log='x',xaxt='n',yaxt='n')
                abline(h=0,col='grey')
                mtext(side=3, "Corrected")
                axis(side=1, padj=-1.5, cex.axis=0.8)
                axis(side=2, padj=1, cex.axis=0.8)
                for(p2 in 1:2){
                    for(SB in 2:6){
                        lines(xval, h2[1,p2,SB,]-h2[2,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                        points(xval, h2[1,p2,SB,]-h2[2,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                    }
                    SB=1;
                    lines(xval, h2[1,p2,SB,]-h2[2,p2,SB,], col="black")
                    points(xval, h2[1,p2,SB,]-h2[2,p2,SB,], col="black")
                }
            }

            {
                for(p1 in 1:2){
                    plot(1,1,type='n',xlim=range(SINGgMAF,GDMAC),ylim=range(cumh2,na.rm=T),log='x',xaxt='n',yaxt='n')
                    abline(h=0,col='grey')
                    if(p1==1){
                        mtext(side=2,expression(paste("cum ", widehat(h^2)[MAF<x])), line=1.1,cex=0.8)
                    }
                    axis(side=1, padj=-1.5, cex.axis=0.8)
                    axis(side=2, padj=1, cex.axis=0.8)
                    for(p2 in 1:2){
                        for(SB in 2:6){
                            lines(xval, cumh2[p1,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]),lty=p1)
                            points(xval, cumh2[p1,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                        }
                        SB=1;
                        lines(xval, cumh2[p1,p2,SB,], col="black",lty=p1)
                        points(xval, cumh2[p1,p2,SB,], col="black")
                    }
                }
                plot(1,1,type='n',xlim=range(SINGgMAF,GDMAC),ylim=range(cumh2,na.rm=T),log='x',xaxt='n',yaxt='n')
                abline(h=0,col='grey')
                axis(side=1, padj=-1.5, cex.axis=0.8)
                axis(side=2, padj=1, cex.axis=0.8)
                for(p2 in 1:2){
                    for(SB in 2:6){
                        lines(xval, cumh2[1,p2,SB,]-cumh2[2,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                        points(xval, cumh2[1,p2,SB,]-cumh2[2,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                    }
                    SB=1;
                    lines(xval, cumh2[1,p2,SB,]-cumh2[2,p2,SB,], col="black",lty=p1)
                    points(xval, cumh2[1,p2,SB,]-cumh2[2,p2,SB,], col="black")
                }
            }

            {
                for(p1 in 1:2){
                    plot(1,1,type='n',xlim=range(SINGgMAF,GDMAC),ylim=range(rcumh2,na.rm=T),log='x',xaxt='n',yaxt='n')
                    abline(h=0,col='grey')
                    if(p1==1){
                        mtext(side=2,expression(paste("rev cum ", widehat(h^2)[MAF>x])), line=1.1,cex=0.8)
                    }
                    mtext(side=1,"MAF (x; log-scale)", line=1.5,cex=0.8)
                    
                    axis(side=1, padj=-1.5, cex.axis=0.8)
                    axis(side=2, padj=1, cex.axis=0.8)
                    for(p2 in 1:2){
                        for(SB in 2:6){
                            lines(xval, rcumh2[p1,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]),lty=p1)
                            points(xval, rcumh2[p1,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                        }
                        SB=1;
                        lines(xval, rcumh2[p1,p2,SB,], col="black",lty=p1)
                        points(xval, rcumh2[p1,p2,SB,], col="black")
                    }
                }
                plot(1,1,type='n',xlim=range(SINGgMAF,GDMAC),ylim=range(rcumh2,na.rm=T),log='x',xaxt='n',yaxt='n')
                abline(h=0,col='grey')
                mtext(side=1,"MAF (x; log-scale)", line=1.5,cex=0.8)
                axis(side=1, padj=-1.5, cex.axis=0.8)
                axis(side=2, padj=1, cex.axis=0.8)
                for(p2 in 1:2){
                    for(SB in 2:6){
                        lines(xval, rcumh2[1,p2,SB,]-rcumh2[2,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                        points(xval, rcumh2[1,p2,SB,]-rcumh2[2,p2,SB,], col=rgb(abs(1-p2),0,2-p2,seq(0.14,1,length=6)[SB]))
                    }
                    SB=1;
                    lines(xval, rcumh2[1,p2,SB,]-rcumh2[2,p2,SB,], col="black",lty=p1)
                    points(xval, rcumh2[1,p2,SB,]-rcumh2[2,p2,SB,], col="black")
                }
            }
        }
        dev.off()
    }
}





########  Generate sub figures for presentations

{ # Figure 2.
    ## This is based on reading *.meanQuant files
    ##cum h2:  dat.0[,1:K]
    ##mean maxMAF per bin: dat.0[,K+(gPC+pPC)+1:K]
    ##Number of SNPs:  dat.0[,2*K+(gPC+pPC)+1:K]
    
    {  
        W="1000000"
        K=20
        gPC=10
        pPC=10
        QN=1
        
        {##get data
            MAF = array(,dim=c(5,20))
            QMAT = array(,dim=c(2,5,20))
            H2MATS = array(,dim=c(5,20))
            MM = c(36,5,2,1);
            for(i in 1:length(MM)){
                ##first read in data h2 (cumsum)
                filename=paste("merge_1000000_K=20_MAC=",MM[i],"_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="");
                foo=read.table(filename);
                MAF[i,1:K] = unlist(foo[1,K+(gPC+pPC)+1:K])
                h2 = foo[1,1:K];
                ql = foo[2,1:K];
                qu = foo[3,1:K];
                ##transform to h2 per bin
                for(j in K:2){
                    h2[j] = h2[j] - h2[j-1];
                    ql[j] = ql[j] - ql[j-1];
                    qu[j] = qu[j] - qu[j-1];
                }
                ##now get permuted values
                filename=paste("merge_1000000_K=20_MAC=",MM[i],"_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="");
                foo=read.table(filename);
                h2p = foo[1,1:K];
                qlp = foo[2,1:K];
                qup = foo[3,1:K];
                ##transform to h2 per bin
                for(j in K:2){
                    h2p[j] = h2p[j] - h2p[j-1];
                    qlp[j] = qlp[j] - qlp[j-1];
                    qup[j] = qup[j] - qup[j-1];
                }
                
                H2MATS[i,1:K] = unlist(h2-h2p)
                QMAT[1,i,1:K] = unlist(ql-qlp) #lower
                QMAT[2,i,1:K] = unlist(qu-qup) #upper
            }
            
            ##now get gnomAD partition data
            foo=read.table(paste("ALLgnomAD_SFS_K=",K,".txt",sep=""));
            MAF[5,1:K] = foo[,1];
            MAF[5,1] = 1e-5
            filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_ALLgnomAD.he.meanQuant";
            foo=read.table(filename);
            h2 = foo[1,1:K];
            ql = foo[2,1:K];
            qu = foo[3,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2[j] = h2[j] - h2[j-1];
                ql[j] = ql[j] - ql[j-1];
                qu[j] = qu[j] - qu[j-1];
            }
            ##now get permuted values
            filename="merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_ALLgnomAD.he.meanQuant";
            foo=read.table(filename);
            h2p = foo[1,1:K];
            qlp = foo[2,1:K];
            qup = foo[3,1:K];
            ##transform to h2 per bin
            for(j in K:2){
                h2p[j] = h2p[j] - h2p[j-1];
                qlp[j] = qlp[j] - qlp[j-1];
                qup[j] = qup[j] - qup[j-1];
            }
            
            H2MATS[5,1:K] = unlist(h2-h2p)
            QMAT[1,5,1:K] = unlist(ql-qlp) #lower
            QMAT[2,5,1:K] = unlist(qu-qup) #upper
        }
        
        { ##plot figure 2A: cum h2 vs MAF across minMAC
            pdf(paste(FIGDIR,"/PresFigs/2A_cum_h2_vs_MAF_by_minMAC_",0,".pdf",sep=""),height=2,width=2.75,pointsize=10)
            par(mar=c(2,2.8,0.3,0.3))
            plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=c(0,0.09),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
            abline(h=seq(0,0.1,by=0.02),col=rgb(0,0,0,0.1),lwd=0.75)
            abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
            box()
            axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
            axis(side=2,padj=1.3,cex.axis=0.7)
            mtext(side=1,"Minor Allele Frequency",line=1)
            mtext(side=2,expression(paste("Cumulative ",{widehat(h^2)},sep="")),line=1.15)
            dev.off()

            for(STEP in 1:5){
                pdf(paste(FIGDIR,"/PresFigs/2A_cum_h2_vs_MAF_by_minMAC_",STEP,".pdf",sep=""),height=2,width=2.75,pointsize=10)
                par(mar=c(2,2.8,0.3,0.3))
                plot(1,1,type='n',xlim=c(1e-5,0.5),ylim=c(0,0.09),log='x',xaxt='n',yaxt='n',xlab="",ylab="")
                abline(h=seq(0,0.1,by=0.02),col=rgb(0,0,0,0.1),lwd=0.75)
                abline(v=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),col=rgb(0,0,0,0.1),lwd=0.75)
                LEG=c("gnomAD",expression("MAC" >= 1),expression("MAC" >= 2),expression("MAC" >= 5),expression("MAF" >= "5%"));
                legend("topleft",LEG[(5-STEP+1):5], col=pcol[STEP:1], lty=1,lwd=2,bg=rgb(1,1,1,0.8),ncol=2,box.col="white",pch=20,cex=0.8)
                box()
                for(i in 1:STEP){
                    pci = pcol[i]
                    rgbi = col2rgb(pci)/255;
                    
                    polygon(x=c(MAF[i,1:K], rev(MAF[i,1:K])), y=c(cumsum(QMAT[1,i,1:K]),rev(cumsum(QMAT[2,i,1:K]))), col=rgb(rgbi[1],rgbi[2],rgbi[3],0.2),border=NA)
                    lines(MAF[i,1:K],cumsum(H2MATS[i,1:K]),col=pcol[i],lwd=0.5)
                    points(MAF[i,1:K],cumsum(H2MATS[i,1:K]),col=pcol[i],pch=20)
                }
                
                axis(side=1,at=c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5),labels=c(paste("*"),"1e-4","1e-3","0.01","0.1","0.5"),padj=-2, cex.axis=0.7)
                axis(side=2,padj=1.3,cex.axis=0.7)
                mtext(side=1,"Minor Allele Frequency",line=1)
                mtext(side=2,expression(paste("Cumulative ",{widehat(h^2)},sep="")),line=1.15)
                dev.off()
            }
        }
    }
}


{ ##playing around...
    ##normal cumulative
    
    pdf(paste(FIGDIR,"/NotUsed/DCF_nSingBins.pdf",sep=""),height=2,width=4.75,pointsize=10)
    par(mfrow=c(1,2),mar=c(2,2,0,0))
    plot(1,1,type='n', ylim=c(0,0.12),xlim=c(1,20),log='x');
    legend("topleft",legend=1:5,lty=1,col=pcol)
    Kt=19;
    for(nSB in 1:3){
        K=Kt+nSB;
        d1=read.table(paste("merge_1000000_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_gnomAD.SINGpart.meanQuant",sep=""))
        d2=read.table(paste("merge_1000000_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_gnomAD.SINGpart.meanQuant",sep=""))
        t1=unlist(d1[1,1:K])
        tb1 = t1;
        t2=unlist(d2[1,1:K])
        tb2 = t2;
        for(i in K:2){
            tb1[i] = t1[i]-t1[i-1];
            tb2[i] = t2[i]-t2[i-1];
        }
        lines(c(1:nSB/nSB+0.8,2:(Kt+1)), (t1-t2),col=pcol[nSB])
    }

    ##reverse cumulative
    plot(1,1,type='n', ylim=c(0,0.12),xlim=c(1,20),log='x');
    legend("topright",legend=1:5,lty=1,col=pcol)
    Kt=19;
    for(nSB in 1:5){
        K=Kt+nSB;
        d1=read.table(paste("merge_1000000_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_gnomAD.SINGpart.meanQuant",sep=""))
        d2=read.table(paste("merge_1000000_K=",K,"_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_gnomAD.SINGpart.meanQuant",sep=""))
        t1=unlist(d1[1,1:K])
        tb1 = t1;
        t2=unlist(d2[1,1:K])
        tb2 = t2;
        for(i in K:2){
            tb1[i] = t1[i]-t1[i-1];
            tb2[i] = t2[i]-t2[i-1];
        }
        lines(c(1:nSB/nSB+0.8,2:(Kt+1)), rev(cumsum(rev(tb1-tb2))),col=pcol[nSB])
    }
    dev.off()
}



{ ##Play with Lawrence data
    Lh2=scan("~/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures5/LawrenceData/transformed.real.data")
    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.h2.gz",sep="");
    F=gzfile(filename);
    h2=read.table(F);
    dim(h2)
    h2 = subset(h2,h2[,1] %in% expGenes[,1])
    h2 = t(apply(h2[,2+1:20],1,function(x){for(j in which(is.na(x))){x[j] = 0}; x}));
    dim(h2)
    
    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.h2.gz",sep="");
    F=gzfile(filename);
    h2p=read.table(F);
    ##close(F);
    dim(h2p)
    h2p = subset(h2p, h2p[,1] %in% expGenes[,1])
    h2p = t(apply(h2p[,2+1:20],1,function(x){for(j in which(is.na(x))){x[j] = 0}; x}));
    dim(h2p);

    he = apply(t(apply(h2,1,cumsum)),2,function(x){mean(x,na.rm=T)})
    hep = apply(t(apply(h2p,1,cumsum)),2,function(x){mean(x,na.rm=T)})
    hec=he-hep;

    pdf(paste(FIGDIR,"/NotUsed/CompLawrence_CDF.pdf",sep=""),height=2,width=2.75,pointsize=10)
    par(mar=c(2,2.8,0.3,0.3))
    plot(tmaf/720,cumsum(Lh2),ylim=c(0,1),log='x',col='red')
    lines(MAF,he/max(he),col='blue',pch=20);
    lines(MAF,hec/max(hec),col=pcol[4],pch=20);
    legend("topleft",c("Lawrence","Raw PCA-corrected","Permutation-Corrected"),col=c("red","blue",pcol[4]),lty=c(0,1,1),pch=c(1,-1,-1),cex=0.8)
    dev.off()

    foo=read.table("~/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures5/LawrenceData/real_data.posterior.txt")
    pdf(paste(FIGDIR,"/NotUsed/posterior_params.pdf",sep=""),height=8,width=2.75,pointsize=10)
    par(mfrow=c(3,1),mar=c(2.1,2.8,0.3,0.3),oma=c(.2,0,0,0))
    plot(density(foo[,1],adjust=.2),xlab="",ylab="",xaxt='n',yaxt='n',main="",xlim=c(0,1))
    abline(v=mean(foo[,1]))
    axis(side=1,padj=-2, cex.axis=0.7)
    axis(side=2,padj=1.3,cex.axis=0.7)
    mtext(side=1,"rho",line=1)
    mtext(side=2,"Posterior",line=1.15)
    plot(density(foo[,2],adjust=.2),xlab="",ylab="",xaxt='n',yaxt='n',main="",xlim=c(0,1))
    abline(v=mean(foo[,2]))
    axis(side=1,padj=-2, cex.axis=0.7)
    axis(side=2,padj=1.3,cex.axis=0.7)
    mtext(side=1,"tau",line=1)
    mtext(side=2,"Posterior",line=1.15)
    plot(density(foo[,3],adjust=.2),xlab="",ylab="",xaxt='n',yaxt='n',main="",xlim=c(0,1))
    abline(v=mean(foo[,3]))
    axis(side=1,padj=-2, cex.axis=0.7)
    axis(side=2,padj=1.3,cex.axis=0.7)
    mtext(side=1,"phi",line=1)
    mtext(side=2,"Posterior",line=1.15)
    dev.off();

    pdf(paste(FIGDIR,"/NotUsed/posterior_paramsJoint_rhotau.pdf",sep=""),height=2.75,width=2.75,pointsize=10)
    par(mar=c(2,2,0.3,0.3))
    plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab="",ylab="",main="")
    axis(side=1,padj=-2, cex.axis=0.7)
    axis(side=2,padj=1.3,cex.axis=0.7)
    mtext(side=1,expression(rho),line=1)
    mtext(side=2,expression(tau),line=1.15)

    apply(foo,1,function(x){points(x[1],x[2],pch=".",col=rgb(0,0,0,0.25))})
    dev.off()

    foo=read.table("~/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures5/LawrenceData/resamp_sims.txt")
    pdf(paste(FIGDIR,"/NotUsed/posterior_draws.pdf",sep=""),height=2,width=2.75,pointsize=10)
    par(mar=c(2,2.8,0.3,0.3))
    plot(1,1,type='n',xlim=c(1e-3,0.5),ylim=c(0,1),xlab="",ylab="",main="",xaxt='n',yaxt='n',log='x')
    abline(v=mean(foo[,1]))
    axis(side=1,padj=-2, cex.axis=0.7)
    axis(side=2,padj=1.3,cex.axis=0.7)
    mtext(side=1,"MAF",line=1)
    mtext(side=2,expression(paste("cum. prop. ",widehat(h^2))),line=1.15)
    for(i in 1:nrow(foo)){
        lines(tmaf/720, cumsum(unlist(foo[i,7:ncol(foo)])), col=rgb(0,0,0,0.1))
    }
    lines(tmaf/720,cumsum(Lh2),col=pcol[4],lwd=2)
    points(tmaf/720,cumsum(Lh2),col=pcol[4],pch=20)
    dev.off()
}

{ #new figs for Lawrence' data
    t = read.table("/Users/rhernandez/Dropbox/Research/Data/GEUVADIS/Scripts/h2Rare/Manuscript/Figures5/LawrenceData/real_infer.5000.txt");
    al<-function(m,s) {
	return( (m^2 - m^3 -m*(s^2))/s^2)
    }
    
    be <-function(m,s) {
	return( (m-1)*(m^2-m+s^2)/(s^2))
    }

    ## sample marginal tau
    marg_tau<-c()
    numE = 0;
    for (i in seq(1,5000)) {
        numE = numE+1;
        ##cat(numE,": ",t$V1[i], t$V2[i], t$V4[i], t$V5[i],"\n");
        marg_tau<-c(marg_tau,rbeta(100,al(t$V2[i],t$V5[i]),be(t$V2[i],t$V5[i])))
    }

    ## samples marginal rho
    marg_rho<-c()
    numS = 0;
    for (i in seq(1,5000)) {
        numS = numS+1;
        marg_rho<-c(marg_rho,rbeta(100,al(t$V1[i],t$V4[i]),be(t$V1[i],t$V4[i])))
    }

    ## samples marginal phi
    marg_phi<-c()
    numS = 0;
    for (i in seq(1,5000)) {
        numS = numS+1;
        marg_phi<-c(marg_phi,rbeta(100,al(t$V3[i],t$V6[i]),be(t$V3[i],t$V6[i])))
    }

    
    ## samples in vicinty of eyre-walker (i.e., rho ~ 1)
    data_tau<-c()
    numE = 0;
    for (i in seq(1,5000)) {
	if(t$V1[i] > 0.9995 ) {
            numE = numE+1;
            cat(numE,": ",t$V1[i], t$V2[i], t$V4[i], t$V5[i],"\n");
	    data_tau<-c(data_tau,rbeta(10000,al(t$V2[i],t$V5[i]),be(t$V2[i],t$V5[i])))
	}
    }

    ## samples in vicinty of simons (i.e., tau ~ 1)
    data_rho<-c()
    numS = 0;
    for (i in seq(1,5000)) {
	if(t$V2[i] > 0.9995 ) {
            numS = numS+1;
	    data_rho<-c(data_rho,rbeta(10000,al(t$V1[i],t$V4[i]),be(t$V1[i],t$V4[i])))
	}
    }

    foot = t[rev(order(t[,1])),]
    foor = t[rev(order(t[,2])),]

    Nlines=9;
    pdf(paste(FIGDIR,"/NotUsed/Eyrewalker_Simons_sort.pdf",sep=""),height=4,width=8,pointsize=10);
    par(mfrow=c(1,2))
    plot(1,1,type='n',xlim=c(0,1),ylim=c(0,5),main="Eyre-Walker Model (rho=1)",ylab="Density",xlab="tau")
    for(i in 1:Nlines){
        data_tau = rbeta(10000,al(foot$V2[i],foot$V5[i]),be(foot$V2[i],foot$V5[i]))
        lines(density(data_tau),col=pcol[i])
    }
    legend("topright",legend=foot[1:Nlines,1],col=c(pcol[1:Nlines]),lty=1)

    plot(1,1,type='n',xlim=c(0,1),ylim=c(0,5),main="Simons et al Model (tau=1)",xlab="rho",ylab="Density")
    for(i in 1:Nlines){
        data_rho = rbeta(10000,al(foor$V1[i],foor$V4[i]),be(foor$V1[i],foor$V4[i]))
        lines(density(data_rho),col=pcol[i])
    }
    legend("topright",legend=foor[1:Nlines,2],col=c(pcol[1:Nlines]),lty=1)
    dev.off()

    { #plot all densities for 3 parameter model
        
        jpeg(paste(FIGDIR,"/NotUsed/3param_all.jpg",sep=""),height=4,width=8,units="in",pointsize=10,res=300);
        par(mfrow=c(1,3),mar=c(3.5,3.5,0,0),oma=c(0,0,0.5,0.5))
        param=c("rho","tau","phi");
        x=seq(1e-5,1-1e-5,by=1e-4);
        v=viridis(nrow(t),option="C")
        for(i in 1:3){ #rho, tau, phi
            cat("working on",param[i],"\n");
            plot(1,1,type='n',xlim=c(0,1),ylim=c(0,10),xlab="",ylab="",main="")
            mtext(side=1,param[i],line=2)
            mtext(side=2,"Density",line=2)
            foo = t[order(t[,i]),]
            d = c();
            for(j in sample(1:nrow(foo))){
                y=dbeta(x,al(foo[j,i],foo[j,i+3]),be(foo[j,i],foo[j,i+3]))
                d=c(d,rbeta(1000,al(foo[j,i],foo[j,i+3]),be(foo[j,i],foo[j,i+3])));
                col=col2rgb(v[j])/255;
                lines(x,y,col=rgb(col[1],col[2],col[3],0.1))
            }
            lines(density(d),lwd=2)
        }
        dev.off()
    }

    { #plot ~joint density for rho-tau from 3 parameter model
        
        d1 = c();
        d2 = c();
        for(j in sample(1:nrow(foo))){
            d1=c(d1,rbeta(1000,al(foo[j,1],foo[j,4]),be(foo[j,1],foo[j,4])));
            d2=c(d2,rbeta(1000,al(foo[j,2],foo[j,5]),be(foo[j,2],foo[j,5])));
        }

        
        BREAKS = 50;
        Bx = seq(0, 1.001, length=BREAKS+1)
        By = seq(0, 1.001, length=BREAKS+1)
        
        JOINT=matrix(,nr=BREAKS,nc=BREAKS);
        for(i in 1:BREAKS){
            cat("pulling out i=",i,"\n");
            for(j in 1:BREAKS){
                JOINT[i,j] = sum(d1>=Bx[i] & d1<Bx[i+1] & d2>=By[j] & d2<By[j+1])
            }
        }

        jpeg(paste(FIGDIR,"/NotUsed/3param_joint_rho-tau_log.jpg",sep=""),height=4,width=8,units="in",pointsize=10,res=300);
        par(mar=c(3.5,3.5,0.5,0.5))
        param=c("rho","tau","phi");
        plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="")
        mtext(side=1,expression(rho),line=2)
        mtext(side=2,expression(tau),line=2)
        ##        points(d1,d2,pch=".")
        ##        dev.off()
        image(Bx[1:BREAKS], By[1:BREAKS], log(JOINT),xlim=c(0,1),ylim=c(0,1),
              col = viridis(BREAKS, option = "C"),
              xlab = "", ylab = "",xaxt = "n", yaxt = "n", bty = "n")
        dev.off();
    }

    { #try it as a ridgeline plot?
        maxy1 = 1;
        maxy2 = 1;
        B2 = 10;
        for(i in 1:B2){
            foo=which(d1<=i/B2 & d1>(i-1)/B2);
            foo2=density(d2[foo],adjust=0.5);
            if(max(foo2$y) > maxy1){
                maxy1=max(foo2$y);
            }
            foo=which(d2<=i/B2 & d2>(i-1)/B2);
            foo2=density(d1[foo],adjust=0.5);
            if(max(foo2$y) > maxy2){
                maxy2=max(foo2$y);
            }
        }
        pdf(paste(FIGDIR,"/NotUsed/3param_ridgeline_rho-tau.pdf",sep=""),height=2,width=4.25,pointsize=10);
        par(mfrow=c(1,2),mar=c(2.5,2.5,.5,.5))
        plot(1,1,type='n',xlim=c(0,1),ylim=c(1,11),xaxt='n',yaxt='n',main="",xlab="",ylab="")
        mtext(side=1,expression(tau),line=1.5)
        mtext(side=2,expression(rho),line=1.5)
        axis(side=1,cex.axis=0.8,padj=-1)
        axis(side=2,at=1:(B2+1),labels=0:B2/B2,cex.axis=0.8,padj=1)
        abline(h=0:(B2+1),col=rgb(0,0,0,0.5));
        vcol=viridis(B2,option="C");
        for(i in B2:1){
            foo=which(d1<=i/B2 & d1>(i-1)/B2);
            foo2=density(d2[foo],adjust=0.5);
            y2=i+foo2$y/maxy1*1.2;
            tc=col2rgb(vcol[i])/255;
            polygon(x=c(foo2$x,foo2$x[1]), y=c(y2,y2[1]), col=rgb(tc[1],tc[2],tc[3],0.2),border=NA)
            lines(foo2$x,y2);
        }
        
        plot(1,1,type='n',xlim=c(0,1),ylim=c(1,11),xaxt='n',yaxt='n',main="",xlab="",ylab="")
        mtext(side=1,expression(rho),line=1.5)
        mtext(side=2,expression(tau),line=1.5)
        axis(side=1,cex.axis=0.8,padj=-1)
        axis(side=2,at=1:(B2+1),labels=0:B2/B2,cex.axis=0.8,padj=1)
        abline(h=0:(B2+1),col=rgb(0,0,0,0.5));
        vcol=viridis(B2,option="C");
        for(i in B2:1){
            foo=which(d2<=i/B2 & d2>(i-1)/B2);
            foo2=density(d1[foo],adjust=0.5);
            y2=i+foo2$y/maxy2*1.2;
            tc=col2rgb(vcol[i])/255;
            polygon(x=c(foo2$x,foo2$x[1]), y=c(y2,y2[1]), col=rgb(tc[1],tc[2],tc[3],0.2),border=NA)
            lines(foo2$x,y2);
        }
        dev.off()
    }

    { #try as ecdfs
        pdf(paste(FIGDIR,"/NotUsed/posterior_cdf_rho-tau.pdf",sep=""),height=2,width=4.5,pointsize=10)
        par(mfrow=c(1,2),mar=c(2,2,0.3,0.3))
        plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab="",ylab="",main="")
        axis(side=1,padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=1,expression(tau),line=1)
        mtext(side=2,"Cumulative Density",line=1.15)
        
        B2 = 10;
        B3 = 100;
        x=seq(1/B3,1,length=B3);
        legend("topleft",rev(c(expression(paste("0.4<",rho,"<=0.5")),expression(paste("0.5<",rho,"<=0.6")),expression(paste("0.6<",rho,"<=0.7")),expression(paste("0.7<",rho,"<=0.8")),expression(paste("0.8<",rho,"<=0.9")),expression(paste("0.9<",rho,"<=1.0")) )), col=vcol[10:5], bg=rgb(1,1,1,0.8),lty=1,box.col=rgb(1,1,1,0.8), cex=0.7)
        legend("bottomright",rev(c(expression(rho <= 0.1),expression(paste("0.1<",rho,"<=0.2")),expression(paste("0.2<",rho,"<=0.3")),expression(paste("0.3<",rho,"<=0.4")))), col=vcol[4:1], bg=rgb(1,1,1,0.8),lty=1,box.col=rgb(1,1,1,0.8), cex=0.7)
        box()
        for(i in 1:B2){
            foo=which(d1<=i/B2 & d1>(i-1)/B2);
            y=c();
            for(j in x){
                y=c(y,mean(d2[foo] <= j));
            }
            lines(x,y,col=vcol[i]);
        }



        plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab="",ylab="",main="")
        axis(side=1,padj=-2, cex.axis=0.7)
        axis(side=2,padj=1.3,cex.axis=0.7)
        mtext(side=1,expression(rho),line=1)
        mtext(side=2,"Cumulative Density",line=1.15)
        
        B2 = 10;
        B3 = 100;
        x=seq(1/B3,1,length=B3);
        legend("topleft",rev(c(expression(paste("0.3<",tau,"<=0.4")),expression(paste("0.4<",tau,"<=0.5")),expression(paste("0.5<",tau,"<=0.6")),expression(paste("0.6<",tau,"<=0.7")),expression(paste("0.7<",tau,"<=0.8")),expression(paste("0.8<",tau,"<=0.9")),expression(paste("0.9<",tau,"<=1.0")) )), col=vcol[10:4], bg=rgb(1,1,1,0.8),lty=1,box.col=rgb(1,1,1,0.8), cex=0.7)
        legend("bottomright",rev(c(expression(tau <= 0.1),expression(paste("0.1<",tau,"<=0.2")),expression(paste("0.2<",tau,"<=0.3")))), col=vcol[3:1], bg=rgb(1,1,1,0.8),lty=1,box.col=rgb(1,1,1,0.8), cex=0.7)
        box()
        for(i in 1:B2){
            foo=which(d2<=i/B2 & d2>(i-1)/B2);
            y=c();
            for(j in x){
                y=c(y,mean(d1[foo] <= j));
            }
            lines(x,y,col=vcol[i]);
        }

        dev.off()
    }

    

}

{## Supp figure comparing 10PCs vs population labels
    H2MATSs = array(,dim=c(2,2,20)) ## [pop, PC; raw, perm; K=1..20]

    filename=paste("merge_1000000_K=20_MAC=1_ngPC=0_npPC=0_CV=12-13-14-15_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="");
    h2POP=read.table(filename);
    H2MATSs[1,1,1:20] = unlist(h2POP[1,1:20])

    filename=paste("merge_1000000_K=20_MAC=1_ngPC=0_npPC=0_CV=12-13-14-15_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="");
    h2POP=read.table(filename);
    H2MATSs[1,2,1:20] = unlist(h2POP[1,1:20])

    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NORESID.he.meanQuant",sep="");
    h2PC=read.table(filename);
    H2MATSs[2,1,1:20] = unlist(h2PC[1,1:20])

    filename=paste("merge_1000000_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_PERM=3_NORESID.he.meanQuant",sep="");
    h2PC=read.table(filename);
    H2MATSs[2,2,1:20] = unlist(h2PC[1,1:20])


    H2c = array(,dim=c(2,20)) ## [pop, PC; K=1..20]
    H2c[1,] = H2MATSs[1,1,]-H2MATSs[1,2,];
    H2c[2,] = H2MATSs[2,1,]-H2MATSs[2,2,];

    h2c = array(,dim=c(2,20))
    h2MATSs = H2MATSs ## not cum [pop, PC; raw, perm; K=1..20]
   
    h2c[1,1] = H2c[1,1]
    h2c[2,1] = H2c[2,1]
    for(i in 2:20){
        h2c[1,i] = H2c[1,i]-H2c[1,i-1]
        h2c[2,i] = H2c[2,i]-H2c[2,i-1]
        for(j in 1:2){
            for(k in 1:2){
                h2MATSs[j,k,i] = H2MATSs[j,k,i]-H2MATSs[j,k,i-1];
            }
        }
    }

    d = c(H2MATSs[1,1,1], H2MATSs[1,2,1], h2c[1,1], H2MATSs[2,1,1], H2MATSs[2,2,1], h2c[2,1])

    pdf(paste(FIGDIR,"SUPP/S24_PopCovar_vs_10PCs.pdf",sep=""),height=2,width=2.25,pointsize=10)
    par(mar=c(5.1,2.8,0.4,0.3))
    plot(c(1:3,1:3),d,ylim=c(0,0.03),pch=20,col=c(1,1,1,2,2,2),xaxt='n',yaxt='n',xlab="",ylab="")

    abline(h=seq(0,0.03,by=0.005),col=rgb(0,0,0,0.1),lwd=0.75)
    axis(side=1,at=1:3,labels=c("Raw","trans perm","corrected"),las=3,cex=0.7)
    axis(side=2,padj=1.3, cex.axis=0.7,at=seq(0,0.03,by=0.005),labels=NA)
    axis(side=2,padj=1.3, cex.axis=0.7,at=c(0,0.01,0.02,0.03),labels=c(0,0.01,0.02,0.03))
        
    mtext(side=2,expression(paste({widehat(h^2)[singletons]},sep="")),line=1.15)

    legend("topright",c("Population covariate","10 PCs"), col=c(1,2), pch=20, bg=rgb(1,1,1,0.8),cex=0.8)
    box()

    dev.off()
}
