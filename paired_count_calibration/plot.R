out <- commandArgs(TRUE)[1]
d <- read.csv(file.path(out,'summary.csv'),check.names=FALSE)
colors <- c('#176B91','#D4663A','#23856D','#885F9C')
modes <- c('shared','shared','design_specific','design_specific')
tests <- c('model_based','HC3','model_based','HC3')
panels <- unique(d[,c('depth','exposure_fraction','rho')])
for(kind in c('FDR','coverage')) {
  pdf(file.path(out,if(kind=='FDR')'paired_FDR.pdf' else 'paired_coverage.pdf'),
      width=9,height=max(4,ceiling(nrow(panels)/2)*2.6+1.4),useDingbats=FALSE)
  par(mfrow=c(ceiling(nrow(panels)/2),2),mar=c(3.6,4.3,2.6,1),oma=c(4.3,0,1,0),
      mgp=c(2.2,.65,0),las=1,cex=.85,bty='l',family='sans')
  key <- if(kind=='FDR')'fdp_planned_050' else 'coverage95'
  for(i in seq_len(nrow(panels))) {
    q <- panels[i,]
    z <- d[d$depth==q$depth&d$exposure_fraction==q$exposure_fraction&d$rho==q$rho,,drop=FALSE]
    ymin <- if(kind=='FDR')0 else max(0,min(c(z[[key]],z[[key]]-1.96*z[[paste0(key,'_mcse')]]),na.rm=TRUE)-.01)
    ymax <- if(kind=='FDR')max(.08,max(z$fdp_planned_050_ci_high,na.rm=TRUE)) else 1
    plot(NA,xlim=range(d$n)+c(-5,5),ylim=c(ymin,ymax),
         xlab=if(i>nrow(panels)-2)'Samples' else '',
         ylab=if(kind=='FDR')'BH FDR (paired null)' else '95% interval coverage',xaxt='n')
    axis(1,at=sort(unique(d$n)))
    abline(h=if(kind=='FDR').05 else .95,col='#7B8087',lty=3)
    for(j in 1:4) {
      a <- z[z$weight_mode==modes[j]&z$test==tests[j],,drop=FALSE];a<-a[order(a$n),]
      lo <- if(kind=='FDR')a$fdp_planned_050_ci_low else pmax(0,a[[key]]-1.96*a$coverage95_mcse)
      hi <- if(kind=='FDR')a$fdp_planned_050_ci_high else pmin(1,a[[key]]+1.96*a$coverage95_mcse)
      xx<-a$n+(j-2.5)*1.4
      if(nrow(a)) {
        segments(xx,lo,xx,hi,col=adjustcolor(colors[j],.55),lwd=1.4)
        lines(xx,a[[key]],col=colors[j],lty=if(j%%2)1 else 2,lwd=1.4)
        points(xx,a[[key]],col=colors[j],pch=c(16,17,15,18)[j],cex=.85)
      }
    }
    title(main=sprintf('%s; exposure genes %d%%; rho %.1f',
      if(q$depth=='equal')'Equal depth' else 'IID depth',100*q$exposure_fraction,q$rho),font.main=1,cex.main=.9)
  }
  par(fig=c(0,1,0,.055),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0))
  plot.new()
  legend('center',legend=c('Shared / model-based','Shared / HC3','Design-specific / model-based',
    'Design-specific / HC3'),col=colors,pch=c(16,17,15,18),lty=c(1,2,1,2),ncol=2,bty='n',cex=.8)
  dev.off()
}
