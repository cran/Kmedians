utils::globalVariables(c("K", "x", "y", "coul", "couleur",'SE','Criterion'))

KmedianK=function(X,K=3,ninit=0,niter=20,init=TRUE)
{
  d=ncol(X)
  n=nrow(X)
  finalcenters=matrix(0,nrow=K,ncol=ncol(X))
  centers=matrix(0,nrow=K,ncol=ncol(X))
  dist=matrix(0,nrow=nrow(X),ncol=K)
  finaldist=dist
  cluster=rep(0,nrow(X))
  finalcluster=rep(0,nrow(X))
  finalniter=0
  finalSE=10^10
  if (init==T)
  {
    dist=matrix(0,nrow=nrow(X),ncol=K)
    Sigma=c()
    cluster=genieclust::genie(X,k=K)
    for (k in 1:K)
    {
      I=which(cluster==k)
      if (length(I)>1)
      {
        centers[k,]=Gmedian::Weiszfeld(X[I,])$median
      }
      if (length(I) == 1)
      {
        centers[k,]=X[I,]
      }
    }
    l=0
    distclust=1
    while (l < niter && distclust > 0)
    {
      cluster0=cluster
      l=l+1
      for (k in 1:K)
      {
        dist[,k]=sqrt(rowSums((X-matrix(rep(centers[k,],n),nrow=n,byrow=T))^2))
      }
      cluster=apply(dist,1,FUN=which.min)
      SE=0
      for (k in 1:K)
      {
        I=which(cluster==k)
        SE=SE+sum(dist[I,k])
        if (length(I)>1)
        {
          centers[k,]=Gmedian::Weiszfeld(X[I,])$median
        }
        if (length(I) == 1)
        {
          centers[k,]=X[I,]
        }
      }
      distclust=sum((cluster-cluster0)^2)
      if (SE < finalSE)
      {
        finalcenters=centers
        finalcluster=cluster
        finaldist=dist
        finalSE=SE
      }
    }

  }
  if (ninit>0)
  {
    for (o in 1:ninit)
    {
      dist=matrix(0,nrow=nrow(X),ncol=K)
      Sigma=c()
      centers=X[sample(1:n,K),]
      l=0
      distclust=1
      while (l < niter && distclust > 0)
      {
        cluster0=cluster
        l=l+1
        for (k in 1:K)
        {
          dist[,k]=sqrt(rowSums((X-matrix(rep(centers[k,],n),nrow=n,byrow=T))^2))
        }
        cluster=apply(dist,1,FUN=which.min)
        distclust=sum((cluster-cluster0)^2)
        SE=0
        for (k in 1:K)
        {
          I=which(cluster==k)
          SE=SE+sum(dist[I,k])
          if (length(I)>1)
          {
            centers[k,]=Weiszfeld(X[I,])$median
          }
          if (length(I) == 1)
          {
            centers[k,]=X[I,]
          }
        }
        if (SE < finalSE)
        {
          finalcenters=centers
          finalcluster=cluster
          finaldist=dist
          finalSE=SE
        }
      }
    }
  }
  resultat=list(cluster=finalcluster,centers=finalcenters,
                SE=finalSE)
  return(resultat)
}


OnlineKmedianK=function(X,K=3,ninit=0,niter=20,init=TRUE)
{
  d=ncol(X)
  n=nrow(X)
  finalcenters=matrix(0,nrow=K,ncol=ncol(X))
  centers=matrix(0,nrow=K,ncol=ncol(X))
  dist=matrix(0,nrow=nrow(X),ncol=K)
  finaldist=dist
  cluster=rep(0,nrow(X))
  finalcluster=rep(0,nrow(X))
  finalniter=0
  finalSE=10^10
  if (init==T)
  {
    dist=matrix(0,nrow=nrow(X),ncol=K)
    Sigma=c()
    cluster=genieclust::genie(X,k=K)
    for (k in 1:K)
    {
      I=which(cluster==k)
      if (length(I)>1)
      {
        centers[k,]=Gmedian::Weiszfeld(X[I,])$median
      }
      if (length(I) == 1)
      {
        centers[k,]=X[I,]
      }
    }
    l=0
    resultat=kGmedian(X,ncenters=centers,nstart=1,nstartkmeans=0,iter.max=0)
    finalcluster=resultat$cluster
    finalcenters=resultat$centers
    SE=0
    for (k in 1:K)
    {
      dist[,k]=sqrt(rowSums((X-matrix(rep(resultat$centers[k,],n),nrow=n,byrow=T))^2))
      I=which(resultat$cluster==k)
      SE=SE+sum(dist[I,k])

    }
    finaldist=dist
    finalSE=SE
  }
  if (ninit>0)
  {
    for (o in 1:ninit)
    {
      dist=matrix(0,nrow=nrow(X),ncol=K)
      Sigma=c()
      centers=X[sample(1:n,K),]
      l=0
      distclust=1
      resultat=kGmedian(X,ncenters=centers,nstart=1,nstartkmeans=0,iter.max=0)
      SE=0
      for (k in 1:K)
      {
        dist[,k]=sqrt(rowSums((X-matrix(rep(resultat$centers[k,],n),nrow=n,byrow=T))^2))
        I=which(resultat$cluster==k)
        SE=SE+sum(dist[I,k])

      }
      if (SE < finalSE)
      {
        finalcluster=resultat$cluster
        finalcenters=resultat$centers
        finaldist=resultat$withinsrs
        finalSE=resultat=sum(resultat$withinsrs)
      }
    }
  }
  resultat=list(cluster=finalcluster,centers=finalcenters,dist=finaldist,
                SE=finalSE)
  return(resultat)
}


KGmedianK=function(X,K=3,ninit=0,niter=20,init=TRUE)
{
  d=ncol(X)
  n=nrow(X)
  finalcenters=matrix(0,nrow=K,ncol=ncol(X))
  centers=matrix(0,nrow=K,ncol=ncol(X))
  dist=matrix(0,nrow=nrow(X),ncol=K)
  finaldist=dist
  cluster=rep(0,nrow(X))
  finalcluster=rep(0,nrow(X))
  finalniter=0
  finalSE=10^10
  if (init==T)
  {
    dist=matrix(0,nrow=nrow(X),ncol=K)
    Sigma=c()
    cluster=genieclust::genie(X,k=K)
    for (k in 1:K)
    {
      I=which(cluster==k)
      if (length(I)>1)
      {
        centers[k,]=Gmedian::Gmedian(X[I,])
      }
      if (length(I) == 1)
      {
        centers[k,]=X[I,]
      }
    }
    l=0
    distclust=1
    while (l < niter && distclust > 0)
    {
      cluster0=cluster
      l=l+1
      for (k in 1:K)
      {
        dist[,k]=sqrt(rowSums((X-matrix(rep(centers[k,],n),nrow=n,byrow=T))^2))
      }
      cluster=apply(dist,1,FUN=which.min)
      SE=0
      for (k in 1:K)
      {
        I=which(cluster==k)
        SE=SE+sum(dist[I,k])
        if (length(I)>1)
        {
          centers[k,]=Gmedian::Gmedian(X[I,])
        }
        if (length(I) == 1)
        {
          centers[k,]=X[I,]
        }
      }
      distclust=sum((cluster-cluster0)^2)
      if (SE < finalSE)
      {
        finalcenters=centers
        finalcluster=cluster
        finaldist=dist
        finalSE=SE
      }
    }

  }
  if (ninit>0)
  {
    for (o in 1:ninit)
    {
      dist=matrix(0,nrow=nrow(X),ncol=K)
      Sigma=c()
      centers=X[sample(1:n,K),]
      l=0
      distclust=1
      while (l < niter && distclust > 0)
      {
        cluster0=cluster
        l=l+1
        for (k in 1:K)
        {
          dist[,k]=sqrt(rowSums((X-matrix(rep(centers[k,],n),nrow=n,byrow=T))^2))
        }
        cluster=apply(dist,1,FUN=which.min)
        distclust=sum((cluster-cluster0)^2)
        SE=0
        for (k in 1:K)
        {
          I=which(cluster==k)
          SE=SE+sum(dist[I,k])
          if (length(I)>1)
          {
            centers[k,]=Gmedian(X[I,])
          }
          if (length(I) == 1)
          {
            centers[k,]=X[I,]
          }
        }
        if (SE < finalSE)
        {
          finalcenters=centers
          finalcluster=cluster
          finaldist=dist
          finalSE=SE
        }
      }
    }
  }
  resultat=list(cluster=finalcluster,centers=finalcenters,
                SE=finalSE)
  return(resultat)
}

Kmed=function(X,K=3,ninit=0,niter=20,init=TRUE,method='Offline')
{
  if (method=='Offline')
  {
    resultat=KmedianK(X,K=K,ninit=ninit,
                      niter=niter,init=init)
  }
  if (method=='Semi-Online')
  {
    resultat=KGmedianK(X,K=K,ninit=ninit,
                       niter=niter,init=init)
  }
  if (method=='Online')
  {
    resultat=OnlineKmedianK(X,K=K,ninit=ninit,
                            niter=niter,init=init)
  }
  return(resultat)
}



Kmedians=function(X,nclust=1:15,ninit=0,niter=20,
                  method='Offline', init=TRUE,par=TRUE)
{
  if (length(nclust)==1)
  {
    resultat=Kmed(X,K=nclust,ninit=ninit,niter=niter,
                  method=method,init=init)
    resultat2=list(bestresult=resultat,allresults=resultat,SE=resultat$SE,data=X,nclust=nclust)
  }
  if (length(nclust)>1)
  {
    if (par ==TRUE)
    {
      numCores = min(detectCores()-2,length(nclust))
      cl = makeCluster(numCores  )
      registerDoParallel(cl)
      resultat=foreach(K=nclust,.multicombine = TRUE ,.inorder = TRUE,
                       .packages = c('genieclust','Gmedian'),
                       .combine='list')  %dopar%
        {
          resultatk=Kmed(X,K=K,ninit=ninit,niter=niter,
                         method=method,init=init)
          return(resultatk)
        }

      stopCluster(cl)
    }

    if (par ==FALSE)
    {
      resultat=foreach(K=nclust,.multicombine = TRUE, .inorder = TRUE,
                       .packages = c('genieclust','Gmedian'),
                       .combine='list')  %do%
        {
          resultatk=Kmed(X,K=K,ninit=ninit,niter=niter,
                         method=method,init=init)
          return(resultatk)
        }

    }
    SE=c()
    for (i in 1:length(nclust))
    {
      SE=c(SE,resultat[[i]]$SE)
    }
    if (length(nclust)<= 10)
    {
      message('Capushe cannot be used, nclust must be of length larger that 10. No model has been selected.')
      resultat2=list(allresults=resultat,SE=SE,data=X,nclust=nclust)
    }
    if (length(nclust) > 10)
    {
      cap=capushe(cbind(nclust,sqrt(nclust),nclust,SE))
      Kopt=as.numeric(cap@DDSE@model)
      I=which(nclust==Kopt)
      bestresult=resultat[[I]]
      resultat2=list(bestresult=bestresult,allresults=resultat,
                     SE=SE,cap=cap,Ksel=Kopt,data=X,nclust=nclust)
    }

  }
  return(resultat2)
}



kmeansmaison=function(X,K=5,ninit=1,niter=20)
{
  rkm=kmeans(X,centers=as.numeric(K),nstart=ninit,
             iter.max=niter)
  resultat=list(cluster=rkm$cluster,centers=rkm$centers,
                SE=rkm$tot.withinss)
  return(resultat)
}

Kmeans=function(X,nclust=1:15,ninit=1,niter=20,par=TRUE)
{
  if (length(nclust)==1)
  {
    resultat=kmeansmaison(X,K=nclust,ninit=ninit,
                          niter=niter)
    resultat2=list(bestresult=resultat,allresults=resultat,SE=resultat$SE,data=X,nclust=nclust)
  }
  if (length(nclust)>1)
  {
    if (par ==TRUE)
    {
      numCores = min(detectCores()-2,length(nclust))
      cl = makeCluster(numCores)
      registerDoParallel(cl)
      resultat=foreach(K=nclust,.multicombine = TRUE ,.inorder = TRUE,
                       .packages = c('genieclust','Gmedian','stats'),
                       .combine='list')  %dopar%
        {
          resultatk=kmeansmaison(X,K=K,ninit=ninit,niter=niter)
          return(resultatk)
        }

      stopCluster(cl)
    }

    if (par ==FALSE)
    {
      resultat=foreach(K=nclust,.multicombine = TRUE, .inorder = TRUE,
                       .packages = c('genieclust','Gmedian'),
                       .combine='list')  %do%
        {
          resultatk=kmeansmaison(X,K=K,ninit=ninit,niter=niter)
          return(resultatk)
        }

    }
    SE=c()
    for (i in 1:length(nclust))
    {
      SE=c(SE,resultat[[i]]$SE)
    }
    if (length(nclust)<= 10)
    {
      message('Capushe cannot be used, nclust must be of length larger that 10. No model has been selected.')
      resultat2=list(allresults=resultat,SE=SE,data=X,nclust=nclust)
    }
    if (length(nclust) > 10)
    {
      cap=capushe(cbind(nclust,sqrt(nclust),nclust,SE))
      Kopt=as.numeric(cap@DDSE@model)
      I=which(nclust==Kopt)
      bestresult=resultat[[I]]
      resultat2=list(bestresult=bestresult,allresults=resultat,
                     SE=SE,cap=cap,Ksel=Kopt,data=X,nclust=nclust)

    }
  }
  return(resultat2)
}


gen_K=function(n=500,d=5,K=3,pcont=0,df=1,cont="Student",min=-5,max=5,radius=5)
{
  Sigma= (diag(d)*2)
  X=c()
  Tclassif=c()
  mean=rep(0,d)
  for (k in 1:K)
  {
    Z=rnorm(d)
    Z=radius*Z/sqrt(sum(Z^2))
    X=rbind(X,mvtnorm::rmvnorm(n ,mean=Z))
    Tclassif=c(Tclassif,rep(k,n))
  }

  if (pcont>0)
  {
    if(cont=='Unif')
    {
      Z=matrix(runif(pcont*n*K*d,min,max),ncol=d)
    }
    if (cont=='Student')
    {
      Z=matrix(rt(pcont*n*K*d,df=df),ncol=d)
    }
    I=sample(1:(K*n),size=pcont*n*K)
    X[I,]=Z
    Tclassif[I]="outliers"
  }


  mel=sample.int(K*n)
  Tclassif=Tclassif[mel]
  X=X[mel,]
  return(list(X=X,cluster=Tclassif))
}

Kplot=function(a,propplot=0.95,graph=c('Two_Dim','Capushe','Profiles','SE','Criterion'),bestresult=TRUE,Ksel=FALSE,bycluster=TRUE)
{
  nclust=a$nclust
  X=a$data
  d=ncol(X)

  if (length(nclust) <=10)
  {
    message('no model has been selected. The length of nclust must be larger than 10')
  }
  if (d==1)
  {
    message('d must be larger than 1')
  }
  if (d >1 & (length(nclust)>10) & bestresult==T)
  {
    cluster=a$bestresult$cluster
    if (d > 2){
      Gmed=Gmedian::WeiszfeldCov(X,scores = 2)
      vec=Gmed$vectors
    }
    if (d==2)
    {
      Gmed=Gmedian::WeiszfeldCov(X,scores = F)
      eig=eigen(Gmed$covmedian)
      vec=eig$vectors
    }
    med0=Gmed$median
    if (length(which(graph=='Two_Dim'))>0 & bycluster==TRUE){
      med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
      Xtrans=matrix(0,nrow=nrow(X),ncol=2)
      for(i in 1:nrow(X))
      {
        Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
        Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
      }
      Xtr=c()
      clusplot=c()
      centers=a$bestresult$centers
      for (k in 1:a$Ksel){
        distk=c()
        medk=c(sum(centers[k,]*vec[,1])-med[1],sum(centers[k,]*vec[,2])-med[2])
        I=which(a$bestresult$cluster==k)
        Xtemp=Xtrans[I,]
        if (length(I)>1){
          distk=rowSums((X[I,]-centers[k,])^2)
          J=order(distk,decreasing = FALSE)[1:(floor((propplot)*length(I)))]
          Xtr=rbind(Xtr,Xtemp[J,])
          clusplot=c(clusplot,rep(k,floor((propplot)*length(I))))
        }
        if (length(I)==1){
          Xtr=rbind(Xtr,Xtemp)
          clusplot=c(clusplot,k)
        }
      }
      dataplot=data.frame(x=Xtr[,1],y=Xtr[,2],K=as.character(clusplot))
      plo=  ggplot(data=dataplot) +
        scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
        scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
        geom_point(aes(x = x, y = y,color = K))
      print(plo)
    }
    if (length(which(graph=='Two_Dim'))>0 & bycluster==FALSE){
      med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
      Xtrans=matrix(0,nrow=nrow(X),ncol=2)
      for(i in 1:nrow(X))
      {
        Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
        Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
      }
      I=order(abs(Xtrans[,1]-med[1]),decreasing=TRUE)[1:(0.5*(1-propplot)*nrow(X))]
      J=order(abs(Xtrans[,2]-med[2]),decreasing=TRUE)[1:(0.5*(1-propplot)*nrow(X))]
      dataplot=data.frame(x=Xtrans[-c(I,J),1],y=Xtrans[-c(I,J),2],K=as.character(cluster[-c(I,J)]))
      plo=  ggplot(data=dataplot) +
        scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
        scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
        geom_point(aes(x = x, y = y,color = K))
      print(plo)
    }
    if (length(which(graph=='Capushe'))>0){
      capushe::plot(a$cap@DDSE,newwindow=FALSE)
    }
    if (length(which(graph=='SE'))>0){
      data_SE=data.frame(K=nclust,SE=a$SE)
      se_plot=  ggplot(data_SE,aes(x=K , y=SE))+
        scale_x_continuous('K') +
        scale_y_continuous('Sum of errors') +
        geom_line(mapping=aes(x=K , y=SE),show.legend = F)
      print(se_plot)
    }
    if (length(which(graph=='Criterion'))>0){
      aopt=a$cap@DDSE@graph$reg$coefficients[2]
      data_crit=data.frame(K=nclust,Criterion=(2*aopt*a$cap@DDSE@graph$pen+ a$SE)/nrow(X))
      crit_plot=  ggplot(data_crit,aes(x=K , y=Criterion))+
        scale_x_continuous('K') +
        geom_line(mapping=aes(x=K , y=Criterion),show.legend = F)
      print(crit_plot)
    }
    if (length(which(graph=='Profiles'))>0 & bycluster==TRUE){
      Xt=c()
      centers=a$bestresult$centers
      clustplot=c()
      for (k in 1:a$Ksel){
        I=which(a$bestresult$cluster==k)
        Xtemp=X[I,]
        dist=c()
        if (length(I)>1){
          dist=rowSums((X[I,]-centers[k,])^2)
          J=order(dist,decreasing = FALSE)[1:(floor((propplot)*length(I)))]
          Xt=rbind(Xt,Xtemp[J,])
          clustplot=c(clustplot,rep(k,floor((propplot)*length(I))))
        }
        if (length(I)==1){
          Xt=rbind(Xt,Xtemp)
          clustplot=c(clustplot,k)
        }
      }
      if (d>2){
        dataclust=data.frame(x=1:d,
                             y=t(Xt))
        datacen=data.frame(xt=1:d,t(a$bestresult$centers))
        ind=c()
        Kclu=c()
        Kcen=c()
        indK=c()
        col=c()
        for (o in 1:nrow(Xt))
        {
          ind=c(ind,rep(o,d))
        }
        for (o in 1:a$Ksel)
        {
          I=(which(clustplot==o))
          Kclu=c(Kclu,as.character(rep(clustplot[I],d)))
          Kcen=c(Kcen,rep(as.character(o),d))
          indK=c(indK,rep((o+max(ind)),d))
        }
        data_new <- melt(dataclust, id = c("x"))
        data_cen <- melt(datacen, id = c("xt"))
        datanew=data.frame( x=c(data_new$x,data_cen$xt),y=c(data_new$value,data_cen$value),
                            couleur=c(ind,indK),Kclu=c(Kclu,Kcen),
                            coul=c(rep('grey',length(data_new$x)),rep('red',length(data_cen$xt))))

        prof_plot=  ggplot(datanew,aes(x=x , y=y,color=coul))+
          scale_color_manual(values=c('grey','red'))+
          scale_x_continuous('') +
          scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
          geom_line(mapping=aes(x=x , y=y , group=couleur),show.legend = F,alpha=0.4)+
          #  geom_line(aes(color=(coul)))+
          facet_wrap(~as.factor(Kclu))
        print(prof_plot)
      }
      if (d==2){
        dataclust=data.frame(x=1:d,
                             y=t(Xt))
        Kclu=as.character(clustplot)
        Kcen=as.character(1:(a$Ksel))
        data_clust=data.frame(x=c(Xt[,1],(a$bestresult$centers)[,1]),y=c(Xt[,2],(a$bestresult$centers)[,2]),
                              K=c(Kclu,Kcen),coul=c(rep('grey',nrow(Xt)),rep('red',a$Ksel)))
        prof_plot=  ggplot(data_clust,aes(x=x , y=y,color=coul))+
          scale_color_manual(values=c('grey','red'))+
          scale_x_continuous('') +
          scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
          geom_point(mapping=aes(x=x , y=y),show.legend = F,alpha=0.4)+
          #  geom_line(aes(color=(coul)))+
          facet_wrap(~as.factor(K))
        print(prof_plot)
      }
    }
    if (length(which(graph=='Profiles'))>0 & bycluster==FALSE){
      dist=c()
      for (i in 1:nrow(X))
      {
        dist=c(dist,sum((X[i,]-Gmed$median)^2))
      }
      ord=order(dist,decreasing = FALSE)
      Xt=X[ord,]
      Xt=Xt[1:floor(propplot*nrow(X)),]
      cluster2=cluster[ord]
      cluster2=cluster2[1:floor(propplot*nrow(X))]
      if (d>2){
        dataclust=data.frame(x=1:d,
                             y=t(Xt))
        datacen=data.frame(xt=1:d,t(a$bestresult$centers))
        ind=c()
        Kclu=c()
        for (o in 1:nrow(Xt))
        {
          ind=c(ind,rep(o,d))
          Kclu=c(Kclu,rep(as.character(cluster2[o]),d))
        }
        Kcen=c()
        indK=c()
        col=c()
        for (o in 1:a$Ksel)
        {
          I=length(which(Kclu==as.character(o)))
          Kcen=c(Kcen,rep(as.character(o),d))
          indK=c(indK,rep((o+max(ind)),d))
        }
        data_new <- melt(dataclust, id = c("x"))
        data_cen <- melt(datacen, id = c("xt"))
        datanew=data.frame( x=c(data_new$x,data_cen$xt),y=c(data_new$value,data_cen$value),
                            couleur=c(ind,indK),Kclu=c(Kclu,Kcen),
                            coul=c(rep('grey',length(data_new$x)),rep('red',length(data_cen$xt))))

        prof_plot=  ggplot(datanew,aes(x=x , y=y,color=coul))+
          scale_color_manual(values=c('grey','red'))+
          scale_x_continuous('') +
          scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
          geom_line(mapping=aes(x=x , y=y , group=couleur),show.legend = F,alpha=0.4)+
          #  geom_line(aes(color=(coul)))+
          facet_wrap(~as.factor(Kclu))
        print(prof_plot)
      }
      if (d==2){
        dataclust=data.frame(x=1:d,
                             y=t(Xt))
        Kclu=c()
        for (o in 1:nrow(Xt))
        {
          Kclu=c(Kclu,as.character(cluster2[o]))
        }
        Kcen=c()
        for (o in 1:a$Ksel)
        {
          Kcen=c(Kcen,as.character(o))
        }
        data_clust=data.frame(x=c(Xt[,1],(a$bestresult$centers)[,1]),y=c(Xt[,2],(a$bestresult$centers)[,2]),
                              K=c(Kclu,Kcen),coul=c(rep('grey',nrow(Xt)),rep('red',a$Ksel)))
        prof_plot=  ggplot(data_clust,aes(x=x , y=y,color=coul))+
          scale_color_manual(values=c('grey','red'))+
          scale_x_continuous('') +
          scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
          geom_point(mapping=aes(x=x , y=y),show.legend = F,alpha=0.4)+
          #  geom_line(aes(color=(coul)))+
          facet_wrap(~as.factor(K))
        print(prof_plot)
      }
    }
  }
  if (d >1 & Ksel!=FALSE)
  {
    if (length(which(nclust==Ksel))==0){
      message('K is not in nclust')
    }
    if (length(which(nclust==Ksel))>0){
      Knclu=which(nclust==Ksel)
      cluster=a$allresults[[Knclu]]$cluster
      if (length(nclust) >1){
        cluster=a$allresults[[Knclu]]$cluster
      }
      if (length(nclust)==1){
        cluster=a$allresults[[Knclu]]$cluster
      }
      if (d > 2){
        Gmed=Gmedian::WeiszfeldCov(X,scores = 2)
        vec=Gmed$vectors
      }
      if (d==2)
      {
        Gmed=Gmedian::WeiszfeldCov(X,scores = F)
        eig=eigen(Gmed$covmedian)
        vec=eig$vectors
      }
      med0=Gmed$median
      if (length(which(graph=='Two_Dim'))>0 & bycluster==TRUE){
        med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
        Xtrans=matrix(0,nrow=nrow(X),ncol=2)
        for(i in 1:nrow(X))
        {
          Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
          Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
        }
        Xtr=c()
        clusplot=c()
        centers=a$allresults[[Knclu]]$centers
        for (k in 1:Ksel){
          distk=c()
          medk=c(sum(centers[k,]*vec[,1])-med[1],sum(centers[k,]*vec[,2])-med[2])
          I=which(a$allresults[[Knclu]]$cluster==k)
          Xtemp=Xtrans[I,]
          if (length(I)>1){
            distk=rowSums((X[I,]-centers[k,])^2)
            J=order(distk,decreasing = FALSE)[1:(floor((propplot)*length(I)))]
            Xtr=rbind(Xtr,Xtemp[J,])
            clusplot=c(clusplot,rep(k,floor((propplot)*length(I))))
          }
          if (length(I)==1){
            Xtr=rbind(Xtr,Xtemp)
            clusplot=c(clusplot,k)
          }

        }
        dataplot=data.frame(x=Xtr[,1],y=Xtr[,2],K=as.character(clusplot))
        plo=  ggplot(data=dataplot) +
          scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
          scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
          geom_point(aes(x = x, y = y,color = K))
        print(plo)
      }
      if (length(which(graph=='Two_Dim'))>0 & bycluster==FALSE){
        med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
        Xtrans=matrix(0,nrow=nrow(X),ncol=2)
        for(i in 1:nrow(X))
        {
          Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
          Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
        }
        I=order(abs(Xtrans[,1]-med[1]),decreasing=TRUE)[1:(0.5*(1-propplot)*nrow(X))]
        J=order(abs(Xtrans[,2]-med[2]),decreasing=TRUE)[1:(0.5*(1-propplot)*nrow(X))]
        dataplot=data.frame(x=Xtrans[-c(I,J),1],y=Xtrans[-c(I,J),2],K=as.character(cluster[-c(I,J)]))
        plo=  ggplot(data=dataplot) +
          scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
          scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
          geom_point(aes(x = x, y = y,color = K))
        print(plo)
      }
      if (length(which(graph=='Profiles'))>0 & bycluster==TRUE){
        Xt=c()
        centers=a$allresults[[Knclu]]$centers
        clustplot=c()
        for (k in 1:Ksel){
          I=which(a$allresults[[Knclu]]$cluster==k)
          Xtemp=X[I,]
          dist=c()
          if (length(I)>1){
            dist=rowSums((X[I,]-centers[k,])^2)
            J=order(dist,decreasing = FALSE)[1:(floor((propplot)*length(I)))]
            Xt=rbind(Xt,Xtemp[J,])
            clustplot=c(clustplot,rep(k,floor((propplot)*length(I))))
          }
          if (length(I)==1){
            Xt=rbind(Xt,Xtemp)
            clustplot=c(clustplot,k)
          }
        }
        if (d>2){
          dataclust=data.frame(x=1:d,
                               y=t(Xt))
          datacen=data.frame(xt=1:d,t(a$allresults[[Knclu]]$centers))
          ind=c()
          Kclu=c()
          Kcen=c()
          indK=c()
          col=c()
          for (o in 1:nrow(Xt))
          {
            ind=c(ind,rep(o,d))
          }
          for (o in 1:Ksel)
          {
            I=(which(clustplot==o))
            Kclu=c(Kclu,as.character(rep(clustplot[I],d)))
            Kcen=c(Kcen,rep(as.character(o),d))
            indK=c(indK,rep((o+max(ind)),d))
          }
          data_new <- melt(dataclust, id = c("x"))
          data_cen <- melt(datacen, id = c("xt"))
          datanew=data.frame( x=c(data_new$x,data_cen$xt),y=c(data_new$value,data_cen$value),
                              couleur=c(ind,indK),Kclu=c(Kclu,Kcen),
                              coul=c(rep('grey',length(data_new$x)),rep('red',length(data_cen$xt))))

          prof_plot=  ggplot(datanew,aes(x=x , y=y,color=coul))+
            scale_color_manual(values=c('grey','red'))+
            scale_x_continuous('') +
            scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
            geom_line(mapping=aes(x=x , y=y , group=couleur),show.legend = F,alpha=0.4)+
            #  geom_line(aes(color=(coul)))+
            facet_wrap(~as.factor(Kclu))
          print(prof_plot)
        }
        if (d==2){
          dataclust=data.frame(x=1:d,
                               y=t(Xt))
          Kclu=as.character(clustplot)
          Kcen=as.character(1:(Ksel))
          data_clust=data.frame(x=c(Xt[,1],(a$allresults[[Knclu]]$centers)[,1]),y=c(Xt[,2],(a$allresults[[Knclu]]$centers)[,2]),
                                K=c(Kclu,Kcen),coul=c(rep('grey',nrow(Xt)),rep('red',Ksel)))
          prof_plot=  ggplot(data_clust,aes(x=x , y=y,color=coul))+
            scale_color_manual(values=c('grey','red'))+
            scale_x_continuous('') +
            scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
            geom_point(mapping=aes(x=x , y=y),show.legend = F,alpha=0.4)+
            #  geom_line(aes(color=(coul)))+
            facet_wrap(~as.factor(K))
          print(prof_plot)
        }
      }
      if (length(which(graph=='Profiles'))>0 & bycluster==FALSE){
        dist=c()
        for (i in 1:nrow(X))
        {
          dist=c(dist,sum((X[i,]-Gmed$median)^2))
        }
        ord=order(dist,decreasing = FALSE)
        Xt=X[ord,]
        Xt=Xt[1:floor(propplot*nrow(X)),]
        cluster2=cluster[ord]
        cluster2=cluster2[1:floor(propplot*nrow(X))]
        if (d>2){
          dataclust=data.frame(x=1:d,
                               y=t(Xt))
          datacen=data.frame(xt=1:d,t(a$allresults[[Knclu]]$centers))
          ind=c()
          Kclu=c()
          for (o in 1:nrow(Xt))
          {
            ind=c(ind,rep(o,d))
            Kclu=c(Kclu,rep(as.character(cluster2[o]),d))
          }
          Kcen=c()
          indK=c()
          col=c()
          for (o in 1:Ksel)
          {
            I=length(which(Kclu==as.character(o)))
            Kcen=c(Kcen,rep(as.character(o),d))
            indK=c(indK,rep((o+max(ind)),d))
          }
          data_new <- melt(dataclust, id = c("x"))
          data_cen <- melt(datacen, id = c("xt"))
          datanew=data.frame( x=c(data_new$x,data_cen$xt),y=c(data_new$value,data_cen$value),
                              couleur=c(ind,indK),Kclu=c(Kclu,Kcen),
                              coul=c(rep('grey',length(data_new$x)),rep('red',length(data_cen$xt))))

          prof_plot=  ggplot(datanew,aes(x=x , y=y,color=coul))+
            scale_color_manual(values=c('grey','red'))+
            scale_x_continuous('') +
            scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
            geom_line(mapping=aes(x=x , y=y , group=couleur),show.legend = F,alpha=0.4)+
            #  geom_line(aes(color=(coul)))+
            facet_wrap(~as.factor(Kclu))
          print(prof_plot)
        }
        if (d==2){
          dataclust=data.frame(x=1:d,
                               y=t(Xt))
          Kclu=c()
          for (o in 1:nrow(Xt))
          {
            Kclu=c(Kclu,as.character(cluster2[o]))
          }
          Kcen=c()
          for (o in 1:Ksel)
          {
            Kcen=c(Kcen,as.character(o))
          }
          data_clust=data.frame(x=c(Xt[,1],(a$allresults[[Knclu]]$centers)[,1]),y=c(Xt[,2],(a$allresults[[Knclu]]$centers)[,2]),
                                K=c(Kclu,Kcen),coul=c(rep('grey',nrow(Xt)),rep('red',Ksel)))
          prof_plot=  ggplot(data_clust,aes(x=x , y=y,color=coul))+
            scale_color_manual(values=c('grey','red'))+
            scale_x_continuous('') +
            scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
            geom_point(mapping=aes(x=x , y=y),show.legend = F,alpha=0.4)+
            #  geom_line(aes(color=(coul)))+
            facet_wrap(~as.factor(K))
          print(prof_plot)
        }
      }
    }
  }
}
