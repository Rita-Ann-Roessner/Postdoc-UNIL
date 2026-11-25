.sourced_for_benchmark <- TRUE; source("ESM_motif_builder.R")
RUN <- "TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles"
EP <- c("A0201_ELAGIGILTV","A0201_GILGFVFTL","A0201_LLWNGPMAV")
eff <- function(p){p<-p[p>0];exp(-sum(p*log(p)))}
# posterior with a given weight vector and shrinkage alpha
post_w <- function(sc,w,col,prior,alpha){
  g<-sc[[col]];k<-!is.na(g);pg<-sum(w)/length(w)
  nt<-tapply(w[k],g[k],sum);na<-table(g[k]);gn<-names(prior)
  ntv<-vapply(gn,function(x) if(x%in%names(nt)) as.numeric(nt[[x]]) else 0,numeric(1))
  nav<-vapply(gn,function(x) if(x%in%names(na)) as.numeric(na[[x]]) else 0,numeric(1))
  p<-alpha*unlist(prior)+pmax(0,ntv-nav*pg); if(sum(p)==0) return(NULL); p/sum(p)}
rows<-list()
for (ep in EP){val<-read.csv(file.path(RUN,ep,"validation.csv"));pos<-val[val$Label==1,]
 for(st in 0:1){d<-file.path(RUN,ep,sprintf("step%d",st));thr<-ESM_THRESHOLDS[st+1]
  for(ch in c("A","B")){cn<-if(ch=="A")"alpha" else "beta"
   sf<-file.path(d,sprintf("output_%s.csv",cn));mf<-meta_csv(d,cn);if(!file.exists(sf))next
   sc<-load_esm_scores(sf,mf);s<-sc[[SCORE_COL]];if(sum(s>=thr,na.rm=TRUE)<10)next
   F<-ecdf(s)(s); q<-mean(s<thr,na.rm=TRUE)   # quantile the fixed threshold sits at
   W<-list(hard = as.numeric(s>=thr),
           hinge= as.numeric(s>=thr)*(1+(s-thr)/(0.5*sd(s,na.rm=TRUE))),
           ecdf = F, ecdf2 = F^2, ecdf4 = F^4, ecdf8 = F^8,
           qgate= pmax(0,(F-q)/(1-q)))
   for(gt in c("V","J")){col<-paste0("TR",ch,gt);if(!col%in%names(pos))next
    tf<-table(pos[[col]]);tf<-tf/sum(tf);pr<-vj_baseline_prior[[ch]][[gt]]
    base<-sum((pr/sum(pr))*tf[names(pr)],na.rm=TRUE)
    ph<-post_w(sc,W$hard,col,pr,VJ_PRIOR_STRENGTH); target<-eff(ph)
    for(nm in names(W)){w<-W[[nm]];w[is.na(w)]<-0
      p1<-post_w(sc,w,col,pr,VJ_PRIOR_STRENGTH)
      # re-shrink so this scheme spreads over as many genes as the hard one
      f<-function(a){p<-post_w(sc,w,col,pr,a); if(is.null(p)) return(1e6); eff(p)-target}
      am<-tryCatch(uniroot(f,c(1e-3,1e6))$root,error=function(e)NA)
      pm<-if(is.na(am)) p1 else post_w(sc,w,col,pr,am)
      rows[[length(rows)+1]]<-data.frame(ep=sub("^A0201_","",ep),st=st,col=col,scheme=nm,
        eff_raw=eff(p1), lift_raw=sum(p1*tf[names(p1)],na.rm=TRUE)/base,
        lift_matched=sum(pm*tf[names(pm)],na.rm=TRUE)/base)}}}}}
r<-do.call(rbind,rows)
h<-r[r$scheme=="hard",c("ep","st","col","lift_raw","lift_matched")]
names(h)[4:5]<-c("h_raw","h_match"); k<-merge(r,h,by=c("ep","st","col"))
cat("\n=== 24 comparisons: raw, and re-shrunk to the same effective gene count ===\n")
a<-aggregate(cbind(lift_raw,lift_matched,eff_raw)~scheme,r,mean)
a$win_raw   <-aggregate(I(lift_raw>h_raw)~scheme,k,sum)[,2]
a$win_match <-aggregate(I(lift_matched>h_match)~scheme,k,sum)[,2]
a$med_match <-aggregate(I(lift_matched/h_match)~scheme,k,median)[,2]
print(a[order(-a$lift_matched),],row.names=FALSE,digits=4)
