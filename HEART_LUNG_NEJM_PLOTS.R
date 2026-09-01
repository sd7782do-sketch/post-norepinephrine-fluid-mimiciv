# Heart & Lung revision - NEJM-style figures
# Uses the definitive output folder produced by MIMIC_IV_Fluid_04_MASTER_REBUILD.R v1.1
# Required packages: data.table, ggplot2, survey, scales

RESULTS_DIR <- ""
FIGURE_DIR <- ""

run_nejm_figures <- function(results_dir = "", figure_dir = "") {
  pkgs <- c("data.table","ggplot2","survey","scales")
  miss <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(miss)) stop("Installare prima: ",paste(miss,collapse=", "))
  if(!nzchar(results_dir)) {
    if(.Platform$OS.type=="windows") results_dir <- utils::choose.dir(caption="Select definitive MIMIC-IV results folder")
    if(is.null(results_dir)||!nzchar(results_dir)) stop("Impostare RESULTS_DIR.")
  }
  results_dir <- normalizePath(results_dir,winslash="/",mustWork=TRUE)
  if(!nzchar(figure_dir)) figure_dir <- file.path(results_dir,"NEJM_Figures")
  dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
  need <- c("02_screening_flow.csv","04_analytic_dataset_weighted.csv",
            "07_primary_and_secondary_models.csv","08_first_stay_sensitivity.csv",
            "09_threshold_sensitivity.csv","12_spline_adjusted_risk_curve.csv")
  if(any(!file.exists(file.path(results_dir,need)))) stop("File mancanti: ",paste(need[!file.exists(file.path(results_dir,need))],collapse=", "))
  rd <- function(x) data.table::fread(file.path(results_dir,x))
  flow<-rd(need[1]); dat<-rd(need[2]); main<-rd(need[3]); first<-rd(need[4]);thr<-rd(need[5]);curve<-rd(need[6])

  blue <- "#006BA6"; orange <- "#D55E00"; dark <- "#222222"; gray <- "#6B6B6B"; light <- "#E8E8E8"
  theme_nejm <- function(base_size=10) ggplot2::theme_classic(base_size=base_size,base_family="Arial")+
    ggplot2::theme(axis.title=ggplot2::element_text(color=dark),axis.text=ggplot2::element_text(color=dark),
      plot.title=ggplot2::element_text(face="bold",size=base_size+1,hjust=0),
      plot.subtitle=ggplot2::element_text(color=gray,size=base_size-1),
      legend.position="top",legend.title=ggplot2::element_blank(),
      plot.margin=ggplot2::margin(8,12,8,8))
  save_tiff <- function(plot,name,w=7,h=5) ggplot2::ggsave(file.path(figure_dir,name),plot,
    device="tiff",width=w,height=h,units="in",dpi=600,compression="lzw",bg="white")

  # Figure 1: source-level flow, with no unnumbered opening box.
  groups <- dat[,.(n=.N),by=restrictive_1000]
  lower_n <- groups[restrictive_1000==1,n]; higher_n <- groups[restrictive_1000==0,n]
  labels <- c(sprintf("All ICU stays\nN = %s",scales::comma(flow$n[1])),
    sprintf("Adult stays with coded sepsis\nand norepinephrine\nN = %s",scales::comma(flow$n[3])),
    sprintf("Norepinephrine initiated\nduring ICU hours 0-24\nN = %s",scales::comma(flow$n[4])),
    sprintf("Complete required pre-index data\nN = %s",scales::comma(flow$n[5])),
    sprintf("Alive and in ICU at 6-hour landmark\nN = %s (%s patients)",scales::comma(flow$n[6]),scales::comma(data.table::uniqueN(dat$subject_id))))
  y <- 5:1
  boxes <- data.frame(x=0,y=y,label=labels)
  excluded <- data.frame(
    x=1.48,y=c(4.5,3.5,2.5,1.5),
    label=c(
      sprintf("Excluded: no coded sepsis\nand/or no norepinephrine\nN = %s",scales::comma(flow$n[1]-flow$n[3])),
      sprintf("Excluded: norepinephrine outside\nICU hours 0-24\nN = %s",scales::comma(flow$n[3]-flow$n[4])),
      sprintf("Excluded: incomplete required\npre-index data\nN = %s",scales::comma(flow$n[4]-flow$n[5])),
      sprintf("Excluded: death or ICU discharge\nbefore the 6-hour landmark\nN = %s",scales::comma(flow$n[5]-flow$n[6]))))
  p1 <- ggplot2::ggplot(boxes)+
    ggplot2::geom_segment(data=data.frame(x=0,xend=0,y=y[-length(y)]-.37,yend=y[-1]+.37),
      ggplot2::aes(x=x,xend=xend,y=y,yend=yend),linewidth=.45,color=gray,
      arrow=ggplot2::arrow(length=grid::unit(.09,"in"),type="closed"))+
    ggplot2::geom_segment(data=data.frame(x=.48,xend=.96,y=c(4.5,3.5,2.5,1.5),yend=c(4.5,3.5,2.5,1.5)),
      ggplot2::aes(x=x,xend=xend,y=y,yend=yend),linewidth=.42,color=gray,
      arrow=ggplot2::arrow(length=grid::unit(.08,"in"),type="closed"))+
    ggplot2::geom_label(ggplot2::aes(x=x,y=y,label=label),family="Arial",size=3.2,
      label.size=.35,label.padding=grid::unit(.18,"lines"),fill="white",color=dark)+
    ggplot2::geom_label(data=excluded,ggplot2::aes(x=x,y=y,label=label),family="Arial",size=2.8,
      label.size=.3,label.padding=grid::unit(.15,"lines"),fill="#F4F4F4",color=dark)+
    ggplot2::annotate("segment",x=0,xend=0,y=.64,yend=.43,color=gray,linewidth=.45)+
    ggplot2::annotate("segment",x=-.70,xend=.70,y=.43,yend=.43,color=gray,linewidth=.45)+
    ggplot2::annotate("segment",x=-.70,xend=-.70,y=.43,yend=.22,color=gray,linewidth=.45,
      arrow=ggplot2::arrow(length=grid::unit(.08,"in"),type="closed"))+
    ggplot2::annotate("segment",x=.70,xend=.70,y=.43,yend=.22,color=gray,linewidth=.45,
      arrow=ggplot2::arrow(length=grid::unit(.08,"in"),type="closed"))+
    ggplot2::annotate("label",x=-.70,y=0,label=sprintf("Lower-volume exposure\n<=1,000 mL\nN = %s",scales::comma(lower_n)),family="Arial",size=3.1,label.size=.35,fill="white")+
    ggplot2::annotate("label",x=.70,y=0,label=sprintf("Higher-volume exposure\n>1,000 mL\nN = %s",scales::comma(higher_n)),family="Arial",size=3.1,label.size=.35,fill="white")+
    ggplot2::coord_cartesian(xlim=c(-1.25,2.2),ylim=c(-.45,5.45),clip="off")+
    ggplot2::theme_void(base_family="Arial")+
    ggplot2::labs(title="Figure 1. Source-level cohort reconstruction")+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=11,hjust=0),plot.margin=ggplot2::margin(10,20,10,20))
  save_tiff(p1,"Figure_1_Cohort_Flow_NEJM.tiff",8.5,7)

  # Refit definitive clustered models to obtain delta-method CIs for marginal risks.
  dat[,log1p_pre_ne_fluid_6h_ml:=log1p(pre_ne_fluid_6h_ml)]
  rhs <- paste(c("scale(age)","male","scale(hours_icu_to_ne)","scale(base_map)","scale(base_hr)",
    "scale(log1p(base_lactate))","scale(log1p(base_creatinine))","baseline_invasive_vent",
    "scale(charlson)","scale(pre_ne_fluid_6h_ml)","scale(log1p_pre_ne_fluid_6h_ml)",
    "scale(log1p(initial_ne_rate))","scale(modified_sofa_5c)","modified_sofa_components",
    "pneumonia","urinary_source","abdominal_source","bloodstream_sepsis_code"),collapse="+")
  marginal <- function(d,outcome) {
    des<-survey::svydesign(ids=~subject_id,weights=~sw_trunc,data=d)
    fit<-survey::svyglm(stats::as.formula(paste(outcome,"~restrictive_1000+",rhs)),design=des,family=stats::quasibinomial())
    V<-stats::vcov(fit); b<-stats::coef(fit)
    one <- function(g) {
      nd<-data.table::copy(d);nd[,restrictive_1000:=g]
      X<-stats::model.matrix(stats::delete.response(stats::terms(fit)),data=nd)
      X<-X[,names(b),drop=FALSE]; pr<-stats::plogis(as.numeric(X%*%b));ww<-nd$sw_trunc
      risk<-sum(ww*pr)/sum(ww);grad<-colSums(X*(ww*pr*(1-pr)))/sum(ww)
      se<-sqrt(as.numeric(t(grad)%*%V%*%grad));sel<-se/(risk*(1-risk))
      c(risk=risk,lo=stats::plogis(stats::qlogis(risk)-1.96*sel),hi=stats::plogis(stats::qlogis(risk)+1.96*sel))
    }
    rbind(one(1),one(0))
  }
  m28<-marginal(dat,"death28_landmark")
  ventdat<-dat[baseline_invasive_vent==0];mv<-marginal(ventdat,"incident_vent_72h")
  riskplot<-data.frame(outcome=rep(c("28-day mortality","Incident invasive ventilation"),each=2),
    exposure=rep(c("<=1,000 mL",">1,000 mL"),2),
    risk=c(m28[,"risk"],mv[,"risk"]),lo=c(m28[,"lo"],mv[,"lo"]),hi=c(m28[,"hi"],mv[,"hi"]))
  riskplot$outcome<-factor(riskplot$outcome,levels=c("28-day mortality","Incident invasive ventilation"))
  p2<-ggplot2::ggplot(riskplot,ggplot2::aes(x=exposure,y=risk,color=exposure))+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lo,ymax=hi),width=.08,linewidth=.55)+
    ggplot2::geom_point(size=2.6)+ggplot2::facet_wrap(~outcome,scales="free_y")+
    ggplot2::scale_color_manual(values=c("<=1,000 mL"=blue,">1,000 mL"=orange))+
    ggplot2::scale_y_continuous(labels=scales::percent_format(accuracy=1),expand=ggplot2::expansion(mult=c(.08,.16)))+
    ggplot2::geom_text(ggplot2::aes(label=scales::percent(risk,accuracy=.1)),vjust=-1.25,size=3,color=dark,show.legend=FALSE)+
    ggplot2::labs(title="Figure 2. Adjusted clinical outcome risks",subtitle="Marginal standardized estimates; error bars indicate 95% confidence intervals",x=NULL,y="Adjusted risk")+
    theme_nejm()+ggplot2::theme(legend.position="top",panel.spacing=grid::unit(1,"lines"))
  save_tiff(p2,"Figure_2_Adjusted_Risks_NEJM.tiff",7,4.2)

  # Figure 3: forest plot. No 'favours' annotations because these are exposure associations.
  forest<-data.table::rbindlist(list(
    main[,.(label=analysis,OR,CI_low,CI_high,p_value)],
    first[,.(label=analysis,OR,CI_low,CI_high,p_value)],
    thr[threshold_ml!=1000,.(label=analysis,OR,CI_low,CI_high,p_value)]),fill=TRUE)
  labels_order<-c("Primary: 28-day mortality","Sensitivity: 30-day mortality","Sensitivity: hospital mortality",
    "Incident invasive ventilation","Sensitivity: first eligible stay","Threshold <=500 mL","Threshold <=1500 mL")
  forest<-forest[match(labels_order,label)];forest[,label:=factor(label,levels=rev(labels_order))]
  p3<-ggplot2::ggplot(forest,ggplot2::aes(y=label,x=OR,xmin=CI_low,xmax=CI_high))+
    ggplot2::geom_vline(xintercept=1,color=gray,linewidth=.45,linetype=2)+
    ggplot2::geom_errorbarh(height=.16,linewidth=.55,color=dark)+ggplot2::geom_point(size=2.4,color=blue)+
    ggplot2::scale_x_log10(breaks=c(.5,.75,1,1.5,2),limits=c(.48,2.05))+
    ggplot2::labs(
      title="Figure 3. Adjusted associations with fluid exposure",
      subtitle=paste0(
        "Odds ratios compare the lower-volume with the higher-volume group;\n",
        "error bars indicate 95% confidence intervals"
      ),
      x="Odds ratio (log scale)",y=NULL
    )+theme_nejm()+
    ggplot2::theme(axis.line.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
  save_tiff(p3,"Figure_3_Forest_NEJM.tiff",7.2,5.2)

  # Figure 4: marginal spline curve, restricted to 1st-99th percentiles by source file.
  p4<-ggplot2::ggplot(curve,ggplot2::aes(x=fluid_ml,y=adjusted_risk))+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=CI_low,ymax=CI_high),fill=blue,alpha=.16,color=NA)+
    ggplot2::geom_line(color=blue,linewidth=.8)+
    ggplot2::geom_vline(xintercept=1000,color=gray,linetype=2,linewidth=.45)+
    ggplot2::scale_x_continuous(labels=scales::label_comma(),breaks=c(0,1000,2000,3000,4000,5000))+
    ggplot2::scale_y_continuous(labels=scales::percent_format(accuracy=1),limits=c(.30,.52))+
    ggplot2::labs(title="Figure 4. Adjusted mortality across post-norepinephrine fluid volume",
      subtitle="Natural spline with patient-clustered covariance; shaded area indicates 95% confidence interval; overall P=0.983",
      x="Fluid volume during the first 6 hours after norepinephrine initiation (mL)",y="Adjusted 28-day mortality risk")+
    theme_nejm()+ggplot2::theme(legend.position="none")
  save_tiff(p4,"Figure_4_Spline_NEJM.tiff",7.2,5)

  message("Figures saved in:\n",normalizePath(figure_dir,winslash="/",mustWork=TRUE))
  invisible(figure_dir)
}

run_nejm_figures(RESULTS_DIR,FIGURE_DIR)
