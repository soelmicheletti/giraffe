tigress <-
  function(expdata, tflist=colnames(expdata), targetlist=colnames(expdata), alpha=0.2, nstepsLARS=5, nsplit=100, normalizeexp=TRUE, scoring="area", allsteps=TRUE, verb=FALSE, usemulticore=FALSE)
  {
    # Check if we can run multicore
    if (usemulticore) {
      require(parallel)
    }
    
    # Gene names
    genenames <- colnames(expdata)
    ngenes <- length(genenames)
      
    # Normalize expression data for each gene
    if (normalizeexp)
      expdata <- scale(expdata)
    
    # Make sure there are no more steps than variables
    if (nstepsLARS>length(tflist)-1){
      nstepsLARS<-length(tflist)-1
      if (nstepsLARS==0){cat('Too few transcription factors! \n',stderr())}
      if (verb){
        cat(paste('Variable nstepsLARS was changed to: ',nstepsLARS,'\n')) }}
    
    # Locate TF in gene list by matching their names
    ntf <- length(tflist)
    tfindices <- match(tflist,genenames)
    if (max(is.na(tfindices))) {
      stop('Error: could not find all TF in the gene list!')
    }
    
    # Locate targets in gene list by matching their names
    ntargets <- length(targetlist)
    targetindices <- match(targetlist,genenames)
    if (max(is.na(targetindices))) {
      stop('Error: could not find all targets in the gene list!')
    }
    print("A")
    # Prepare scoring matrix
    if (allsteps) {
      scorestokeep <- nstepsLARS
    } else {
      scorestokeep <- 1
    }
    score <- list()
    
    # A small function to score the regulators of a single gene
    stabselonegene <- function(itarget) {
      print(itarget)
      if (verb) {
        cat('.')
      }
      
      # Name of the target gene
      targetname <- targetlist[itarget]
      # Find the TF to be used for prediction (all TF except the target if the target is itself a TF)
      predTF <- !match(tflist,targetname,nomatch=0)
      r <- stabilityselection(as.matrix(expdata[,tfindices[predTF]]), as.matrix(expdata[,targetname]), nsplit=nsplit, nsteps=nstepsLARS, alpha=alpha)
      sc <- array(0,dim=c(ntf,scorestokeep),dimnames = list(tflist,seq(scorestokeep)))
      if (allsteps) {
        sc[predTF,] <- t(r)
      } else {
        sc[predTF,] <- t(r[nstepsLARS,])
      }
      invisible(sc)
    }
    
    # Treat target genes one by one
    if (usemulticore) {
      if (requireNamespace("foreach") && foreach::getDoParRegistered()) {
        `%dopar%` = foreach::`%dopar%`
        score <- foreach::foreach(i=seq(ntargets), .packages=c("tigress", "lars"),
                                  .export=c("targetlist", "tflist", "tfindices", "expdata",
                                            "verb", "nsplit", "nstepsLARS", "alpha", "ntf",
                                            "scorestokeep", "allsteps")) %dopar% stabselonegene(i)
      } else {
        score <- mclapply(seq(ntargets),stabselonegene,mc.cores=detectCores()-1)
      }
    } else {
      score <- lapply(seq(ntargets),stabselonegene)
    }
    print("B")
    # Rank scores
    edgepred <- list()
    for (i in seq(scorestokeep)) {
      # Combine all scores in a single vectors
      edgepred[[i]] <- matrix(unlist(lapply(score,function(x) x[,i,drop=FALSE])), nrow=ntf)
      rownames(edgepred[[i]]) <- tflist
      colnames(edgepred[[i]]) <- targetlist
    }
    
    # Return the result
    if (allsteps) {
      return(edgepred)
    } else {
      return(edgepred[[1]])
    }
  }

stabilityselection <-
  function(x,y,nsplit=100,nstepsLARS=20,alpha=0.2,scoring="area")
  {
    if (!is.numeric(y) || sd(y)==0) stop("y should be a vector of scalars not constant.")
    n <- nrow(x)
    p <- ncol(x)
    halfsize <- as.integer(n/2)
    freq <- matrix(0,nstepsLARS,p)
    
    i <- 0
    while (i < 2*nsplit) {
      # Randomly reweight each variable
      xs <- t(t(x)*runif(p,alpha,1))
      
      # Ramdomly split the sample in two sets
      badsplit <- TRUE
      while (badsplit) {
        perm <- sample(n)
        i1 <- perm[1:halfsize]
        i2 <- perm[(halfsize+1):n]
        if (max(sd(y[i1]),sd(y[i2]))>0) {badsplit=FALSE}
      }
      
      # run LARS on each randomized, sample and check which variables are selected
      if (sd(y[i1]>0)) {
        r <- lars(xs[i1,],y[i1],max.steps=nstepsLARS,normalize=FALSE,trace=FALSE)
        freq<-freq + abs(sign(r$beta[2:(nstepsLARS+1),]))
        i <- i+1
      }
      if (sd(y[i2]>0)) {
        r <- lars(xs[i2,],y[i2],max.steps=nstepsLARS,normalize=FALSE,trace=FALSE)
        freq<-freq + abs(sign(r$beta[2:(nstepsLARS+1),]))
        i <- i+1
      }
    }
    
    # normalize frequence in [0,1] to get the stability curves
    freq <- freq/i
    
    # Compute normalized area under the stability curve
    if (scoring=="area")
      score <- apply(freq,2,cumsum)/seq(nstepsLARS)
    else
      score <- apply(freq, 2, cummax)
    
    invisible(score)
  }

setwd("/home/soel/giraffe/notebooks/data/prostate/")
expr <- read.csv("expr_t.csv")
rownames(expr)=expr[,1]

# Obtained in Python through
#translate = pd.read_csv("data/breast/raw/gen_v26_mapping.csv")
#tf = []
#for t in motif.columns:
#  translations = list(translate[translate['gene_name'] == t]['gene_id'])
#for translation in translations:
#  if translation[0:15] in motif.index:
#  tf.append(translation[0:15])

tf_comp <- c('ENSG00000169297',
        'ENSG00000165606',
        'ENSG00000150907',
        'ENSG00000165495',
        'ENSG00000136367',
        'ENSG00000164900',
        'ENSG00000143355',
        'ENSG00000012504',
        'ENSG00000137203',
        'ENSG00000166211',
        'ENSG00000107249',
        'ENSG00000166823',
        'ENSG00000138336',
        'ENSG00000010030',
        'ENSG00000108813',
        'ENSG00000179348',
        'ENSG00000148516',
        'ENSG00000215612',
        'ENSG00000095794',
        'ENSG00000074800',
        'ENSG00000167081',
        'ENSG00000116017',
        'ENSG00000162924',
        'ENSG00000136630',
        'ENSG00000131759',
        'ENSG00000123268',
        'ENSG00000119866',
        'ENSG00000179111',
        'ENSG00000073861',
        'ENSG00000151090',
        'ENSG00000187079',
        'ENSG00000176083',
        'ENSG00000188816',
        'ENSG00000119508',
        'ENSG00000159216',
        'ENSG00000169840',
        'ENSG00000133937',
        'ENSG00000162613',
        'ENSG00000204231',
        'ENSG00000103241',
        'ENSG00000170345',
        'ENSG00000106546',
        'ENSG00000157554',
        'ENSG00000196812',
        'ENSG00000082175',
        'ENSG00000203883',
        'ENSG00000185697',
        'ENSG00000184895',
        'ENSG00000100105',
        'ENSG00000198815',
        'ENSG00000163132',
        'ENSG00000152284',
        'ENSG00000089225',
        'ENSG00000185551',
        'ENSG00000106511',
        'ENSG00000157613',
        'ENSG00000173153',
        'ENSG00000095951',
        'ENSG00000101883',
        'ENSG00000143032',
        'ENSG00000185022',
        'ENSG00000172379',
        'ENSG00000139352',
        'ENSG00000116035',
        'ENSG00000105419',
        'ENSG00000101057',
        'ENSG00000166888',
        'ENSG00000172273',
        'ENSG00000115816',
        'ENSG00000178951',
        'ENSG00000101076',
        'ENSG00000006377',
        'ENSG00000213928',
        'ENSG00000007372',
        'ENSG00000163884',
        'ENSG00000066336',
        'ENSG00000185668',
        'ENSG00000120075',
        'ENSG00000128713',
        'ENSG00000174332',
        'ENSG00000071794',
        'ENSG00000126368',
        'ENSG00000177045',
        'ENSG00000170653',
        'ENSG00000177374',
        'ENSG00000179388',
        'ENSG00000102554',
        'ENSG00000119547',
        'ENSG00000064218',
        'ENSG00000166261',
        'ENSG00000129514',
        'ENSG00000146592',
        'ENSG00000168875',
        'ENSG00000148200',
        'ENSG00000151650',
        'ENSG00000156127',
        'ENSG00000185155',
        'ENSG00000101412',
        'ENSG00000130816',
        'ENSG00000168269',
        'ENSG00000134595',
        'ENSG00000213999',
        'ENSG00000109705',
        'ENSG00000136352',
        'ENSG00000110851',
        'ENSG00000004848',
        'ENSG00000071564',
        'ENSG00000001167',
        'ENSG00000069667',
        'ENSG00000245848',
        'ENSG00000162772',
        'ENSG00000116833',
        'ENSG00000197757',
        'ENSG00000171056',
        'ENSG00000166478',
        'ENSG00000099326',
        'ENSG00000075891',
        'ENSG00000204531',
        'ENSG00000197063',
        'ENSG00000065970',
        'ENSG00000186951',
        'ENSG00000084093',
        'ENSG00000162419',
        'ENSG00000118260',
        'ENSG00000125740',
        'ENSG00000105610',
        'ENSG00000075426',
        'ENSG00000135457',
        'ENSG00000105722',
        'ENSG00000198807',
        'ENSG00000106038',
        'ENSG00000260027',
        'ENSG00000134107',
        'ENSG00000167034',
        'ENSG00000160199',
        'ENSG00000100811',
        'ENSG00000179922',
        'ENSG00000163623',
        'ENSG00000037965',
        'ENSG00000129194',
        'ENSG00000251493',
        'ENSG00000152977',
        'ENSG00000143390',
        'ENSG00000176842',
        'ENSG00000214575',
        'ENSG00000186350',
        'ENSG00000147862',
        'ENSG00000184271',
        'ENSG00000072310',
        'ENSG00000068305',
        'ENSG00000187098',
        'ENSG00000102935',
        'ENSG00000164093',
        'ENSG00000168610',
        'ENSG00000179528',
        'ENSG00000189079',
        'ENSG00000162367',
        'ENSG00000068323',
        'ENSG00000189298',
        'ENSG00000197921',
        'ENSG00000171169',
        'ENSG00000163064',
        'ENSG00000205927',
        'ENSG00000175745',
        'ENSG00000117595',
        'ENSG00000174279',
        'ENSG00000120738',
        'ENSG00000169016',
        'ENSG00000105672',
        'ENSG00000119950',
        'ENSG00000163909',
        'ENSG00000180818',
        'ENSG00000140968',
        'ENSG00000188290',
        'ENSG00000077150',
        'ENSG00000184481',
        'ENSG00000112033',
        'ENSG00000124827',
        'ENSG00000173473',
        'ENSG00000198176',
        'ENSG00000136327',
        'ENSG00000185960',
        'ENSG00000185960',
        'ENSG00000008196',
        'ENSG00000173404',
        'ENSG00000112182',
        'ENSG00000164853',
        'ENSG00000159184',
        'ENSG00000140009',
        'ENSG00000139651',
        'ENSG00000111046',
        'ENSG00000111049',
        'ENSG00000185129',
        'ENSG00000126456',
        'ENSG00000167182',
        'ENSG00000162676',
        'ENSG00000196767',
        'ENSG00000186564',
        'ENSG00000196628',
        'ENSG00000162086',
        'ENSG00000188620',
        'ENSG00000261787',
        'ENSG00000096401',
        'ENSG00000140044',
        'ENSG00000169953',
        'ENSG00000107859',
        'ENSG00000197905',
        'ENSG00000106459',
        'ENSG00000124664',
        'ENSG00000164438',
        'ENSG00000212993',
        'ENSG00000125817',
        'ENSG00000138136',
        'ENSG00000160685',
        'ENSG00000143190',
        'ENSG00000106410',
        'ENSG00000169083',
        'ENSG00000157557',
        'ENSG00000244405',
        'ENSG00000120149',
        'ENSG00000136997',
        'ENSG00000168267',
        'ENSG00000134852',
        'ENSG00000179456',
        'ENSG00000009950',
        'ENSG00000089094',
        'ENSG00000164778',
        'ENSG00000152192',
        'ENSG00000175879',
        'ENSG00000135903',
        'ENSG00000130700',
        'ENSG00000168826',
        'ENSG00000102974',
        'ENSG00000173917',
        'ENSG00000100146',
        'ENSG00000147421',
        'ENSG00000107175',
        'ENSG00000198517',
        'ENSG00000168214',
        'ENSG00000111087',
        'ENSG00000130675',
        'ENSG00000172818',
        'ENSG00000164330',
        'ENSG00000175387',
        'ENSG00000165462',
        'ENSG00000125798',
        'ENSG00000131931',
        'ENSG00000116044',
        'ENSG00000168874',
        'ENSG00000185122',
        'ENSG00000128272',
        'ENSG00000171606',
        'ENSG00000163497',
        'ENSG00000132604',
        'ENSG00000064195',
        'ENSG00000187140',
        'ENSG00000123358',
        'ENSG00000185811',
        'ENSG00000160224',
        'ENSG00000073282',
        'ENSG00000006468',
        'ENSG00000109787',
        'ENSG00000085276',
        'ENSG00000175832',
        'ENSG00000114315',
        'ENSG00000025156',
        'ENSG00000105856',
        'ENSG00000129173',
        'ENSG00000010818',
        'ENSG00000122691',
        'ENSG00000052850',
        'ENSG00000135373',
        'ENSG00000113196',
        'ENSG00000020633',
        'ENSG00000165156',
        'ENSG00000139083',
        'ENSG00000164299',
        'ENSG00000174282',
        'ENSG00000008441',
        'ENSG00000135547',
        'ENSG00000185610',
        'ENSG00000107807',
        'ENSG00000185591',
        'ENSG00000081059',
        'ENSG00000115112',
        'ENSG00000107485',
        'ENSG00000153234',
        'ENSG00000100987',
        'ENSG00000006194',
        'ENSG00000118513',
        'ENSG00000169740',
        'ENSG00000109320',
        'ENSG00000170549',
        'ENSG00000161940',
        'ENSG00000124813',
        'ENSG00000111249',
        'ENSG00000130522',
        'ENSG00000176165',
        'ENSG00000182759',
        'ENSG00000106689',
        'ENSG00000150347',
        'ENSG00000143437',
        'ENSG00000148737',
        'ENSG00000135374',
        'ENSG00000160113',
        'ENSG00000112561',
        'ENSG00000012048',
        'ENSG00000122592',
        'ENSG00000167074',
        'ENSG00000175197',
        'ENSG00000184828',
        'ENSG00000130182',
        'ENSG00000105516',
        'ENSG00000128604',
        'ENSG00000101096',
        'ENSG00000101190',
        'ENSG00000090447',
        'ENSG00000176678',
        'ENSG00000129911',
        'ENSG00000176182',
        'ENSG00000126351',
        'ENSG00000119919',
        'ENSG00000215397',
        'ENSG00000149948',
        'ENSG00000092607',
        'ENSG00000111206',
        'ENSG00000172819',
        'ENSG00000164749',
        'ENSG00000172201',
        'ENSG00000063515',
        'ENSG00000049768',
        'ENSG00000043355',
        'ENSG00000131408',
        'ENSG00000135363',
        'ENSG00000158773',
        'ENSG00000110693',
        'ENSG00000139613',
        'ENSG00000169926',
        'ENSG00000128710',
        'ENSG00000109072',
        'ENSG00000189308',
        'ENSG00000188909',
        'ENSG00000139515',
        'ENSG00000115507',
        'ENSG00000101216',
        'ENSG00000120068',
        'ENSG00000198914',
        'ENSG00000170370',
        'ENSG00000106261',
        'ENSG00000126603',
        'ENSG00000136944',
        'ENSG00000121068',
        'ENSG00000008197',
        'ENSG00000126767',
        'ENSG00000162761',
        'ENSG00000153879',
        'ENSG00000125347',
        'ENSG00000125820',
        'ENSG00000069011',
        'ENSG00000165588',
        'ENSG00000091831',
        'ENSG00000072736',
        'ENSG00000185630',
        'ENSG00000064961',
        'ENSG00000172845',
        'ENSG00000204595',
        'ENSG00000162992',
        'ENSG00000171956',
        'ENSG00000166949',
        'ENSG00000108511',
        'ENSG00000133794',
        'ENSG00000175325',
        'ENSG00000164736',
        'ENSG00000105866',
        'ENSG00000112333',
        'ENSG00000112592',
        'ENSG00000105997',
        'ENSG00000114439',
        'ENSG00000169635',
        'ENSG00000165891',
        'ENSG00000100393',
        'ENSG00000087510',
        'ENSG00000143171',
        'ENSG00000124766',
        'ENSG00000143178',
        'ENSG00000143867',
        'ENSG00000204103',
        'ENSG00000105996',
        'ENSG00000151615',
        'ENSG00000082641',
        'ENSG00000116132',
        'ENSG00000106571',
        'ENSG00000169981',
        'ENSG00000135638',
        'ENSG00000078900',
        'ENSG00000111424',
        'ENSG00000123685',
        'ENSG00000164107',
        'ENSG00000119715',
        'ENSG00000160973',
        'ENSG00000115297',
        'ENSG00000115415',
        'ENSG00000102349',
        'ENSG00000170485',
        'ENSG00000119614',
        'ENSG00000136574',
        'ENSG00000204304',
        'ENSG00000197587',
        'ENSG00000168505',
        'ENSG00000116016',
        'ENSG00000198081',
        'ENSG00000142700',
        'ENSG00000143257',
        'ENSG00000135111',
        'ENSG00000118058',
        'ENSG00000126804',
        'ENSG00000115966',
        'ENSG00000113722',
        'ENSG00000177606',
        'ENSG00000078399',
        'ENSG00000177426',
        'ENSG00000185507',
        'ENSG00000175592',
        'ENSG00000109381',
        'ENSG00000181449',
        'ENSG00000144852',
        'ENSG00000170178',
        'ENSG00000134138',
        'ENSG00000171532',
        'ENSG00000164651',
        'ENSG00000171223',
        'ENSG00000132170',
        'ENSG00000180613',
        'ENSG00000111704',
        'ENSG00000163435',
        'ENSG00000109132',
        'ENSG00000142025',
        'ENSG00000113916',
        'ENSG00000197579',
        'ENSG00000170365',
        'ENSG00000136826',
        'ENSG00000108788',
        'ENSG00000163508',
        'ENSG00000134317',
        'ENSG00000132005',
        'ENSG00000204060',
        'ENSG00000077092',
        'ENSG00000107187',
        'ENSG00000125813',
        'ENSG00000127528',
        'ENSG00000165702',
        'ENSG00000120094',
        'ENSG00000163848',
        'ENSG00000173253',
        'ENSG00000081189',
        'ENSG00000102145',
        'ENSG00000183072',
        'ENSG00000177463',
        'ENSG00000186130',
        'ENSG00000105880',
        'ENSG00000005102',
        'ENSG00000180318',
        'ENSG00000116604',
        'ENSG00000108924',
        'ENSG00000184486',
        'ENSG00000137090',
        'ENSG00000114861',
        'ENSG00000164458',
        'ENSG00000120837',
        'ENSG00000170581',
        'ENSG00000126561',
        'ENSG00000181690',
        'ENSG00000070444',
        'ENSG00000177551',
        'ENSG00000151702',
        'ENSG00000151623',
        'ENSG00000135625',
        'ENSG00000123095',
        'ENSG00000138083',
        'ENSG00000128610',
        'ENSG00000115844',
        'ENSG00000257923',
        'ENSG00000102034',
        'ENSG00000134438',
        'ENSG00000171786',
        'ENSG00000043039',
        'ENSG00000120093',
        'ENSG00000009709',
        'ENSG00000141646',
        'ENSG00000122877',
        'ENSG00000137709',
        'ENSG00000113580',
        'ENSG00000134532',
        'ENSG00000214717',
        'ENSG00000214717',
        'ENSG00000178860',
        'ENSG00000128709',
        'ENSG00000122180',
        'ENSG00000118495',
        'ENSG00000185670',
        'ENSG00000184221',
        'ENSG00000106004',
        'ENSG00000180535',
        'ENSG00000164379',
        'ENSG00000173976',
        'ENSG00000128573',
        'ENSG00000137166',
        'ENSG00000057657',
        'ENSG00000172216',
        'ENSG00000067955',
        'ENSG00000184058',
        'ENSG00000168779',
        'ENSG00000112658',
        'ENSG00000164532',
        'ENSG00000137265',
        'ENSG00000080298',
        'ENSG00000104856',
        'ENSG00000143995',
        'ENSG00000168310',
        'ENSG00000176692',
        'ENSG00000138795',
        'ENSG00000164048',
        'ENSG00000129152',
        'ENSG00000123405',
        'ENSG00000131196',
        'ENSG00000054598',
        'ENSG00000173757',
        'ENSG00000174197',
        'ENSG00000136535',
        'ENSG00000136931',
        'ENSG00000105967',
        'ENSG00000138378',
        'ENSG00000156273',
        'ENSG00000162599',
        'ENSG00000123388',
        'ENSG00000114853',
        'ENSG00000066136',
        'ENSG00000120690',
        'ENSG00000167967',
        'ENSG00000205922',
        'ENSG00000137270',
        'ENSG00000196092',
        'ENSG00000140262',
        'ENSG00000118217',
        'ENSG00000074047',
        'ENSG00000100219',
        'ENSG00000118707',
        'ENSG00000102908',
        'ENSG00000140836',
        'ENSG00000123576',
        'ENSG00000105698',
        'ENSG00000111783',
        'ENSG00000113658',
        'ENSG00000177485',
        'ENSG00000100644',
        'ENSG00000005889',
        'ENSG00000178403',
        'ENSG00000174576',
        'ENSG00000125285',
        'ENSG00000007866',
        'ENSG00000106536',
        'ENSG00000123364',
        'ENSG00000133740',
        'ENSG00000087903',
        'ENSG00000175329',
        'ENSG00000144355',
        'ENSG00000158711',
        'ENSG00000167840',
        'ENSG00000134046',
        'ENSG00000165556',
        'ENSG00000128714',
        'ENSG00000160961',
        'ENSG00000177732',
        'ENSG00000125398',
        'ENSG00000148704',
        'ENSG00000178665',
        'ENSG00000079432',
        'ENSG00000156925',
        'ENSG00000129535',
        'ENSG00000028277',
        'ENSG00000137273',
        'ENSG00000164916',
        'ENSG00000198911',
        'ENSG00000170561',
        'ENSG00000153779',
        'ENSG00000182158',
        'ENSG00000131668',
        'ENSG00000173039',
        'ENSG00000105991',
        'ENSG00000117036',
        'ENSG00000117707',
        'ENSG00000121075',
        'ENSG00000221869',
        'ENSG00000154727',
        'ENSG00000143365',
        'ENSG00000039600',
        'ENSG00000169136',
        'ENSG00000165804',
        'ENSG00000180532',
        'ENSG00000103495',
        'ENSG00000124092',
        'ENSG00000106331',
        'ENSG00000102878',
        'ENSG00000162702',
        'ENSG00000125618',
        'ENSG00000124782',
        'ENSG00000092067',
        'ENSG00000205250',
        'ENSG00000106852',
        'ENSG00000134954',
        'ENSG00000120798',
        'ENSG00000125952',
        'ENSG00000170265',
        'ENSG00000169856',
        'ENSG00000171634',
        'ENSG00000007968',
        'ENSG00000170608',
        'ENSG00000164683',
        'ENSG00000126746',
        'ENSG00000005513',
        'ENSG00000188786',
        'ENSG00000112242',
        'ENSG00000106031',
        'ENSG00000111145',
        'ENSG00000156150',
        'ENSG00000141510',
        'ENSG00000100968',
        'ENSG00000118689',
        'ENSG00000152804',
        'ENSG00000064835',
        'ENSG00000215271',
        'ENSG00000184937',
        'ENSG00000135100',
        'ENSG00000123407',
        'ENSG00000256683',
        'ENSG00000174963',
        'ENSG00000163666',
        'ENSG00000180828',
        'ENSG00000196482',
        'ENSG00000065978')
library(lars)
g <- 5000
tf_r <- tf_comp[which(tf_comp %in% colnames(expr)[seq(2, g, 1)])]
R <- tigress(expr[,c(seq(2, g, 1))], tflist = tf_r, nstepsLARS=2, nsplit=10, allsteps=FALSE)
write.csv(R, "R_tigress.csv")