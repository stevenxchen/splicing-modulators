#Junction counter adapted from rMATS "process unique sam"

### import necessary libraries
import re,os,sys,logging,time,datetime,commands,argparse;
import scipy,math,pickle,pysam;
from scipy import stats;

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;

######### basic functions #############
#
#Input:
#insertsize_inclusion: a numerical variable for the paired-end read inclusion form insert size
#insertsize_skipping: a numerical variable for the paired-end read skipping form insert size
#insertsize_mean: a numerical variable for the mean of the insert size
#insertsize_var: a numerical variable for the variance of the insert size
#Output:
#A 2-element vector for the [fraction of the inclusion form, fraction of the skipping form]
def PE_fraction(insertsize_inclusion, insertsize_skipping, insertsize_mean, insertsize_var):
  p1=stats.norm.cdf(insertsize_inclusion+0.5,insertsize_mean,math.sqrt(insertsize_var));
  p1=p1-stats.norm.cdf(insertsize_inclusion-0.5,insertsize_mean,math.sqrt(insertsize_var));
  p2=stats.norm.cdf(insertsize_skipping+0.5,insertsize_mean,math.sqrt(insertsize_var));
  p2=p2-stats.norm.cdf(insertsize_skipping-0.5,insertsize_mean,math.sqrt(insertsize_var));
  return([p1/(p1+p2),p2/(p1+p2)]);

def getInitialCounts(): ## getting initial counts for each AS event
  #rValue = [[[[],[]]],[[[],[]]],[[[],[]]]]; ## count type 1, 2, and 3, Include or Skip
  rValue = [[[],[]],[[],[]],[[],[]]];
  rValue[CT1][I]=0;
  rValue[CT1][S]=0;
  rValue[CT2][I]=0;
  rValue[CT2][S]=0;
  rValue[CT3][I]=0;
  rValue[CT3][S]=0;
  return rValue;
#### end of getInitialCounts()

###############################
## SAMPLE COUNTING FUNCTIONS ##
###############################
def is_unique(r):
  # check if read is uniquely mapped
  for tag in r.tags:
    if tag[0]=='NH': # NH is tag for number of hits
      if int(tag[1])==1: # uniquely mapped if NH=1
        if dataType=='single': ## single end, sufficient
          return True;
        elif r.is_proper_pair:
          return True; 
  return False

#
def processSample(sample): ## call it with processSample(sample_1, S1) something like this
    global anchorLength;  
    readCount=0;
  ### process the given sample ###
  #for s1 in sample: ## for each sam file 
  #  rep = sample.index(s1);
  #  if len(s1.strip())<1: ## incorrect split. a user might accidently put comma at the end of input sam file list
  #    continue; ### just skip this entry, probably the last one though
    s1=sample;
    sFile = pysam.Samfile(s1.strip(),'rb'); ## open bam file
    e1 = {}; ## edge count here
    for read in sFile.fetch():
      readCount=readCount+1;
      if (readCount % 1000000) == 0:
        print(str(readCount)+" reads processed...");

      if not is_unique(read):
        continue;

      rL = readLength;  
      chr = sFile.getrname(read.tid)
      if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
        chr = 'chr'+chr;
      mc = read.pos+1
      mString = read.cigarstring
      #rID = read.qname

      group = mc/chunk; ## group does not change, it's okay to check only one group for a junction read
##      if showMe: print 'mc:'+ str(group);
      if 'D' in mString or 'I' in mString or 'S' in mString or 'H' in mString or 'P' in mString or 'X' in mString or '=' in mString: ## skip
        continue; ## go to next line

      ### check to see if the line is either exonic read or junction read
      split_mString = mString.split('M');
      tor = 0; ## type of read, 0 nothing, 1 exonic read, 2 junction read
      if len(split_mString)==2:
        tor = 1; ############ exonic read ######
        rL = int(split_mString[0]); ## read length specified in this mapping string
        mec = mc+rL-1; ## mapping end coord 
        egroup = mec/chunk; ## group for the mapping end coord of genomic reads

        if rL != readLength:
            continue

        ## SE ###
        if chr in se: ## this chromosome has se event(s)
          if group in se[chr]: ## this group has skipped exon event(s)
            for c in se[chr][group]: ## for each skipped exon event in this group
              if (mc>se[chr][group][c][0] and mec<=se[chr][group][c][1]): ## read on the target
                c_se[c][CT2][I]+=1; 
        ### end of SE ###

        ### MXE ####
        if chr in mxe: ## this chromosome has mxe event(s)
          if group in mxe[chr]: ## this group has mxe event(s)
            for c in mxe[chr][group]: ## for each mxe event in this group
              if (mc>mxe[chr][group][c][0] and mec<=mxe[chr][group][c][1]): ## read on the target exon
                c_mxe[c][CT2][I]+=1;
              elif (mc>mxe[chr][group][c][2] and mec<=mxe[chr][group][c][3]): ## read on the second exon
                c_mxe[c][CT2][S]+=1;
        ## end of MXE ###

        ## A5SS ##
        if chr in a5ss: ## this chromosome has a5ss event(s)
          if egroup in a5ss[chr]: ## this egroup has a5ss event(s), this group is valid for positive strand only
            for c in a5ss[chr][egroup]: ## for each a5ss event in this group

              if a5ss[chr][egroup][c][4]>a5ss[chr][egroup][c][1]: ## positive strand
                if (mc<=(a5ss[chr][egroup][c][3]-(rL-junctionLength/2)+1) and mec<=a5ss[chr][egroup][c][1] and mec>=(a5ss[chr][egroup][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a5ss[c][CT1][I]+=1;
                  c_a5ss[c][CT2][I]+=1;
                if (mc>a5ss[chr][egroup][c][3] and mec<=a5ss[chr][egroup][c][1]): ## exon read supporting target
                  c_a5ss[c][CT2][I]+=1;

          if group in a5ss[chr]: ## this group has a5ss event(s), this group is valid for negative strand only
            for c in a5ss[chr][group]: ## for each a5ss event in this group
                
              if a5ss[chr][group][c][4]<a5ss[chr][group][c][1]: ## negative strand
                if (mc>a5ss[chr][group][c][0] and mc<=(a5ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec>=(a5ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a5ss[c][CT1][I]+=1;
                  c_a5ss[c][CT2][I]+=1;
                if (mc>a5ss[chr][group][c][0] and mec<=a5ss[chr][group][c][2]): ## exon read supporting target
                  c_a5ss[c][CT2][I]+=1;

        ## end of A5SS ###

        ## A3SS ##
        if chr in a3ss: ## this chromosome has a3ss event(s)
          if egroup in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][egroup]: ## for each a3ss event in this group

              if a3ss[chr][egroup][c][4]>a3ss[chr][egroup][c][1]: ## negative strand
                if (mc<=(a3ss[chr][egroup][c][3]-(rL-junctionLength/2)+1) and mec<=a3ss[chr][egroup][c][1] and mec>=(a3ss[chr][egroup][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a3ss[c][CT1][I]+=1;
                  c_a3ss[c][CT2][I]+=1;
                if (mc>a3ss[chr][egroup][c][3] and mec<=a3ss[chr][egroup][c][1]): ## exon read supporting target
                  c_a3ss[c][CT2][I]+=1;

          if group in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][group]: ## for each a3ss event in this group

              if a3ss[chr][group][c][4]<a3ss[chr][group][c][1]: ## positive strand
                if (mc>a3ss[chr][group][c][0] and mc<=(a3ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec>=(a3ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a3ss[c][CT1][I]+=1;
                  c_a3ss[c][CT2][I]+=1;
                if (mc>a3ss[chr][group][c][0] and mec<=a3ss[chr][group][c][2]): ## exon read supporting target
                  c_a3ss[c][CT2][I]+=1;

        ## end of A3SS ###

        ## RI ##
        if chr in ri: ## this chromosome has ri event(s)
          groupsToExamine = list(set([group,egroup])); ## to examine multiple groups 
          tempProcessedRI_id={}; ## to examine multiple groups
          for ggg in groupsToExamine: ## for each group
            if ggg in ri[chr]: ## this group has ri event(s)
              for c in ri[chr][ggg]: ## for each ri event in this group, strand does not matter for ri events
                if c in tempProcessedRI_id: ## already processed this event, skip it 
                  continue; ## next c please
                else: ## new c here
                  tempProcessedRI_id[c]='1';
                  if (mc<=(ri[chr][ggg][c][3]-(rL-junctionLength/2)+1) and mec>=(ri[chr][ggg][c][3]+(rL-junctionLength/2))) or (mc<=(ri[chr][ggg][c][4]-(rL-junctionLength/2)+1) and mec>=(ri[chr][ggg][c][4]+(rL-junctionLength/2))): ## multi-exon read supporting target
                    c_ri[c][CT1][I]+=1;
                    c_ri[c][CT2][I]+=1;
                  if (mc>ri[chr][ggg][c][3] and mec<=ri[chr][ggg][c][4]): ## exon read supporting target
                    c_ri[c][CT2][I]+=1;
        ## end of RI ###

    ####### now take care of junction counts from junction reads ######    

      elif len(split_mString)>=3: ###### junction read ###########
        beginning=int(split_mString[0].split('N')[-1]);
        ending = int(split_mString[-2].split('N')[-1]);
        jS=mc; jE=mc-1;
        usedE={};## to avoid using the same read more than one time for the same event
        #usedG={}; ## used group, to avoid examining the same group more than once
        prevN=-1; nextN=-1; ## for previous junction and next junction
        for ec in range(0,len(split_mString)-2): ## for each coordinate
          secondNumber = int(split_mString[ec].split('N')[-1]);
          jumpNumber = int(split_mString[ec+1].split('N')[0]);
          lastNumber = int(split_mString[ec+1].split('N')[-1]);
          if ec>0: ## there is prevN
            prevN = int(split_mString[ec].split('N')[0]);
          if (ec>=0 and ec<len(split_mString)-3): ## there is nextN
            nextN = int(split_mString[ec+2].split('N')[0]);
          if (ec>=0 and ec<len(split_mString)-2): ## there iS nextN
            if (ec==len(split_mString)-3):
              nextN=-1;
            else:
              nextN = int(split_mString[ec+2].split('N')[0]);
          jS = jE+secondNumber; ## 1-base
          jE = jS+jumpNumber; ## 0-base
          key = chr+'_'+str(jS)+'_'+str(jE);
          cStart=jS-secondNumber+1; cEnd=jS;
          cStart_last=cStart;
          cEnd_last=cEnd;

          minAnchor = min(int(split_mString[0]), int(split_mString[-2].split('N')[1])); ## min nts going across junction
          if minAnchor<anchorLength: ## not a valid junction read, do not count it
            continue;        

          if ec==0: ## first junction, check the first number
            if beginning<anchorLength: ## not a valid junction read, do not count it
              continue; ## next junction
          if ec==(len(split_mString)-3): ## check the last segment length
            if ending<anchorLength: ## not a valid junction read, do not count it
              continue;

          ## edge counts
          if key in e1: ## exist!
            e1[key] = e1[key]+1;
          else: ## new junction
            e1[key] = 1;

          group = jS/chunk;  ## changing group here
          #if group in usedG: ## already visited this group
          #  continue; ## next junction coord
          #else: ## update usedG dictioanry
          #  usedG[group]='1';

          ## SE ###
          if chr in se: ## this chromosome has se event(s)
            if group in se[chr]: ## this group has skipped exon event(s)
              for c in se[chr][group]: ## for each skipped exon event in this group, examine if the given junction is part of it
                seLength=se[chr][group][c][1]-se[chr][group][c][0];
                upjLength=se[chr][group][c][0]-se[chr][group][c][3];
                dnjLength=se[chr][group][c][4]-se[chr][group][c][1];
                if (jS==se[chr][group][c][3] and jE==se[chr][group][c][0] and ((lastNumber<=seLength and nextN==-1) or (lastNumber==seLength and nextN==dnjLength)))  or (jS==se[chr][group][c][1] and jE==se[chr][group][c][4] and ((secondNumber<=seLength and prevN==-1) or (secondNumber==seLength and prevN==upjLength))): ## IJC
                  key=':'.join(["SE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_se[c][CT1][I]+=1; 
                  c_se[c][CT2][I]+=1; 
                  usedE[key]=1; ## okay to overwrite     
                elif jS==se[chr][group][c][3] and jE==se[chr][group][c][4]: ## SJC
                  key=':'.join(["SE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_se[c][CT1][S]+=1; 
                  c_se[c][CT2][S]+=1; 
                  usedE[key]=1; ## okay to overwrite     
          ### end of SE ###  

          ## MXE ###
          if chr in mxe: ## this chromosome has mxe event(s)
            if group in mxe[chr]: ## this group has mxe event(s)
              for c in mxe[chr][group]: ## for each mxe event in this group, examine if the given junction is part of it
                mxeLen1=mxe[chr][group][c][1]-mxe[chr][group][c][0];
                mxeLen2=mxe[chr][group][c][3]-mxe[chr][group][c][2];
                upj1Len=mxe[chr][group][c][0]-mxe[chr][group][c][5];
                dnj1Len=mxe[chr][group][c][6]-mxe[chr][group][c][1];
                upj2Len=mxe[chr][group][c][2]-mxe[chr][group][c][5];
                dnj2Len=mxe[chr][group][c][6]-mxe[chr][group][c][3];


                if (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][0] and ((lastNumber<=mxeLen1 and nextN==-1) or (lastNumber==mxeLen1 and nextN==dnj1Len))) or (jS==mxe[chr][group][c][1] and jE==mxe[chr][group][c][6] and ((secondNumber<=mxeLen1 and prevN==-1) or (secondNumber==mxeLen1 and prevN==upj1Len))): ## IJC
                  key=':'.join(["MXE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_mxe[c][CT1][I]+=1;
                  c_mxe[c][CT2][I]+=1;
                  usedE[key]=1; ## okay to overwrite     
                elif (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][2] and ((lastNumber<=mxeLen2 and nextN==-1) or (lastNumber==mxeLen2 and nextN==dnj2Len))) or (jS==mxe[chr][group][c][3] and jE==mxe[chr][group][c][6] and ((secondNumber<=mxeLen2 and prevN==-1) or (secondNumber==mxeLen2 and prevN==upj2Len))): ## SJC
                  key=':'.join(["MXE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_mxe[c][CT1][S]+=1;
                  c_mxe[c][CT2][S]+=1;
                  usedE[key]=1; ## okay to overwrite     
          ### end of MXE ###  

          ## A5SS ###
          if chr in a5ss: ## this chromosome has a5ss event(s)
            if group in a5ss[chr]: ## this group has a5ss event(s)
              for c in a5ss[chr][group]: ## for each a5ss event in this group, examine if the given junction is part of it
                if a5ss[chr][group][c][4]>a5ss[chr][group][c][1]: ## positive strand
                  if ( jS==a5ss[chr][group][c][1] and jE==a5ss[chr][group][c][4] ): ## IJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][I]+=1;
                    c_a5ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a5ss[chr][group][c][3] and jE==a5ss[chr][group][c][4] ): ## SJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][S]+=1;
                    c_a5ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif nextN == -1:
                    cStart_last = jE+1; cEnd_last=jE+lastNumber;
                    if (cStart_last<=(a5ss[chr][group][c][3]-(rL-junctionLength/2)+1) and cEnd_last<=a5ss[chr][group][c][1] and cEnd_last>=(a5ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A5SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a5ss[c][CT1][I]+=1;
                        c_a5ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
                else: ## negative strand
                  if ( jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][0]  ): ## IJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][I]+=1;
                    c_a5ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][2] ): ## SJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][S]+=1;
                    c_a5ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif prevN == -1:
                    if (cStart>a5ss[chr][group][c][0] and cStart<=(a5ss[chr][group][c][2]-(rL-junctionLength/2)+1) and cEnd>=(a5ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A5SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a5ss[c][CT1][I]+=1;
                        c_a5ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
          ### end of A5SS ###

          ## A3SS ###
          if chr in a3ss: ## this chromosome has a3ss event(s)
            if group in a3ss[chr]: ## this group has a3ss event(s)
              for c in a3ss[chr][group]: ## for each a3ss event in this group, examine if the given junction is part of it
                if a3ss[chr][group][c][4]>a3ss[chr][group][c][1]: ## negative strand
                  if ( jS==a3ss[chr][group][c][1] and jE==a3ss[chr][group][c][4]  ): ## IJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][I]+=1;
                    c_a3ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a3ss[chr][group][c][3] and jE==a3ss[chr][group][c][4] ): ## SJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][S]+=1;
                    c_a3ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif nextN == -1:
                    cStart_last = jE+1; cEnd_last=jE+lastNumber;
                    if (cStart_last<=(a3ss[chr][group][c][3]-(rL-junctionLength/2)+1) and cEnd_last<=a3ss[chr][group][c][1] and cEnd_last>=(a3ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A3SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a3ss[c][CT1][I]+=1;
                        c_a3ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
                else: ## positive strand
                  if ( jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][0] ): ## IJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][I]+=1;
                    c_a3ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][2] ): ## SJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][S]+=1;
                    c_a3ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif prevN == -1:
                    if (cStart>a3ss[chr][group][c][0] and cStart<=(a3ss[chr][group][c][2]-(rL-junctionLength/2)+1) and cEnd>=(a3ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A3SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a3ss[c][CT1][I]+=1;
                        c_a3ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
          ### end of A3SS ###

          ## RI ###
          if chr in ri: ## this chromosome has ri event(s)
            if group in ri[chr]: ## this group has ri event(s)
              for c in ri[chr][group]: ## for each ri event in this group, examine if the given junction is part of it
                if jS==ri[chr][group][c][3] and jE==ri[chr][group][c][4] : ## SJC
                  key=':'.join(["RI",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_ri[c][CT1][S]+=1;
                  c_ri[c][CT2][S]+=1;
                  usedE[key]=1; ## okay to overwrite     
                if prevN == -1:
                  if cStart<=(ri[chr][group][c][4]-(rL-junctionLength/2)+1) and cEnd>=(ri[chr][group][c][4]+(rL-junctionLength/2)): ## multi-exon read supporting target
                    key=':'.join(["RI",str(c)]);
                    if key not in usedE: ## read didn't cover this event before, add it
                      c_ri[c][CT1][I]+=1;
                      c_ri[c][CT2][I]+=1;
                      usedE[key]=1; ## okay to overwrite
                if nextN == -1:
                  cStart_last = jE+1; cEnd_last=jE+lastNumber; ## just change cStart and cEnd then do the same process
                  if cStart_last<=(ri[chr][group][c][3]-(rL-junctionLength/2)+1) and cEnd_last>=(ri[chr][group][c][3]+(rL-junctionLength/2)):
                    key=':'.join(["RI",str(c)]);
                    if key not in usedE: ## read didn't cover this event before, add it
                      c_ri[c][CT1][I]+=1;
                      c_ri[c][CT2][I]+=1;
                      usedE[key]=1; ## okay to overwrite
          ### end of RI ###

          #### end of for ec #####

      else: ## it is not exonic nor junction read. proceed to the next line
        continue;

    sFile.close();  
    #edgeFile = open(s1.strip()+'.edgeCount', 'w');
    #for k in e1:
    #  edgeFile.write(k+'\t'+str(e1[k])+'\n');
    #edgeFile.close();
##### end of processSample ######


def processSample_stranded(sample,dt,lt): ## call it with processSample_firststrand(sample_1, S1,dataType,libraryType) something like this
    global anchorLength;  

  ### process the given sample ###
  #for s1 in sample: ## for each sam file 
  #  rep = sample.index(s1);
  #  if len(s1.strip())<1: ## incorrect split. a user might accidently put comma at the end of input sam file list
  #    continue; ### just skip this entry, probably the last one though
    s1=sample;
    sFile = pysam.Samfile(s1.strip(),'rb'); ## open bam file
    e1 = {}; ## edge count here
    for read in sFile.fetch():

      if not is_unique(read):
        continue;

      rL = readLength;
      chr = sFile.getrname(read.tid)
      if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
        chr = 'chr'+chr;
      mc = read.pos+1
      mString = read.cigarstring
      #rID = read.qname
      group = mc/chunk; ## group does not change, it's okay to check only one group for a junction read
      if 'D' in mString or 'I' in mString or 'S' in mString or 'H' in mString or 'P' in mString or 'X' in mString or '=' in mString: ## skip
        continue; ## go to next line


      ### check to see if the line is either exonic read or junction read
      split_mString = mString.split('M');
      tor = 0; ## type of read, 0 nothing, 1 exonic read, 2 junction read
      if len(split_mString)==2:
        tor = 1; ############ exonic read ######
        rL = int(split_mString[0]); ## read length specified in this mapping string
        mec = mc+rL-1; ## mapping end coord 
        egroup = mec/chunk;
        myStrand='';
        sFlag=read.flag;
        if dt=='paired' and lt=='first': ## paired-read, firststrand
          if (sFlag&16 and sFlag&64) or (sFlag&32 and sFlag&128): ##'+' strand
            myStrand='+';
          elif (sFlag&32 and sFlag&64) or (sFlag&16 and sFlag&128): ##'-' strand
            myStrand='-';
        elif dt=='paired' and lt=='second': ## paired-read, secondstrand
          if (sFlag&32 and sFlag&64) or (sFlag&16 and sFlag&128): ##'+' strand
            myStrand='+';
          elif (sFlag&16 and sFlag&64) or (sFlag&32 and sFlag&128): ##'-' strand
            myStrand='-';
        elif dt=='single' and lt=='first': ##single-end, firststrand
          if sFlag&16 : ##'+' strand
            myStrand='+';
          elif sFlag&32: ##'-' strand
            myStrand='-';
        elif dt=='single' and lt=='second': ##single-end, secondstrand
          if sFlag&32 : ##'+' strand
            myStrand='+';
          elif sFlag&16: ##'-' strand
            myStrand='-';

        if len(myStrand)==0 or rL != readLength: ## skip this read, cannot determine the strand
          continue;

        ## SE ###
        if chr in se: ## this chromosome has se event(s)
          if group in se[chr]: ## this group has skipped exon event(s)
            for c in se[chr][group]: ## for each skipped exon event in this group
              if se[chr][group][c][-1]==myStrand: ## same strand
                if (mc>se[chr][group][c][0] and mec<=se[chr][group][c][1]): ## read on the target
                  c_se[c][CT2][I]+=1; 
        ### end of SE ###

        ### MXE ####
        if chr in mxe: ## this chromosome has mxe event(s)
          if group in mxe[chr]: ## this group has mxe event(s)
            for c in mxe[chr][group]: ## for each mxe event in this group
              if mxe[chr][group][c][-1]==myStrand: ## same strand
                if (mc>mxe[chr][group][c][0] and mec<=mxe[chr][group][c][1]): ## read on the target exon
                  c_mxe[c][CT2][I]+=1;
                elif (mc>mxe[chr][group][c][2] and mec<=mxe[chr][group][c][3]): ## read on the second exon
                  c_mxe[c][CT2][S]+=1;
        ## end of MXE ###

        ## A5SS ##
        if chr in a5ss: ## this chromosome has a5ss event(s)
          if egroup in a5ss[chr]: ## this egroup has a5ss event(s), valid only for positive strand only
            for c in a5ss[chr][egroup]: ## for each a5ss event in this group
              if a5ss[chr][egroup][c][-1]==myStrand: ## same strand
                if a5ss[chr][egroup][c][4]>a5ss[chr][egroup][c][1]: ## positive strand
                  if (mc<=(a5ss[chr][egroup][c][3]-(rL-junctionLength/2)+1) and mec<=a5ss[chr][egroup][c][1] and mec>=(a5ss[chr][egroup][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                    c_a5ss[c][CT1][I]+=1;
                    c_a5ss[c][CT2][I]+=1;
                  if (mc>a5ss[chr][egroup][c][3] and mec<=a5ss[chr][egroup][c][1]): ## exon read supporting target
                    c_a5ss[c][CT2][I]+=1;

          if group in a5ss[chr]: ## this group has a5ss event(s), for negative strand ones only
            for c in a5ss[chr][group]: ## for each a5ss event in this group
              if a5ss[chr][group][c][-1]==myStrand: ## same strand
                  
                if a5ss[chr][group][c][4]<a5ss[chr][group][c][1]: ## negative strand
                  if (mc>a5ss[chr][group][c][0] and mc<=(a5ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec>=(a5ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                    c_a5ss[c][CT1][I]+=1;
                    c_a5ss[c][CT2][I]+=1;
                  if (mc>a5ss[chr][group][c][0] and mec<=a5ss[chr][group][c][2]): ## exon read supporting target
                    c_a5ss[c][CT2][I]+=1;
        ## end of A5SS ###
       
        ## A3SS ##
        if chr in a3ss: ## this chromosome has a3ss event(s)
          if egroup in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][egroup]: ## for each a3ss event in this group
              if a3ss[chr][egroup][c][-1]==myStrand: ## same strand
                if a3ss[chr][egroup][c][4]>a3ss[chr][egroup][c][1]: ## negative strand
                  if (mc<=(a3ss[chr][egroup][c][3]-(rL-junctionLength/2)+1) and mec<=a3ss[chr][egroup][c][1] and mec>=(a3ss[chr][egroup][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                    c_a3ss[c][CT1][I]+=1;
                    c_a3ss[c][CT2][I]+=1;
                  if (mc>a3ss[chr][egroup][c][3] and mec<=a3ss[chr][egroup][c][1]): ## exon read supporting target
                    c_a3ss[c][CT2][I]+=1;
 
          if group in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][group]: ## for each a3ss event in this group
              if a3ss[chr][group][c][-1]==myStrand: ## same strand
                if a3ss[chr][group][c][4]<a3ss[chr][group][c][1]: ## positive strand
                  if (mc>a3ss[chr][group][c][0] and mc<=(a3ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec>=(a3ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                    c_a3ss[c][CT1][I]+=1;
                    c_a3ss[c][CT2][I]+=1;
                  if (mc>a3ss[chr][group][c][0] and mec<=a3ss[chr][group][c][2]): ## exon read supporting target
                    c_a3ss[c][CT2][I]+=1;
        ## end of A3SS ###

        ## RI ##
        if chr in ri: ## this chromosome has ri event(s)
          groupsToExamine = list(set([group,egroup])); ## to examine multiple groups 
          tempProcessedRI_id={}; ## to examine multiple groups
          for ggg in groupsToExamine: ## for each group
            if ggg in ri[chr]: ## this group has ri event(s)
              for c in ri[chr][ggg]: ## for each ri event in this group, strand does not matter for ri events

                if c in tempProcessedRI_id: ## already processed this event, skip it 
                  continue; ## next c please
                else: ## new c here
                  if ri[chr][ggg][c][-1]==myStrand: ## same strand
                    if (mc<=(ri[chr][ggg][c][3]-(rL-junctionLength/2)+1) and mec>=(ri[chr][ggg][c][3]+(rL-junctionLength/2))) or (mc<=(ri[chr][ggg][c][4]-(rL-junctionLength/2)+1) and mec>=(ri[chr][ggg][c][4]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      c_ri[c][CT1][I]+=1;
                      c_ri[c][CT2][I]+=1;
                    if (mc>ri[chr][ggg][c][3] and mec<=ri[chr][ggg][c][4]): ## exon read supporting target
                      c_ri[c][CT2][I]+=1;
        ## end of RI ###

    ####### now take care of junction counts from junction reads ######    

      elif len(split_mString)>=3: ###### junction read ###########
        beginning=int(split_mString[0].split('N')[-1]);
        ending = int(split_mString[-2].split('N')[-1]);
        jS=mc; jE=mc-1;
        usedE={};## to avoid using the same read more than one time for the same event
        #usedG={}; ## used group, to avoid examining the same group more than once
        prevN=-1; nextN=-1; ## for previous junction and next junction
        for ec in range(0,len(split_mString)-2): ## for each coordinate
          secondNumber = int(split_mString[ec].split('N')[-1]);
          jumpNumber = int(split_mString[ec+1].split('N')[0]);
          lastNumber = int(split_mString[ec+1].split('N')[-1]);
          if ec>0: ## there is prevN
            prevN = int(split_mString[ec].split('N')[0]);
          if (ec>=0 and ec<len(split_mString)-3): ## there is nextN
            nextN = int(split_mString[ec+2].split('N')[0]);
          if (ec>=0 and ec<len(split_mString)-2): ## there iS nextN
            if (ec==len(split_mString)-3):
              nextN=-1;
            else:
              nextN = int(split_mString[ec+2].split('N')[0]);
          jS = jE+secondNumber; ## 1-base
          jE = jS+jumpNumber; ## 0-base
          key = chr+'_'+str(jS)+'_'+str(jE);
          cStart=jS-secondNumber+1; cEnd=jS;

          minAnchor = min(int(split_mString[0]), int(split_mString[-2].split('N')[1])); ## min nts going across junction
          if minAnchor<anchorLength: ## not a valid junction read, do not count it
            continue;        

          if ec==0: ## first junction, check the first number
            if beginning<anchorLength: ## not a valid junction read, do not count it
              continue; ## next junction
          if ec==(len(split_mString)-3): ## check the last segment length
            if ending<anchorLength: ## not a valid junction read, do not count it
              continue;

          ## edge counts
          if key in e1: ## exist!
            e1[key] = e1[key]+1;
          else: ## new junction
            e1[key] = 1;

          group = jS/chunk;  ## changing group here
          #if group in usedG: ## already visited this group
          #  continue; ## next junction coord
          #else: ## update usedG dictioanry
          #  usedG[group]='1';

          ## SE ###
          if chr in se: ## this chromosome has se event(s)
            if group in se[chr]: ## this group has skipped exon event(s)
              for c in se[chr][group]: ## for each skipped exon event in this group, examine if the given junction is part of it
                seLength=se[chr][group][c][1]-se[chr][group][c][0];
                upjLength=se[chr][group][c][0]-se[chr][group][c][3];
                dnjLength=se[chr][group][c][4]-se[chr][group][c][1];
                if (jS==se[chr][group][c][3] and jE==se[chr][group][c][0] and ((lastNumber<=seLength and nextN==-1) or (lastNumber==seLength and nextN==dnjLength)))  or (jS==se[chr][group][c][1] and jE==se[chr][group][c][4] and ((secondNumber<=seLength and prevN==-1) or (secondNumber==seLength and prevN==upjLength))): ## IJC
                  key=':'.join(["SE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_se[c][CT1][I]+=1; 
                  c_se[c][CT2][I]+=1; 
                  usedE[key]=1; ## okay to overwrite     
                elif jS==se[chr][group][c][3] and jE==se[chr][group][c][4]: ## SJC
                  key=':'.join(["SE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_se[c][CT1][S]+=1; 
                  c_se[c][CT2][S]+=1; 
                  usedE[key]=1; ## okay to overwrite     
          ### end of SE ###  

          ## MXE ###
          if chr in mxe: ## this chromosome has mxe event(s)
            if group in mxe[chr]: ## this group has mxe event(s)
              for c in mxe[chr][group]: ## for each mxe event in this group, examine if the given junction is part of it
                mxeLen1=mxe[chr][group][c][1]-mxe[chr][group][c][0];
                mxeLen2=mxe[chr][group][c][3]-mxe[chr][group][c][2];
                upj1Len=mxe[chr][group][c][0]-mxe[chr][group][c][5];
                dnj1Len=mxe[chr][group][c][6]-mxe[chr][group][c][1];
                upj2Len=mxe[chr][group][c][2]-mxe[chr][group][c][5];
                dnj2Len=mxe[chr][group][c][6]-mxe[chr][group][c][3];


                if (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][0] and ((lastNumber<=mxeLen1 and nextN==-1) or (lastNumber==mxeLen1 and nextN==dnj1Len))) or (jS==mxe[chr][group][c][1] and jE==mxe[chr][group][c][6] and ((secondNumber<=mxeLen1 and prevN==-1) or (secondNumber==mxeLen1 and prevN==upj1Len))): ## IJC
                  key=':'.join(["MXE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_mxe[c][CT1][I]+=1;
                  c_mxe[c][CT2][I]+=1;
                  usedE[key]=1; ## okay to overwrite     
                elif (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][2] and ((lastNumber<=mxeLen2 and nextN==-1) or (lastNumber==mxeLen2 and nextN==dnj2Len))) or (jS==mxe[chr][group][c][3] and jE==mxe[chr][group][c][6] and ((secondNumber<=mxeLen2 and prevN==-1) or (secondNumber==mxeLen2 and prevN==upj2Len))): ## SJC
                  key=':'.join(["MXE",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_mxe[c][CT1][S]+=1;
                  c_mxe[c][CT2][S]+=1;
                  usedE[key]=1; ## okay to overwrite     
          ### end of MXE ###  

          ## A5SS ###
          if chr in a5ss: ## this chromosome has a5ss event(s)
            if group in a5ss[chr]: ## this group has a5ss event(s)
              for c in a5ss[chr][group]: ## for each a5ss event in this group, examine if the given junction is part of it
                if a5ss[chr][group][c][4]>a5ss[chr][group][c][1]: ## positive strand
                  if ( jS==a5ss[chr][group][c][1] and jE==a5ss[chr][group][c][4] ): ## IJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][I]+=1;
                    c_a5ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a5ss[chr][group][c][3] and jE==a5ss[chr][group][c][4] ): ## SJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][S]+=1;
                    c_a5ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif nextN == -1:
                    cStart = jE+1; cEnd=jE+lastNumber;
                    if (cStart<=(a5ss[chr][group][c][3]-(rL-junctionLength/2)+1) and cEnd<=a5ss[chr][group][c][1] and cEnd>=(a5ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A5SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a5ss[c][CT1][I]+=1;
                        c_a5ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
                else: ## negative strand
                  if ( jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][0]  ): ## IJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][I]+=1;
                    c_a5ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][2] ): ## SJC
                    key=':'.join(["A5SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a5ss[c][CT1][S]+=1;
                    c_a5ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif prevN == -1:
                    cStart = jE+1; cEnd=jE+lastNumber;
                    if (cStart>a5ss[chr][group][c][0] and cStart<=(a5ss[chr][group][c][2]-(rL-junctionLength/2)+1) and cEnd>=(a5ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A5SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a5ss[c][CT1][I]+=1;
                        c_a5ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
          ### end of A5SS ###

          ## A3SS ###
          if chr in a3ss: ## this chromosome has a3ss event(s)
            if group in a3ss[chr]: ## this group has a3ss event(s)
              for c in a3ss[chr][group]: ## for each a3ss event in this group, examine if the given junction is part of it
                if a3ss[chr][group][c][4]>a3ss[chr][group][c][1]: ## negative strand
                  if ( jS==a3ss[chr][group][c][1] and jE==a3ss[chr][group][c][4]  ): ## IJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][I]+=1;
                    c_a3ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a3ss[chr][group][c][3] and jE==a3ss[chr][group][c][4] ): ## SJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][S]+=1;
                    c_a3ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif nextN == -1:
                    cStart = jE+1; cEnd=jE+lastNumber;
                    if (cStart<=(a3ss[chr][group][c][3]-(rL-junctionLength/2)+1) and cEnd<=a3ss[chr][group][c][1] and cEnd>=(a3ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A5SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a3ss[c][CT1][I]+=1;
                        c_a3ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
                else: ## positive strand
                  if ( jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][0] ): ## IJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][I]+=1;
                    c_a3ss[c][CT2][I]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif ( jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][2] ): ## SJC
                    key=':'.join(["A3SS",str(c)]);
                    if key in usedE: ## already used;
                      continue; ## skip and go to next c
                    c_a3ss[c][CT1][S]+=1;
                    c_a3ss[c][CT2][S]+=1;
                    usedE[key]=1; ## okay to overwrite     
                  elif prevN == -1:
                    cStart = jE+1; cEnd=jE+lastNumber;
                    if (cStart>a3ss[chr][group][c][0] and cStart<=(a3ss[chr][group][c][2]-(rL-junctionLength/2)+1) and cEnd>=(a3ss[chr][group][c][2]+(rL-junctionLength/2))): ## multi-exon read supporting target
                      key=':'.join(["A5SS",str(c)]);
                      if key not in usedE: ## read didn't cover this event before, add it
                        c_a3ss[c][CT1][I]+=1;
                        c_a3ss[c][CT2][I]+=1;
                        usedE[key]=1; ## okay to overwrite
          ### end of A3SS ###

          ## RI ###
          if chr in ri: ## this chromosome has ri event(s)
            if group in ri[chr]: ## this group has ri event(s)
              for c in ri[chr][group]: ## for each ri event in this group, examine if the given junction is part of it
                if jS==ri[chr][group][c][3] and jE==ri[chr][group][c][4] : ## SJC
                  key=':'.join(["RI",str(c)]);
                  if key in usedE: ## already used;
                    continue; ## skip and go to next c
                  c_ri[c][CT1][S]+=1;
                  c_ri[c][CT2][S]+=1;
                  usedE[key]=1; ## okay to overwrite     
                if prevN == -1:
                  if cStart<=(ri[chr][group][c][4]-(rL-junctionLength/2)+1) and cEnd>=(ri[chr][group][c][4]+(rL-junctionLength/2)): ## multi-exon read supporting target
                    key=':'.join(["RI",str(c)]);
                    if key not in usedE: ## read didn't cover this event before, add it
                      c_ri[c][CT1][I]+=1;
                      c_ri[c][CT2][I]+=1;
                      usedE[key]=1; ## okay to overwrite
                if nextN == -1:
                  cStart = jE+1; cEnd=jE+lastNumber; ## just change cStart and cEnd then do the same process
                  if cStart<=(ri[chr][group][c][3]-(rL-junctionLength/2)+1) and cEnd>=(ri[chr][group][c][3]+(rL-junctionLength/2)):
                    key=':'.join(["RI",str(c)]);
                    if key not in usedE: ## read didn't cover this event before, add it
                      c_ri[c][CT1][I]+=1;
                      c_ri[c][CT2][I]+=1;
                      usedE[key]=1; ## okay to overwrite
          ### end of RI ###

          #### end of for ec #####

      else: ## it is not exonic nor junction read. proceed to the next line
        continue;

    sFile.close();  
    #edgeFile = open(s1.strip()+'.edgeCount', 'w');
    #for k in e1:
    #  edgeFile.write(k+'\t'+str(e1[k])+'\n');
    #edgeFile.close();
##### end of processSample_stranded ######

#def getId(ids,key):
#    for chr in ids: ## this chromosome has se event(s)
#        for group in ids[chr]: ## this group has skipped exon event(s)
#            if key in ids[chr][group]: 
#               event_location=ids[chr][group][key];
               #[tS,tE,uS,uE,dS,dE,strand]
#               return [chr,event_location[-1],":".join(map(str,event_location[:-1]))];


def writeInputFile2(h1,f1,id_index,cnt,sup): ## header 1,2,3, file 1,2,3, count dict, supple dict)
  ## print header first
  f1.write(h1+'\n');

  for k in sorted(sup.keys()):
    f1.write(str(id_index[k])+'\t'+str(cnt[k][CT1][I])+'\t'+str(cnt[k][CT1][S])+'\t'+str(sup[k][CT1][0])+'\t'+str(sup[k][CT1][1])+'\n');
    #f1.write(':'.join(map(str,getId(id_index,k)))+'\t'+str(cnt[k][CT1][I])+'\t'+str(cnt[k][CT1][S])+'\t'+str(sup[k][CT1][0])+'\t'+str(sup[k][CT1][1])+'\n');

##### end of writeInputFile function #######

def main():
    parser = argparse.ArgumentParser(description='''
            Given a BAM file and rMATS event files, 
            count the reads supporting the events.''')
    parser.add_argument('SE',
            help='skipped event file, generated by rMATS')
    parser.add_argument('MXE',
            help='mutually exclusive event file, generated by rMATS')
    parser.add_argument('A5SS',
            help='alternative 5 prime event file, generated by rMATS')
    parser.add_argument('A3SS',
            help='alternative 3 prime event file, generated by rMATS')
    parser.add_argument('RI',
            help='retianed intron event file, generated by rMATS')
    parser.add_argument('input_file',
            help='BAM file, previously aligned')
    parser.add_argument('outDir',
            help='output directory.')
    parser.add_argument('-rl', '--read_length', type=int, default=50,
            help='read length. Default: 50.')
    parser.add_argument('-jl', '--junction_length', type=int, default=84,
            help='junction length. Default: 84.')
    #parser.add_argument('-p', '--processor', type=int, default=1,
    #        help='number of processors to use. Default: 1.')
    print("initializing variables...")
    #### global variables  ################
    #### default configuration values #####
    #readLength=50;
    #junctionLength=84;
    global SE
    global MXE
    global A5SS
    global A3SS
    #AFE='AFE';
    #ALE='ALE';
    global RI
    global experiment 
    experiment = 'custom';
    global base
    base ='base';
    global dataType 
    dataType = 'paired'; ## either single or paired
    global samDir 
    samDir = '.';
    #input_file='.';
    #outDir = 'myOutput';
    email = 'yourID@your.domain';
    libType='fr-unstranded';
    #######################################
    #import readLength,junctionLength,outDir,input_file
    args = parser.parse_args();
    #events = parse_events(args.efname, args.column)
    #sys.stderr.write('Total number of events: {}.\n'.format(len(events)))
    #sample_info = parse_sample(args.sfname)
    #sys.stderr.write('Total number of samples: {}.\n'.format(len(sample_info['sample'])))
    #psi_diff(events, sample_info, args.ofname, args.iteration, args.seed, args.group_size, args.processor)
    global readLength; 
    readLength = args.read_length;
    global junctionLength; 
    junctionLength = args.junction_length;
    global input_file; 
    input_file = args.input_file;
    global outDir; 
    outDir = args.outDir;
    SE=args.SE;
    MXE=args.MXE;
    A5SS=args.A5SS;
    A3SS=args.A3SS;
    RI=args.RI;
    #
    print("Input File: "+input_file)
    print("Output Folder: "+outDir)
    print("Read Length: "+str(readLength))
    print("Junction Length: "+str(junctionLength))
    #
    commands.getstatusoutput('mkdir '+outDir);
    #
    ejLength = junctionLength-readLength+1; ## effective junction length
    global anchorLength; 
    anchorLength = readLength-junctionLength/2; ## anchor length
    global CT1,CT2,CT3;
    CT1,CT2,CT3=0,1,2; ### count type 1,2, or 3
    #S1,S2=0,1; ## sample 1 or sample 2
    global I,S;
    I,S=0,1; ## inclusion isoform or skipping form
    #
    #pPath = outDir+'/junctions.per.sample.pickle';
    #junctions = pickle.load(open(pPath)); ## junctions per sample
    #

    prefix = experiment;
    #
    ### open AS event files..
    #
    seFile = open(SE); ## skipped exon event file
    mxeFile = open(MXE); ## mxe event file
    a5ssFile = open(A5SS); ## a5ss event file
    a3ssFile = open(A3SS); ## a3ss event file
    riFile = open(RI); ## ri event file
    #
    ### open output files here...
    #
    JC_seFile = open(outDir+'/JC.'+prefix+'.SE.MATS.input.txt', 'w');
    JC_mxeFile = open(outDir+'/JC.'+prefix+'.MXE.MATS.input.txt', 'w');
    JC_a5ssFile = open(outDir+'/JC.'+prefix+'.A5SS.MATS.input.txt', 'w');
    JC_a3ssFile = open(outDir+'/JC.'+prefix+'.A3SS.MATS.input.txt', 'w');
    JC_riFile = open(outDir+'/JC.'+prefix+'.RI.MATS.input.txt', 'w');

    global chunk;
    chunk=1000; ## to speed up the sam file processing

    ##############################
    ## BUILD EVENT DICTIONARIES ##
    ##############################
    print("Building event dictionaries from annotation...")
    global se,mxe,a5ss,a3ss,afe,ale,ri
    se={};mxe={};a5ss={};a3ss={};afe={};ale={};ri={}; ## 7 dictionaries
    global id_se,id_mxe,id_a5ss,id_a3ss,id_afe,id_ale,id_ri
    id_se={};id_mxe={};id_a5ss={};id_a3ss={};id_afe={};id_ale={};id_ri={};
    global e_se,e_mxe,e_a5ss,e_a3ss,e_afe,e_ale,e_ri
    e_se={};e_mxe={};e_a5ss={};e_a3ss={};e_afe={};e_ale={};e_ri={}; ## exons dictionaries
    global c_se,c_mxe,c_a5ss,c_a3ss,c_afe,c_ale,c_ri
    c_se={};c_mxe={};c_a5ss={};c_a3ss={};c_afe={};c_ale={};c_ri={}; ## count dictionaries
    global s_se,s_mxe,s_a5ss,s_a3ss,s_afe,s_ale,s_ri
    s_se={};s_mxe={};s_a5ss={};s_a3ss={};s_afe={};s_ale={};s_ri={}; ## supple dict. effective length of inclusion or skipping form for CT1,CT2,CT3
    #
    #### SE ######
    c=0;     ## count
    numSE=0; ## number of SE
    numSEDup=0; ## duplicate SE id 
    line=seFile.readline(); ## skipping header
    for line in seFile: ## process skipped exon events file
      c+=1;
      ele = line.strip().split('\t');
      id = int(ele[0]);
      chr = ele[3]; 
      if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
        chr = 'chr'+chr;
      strand = ele[4];
      tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
      uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
      dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord

      tLen=tE-tS; ## target exon length
      uLen=min(uE-uS,junctionLength/2); ## upstream exon length
      dLen=min(dE-dS,junctionLength/2); ## downstream exon length
      id_se[id]=":".join(map(str,[chr,strand,tS,tE,uS,uE,dS,dE]));
      e_se[id] = [tS,tE,uS,uE,dS,dE];
      c_se[id] = getInitialCounts(); 
      I_0= min(tLen,readLength)+readLength-4*anchorLength+2; ## effective inclusion form length for JC
      S_0= readLength-2*anchorLength+1; ## effective skipping form length for JC
      I_1= readLength-2*anchorLength+tLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target
      I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet 
  
      s_se[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3 

      group = range(uS/chunk, uE/chunk+1) +  range(tS/chunk, tE/chunk+1) +  range(dS/chunk, dE/chunk+1); ## groups this event could belong 
      group = list(set(group));  ## remove duplicate groups

      if chr in se: ## already processed this chromosome
        for i in group: ## for each possible group
          if i in se[chr]: ## this group is already there
            if id in se[chr][i]: ## not likely but this group already has the id
              numSEDup+=1;
            else: ## new SE ID
              se[chr][i][id] = [tS,tE,uS,uE,dS,dE,strand]; ## skipping event with coords
              numSE+=1;
          else: ## new group to this chromosome
            se[chr][i]={};
            se[chr][i][id] = [tS,tE,uS,uE,dS,dE,strand]; ## skipping event with coords
            numSE+=1;
      else: ## first time accesing this chromosome
        se[chr]={};
        for i in group: ## for each possible group
          se[chr][i]={};
          se[chr][i][id] = [tS,tE,uS,uE,dS,dE,strand]; ## skipping event with coords
          numSE+=1;
    #
    #### MXE ####
    c=0;     ## count
    numMXE=0; ## number of MXE
    numMXEDup=0; ## duplicate MXE id
    line=mxeFile.readline(); ## mxe header
    for line in mxeFile: ## process mxe events file
      c+=1;
      ele = line.strip().split('\t');
      id = int(ele[0]);
      chr = ele[3];
      if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
        chr = 'chr'+chr;
      strand = ele[4]; ## '+' or '-'
      tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
      sS = int(ele[7]); sE = int(ele[8]); ## second exon coord
      if strand=='-': ## negative strand, switch target exon and second exon
        tS = int(ele[7]); tE = int(ele[8]); ## target exon coord
        sS = int(ele[5]); sE = int(ele[6]); ## second exon coord    
      uS = int(ele[9]); uE = int(ele[10]); ## upstream exon coord (samller coord)
      dS = int(ele[11]); dE = int(ele[12]); ## downstream exon coord (bigger coord)

      tLen=tE-tS; ## target exon length
      sLen=sE-sS; ## second exon length
      uLen=min(uE-uS,junctionLength/2); ## upstream exon length
      dLen=min(dE-dS,junctionLength/2); ## downstream exon length
      id_mxe[id]=":".join(map(str,[chr,strand,tS,tE,sS,sE,uS,uE,dS,dE]));
      e_mxe[id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## target, second, up, down (This is different from the input file)
      c_mxe[id] = getInitialCounts();
      I_0= min(tLen,readLength)+readLength-4*anchorLength+2; ## effective inclusion form length for JC
      S_0= min(sLen,readLength)+readLength-4*anchorLength+2; ## effective skipping form length for JC
      I_1= readLength-2*anchorLength+tLen+1; ## effective inclusion form length for JC+reads on target
      S_1= readLength-2*anchorLength+sLen+1; ## effective inclusion form length for JC+reads on target
      I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

      s_mxe[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

      group = range(uS/chunk,uE/chunk+1)+range(tS/chunk,tE/chunk+1);
      group = group + range(sS/chunk,sE/chunk+1)+range(dS/chunk,dE/chunk+1);## groups this event could belong
      group = list(set(group));  ## remove duplicate groups

      if chr in mxe: ## already processed this chromosome
        for i in group: ## for each possible group
          if i in mxe[chr]: ## this group is already there
            if id in mxe[chr][i]: ## not likely but this group already has the id
              numMXEDup+=1;
            else: ## new MXE ID
              mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE,strand]; ## mxe event with coords
              numMXE += 1;
          else: ## new group to this chromosome
            mxe[chr][i]={};
            mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE,strand]; ## mxe event with coords
            numMXE += 1;
      else: ## first time accesing this chromosome
        mxe[chr]={};
        for i in group: ## for each possible group
          mxe[chr][i]={};
          mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE,strand]; ## mxe event with coords
          numMXE+=1;
    #
    #### A5SS ####
    c=0;     ## count
    numA5SS=0; ## number of A5SS
    numA5SSDup=0; ## duplicate A5SS id
    line=a5ssFile.readline(); ## skipping header
    for line in a5ssFile: ## process a5ss events file
      c+=1;
      ele = line.strip().split('\t');
      id = int(ele[0]);
      chr = ele[3];
      if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
        chr = 'chr'+chr;
      strand = ele[4];
      lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
      sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
      fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord

      lLen=min(lE-lS,junctionLength/2); ## long exon length
      sLen=min(sE-sS,junctionLength/2); ## short exon length
      fLen=min(fE-fS,junctionLength/2); ## flanking exon length
      aLen=lE-lS-(sE-sS); ## alternative SS region length
      id_a5ss[id]=":".join(map(str,[chr,strand,lS,lE,sS,sE,fS,fE]));
      e_a5ss[id] = [lS,lE,sS,sE,fS,fE];
      c_a5ss[id] = getInitialCounts();
      I_0= min(aLen,readLength)+readLength-4*anchorLength+2; ## effective inclusion form length for JC
      S_0= readLength-2*anchorLength+1; ## effective skipping form length for JC
      I_1= readLength-2*anchorLength+aLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target
      I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

      s_a5ss[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

      group = range(lS/chunk, lE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
      group = list(set(group));  ## remove duplicate groups

      if chr in a5ss: ## already processed this chromosome
        for i in group: ## for each possible group
          if i in a5ss[chr]: ## this group is already there
            if id in a5ss[chr][i]: ## not likely but this group already has the id
              numA5SSDup+=1;
            else: ## new A5SS ID
              a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE,strand]; ## a5ss event with coords
              numA5SS+=1;
          else: ## new group to this chromosome
            a5ss[chr][i]={};
            a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE,strand]; ## a5ss event with coords
            numA5SS+=1;
      else: ## first time accesing this chromosome
        a5ss[chr]={};
        for i in group: ## for each possible group
          a5ss[chr][i]={};
          a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE,strand]; ## a5ss event with coords
          numA5SS+=1;
    #
    ### A3SS ####
    c=0;     ## count
    numA3SS=0; ## number of A3SS
    numA3SSDup=0; ## duplicate A3SS id
    line=a3ssFile.readline(); ## skipping header
    for line in a3ssFile: ## process a3ss events file
      c+=1;
      ele = line.strip().split('\t');
      id = int(ele[0]);
      chr = ele[3];
      if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
        chr = 'chr'+chr;
      strand = ele[4];
      lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
      sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
      fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord

      lLen=min(lE-lS,junctionLength/2); ## long exon length
      sLen=min(sE-sS,junctionLength/2); ## short exon length
      fLen=min(fE-fS,junctionLength/2); ## flanking exon length
      aLen=lE-lS-(sE-sS); ## alternative SS region length
      id_a3ss[id]=":".join(map(str,[chr,strand,lS,lE,sS,sE,fS,fE]));
      e_a3ss[id] = [lS,lE,sS,sE,fS,fE];
      c_a3ss[id] = getInitialCounts();
      I_0= min(aLen,readLength)+readLength-4*anchorLength+2; ## effective inclusion form length for JC
      S_0= readLength-2*anchorLength+1; ## effective skipping form length for JC
      I_1= readLength-2*anchorLength+aLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target
      I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

      s_a3ss[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

      group = range(lS/chunk, lE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
      group = list(set(group));  ## remove duplicate groups

      if chr in a3ss: ## already processed this chromosome
        for i in group: ## for each possible group
          if i in a3ss[chr]: ## this group is already there
            if id in a3ss[chr][i]: ## not likely but this group already has the id
              numA3SSDup+=1;
            else: ## new A3SS ID
              a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE,strand]; ## a3ss event with coords
              numA3SS+=1;
          else: ## new group to this chromosome
            a3ss[chr][i]={};
            a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE,strand]; ## a3ss event with coords
            numA3SS+=1;
      else: ## first time accesing this chromosome
        a3ss[chr]={};
        for i in group: ## for each possible group
          a3ss[chr][i]={};
          a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE,strand]; ## a3ss event with coords
          numA3SS+=1;
    #
    #### RI ####
    c=0;     ## count
    numRI=0; ## number of RI
    numRIDup=0; ## duplicate RI id
    line=riFile.readline(); ## skipping header
    for line in riFile: ## process ri events file
      c+=1;
      ele = line.strip().split('\t');
      id = int(ele[0]);
      chr = ele[3];
      if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
        chr = 'chr'+chr;
      strand = ele[4];
      rS = int(ele[5]); rE = int(ele[6]); ## ri exon coord (including up- and down-stream exons)
      uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
      dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord

      uLen=min(uE-uS,junctionLength/2); ## upstream exon length
      dLen=min(dE-dS,junctionLength/2); ## downstream exon length
      riLen=rE-rS-(uE-uS)-(dE-dS); ## retained exon length
      id_ri[id]=":".join(map(str,[chr,strand,rS,rE,uS,uE,dS,dE]));
      e_ri[id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
      c_ri[id] = getInitialCounts();
      I_0= min(riLen,readLength)+readLength-4*anchorLength+2; ## effective inclusion form length for JC
      S_0= readLength-2*anchorLength+1; ## effective skipping form length for JC
      I_1= readLength-2*anchorLength+riLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target
      I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

      s_ri[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

      group = range(rS/chunk, rE/chunk+1); ## groups this event could belong
      group = list(set(group));  ## remove duplicate groups

      if chr in ri: ## already processed this chromosome
        for i in group: ## for each possible group
          if i in ri[chr]: ## this group is already there
            if id in ri[chr][i]: ## not likely but this group already has the id
              numRIDup+=1;
            else: ## new RI ID
              ri[chr][i][id] = [rS,rE,uS,uE,dS,dE,strand]; ## ri event with coord
              numRI+=1;
          else: ## new group to this chromosome
            ri[chr][i]={};
            ri[chr][i][id] = [rS,rE,uS,uE,dS,dE,strand]; ## ri event with coord
            numRI+=1;
      else: ## first time accesing this chromosome
        ri[chr]={};
        for i in group: ## for each possible group
          ri[chr][i]={};
          ri[chr][i][id] = [rS,rE,uS,uE,dS,dE,strand]; ## ri event with coord
          numRI+=1;

    print("Counting Reads...")
    ####################################
    ## Running the counting functions ##
    ####################################
    #if dataType=="skipThisForNow": ## call paired..
      #logging.debug("Start processing sample_1: %s" % base_1);
      #processSample_PE(input_file, 72, 60);
    if libType=='fr-unstranded': ## unstranded..
      processSample(input_file);
    elif libType=='fr-firststrand': ## first strand..
      processSample_stranded(input_file,dataType,'first')
    elif libType=='fr-secondstrand': ## second strand..
      processSample_stranded(input_file,dataType,'second')
    ###

    print("Writing Output...")
    ######################################
    ## Writing out the rMATS input file ##
    ######################################
    header = 'ID(chr:strand:SE_s:SE_e:UP_s:UP_e:DN_s:DN_e)\tIC_'+base+'\tSC_'+base+'\tIncFormLen\tSkipFormLen';
    #
    ## SE ##
    writeInputFile2(header,JC_seFile,id_se,c_se,s_se);
    ## MXE ##
    writeInputFile2(header,JC_mxeFile,id_mxe,c_mxe,s_mxe);
    ## A5SS ##...
    writeInputFile2(header,JC_a5ssFile,id_a5ss,c_a5ss,s_a5ss);
    ## A3SS ##...
    writeInputFile2(header,JC_a3ssFile,id_a3ss,c_a3ss,s_a3ss);
    ## RI ##
    writeInputFile2(header,JC_riFile,id_ri,c_ri,s_ri);
    #
    print("Done!")
    #### close all files here ##########
    seFile.close()
    mxeFile.close()
    a5ssFile.close()
    a3ssFile.close()
    riFile.close()
    #
    JC_seFile.close()
    JC_mxeFile.close()
    JC_a5ssFile.close()
    JC_a3ssFile.close()
    JC_riFile.close()
    #
    sys.exit(0);
    #


if __name__ == '__main__':
   main()
