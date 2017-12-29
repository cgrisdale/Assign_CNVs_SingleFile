#!/usr/bin/python
#Cameron Grisdale
#Dec 14, 2017

import sys
import re
from collections import Counter
import resource

#Compare output of CNV software control-freec against gene annotations

def Read_in_cnv(File):
  '''Read in CNV file, return list'''

  cnvl,cnvd,chrm,start,end,copyn,typea,tml,chrmd=[],{},'','','','','',[],{}

  with open(File, 'r') as f:

    fname=str(File)

    for line in f:
      line=line.strip()
      chrm,start,end,copyn,typea=line.split('\t')

      #Define ranges #ranges end at -1 ie. (5,8) is 5,6,7
      cnvrange=xrange(int(start),int(end))

      tml=[chrm,cnvrange,copyn,typea,fname]

      #Load biglist
      cnvl.append(tml)

      #Use range as key in internal dictionary
      #if cnvd[chrm]:
      #  cnvd[chrm][cnvrange]=tml

      #else: #initialize internal dict
      #  cnvd[chrm]={}
      #  cnvd[chrm][cnvrange]=tml

    #print len(cnvd[fname])

  return cnvl

def Read_in_GTF(gtf):
  tmp,gtfdict,gtfcount,dcount=[],{},0,0

  with open(gtf, 'r') as f:

    for line in f:
      gtfcount+=1
      chrm,start,end,strand,gn,ensid=line.strip().split()
      tmp=[chrm,xrange(int(start),int(end)),strand,gn,ensid]

      try:
        if tmp[1][0]>tmp[1][-1] or tmp[1][0]==tmp[1][-1]:
          print "Exiting: end pos is NOT higher than start pos in GTF",start,end,gn,ensid
          sys.exit(0)
      except IndexError:
        print line,tmp
        sys.exit(0)

      if chrm in gtfdict:
        gtfdict[chrm][tmp[1]]=tmp

      else:#initialize internal dict
        gtfdict[chrm]={}
        gtfdict[chrm][tmp[1]]=tmp

  for x in gtfdict:
    num=len(gtfdict[x])
    dcount+=num

  if gtfcount!=dcount:
    print "dict length doesn't match num of lines:",gtfcount,dcount
    sys.exit(0)

  return gtfdict

def Compare_ranges(cnvlist,gtfd):
  '''Go through list of CNV dict, compare to gene coordinate (dict of dict's), count frequency of recurrence'''
  cnvd,genes,types,samples,geneslist={},{},[],[],[]
  tmp=[0,xrange(1),1,'gain','test','gene1','ENSG00000'] #need initial values to avoid IndexError on first round through loop

  #Go through all CNV lines
  for x in cnvlist:
    chrm,cnvr,cpn,etype,filen=x[0:]
    generd={} #make blank gene range dict for each cnv line iteration

    #Go through gtfd for current chromosome, key=xrange
    for i in gtfd[chrm]:

      if range(max(cnvr[0],i[0]),min(cnvr[-1],i[-1])+1): #range matches range key
        #ranges overlap, make note of gene matching cnv region
        genen,ensid=gtfd[chrm][i][3],gtfd[chrm][i][4]
        #print x,i,cnvd

        if cnvr!=tmp[1]: #first match for this CNV
          types.append(etype) #gain or loss
          #samples.append(filen)

        #Make dict for genes with CNV-Gene overlap size so they can be checked when there are duplicates etc.

        generange=xrange(max(cnvr[0],i[0]),min(cnvr[-1],i[-1])+1) #CNV-Gene(current GTF line) overlap
        generange1=len(generange)
        tmp=[chrm,cnvr,cpn,etype,filen,genen,ensid]

        if chrm in cnvd:

          if cnvr in cnvd[chrm]: #cur cnvr in dict (ie. overlapping gene(s) already found)
            print "cnvr in dict"
            generd=cnvd[chrm][cnvr][4]
            generd[genen]=generange1
            tmpl=[chrm,cnvr,cpn,etype,generd]
            cnvd[chrm][cnvr]=tmpl

          else: #current cnvr not in dict, other cnvr's could contain gene
            print "cnvr not in dict"
            generd[genen]=generange1
            tmpl=[chrm,cnvr,cpn,etype,generd]
            cnvd[chrm][cnvr]=tmpl

        else: #Chrm being added, gene can not be in dict yet so no need to check ranges
          print "chrm not in cnvd yet, initialize and add it"
          cnvd[chrm]={}
          generd[genen]=generange1
          tmpl=[chrm,cnvr,cpn,etype,generd]
          cnvd[chrm][cnvr]=tmpl

        #Use dict for gene counts; don't need counts because gene can only be mapped once to CNV range from a single file
        if genen in genes: #
          if filen in genes[genen]:
            pass #gene already exists with current sample name
          else:
            genes[genen].append(filen)
        else: #First time genename is added to dict, along with current sample name
          genes[genen]=[filen]

  print '\t'.join([k+': '+str(v) for k,v in Counter(types).items()]),'\n'

  return cnvd

def Check_overlaps(adict,achrm,agene,arange,grd,cnvrange,copyn,cptype): #(cnvd,chrm,genen,generange1,generd,cnvr,cpn,etype)
  '''This is run every time the ranges overlap between cnv list and gtf dict'''
  delg,delchr,delrange='',0,xrange(0)
  for k,v in adict[achrm].items(): #key=range, value=list including genedict
    if agene in v[4]: #if current gene is in dict of genes
      print "Gene exists",v[4][agene],agene,arange
      #for y,z in v[4].items(): #gene, length; this dictionary should be small (most <20, few ~1000)
      if v[4][agene]<arange: #current gene-GTF overlap is larger than the one in cnvdict
        #delg,delchr,delrange=agene,achrm,cnvrange #del v[4][agene]
        del v[4][agene]
        grd[agene]=arange #genename:generange1
        templist=[achrm,cnvrange,copyn,cptype,grd]
        adict[achrm][cnvrange]=templist #add current one, break, then delete duplicate
        #break
        return adict
      else:
        pass

    else: #gene doesn't exist in any ranges genedict, add it
      grd[agene]=arange
      templist=[achrm,cnvrange,copyn,cptype,grd]
      adict[achrm][cnvrange]=templist
      return adict

  #Found duplicate to delete
  #if delg in adict[delchr][delrange][4]:
  #  del adict[delchr][delrange][4][delg]
  #else:
  #  sys.exit("Key/genename for deletion was not found at deletion step")
  #adict[achrm][cnvrange]=templist
  #return adict
        

def Get_stats(cnvs):
  '''Go through CNVs to get average and cumulative size per gain/loss, for plotting'''
  #cnvs=list of lists: tml=[chrm,cnvrange,copyn,typea,fname]
  altd={}
  altd['gain']={}
  altd['loss']={}
  altd['gain']['arange'],altd['loss']['arange']=[],[]
  altd['gain']['copy'],altd['loss']['copy']=[],[]
  for x in cnvs:
    rangesize,copyn=x[1][-1]-x[1][0],x[2]

    if x[3]=='gain':
      altd['gain']['arange'].append(rangesize)
      altd['gain']['copy'].append(copyn)

    elif x[3]=='loss':
      altd['loss']['arange'].append(rangesize)
      altd['loss']['copy'].append(copyn)

    else:
      sys.exit("Not gain or loss, exiting")

  grsum,lrsum=sum(altd['gain']['arange']),sum(altd['loss']['arange'])
  gcsum,lcsum=sum([int(q) for q in altd['gain']['copy']]),sum([int(q) for q in altd['loss']['copy']])
  grlen,lrlen,gclen,lclen=len(altd['gain']['arange']),len(altd['loss']['arange']),len(altd['gain']['copy']),len(altd['loss']['copy'])
  if grsum==0 or lrsum==0 or gcsum==0 or lcsum==0 or grlen==0 or lrlen==0 or gclen==0 or lclen==0:
    print grsum,lrsum,gcsum,lcsum,grlen,lrlen,gclen,lclen
    sys.exit("Exiting: one of the sums of gain or loss events is equal to zero")
  gravg,lravg,gcavg,lcavg=grsum/grlen,lrsum/lrlen,gcsum/gclen,lcsum/lclen
  print "Cumulative gain:",grsum,"Cumulative loss:",lrsum,'\n'
  print "Avg gain:",gravg,"Avg loss:",lravg,"Avg copy gain:",gcavg,"Avg copy loss:",lcavg,'\n'
  finallist=[grsum,lrsum,gravg,lravg,gcavg,lcavg]
  return finallist


if __name__ == "__main__":

  if len(sys.argv)<4:
    GTF=sys.argv[1]
    Files=sys.argv[2] #

  else:
    sys.exit("Script works for one GTF and one CNV file only; exiting")

  CV,QL=[],[]

  filename=str(Files)

  gtfd=Read_in_GTF(GTF)

  #for i in range(len(Files)):
  cv=Read_in_cnv(Files)
  print "Number of CNVs examined: ",len(cv),'\n'
  CVstats=Get_stats(cv)
  cvoutstats='\t'.join(str(z) for z in CVstats)+'\n'
  header="G_range_sum"+'\t'+"L_range_sum"+'\t'+"G_range_avg"+'\t'+"L_range_avg"+'\t'+"G_copy_avg"+'\t'+"L_copy_avg"+'\n'
  outstats=open(filename+'.stats.tsv', 'w')
  outstats.write(header)
  outstats.write(cvoutstats)
  outstats.close()

  #Flatten list of lists of lists down one level to list of lists
  #QL=[item for sublist in CV for item in sublist]

  c=Compare_ranges(cv,gtfd)

  outf=open('CNV.outfile.'+filename+'.tsv', 'w')
  y=0
  for l,m in c.items():

    for k,v in m.items():
      y+=1
      print k,v
      if y>10:
        sys.exit(0)
      tmpx=[v[0],v[1][0],v[1][-1],v[2],v[3]]
      genec=len(v[4])
      genenm=';'.join(str(z) for z in v[4])
      tmpx.append(genec) #append gene count and string of gene names separated by ;
      tmpx.append(genenm)
      myline='\t'.join(str(y) for y in tmpx)
      outf.write(myline)
      outf.write('\n')
  outf.close()



'''
  #print Counter(samples).most_common(10),'\n'
  #for k,v in genes.items(): #Go through dict of genes with samples as values
  #  for i in range(len(v)): #For each sample in value, add genename to list for counting
  #    geneslist.append(k)
  #print ' '.join([k+':'+str(v) for k,v in Counter(geneslist).most_common(10)]),'\n'
'''



