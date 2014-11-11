#!/share/apps/bin/python2.7
import string, os, sys, time, glob
"""
FusionGeneID    Hgene   Hchr    Hstrand Hbpt    Tgene   Tchr    Tstrand Tbp     Hg-Tg   # of seed read  BPseq
EML4-ALK_42522656-29446394      EML4        chr2    +       42522656        ALK     chr2    -       29446394    EML4-ALK            14  TAGAGCCCACACCTGGGAAAGGACCTAAAGtgtaccgccggaagcaccaggagctgcaag
"""
argvs=sys.argv

fname=argvs[1]
bamfile=argvs[2]

inf = open(fname,'r')
inf2 =open('refGene.txt','r')
gpdic={}
gpdic2={}
glist=[]
gplist=[]
gp_pat_bp_dic={}
gp_chr_dic={}

while 1:
	line=inf.readline()
	if not line: break
	lines=line.split('\t')
	if lines[0]!='FusionGeneID':
		gp=lines[9]
		hg=lines[1]
		tg=lines[5]
		hchr=lines[2]
		tchr=lines[6]
		hbp=lines[4]
		tbp=lines[8]
		gps = gp.split('-')
		if gp not in gplist:
			gplist.append(gp)
		gp_pat_bp_dic[hg]=hbp
		gp_pat_bp_dic[tg]=tbp
		gp_chr_dic[hg]=hchr
		gp_chr_dic[tg]=tchr
		if hg not in glist:
			glist.append(hg)
		if tg not in glist:
			glist.append(tg)


gene_str_dic={}
while 1:
	line= inf2.readline()
	if not line: break
        lines = line[:-1].split('\t')

        ref=lines[1]
	chr=lines[2]
	strand=lines[3]
	exstart_list=lines[9].split(',')[:-1]
	exend_list=lines[10].split(',')[:-1]
	gene=lines[12]
	gene_str_dic[gene]=strand
	if gene in glist:
		if gpdic.has_key(gene):
			list=gpdic[gene]
			list2=gpdic2[gene]
                	for ex in range(len(exstart_list)):
                        	EX=int(exstart_list[ex])
				EX2=int(exend_list[ex])
				if EX not in list:
					list.append(EX)
				if EX2 not in list2:
					list2.append(EX2)
			gpdic[gene]=list
			gpdic2[gene]=list2
		else:
			list=[]
			list2=[]
			for ex in range(len(exstart_list)):
                                EX=int(exstart_list[ex])
				EX2=int(exend_list[ex])
                                if EX not in list:
                                        list.append(EX)
				if EX2 not in list2:
					list2.append(EX2)
                        gpdic[gene]=list
			gpdic2[gene]=list2

for gp in gplist:
	gps=gp.split('-')
	hg=gps[0]
	tg=gps[1]
	hchr=gp_chr_dic[hg]
	tchr=gp_chr_dic[tg]
	hlist=gpdic[hg]
	hlist2=gpdic2[hg]
	tlist=gpdic[tg]
	tlist2=gpdic2[tg]
	hlist.sort()
	hlist2.sort()
	tlist.sort()
	tlist2.sort()
	hmin=hlist[0]
	hmax=hlist2[-1]
	tmin=tlist[0]
	tmax=tlist2[-1]

	
	print 'coverageBed -d -hist -abam '+bamfile+' -b refGenes_bed_files/'+hg+'_exSTR.bed > depth_files/'+hg+'_depth'
	os.system('coverageBed -d -hist -abam '+bamfile+' -b refGenes_bed_files/'+hg+'_exSTR.bed > depth_files/'+hg+'_depth')
	time.sleep(0.1)
	print 'coverageBed -d -hist -abam '+bamfile+' -b refGenes_bed_files/'+tg+'_exSTR.bed > depth_files/'+tg+'_depth'
	os.system('coverageBed -d -hist -abam '+bamfile+' -b refGenes_bed_files/'+tg+'_exSTR.bed > depth_files/'+tg+'_depth')
	time.sleep(0.1)

	print 'sort -k 2,2n -k 3,3n -k 4,4n depth_files/'+hg+'_depth > depth_files/'+hg+'_depth.sorted'
	os.system('sort -k 2,2n -k 3,3n -k 4,4n depth_files/'+hg+'_depth > depth_files/'+hg+'_depth.sorted')
	time.sleep(0.1)
	print 'sort -k 2,2n -k 3,3n -k 4,4n depth_files/'+tg+'_depth > depth_files/'+tg+'_depth.sorted'
	os.system('sort -k 2,2n -k 3,3n -k 4,4n depth_files/'+tg+'_depth > depth_files/'+tg+'_depth.sorted')
	time.sleep(0.1)
	
	hstrand=gene_str_dic[hg]
	tstrand=gene_str_dic[tg]

	hbp=int(gp_pat_bp_dic[hg])
	tbp=int(gp_pat_bp_dic[tg])

	hbploci=0
	tbploci=0
	hflag='f'
	tflag='f'

	Hmin=10000000000
	Hmax=0
	gpout=open(hg+'_graph.txt','w')
	cc=0
	all_xlist1=[]
	list_for_1=[]
	hTlist=[]
	hNlist=[]
	tTlist=[]
	tNlist=[]
	hT=open('depth_files/'+hg+'_depth.sorted','r')
	while 1:
		cc+=1
		all_xlist1.append(str(cc))
		line=hT.readline()
		if not line: break
		lines=line[:-1].split('\t')
		exon_nt=lines[4]
		exbd=int(lines[1])
		exbd2=int(lines[2])
		int_exon_nt=int(exon_nt)
		hTlist.append(exon_nt)
		if hflag=='f' and (hbp-6<=exbd) and (exbd<=hbp+6):
			hbploci=cc
			hflag='t'
		if hflag=='f' and (hbp-6<=exbd2) and (exbd2<=hbp+6):
			hbploci=cc

		if int_exon_nt<Hmin:
			Hmin=int_exon_nt
		if Hmax<int_exon_nt:
			Hmax=int_exon_nt
		nu=lines[2]
		if nu=='1':
			list_for_1.append(cc)

	for h in range(len(hTlist)):
		gpout.write(str(h)+'\t'+hTlist[h]+'\n')
		gpout.flush()
		h+=1
	gpout.close()



	Tmin=10000000000
	Tmax=0
	gpout2=open(tg+'_graph.txt','w')
	cc=0
	all_xlist2=[]
	list_for_2=[]
	tT=open('depth_files/'+tg+'_depth.sorted','r')
	while 1:
		cc+=1
		all_xlist2.append(str(cc))
		line=tT.readline()
		if not line: break
		lines=line[:-1].split('\t')
		exon_nt=lines[4]
		exbd=int(lines[1])
		exbd2=int(lines[2])
		int_exon_nt=int(exon_nt)
		tTlist.append(exon_nt)
		if tflag=='f' and (tbp-6<=exbd) and (exbd<=tbp+6):
			tbploci=cc
			tflag='t'
		if tflag=='f' and (tbp-6<=exbd2) and (exbd2<=tbp+6):
			tbploci=cc


		if int_exon_nt<Tmin:
			Tmin=int_exon_nt
		if Tmax<int_exon_nt:
			Tmax=int_exon_nt
		nu=lines[2]
		if nu=='1':
			list_for_2.append(cc)

	for t in range(len(tTlist)):
		gpout2.write(str(t)+'\t'+tTlist[t]+'\n')
		gpout2.flush()
		t+=1


	all_xlist1=all_xlist1[:-1]
	all_xlist2=all_xlist2[:-1]

	tout=open('temp.py','w')
	totalw='import matplotlib\n'
	totalw+='matplotlib.use(\"Agg\")\n'
	totalw+='import pylab\n'
	totalw+='import numpy as n\n'
	totalw+='import matplotlib.mlab as mlab\n'
	totalw+='import matplotlib.pyplot as plt\n'
	totalw+='import matplotlib.legend as legend\n'
	totalw+='fig = plt.figure()\n'
	totalw+='fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)\n'
	totalw+='ax1 = fig.add_subplot(111)\n'

	totalw+='index = '+str(all_xlist1)+'\n'
	tout.write(totalw)
	tout.flush()

	list_for_1=list_for_1[1:]
	list_for_2=list_for_2[1:]
	list_for_1_=[]
	totalw=''
	if hstrand=='-':
		hTlist.reverse()
		for fo in list_for_1:
			list_for_1_.append(len(all_xlist1)-fo)
		totalw+='y1 = '+str(hTlist)+'\n'
	else:
		list_for_1_=list_for_1
		totalw+='y1 = '+str(hTlist)+'\n'

	totalw+='ax1.plot(index,y1,color=\"red\")\n'
	for f in list_for_1_:
		totalw+='ax1.axvline(x='+str(f)+'., color=\"c\",ls=\"dashed\")\n'

	totalw+='ax1.set_title(\"Coverage Plot for '+hg+'\")\n'
	totalw+='ax1.set_xlabel(\"mRNA nucleotide\")\n'#\\n(Dash line means exon boundary)\")\n'
	totalw+='ax1.set_ylabel(\"Depth\")\n'
	Hlen=len(all_xlist1)/10
	Hgap=Hmax-Hmin
	Hp=Hgap/5
	Hp2=Hgap/4
	#totalw+='ax1.annotate(\"Tumor\", xy=('+str(Hlen)+','+str(Hmax-Hp)+'), color=\"red\", weight=\"bold\")\n'

	if hstrand=='-': 
		hbploci2=len(all_xlist1)-hbploci
	else:
		hbploci2=hbploci
	totalw+='ax1.annotate(\"BP\", xy=('+str(hbploci2)+',0), xytext=('+str(hbploci2)+','+str(Hp)+'), arrowprops=dict(facecolor=\"yellow\", shrink=0.02), horizontalalignment=\"center\")\n'
	totalw+='fig.savefig(\"'+hg+'.png\")\n'


	tout.write(totalw)
	tout.flush()
	os.system('/share/apps/bin/python2.7 temp.py')

	time.sleep(0.1)


	totalw=''
	tout=open('temp.py','w')
	totalw='import matplotlib\n'
	totalw+='matplotlib.use(\"Agg\")\n'
	totalw+='import pylab\n'
	totalw+='import numpy as n\n'
	totalw+='import matplotlib.mlab as mlab\n'
	totalw+='import matplotlib.pyplot as plt\n'
	totalw+='import matplotlib.legend as legend\n'
	totalw+='fig = plt.figure()\n'
	totalw+='fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)\n'
	totalw+='ax1 = fig.add_subplot(111)\n'

	totalw+='index = '+str(all_xlist2)+'\n'
	tout.write(totalw)
	tout.flush()

	list_for_2_=[]
	totalw=''
	if tstrand=='-':
		tTlist.reverse()
		for fo in list_for_2:
			list_for_2_.append(len(all_xlist2)-fo)
		totalw+='y1 = '+str(tTlist)+'\n'
	else:
		list_for_2_=list_for_2
		totalw+='y1 = '+str(tTlist)+'\n'

	totalw+='ax1.plot(index,y1,color=\"red\")\n'
	for f in list_for_2_:
		totalw+='ax1.axvline(x='+str(f)+'., color=\"c\",ls=\"dashed\")\n'

	totalw+='ax1.set_title(\"Coverage Plot for '+tg+' in \")\n'
	totalw+='ax1.set_xlabel(\"mRNA nucleotide\")\n'#\\n(Dash line means exon boundary)\")\n'
	totalw+='ax1.set_ylabel(\"Depth\")\n'
	Tlen=len(all_xlist2)/10
	Tgap=Tmax-Tmin
	Tp=Tgap/5
	Tp2=Tgap/4
	#totalw+='ax1.annotate(\"Tumor\", xy=('+str(Tlen)+','+str(Tmax-Tp)+'), color=\"red\", weight=\"bold\")\n'

	if tstrand=='-':
		tbploci2=len(all_xlist2)-tbploci
	else:
		tbploci2=tbploci
	totalw+='ax1.annotate(\"BP\", xy=('+str(tbploci2)+',0), xytext=('+str(tbploci2)+','+str(Tp)+'), arrowprops=dict(facecolor=\"yellow\", shrink=0.02), horizontalalignment=\"center\")\n'

	totalw+='fig.savefig(\"'+tg+'.png\")\n'


	tout.write(totalw)
	tout.flush()
	os.system('/share/apps/bin/python2.7 temp.py')	

	print 'coverage plot creation finished for '+gp
