
import string, os, sys, glob

dirname=sys.argv[1]
dirname2=sys.argv[2]

inf2=open('known_kinase_fusion_list.txt','r')
inf0=open('refGene_sorted.txt')

ref_gsymbol_dic={}
ref_loci_dic={}
reflist=[]
while 1:
        line=inf0.readline()
        if not line: break
        lines=line[:-1].split('\t')
        ref=lines[1]
        chr=lines[2]
        tstart=lines[4]
        tend=lines[5]
        gsymbol=lines[12]
        ref_gsymbol_dic[ref]=gsymbol
        ref_loci_dic[ref]=chr+':'+tstart+':'+tend
        if ref not in reflist:
                reflist.append(ref)

print dirname
print dirname2

pslx_list=glob.glob(dirname+'/*.pslx.txt')
fa_list=glob.glob(dirname2+'/*.fa')

kin_fusion={}
kin_bp={}
bp_kin={}
bplist=[]
while 1:
	line=inf2.readline()
	if not line: break

	lines=line[:-1].split('\t')
	kin=lines[0]
	fusion=lines[1]
	BP=lines[2]
	kin_fusion[kin]=fusion
	kin_bp[kin]=BP
	bps=BP.split('or')
	for bp in bps:
		bp=bp.strip()
		bplist.append(bp)
		bp_kin[bp]=kin

print bplist
outf=open('kinase_fusion_scan_result_for_PatAllT.txt','w')

print bp_kin
for fn in pslx_list:
#fn="/storage/home/iamlife/nsclc100/pslx_dir_newoption/Pat05T_1.pslx.txt"
#if fn=="/storage/home/iamlife/nsclc100/pslx_dir_newoption/Pat05T_1.pslx.txt":


	Tdic={}
	Tlist=[]

	print '----------'+fn

	for bp in bplist:
		bps=bp.split(':')
		Chr=bps[0]
		loci=int(bps[1])
		kin=bp_kin[bp]
		if kin=='RET':
			print bp
			print loci-6
			print loci
			print loci+6
			
		templine=''
		qname0=''
		qname1=''
		c=0
		d=0
		flag='false'
		fin=open(fn,'r')
		while 1:

			line= fin.readline()
			if not line: break
			
			if line[:4]=='psl:':
				qname=line.split('\t')[9]
				if c==0:
					d+=1
					qname0=qname
					templine=line
					c=1
				else:
					if qname0==qname:
						qname0=qname
						templine=templine+line	
						continue
					else:
						tlines=templine.split('\n')[:-1]
						tlen=len(tlines)
						if (1<tlen) and (tlen<12):
							for tline in tlines:
		                                        	ttlines=tline.split('\t')
		                                        	qname_=ttlines[9]
		                                        	tname=ttlines[13]
		                                        	tstart=int(ttlines[15])
		                                        	tend=int(ttlines[16])
								if tname==Chr:

									if qname0=='8371-12':
										print templine
										print tstart
										print tend
									if ((((loci-6)<=tstart) and (tstart<=(loci+6))) or (((loci-6)<=tend) and (tend<=(loci+6)))):
										#flag='true'
										Tdic[bp_kin[bp]]=templine
										#print templine
						
						qname0=qname
						templine=line


		

	Tlist=Tdic.keys()
	for T in Tlist:
		glist=[]
		TT=Tdic[T]
		Rtempline=TT.split('\n')[:-1]
		for Rline in Rtempline:
			Rlines=Rline.split('\t')
			Rlineqname=Rlines[9]
			Rname=Rlines[13]
			Rstart=int(Rlines[15])
			Rend=int(Rlines[16])


			for ref in reflist:
				locis=ref_loci_dic[ref].split(':')
				chr=locis[0]
				start=int(locis[1])
				end=int(locis[2])
				if chr==Rname:
					if (((start-6<=Rstart) and (Rstart<=end+6)) or ((start-6<=Rend) and (Rend<=end+6))):
						g=ref_gsymbol_dic[ref]
						glist.append(g)

		setg=set(glist)
		listg=list(setg)

		if len(listg)<len(glist):
			#print glist
			#print listg
			print 'pass'
		else:
			outf.write('>'+T+'\t'+str(listg)+'\n'+TT)
			print '>'+T+'\t'+str(listg)+'\n'+TT
			#print listg





