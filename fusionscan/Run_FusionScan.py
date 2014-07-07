#!/usr/bin/python
########################################################################
# Title: Run_FusionScan.py
# Description: pre-procesing of	FusionScan
# Copyright	: Copyright	(c)	2013
# Company: Ewha	Womans University
# @author: Pora	Kim	(iamlife@hanmail.net)
# @version 0.9
########################################################################

import sys, os, glob, string, re
import time, thread, datetime
import os.path
from subprocess import call

argvs=sys.argv

usage = 'Usage:'
usage = usage + '\n\tpython Run_FusionScan.py  <reads1[,reads2]>  <output-dir>  <read-length>  <gfserver host:port>  [options]'
usage = usage + '\nOptions:'
usage = usage + '\n\t --phred33                            qualities are Phred+33                                 [ default           ]'
usage = usage + '\n\t --phred64                            qualities are Phred+64'
usage = usage + '\n\t -P/--threads                         number of threads                          <int>       [ default:  8       ]'
usage = usage + '\n\t -ms/--min-seed                       minimum number of seed reads               <int>       [ default:  2       ]'
usage = usage + '\n\t -md/--min-distance                   minimum distance between two genes         <int>       [ default:  50000  ]'
usage = usage + '\n\t -kfgs/--known-fusion-gene-search     searching known kinase fusion genes'
usage = usage + '\nExample:'
usage = usage + '\n\t python Run_FusionScan.py Test.fa Test 75 00.00.00.00:810 --phred33 -P 1 -ms 2'
dir_name=''
thread_num=1
file_format=''
score=''
known_fusion_search=''
fq_trim="true"
fqlist=['fq','fastq']
falist=['fa','fasta']
flist=[]
trim_ok=''
min_seed=0
time0=datetime.datetime.now()
timea='['+str(time0).split('.')[0]+']'
print timea + ' FusionScan (v0.9)'



########################################################################
# rescue processing
########################################################################
gp_ct_dict={}
def do_rescue_pr(pslx_file_name, align_dir, loc_hash, lockeys):
	t_hash={}
	llist=[]
	bipslx1=open(pslx_file_name,'r')
	temp=[]
	tkeys=[]
	qnn=0
	while 1:
		line = bipslx1.readline()
		if not line:break
		lines4 = line.split('\t')
		if len(lines4)>10:
			if line[:3]=='psl':
				mat=lines4[0].split('psl: ')[1]
				match=int(mat)
				tstart=int(lines4[15])
				tend=int(lines4[16])
				locc=lines4[15]	+ '_' +	lines4[16];
				if locc	not	in llist:
					llist.append(locc)
					for ff in lockeys:
						locs=loc_hash[ff].split('_')
						start=int(locs[0])
						bp=int(locs[1])
						end=int(locs[2])

						if (start-6<tstart) and (tstart<bp) and (bp<tend) and (tend<end+6):
							qnn+=1
							total=str(qnn)+'   ' + str(match) + ' ' + str(start) + ' ' + str(tstart) + ' ' + str(bp) + ' ' + str(tend) + ' ' + str(end)
							if t_hash.has_key(ff):
								temp=t_hash[ff]
								temp.append(total)
								t_hash[ff]=temp
							else:
								if ff not in tkeys:
									tkeys.append(ff)
								temp=[]
								temp.append(total)
								t_hash[ff]=temp
						elif (start-6<tend) and (tend<bp) and (bp<tstart) and (tstart<end+6):
							qnn+=1
							total=str(qnn)+'   ' + str(match) + ' ' + str(start) + ' ' + str(tstart) + ' ' + str(bp) + ' ' + str(tend) + ' ' + str(end)
							if t_hash.has_key(ff):
								temp=t_hash[ff]
								temp.append(total)
								t_hash[ff]=temp
							else:
								if ff not in tkeys:
									tkeys.append(ff)
								temp=[]
								temp.append(total)
								t_hash[ff]=temp

	fr='_1_info'
	for t in tkeys:
		tlist=t_hash[t]
		tfile=open(dir_name+'/'+project_name+'/alignment_dir/'+t+fr,'w')		
		for r in tlist:
			tfile.write(r+'\n')

	print 'finished do_rescue'



########################################################################
# blat result filtering
########################################################################
def filtering3(psl_results):
	chrlist=[]
	remain_chrlist=[]
	remain_result_list=[]
	repeat_flag='true'
	result_list=[]
	results={}
	mark=''
	empty_result_list=[]
	qstartendlist=[]
	matchlist=[]
	matchlistsort=[]
	removelist=[]
	removechrlist=[]
	tst_hash={}
	for u in range(len(psl_results)):
		one_psl_results=psl_results[u]
		tstart=int(one_psl_results['tstart'])
		qstart=one_psl_results['qstart']
		qend=one_psl_results['qend']
		chr=one_psl_results['tname']

		if tst_hash.has_key(chr):
			temp=tst_hash[chr]
			temp.append(tstart)
			tst_hash[chr]=temp
		else:
			temp=[]
			temp.append(tstart)
			tst_hash[chr]=temp

		total =	chr + '_' + str(qstart)	+ ',' +	str(qend)
		qstartendlist.append(total)
		
		if chr not in chrlist:
			chrlist.append(chr)

	print qstartendlist
	cc=0
	for w in range(len(qstartendlist)):
		if repeat_flag=='true':
			aa0=qstartendlist[w].split('_')
			chra=aa0[0]
			aa2=aa0[1].split(',')
			starta=int(aa2[0])
			enda=int(aa2[1])

			for w2 in range(len(qstartendlist)):
				if w!=w2:
					bb0=qstartendlist[w2].split('_')	
					chrb=bb0[0]
					bb2=bb0[1].split(',')
					startb=int(bb2[0])
					endb=int(bb2[1])

					if (startb==starta) and (endb==enda):
						if w in	removelist:
							cc+=1
						if w not in removelist:
							removelist.append(w)
						if w2 not in removelist:
							removelist.append(w2)

			if cc==5:
				repeat_flag='false'
				break
		else:
			break
		cc=0
	return repeat_flag



########################################################################
# making read alignment image
########################################################################
def make_align_image(sequ, dir_name, src_dir, temp_dir, align_dir, info_loc_hash, project_name, host, port, out_file2, loc_hash, gplist_, read_length):
	gp_ct_dict={}
	tempdir=dir_name+'/'+project_name+'/temp_dir/'
	result=[]
	seedseq_hash={}
	for gp2 in gplist_:
		gp_fasta=''
		fff=align_dir+gp2+'_1_info'
		if os.path.isfile(fff):
			qnlist=[]
			loclist=[]
			total_stli=[]
			stlist=[]
			enlist=[]
			lll=loc_hash[gp2]
			llls=lll.split('_')
			ls=int(llls[0])
			bp=int(llls[1])
			le=int(llls[2])
			bp2=bp - ls
			tstring=''
			hereseq=sequ[ls:le]

			basicseq=hereseq[0:bp2]+'BP'+hereseq[bp2:]
			stli=[]
			enli=[]
			qn_loc_hash={}
			st_qn_hash={}
			biff=open(fff,'r')
			st_temp_hash={}
			en_temp_hash={}
			while 1:
				line=biff.readline()
				if not line: break
				qn2 = line.split('   ')[0]
				qn = qn2.split('#')[0]
				lines2 = line.split('   ')[1]
				lines =	lines2.split(' ')
				if len(lines)==6:
					match=int(lines[0])
					start=int(lines[1])
					tstart=int(lines[2])
					tend=int(lines[4])
					end=int(lines[5])
					st=tstart - start
					en=tend	- start
					T=str(st)+'_'+str(en)

					if st_temp_hash.has_key(tstart):
						st_temp=st_temp_hash[tstart]
						st_temp+=1
						st_temp_hash[tstart]=st_temp
					else:
						st_temp_hash[tstart]=1

					if en_temp_hash.has_key(tend):
						en_temp=en_temp_hash[tend]
						en_temp+=1
						en_temp_hash[tend]=en_temp
					else:
						en_temp_hash[tend]=1


					if (st not in total_stli) and (T not in	loclist) and (match >= 45):
						total_stli.append(st)
						qn_loc_hash[qn]=T
						qnlist.append(qn)
						stlist.append(st)
						enlist.append(en)
						loclist.append(T)
						stli.append(st)
						enli.append(en)

						if st_qn_hash.has_key(st):
							temp=st_qn_hash[st]
							temp.append(qn)
							st_qn_hash[st]=temp
						else:
							temp=[]
							temp.append(qn)
							st_qn_hash[st]=temp

			
			stli.sort()
			enli.sort()
			total_stli.sort()
			seedc=0
			resc=0
			seedseq=''
			
			if len(stlist)!=0 and len(enlist)!=0:
				minst=0
				if stli[0]<0:
					minst=0
				else:
					minst=stli[0]

				maxen=0
				if enli[-1]<len(hereseq):
					maxen=enli[-1]
				else:
					maxen=len(hereseq)
				
				for st_ in total_stli:
					if st_qn_hash.has_key(st_):
						qns2=st_qn_hash[st_]
						for q in qns2:
							qn=q.split('#')[0]
							L=''
							st__ = 0
							en__ = 0;
							st2	= 0;
							en2=0
							ss=''
							if(qn_loc_hash.has_key(qn)):
								L=qn_loc_hash[qn]
								Ls=L.split('_')
								st__=int(Ls[0])
								en__=int(Ls[1])
								mst=0
								men=0
								if (bp2<en__) and (st__<bp2):
									stbp = bp2 - st__
									bpen = en__- bp2
									if maxen <=en__:
										men= maxen
									else:
										men= en__
									if st__	< 0:
										mst= 0
									else:
										mst= st__

									if (stbp >=20) and (bpen >= 20):
										seedc += 1
										seedseq	= hereseq[mst:men]
										seedseq_hash[gp2]=seedseq
										gp_fasta = gp_fasta+ '>' +gp2+ '_' +qn + '\n' +	seedseq	+ '\n'
									else:
										if ((stbp >= 7)	and(stbp <20)) or((bpen >= 7) and (bpen < 20)):
											resc +=	1

									s2 = st__ -minst

									for ss_ in range(s2):
										ss=ss+' '
									tstring	= tstring + ss + hereseq[mst:bp2] +' '	+ hereseq[bp2:men] + '\n'


				if info_loc_hash.has_key(gp2):
					writer10=open(tempdir+gp2 +'.fa','w')
					writer10.write('>' + gp2 + '\n' + info_loc_hash[gp2].split('\t')[8])
					writer10.flush()
					writer10.close()

					os.system(dir_name + '/src_dir/gfClient ' +host + ' ' + port + ' -minScore=18 -minIdentity=90 -nohead . ' +tempdir + gp2 + '.fa ' + tempdir + gp2 + '.psl')
					#print dir_name + '/src_dir/gfClient ' +host + ' ' + port + ' -minScore=18 -minIdentity=90 -nohead ' + dir_name+ '/src_dir/ ' +tempdir	+ gp2 + '.fa ' + tempdir + gp2 + '.psl'
					time.sleep(0.1)
					print gp2

					brg=open(tempdir + gp2 + '.psl','r')
					onealign=0
					qName=''
					tempflag='true'
					hashtemplist=[]
					blatpsl_hashtable={}
					while 1:
						lineb=brg.readline()
						if not lineb: break	
						temp = lineb.split('\t')
						matches=int(temp[0])
						qName=temp[9]
						qSize=int(temp[10])
						qStart=int(temp[11])
						qEnd=int(temp[12])
						tStart=int(temp[15])
						tEnd=int(temp[16])
						tName=temp[13]
						chrs=tName.split('_')
						if len(chrs)==1:
							temphash={}
							temphash['matches']=matches
							temphash['tstart']=tStart
							temphash['tend']=tEnd
							temphash['qstart']=qStart
							temphash['qend']=qEnd
							temphash['tname']=tName
							pct=float(float(matches)/float(qSize)*100)
							if (pct	== 100.0) and (matches >= (read_length - 1)) and (tEnd - tStart < 50000):
								onealign = 10
								break
							if (qStart <= 1) and (qEnd >= read_length - 1) and (pct >= 95.0) and (matches >= (read_length - 1)):
								onealign += 1
								if tEnd	- tStart < 50000:
									tempflag = 'false'
									break
							if (onealign < 2) and (tempflag	== 'true'):
								if blatpsl_hashtable.has_key(qName):
									hashtemplist=blatpsl_hashtable[qName]
									hashtemplist.append(temphash)
									blatpsl_hashtable[qName]=hashtemplist
								else:
									hashtemplist=[]
									hashtemplist.append(temphash)
									blatpsl_hashtable[qName]=hashtemplist

					brg.close()

					tempstr=gp2 + ' ' + str(seedc) + ' ' + str(resc)
					if (onealign < 2) and (tempflag == 'true') and (tempstr not in result):
						qnpsl_list=[]
						if blatpsl_hashtable.has_key(qName):
							qnpsl_list=blatpsl_hashtable[qName]

						blat_psl_flag='false'

						if len(qnpsl_list) == 1:
							blatpsl=qnpsl_list[0]
							blat_tstart=blatpsl['tstart']
							blat_tend=blatpsl['tend']
							distance_btw_two = abs(blat_tend - blat_tstart)	
							if distance_btw_two>= 50000:
								blat_psl_flag='true'
						elif len(qnpsl_list)> 1:
							blat_psl_flag='true'

						print blat_psl_flag
						if blat_psl_flag=='true':
							#goo=filtering3(blatpsl_hashtable[qName])
							#print 'goo= '+goo
							goo='true'
							if (goo == 'true') and (tstring!='') and (tstring!=basicseq+ '\n') and (minst < bp2) and (bp2 < maxen) and (seedc >= 1):
								result.append(tempstr)
								linei=info_loc_hash[gp2]
								lines2 = linei.split('\t')
								hg = lines2[0]
								hchr = lines2[1]
								hstr = lines2[2]
								hbp = lines2[3]
								tg = lines2[4]
								tchr = lines2[5]
								tstr = lines2[6]
								tbp = lines2[7]
								BPseq =	lines2[8]
								out_file2.write(gp2+'\t'+hg + '\t' + hchr + '\t' + hstr + '\t' + hbp + '\t' + tg + '\t'	+ tchr + '\t' +	tstr + '\t' + tbp + '\t' + hg + '-' + tg + '\t' + str(seedc) + '\t' + BPseq	+ '\n')
								out_file2.flush()
								out_file=open(align_dir + gp2 + '.RNAalign2','w')
								out_file.write('### '+gp2 + '(seed read count: ' + str(seedc) + ', rescued read count: ' + str(resc) + ')\n')
								out_file.flush()
								out_file.write(hereseq[minst:bp2] + ' ' + hereseq[bp2:maxen] + '\n')
                                                                out_file.flush()
								out_file.write(tstring)
                                                                out_file.flush()
								out_file.write(gp_fasta)
                                                                out_file.flush()
								gp_ct_dict[gp2]=seedc

					else:
						continue

				tstring=''

#	print 'finished	make_alignment_image'
	return gp_ct_dict



########################################################################
# control for parameters 
########################################################################
if len(argvs)>=5:

	flist = argvs[1].split(',') # fq file
	project_name = argvs[2] # project name
	read_length = int(argvs[3]) # thread num by	user

	p1=re.compile('.*\:.*')
	m1=p1.match(argvs[4])
	if m1:
		host_ports = argvs[4].split(':')
		host=host_ports[0]
		port=host_ports[1]
	else:
		print 'Input host:port information.'
		print usage
		sys.exit(1)

	fformat=flist[0].split('.')[1]
	if fformat in fqlist:
		file_format='fastq'
	elif fformat in	falist:
		file_format='fasta'		
		fq_trim='false'
	line=''
	for a in argvs:
		line=line+a+' '

	pwdline=os.popen('pwd')
	for line2 in pwdline.xreadlines():
		if line2[0]=='/':
			dir_name=line2[:-1]

	p2=re.compile('.*--phred33.*')
	m2=p2.match(line)
	if m2:
		score='-Q33'

	p2=re.compile('.*--phred64.*')
	m2=p2.match(line)
	if m2:
		score=''
	else:
		score='-Q33'

	p52=re.compile('.*--threads.*')
	m52=p52.match(line)
	if m52:
		thread_num=int((line.split('--threads')[1].split('-')[0].strip()))

	p5=re.compile('.*-P.*')
	m5=p5.match(line)
	if m5:
		thread_num=int((line.split('-P')[1].split('-')[0].strip()))

	p4=re.compile('.*-ms.*')
	m4=p4.match(line)
	if m4:
		min_seed=int((line.split('-ms')[1].split('-')[0].strip()))

	p42=re.compile('.*--min-seed.*')
	m42=p42.match(line)
	if m42:
		min_seed=int((line.split('--min-seed')[1].split('-')[0].strip()))

	p3=re.compile('.*-md.*')
	m3=p3.match(line)
	if m3:
		min_dis=int((line.split('-md')[1].split('-')[0].strip()))

	p32=re.compile('.*--min-distance.*')
	m32=p32.match(line)
	if m32:
		min_dis=int((line.split('--min-distance')[1].split('-')[0].strip()))


	p8=re.compile('.*-kfgs.*')
	m8=p8.match(line)
	if m8:
		known_fusion_search='go'

	p82=re.compile('.*--known-fusion-gene-search.*')
	m82=p82.match(line)
	if m82:
		known_fusion_search='go'


	if len(flist)==0:
		print 'FusionScan error	: Input	read file is needed.'
		print usage
		sys.exit(1)
	if read_length==0:
		print 'FusionScan error	: Read length is not specified.'
		print usage
		sys.exit(1)
	if project_name=='':
		print 'FusionScan error	: out-dir-name should be specified.'
		print usage
		sys.exit(1)

	if min_seed==0:
		min_seed=2

        ########################################################################
        # preparation for the blat server
        ########################################################################

		
	os.system('mkdir -p '+dir_name+'/'+project_name	+ ' &')
	time.sleep(1)
	time0=datetime.datetime.now()
	timea='['+str(time0).split('.')[0]+']'
	
	print timea + ' 0. Checking gfServer'
	gfserv=	os.popen('ps ux	| grep \"'+host+' ' + port+'\"')
	for lineg in gfserv.xreadlines():
		if lineg:
			len_lineg=len(lineg.split('src_dir'))
			if len_lineg>=2:
				blat_flag='yes'
				break
			else:
				blat_flag='no'

	if blat_flag=='yes':
		time0=datetime.datetime.now()
		timea='['+str(time0).split('.')[0]+']'
		print timea + '	gfServer olready exists.'
	else:
		print '     making gfServer'
		print '     src_dir/gfServer start '+host + ' ' + port + ' -stepSize=5 -log=src_dir/untrans.log src_dir/hg19.2bit &'
		os.system('src_dir/gfServer start '+host + ' ' + port + ' -stepSize=5 -log=src_dir/untrans.log src_dir/hg19.2bit &')
		time.sleep(150)

	

	########################################################################
	# quality trimming
	########################################################################
	
	if fq_trim=='true':
		if len(flist)==1:
			time0=datetime.datetime.now()
			timea='['+str(time0).split('.')[0]+']'
			print timea + ' 1. Quality trimming'
			os.system(dir_name+'/src_dir/fastq_quality_trimmer '+score+' -t	20 -l 40 -i ' + flist[0] +' -o ' + dir_name+'/'+project_name + '/'+project_name + '_trimmed.fq')
			print '         '+dir_name+'/'+project_name + '/'+project_name + '_trimmed.fq'
			
		elif len(flist)==2:
			if thread_num>1:
				ThreadsLeft = 2
				lock = thread.allocate_lock()
		
				def threadexit(id):
					global ThreadsLeft
					#print ' thread	%d is quitting'	% id
					lock.acquire()
					ThreadsLeft-= 1
					lock.release()
		
				def trimming(id,fname):
					os.system(dir_name+'/src_dir/fastq_quality_trimmer '+score+' -t	20 -l 40 -i ' + fname+' -o ' + dir_name+'/'+project_name + '/'+project_name + '_'+str(id)+'_trimmed.fq')
					threadexit(id)
	
				for i in range(2):
					fname =	flist[i]
					thread.start_new_thread(trimming, (i+1,fname))
	
				while ThreadsLeft:
					time.sleep(0.1)
				os.system('cat '+dir_name+'/'+project_name + '/'+project_name + '_1_trimmed.fq ' + dir_name+'/'+project_name + '/'+project_name + '_2_trimmed.fq > ' + dir_name+'/'+project_name + '/'+project_name	+ '_trimmed.fq')
				time.sleep(1)
				os.system('rm '+dir_name+'/'+project_name + '/'+project_name + '_1_trimmed.fq ' + dir_name+'/'+project_name + '/'+project_name + '_2_trimmed.fq')
				time.sleep(1)
				print '                         '+dir_name+'/'+project_name + '/'+project_name+'_trimmed.fq'

			elif thread_num==1:
				time0=datetime.datetime.now()
				timea='['+str(time0).split('.')[0]+']'
				print timea + ' 1. Quality trimming'
				os.system(dir_name+'/src_dir/fastq_quality_trimmer '+score+' -t 20 -l 40 -i ' + flist[0]+' -o ' + dir_name+'/'+project_name + '/'+project_name + '_1_trimmed.fq')
				time.sleep(1)
				os.system(dir_name+'/src_dir/fastq_quality_trimmer '+score+' -t 20 -l 40 -i ' + flist[1]+' -o '+dir_name+'/'+project_name + '/'+project_name + '_2_trimmed.fq')
				time.sleep(1)
				os.system('cat '+dir_name+'/'+project_name + '/'+project_name +	'_1_trimmed.fq ' + dir_name+'/'+project_name + '/'+project_name	+ '_2_trimmed.fq > ' + dir_name+'/'+project_name + '/'+project_name	+ '_trimmed.fq')
				os.system('rm '+dir_name+'/'+project_name +	'/'+project_name + '_1_trimmed.fq ' + dir_name+'/'+project_name	+ '/'+project_name + '_2_trimmed.fq')
				print '                         '+dir_name+'/'+project_name + '/'+project_name+'_trimmed.fq'

	else:
		if len(flist)==1:
			os.system('cp '	+ flist[0] + ' '+dir_name+'/'+project_name + '/'+project_name+'_trimmed.fq')
			time.sleep(1)
			time0=datetime.datetime.now()
			timea='['+str(time0).split('.')[0]+']'
			print timea + ' 1. Skipped Trimmming Step. - It is a FASTA format file.'
			print '                         '+dir_name+'/'+project_name + '/'+project_name+'_trimmed.fq'
		if len(flist)==2:
			os.system('cat '+flist[0]+' ' +flist[1]+ ' > '+dir_name+'/'+project_name + '/'+project_name+'_trimmed.fq')
			time0=datetime.datetime.now()
			timea='['+str(time0).split('.')[0]+']'
			print timea + ' 1. Skipped Trimmming Step. - It is a FASTA format file.'
			print '                         '+dir_name+'/'+project_name + '/'+project_name+'_trimmed.fq'



        ########################################################################
        # getting unmapped reads for refMrna
        ########################################################################
	time11=datetime.datetime.now()
	timea='['+str(time11).split('.')[0]+']'
	print timea + ' 2. Selecting Unmapped Reads'
	bwf=''
	if file_format=='fastq':
		bwf=''
	elif file_format=='fasta':
		bwf='-f'
	print ' 2-1. unmapped reads on refMrna'
	os.system(dir_name+'/src_dir/bowtie2  '+ bwf +' -p ' + str(thread_num) + ' ' + dir_name+'/src_dir/bowtie2_index_refMrna/refMrna	' + dir_name+'/'+project_name +	'/'+project_name+'_trimmed.fq --un ' + dir_name+'/'+project_name + '/'+project_name + '_trimmed_refMrnaUn.fq > '+dir_name+'/'+project_name + '/temp.txt')
        print '                         '+dir_name+'/'+project_name + '/'+project_name+'_trimmed_refMrnaUn.fq'



	########################################################################
	# getting unmapped reads for hg19
	########################################################################
        print ' 2-2. unmapped reads on hg19'
	os.system(dir_name+'/src_dir/bowtie2 '+ bwf +' -p ' + str(thread_num) + ' ' + dir_name+'/src_dir/bowtie2_index_hg19/hg19 ' + dir_name+'/'+project_name + '/'+project_name+'_trimmed_refMrnaUn.fq --un ' + dir_name+'/'+project_name + '/'+project_name + '_trimmed_refMrnaUn_hg19Un.fq > '+dir_name+'/'+project_name + '/temp.txt')

	print '                         '+dir_name+'/'+project_name + '/'+project_name+'_trimmed_refMrnaUn_hg19Un.fq'
	os.system('rm '+dir_name+'/'+project_name + '/temp.txt')




        ########################################################################
        # artifacts removing
        ########################################################################
        time1=datetime.datetime.now()
        timea='['+str(time0).split('.')[0]+']'
        print timea + ' 3. Artifacts Filtering Step'
        os.system(dir_name+'/src_dir/fastx_artifacts_filter -i ' + dir_name+'/'+project_name + '/'+project_name+'_trimmed_refMrnaUn_hg19Un.fq' +' -o ' + dir_name+'/'+project_name + '/'+project_name +  '_trimmed_refMrnaUn_hg19Un_artifactRemoved.fq ' + score)
        print '                         '+dir_name+'/'+project_name + '/'+project_name + '_trimmed_refMrnaUn_hg19Un_artifactRemoved.fq'



	########################################################################
	# getting unique reads
	########################################################################
	time11=datetime.datetime.now()
	timea='['+str(time11).split('.')[0]+']'
	print timea + ' 4. Selecting Unique Reads'
	os.system(dir_name+'/src_dir/fastx_collapser '+score+' -i ' + dir_name+'/'+project_name	+ '/'+project_name + '_trimmed_refMrnaUn_hg19Un_artifactRemoved.fq -o '	+  dir_name+'/'+project_name + '/'+project_name	+ '_trimmed_refMrnaUn_hg19Un_artifactRemoved_collapsed.fa')
	print '                         '+dir_name+'/'+project_name + '/'+project_name + '_trimmed_refMrnaUn_hg19Un_artifactRemoved_collapsed.fa'
	
	
	if (thread_num>1):
		########################################################################
		# split read file per thread number
		########################################################################
		time11=datetime.datetime.now()
		timea='['+str(time11).split('.')[0]+']'
		print timea + ' 5. Split Input Fasta File by Thread Number'
		total_n= os.popen('tail	' + dir_name+'/'+project_name+'/'+project_name + '_trimmed_refMrnaUn_hg19Un_artifactRemoved_collapsed.fa -n 2')	
		total_num=0
		for linex in total_n.xreadlines():
			if linex[0]=='>':
				den = linex[1:].split('-')[0]
				total_num = int(den)
		
		per_num	= int(total_num/thread_num) + 1
		os.system('mkdir -p '+dir_name+'/'+project_name+'/fa_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/pslx_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/log_dir')
		os.system('python '+ dir_name+'/src_dir/split_fasta.py ' + dir_name+'/'+project_name + '/'+project_name + '_trimmed_refMrnaUn_hg19Un_artifactRemoved_collapsed.fa ' + project_name + '	' + str(thread_num) + '	' + str(per_num)+' > ' +dir_name+'/'+project_name + '/log_dir/ddd')
		fa_list=glob.glob(dir_name+'/'+project_name+'/fa_dir/'+project_name+'_[0-9].fa')
		print fa_list
		time5=datetime.datetime.now()
		os.system('rm '+dir_name+'/'+project_name +'/log_dir/ddd')
		print '                         '+dir_name+'/'+project_name+'/fa_dir/'

		
		########################################################################
		# align reads to hg19
		########################################################################
		
		time11=datetime.datetime.now()
		timea='['+str(time11).split('.')[0]+']'
		print timea + ' 6. Align to hg19'
		ThreadsLeft= thread_num
		lock = thread.allocate_lock()
	
		def threadexit(id):
			global ThreadsLeft
			#print '                     thread	%d is quitting'	% id
			lock.acquire()
			ThreadsLeft-= 1
			lock.release()
	
		def align_to_hg19(id,fname):
			os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'_'+str(id)+'/rescue_dir')
			os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'_'+str(id)+'/genefa_dir')
			os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'_'+str(id)+'/genefa_dir2')
			os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'_'+str(id)+'/readfa_dir')
			os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'_'+str(id)+'/info_dir')
			os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'_'+str(id)+'/alignment_dir')
			os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'_'+str(id)+'/temp_dir')
			os.system(dir_name+'/src_dir/ssaha2 -solexa -skip 6 -best 5 -output pslx -save ' + dir_name+'/src_dir/ssaha2_index_hg19/hg19 ' + dir_name+'/'+project_name+'/fa_dir/'+fname + ' > ' + dir_name+'/'+project_name+'/pslx_dir/'+fname.split('.')[0] + '.pslx.txt')
			threadexit(id)
	
                print '----------'
		print fa_list
		print thread_num
		print '----------'

		for i in range(thread_num):
			print '                         '+fa_list[i]
			fname =	fa_list[i].split('/')[-1]
			thread.start_new_thread(align_to_hg19, (i+1,fname))
	
		while ThreadsLeft: 
			time.sleep(0.1)
		print '                     SSAHA2 alignment done.'
		print '                         '+dir_name+'/'+project_name+'/pslx_dir/'

		
		os.system('mkdir -p ' +	dir_name+'/'+project_name+'/genefa_dir')
		os.system('mkdir -p ' +	dir_name+'/'+project_name+'/rescue_dir')
		os.system('mkdir -p ' +	dir_name+'/'+project_name+'/alignment_dir')
		os.system('mkdir -p ' +	dir_name+'/'+project_name+'/temp_dir')
		

		pslx_flist=glob.glob(dir_name+'/'+project_name+'/pslx_dir/'+project_name+'_*.pslx.txt')
		rescue_dir=dir_name+'/'+project_name+'/rescue_dir/'
		align_dir=dir_name+'/'+project_name+'/alignment_dir/'
		temp_dir=dir_name+'/'+project_name+'/temp_dir/'
		src_dir=dir_name+'/src_dir/'

		########################################################################
		# run main algorithm
		########################################################################		
		
		time11=datetime.datetime.now()
		timea='['+str(time11).split('.')[0]+']'
		print timea + ' 7. run FusionRNAscan program * thread_num'
		ThreadsLeft = thread_num
		lock = thread.allocate_lock()

		def threadexit2(id):
			global ThreadsLeft
			#print '                     thread	%d is quitting'	% id
			lock.acquire()
			ThreadsLeft -= 1
			lock.release()

		ofile3 = open(dir_name+'/'+project_name	+ '/log_dir/log_count_'+project_name,'w')

		def run_FusionRNAscan(id,fname):
			fname2=fname.split('_')[1]
			fnum=fname2.split('.')[0]
			print '                         java -jar ' + dir_name + '/src_dir/FusionScan.jar ' +dir_name + ' ' + fname + ' ' +  dir_name+'/'+project_name + '/'+project_name + '_trimmed_artifactRemoved.fq' + ' ' + project_name	+ ' ' +project_name + '_'+str(id)+' ' + str(read_length) + ' ' + host +	' ' + port + ' > ' +dir_name+'/'+project_name + '/log_dir/log_'+project_name + '_'+str(id)
			os.system('java -jar ' + dir_name + '/src_dir/FusionScan.jar ' +dir_name + ' ' + fname + ' ' +  dir_name+'/'+project_name + '/'+project_name + '_trimmed_artifactRemoved.fq' + ' ' + project_name  + ' ' +project_name + '_'+str(id)+' ' + str(read_length) + ' ' + host + ' ' + port + ' > ' +dir_name+'/'+project_name + '/log_dir/log_'+project_name + '_'+str(id))
			#os.system('java -cp ' + dir_name + ' FusionScan.Main ' +dir_name + ' ' + fname +' '  +  dir_name+'/'+project_name + '/'+project_name + '_trimmed_artifactRemoved.fq' + ' ' + project_name + ' ' +project_name + '_'+str(id)+' ' + str(read_length) + ' '+ host + ' ' + port  + ' > ' +dir_name+'/'+project_name + '/log_dir/log_'+project_name +'_'+str(id))
			threadexit2(id)

		for f in range(thread_num):
			fname =	pslx_flist[f]
			thread.start_new_thread(run_FusionRNAscan, (f+1,fname))

		while ThreadsLeft:
			time.sleep(0.1)
		
		rescue_dic={}
		cout=open(rescue_dir+'rescue_all.fa','w')
		catlist=glob.glob(dir_name+'/'+project_name+'/'+project_name+'_*/rescue_dir/rescue_all.fa')
		for cfile in catlist:
			cin=open(cfile,'r')
			ful=''
			while 1:
				cline=cin.readline()
				if not cline: 
					cout.write(ful)
					cout.flush()
					break
				if cline[0]=='>':
					cout.write(ful)
					cout.flush()
					gpct=0
					gp=cline[1:].split('-')[0]
					if rescue_dic.has_key(gp):
						gpct=rescue_dic[gp]
						gpct+=1
						rescue_dic[gp]=gpct
					else:
						rescue_dic[gp]=1
						gpct=1
					ful='>'+gp+'-'+str(gpct)+'\n'
					
				else:
					ful=ful+cline
			cin.close()
		cout.close()

		# rescue step to here.
		gp_chr_hash={}
		info_loc_hash={}
	
		for f in range(thread_num):
			os.system('cp ' +dir_name+'/'+project_name + '/'+project_name + '_'+str(f+1)+'/genefa_dir/* '+ dir_name+'/'+project_name+'/genefa_dir/')
			infoi=open(dir_name+'/'+project_name + '/'+project_name	+'_'+ str(f+1)+'/info_dir/info_'+project_name+'.txt','r')
			while 1:
				readi=infoi.readline()
				if not readi: break
				lines2=readi[:-1].split('\t')
				gp = lines2[0]
				hg = lines2[1]
				hchr = lines2[2]
				hstr = lines2[3]
				hbp= lines2[4]
				tg = lines2[5]
				tchr = lines2[6]
				tstr = lines2[7]
				tbp= lines2[8]
				BPseq =	lines2[9]
				info_loc_hash[gp]=hg + "\t"+ hchr + "\t" +hstr + "\t"+ hbp +"\t" + tg + "\t" + tchr+ "\t" + tstr +"\t" + tbp + "\t" +BPseq

		sequ = ''
		tt = ''
		tlen = 0
		writerg	= open(rescue_dir+'pseudo_list.txt','w')
		writergg = open(rescue_dir+'pseudo.fa','w')

		
		listofFiles=glob.glob(dir_name+'/'+project_name+'/genefa_dir/*.fa')
		fcount=0
		for h in range(len(listofFiles)):

			ff=listofFiles[h]
			gps2=ff.split('/')[-1]
			gps3=gps2.split('_')
			gps=gps2.split('.fa')
			fcount+=1
			ddd=ff[len(ff)-2:]

			if ddd=="fa" and len(gps3)==3:
				gp=gps[0]
				ingf=open(ff,'r')
				cc=0
				bplen=0
				seqq=''
				while 1:
					linegf=ingf.readline()
					if not linegf: break
					if cc==0:
						linegf
						linegfs=linegf.split('>')[1]
						linegfss=linegfs.split('_')
						bplen=int(linegfss[1])
						gp_chr_hash[gp]=linegfss[2]+'_'+linegfss[3]
						cc+=1
					else:
						seqq=linegf[:-1]
						sequ=sequ+seqq
						tt = tt	+ gp + '\t'+ str(tlen)	+ '\t' + str(bplen + tlen) + '\t' +str(len(seqq)+tlen)+ '\n'
						tlen=tlen+len(seqq)		
			
		writerg.write(tt)
		writerg.flush()
		writerg.close()
		writergg.write('>pseudo\n' + sequ+'\n')
		writergg.flush()
		writergg.close()
		
		os.system(src_dir+'ssaha2Build -solexa -save ' + rescue_dir+ 'pseudo '+ rescue_dir+'pseudo.fa')
		print '                         '+src_dir+'ssaha2Build -solexa -save ' + rescue_dir + 'pseudo '+ rescue_dir+'pseudo.fa'
		time.sleep(0.1)
		os.system(src_dir+'ssaha2 -solexa -best	5 -output pslx -save ' +rescue_dir+'pseudo ' + rescue_dir+'rescue_all.fa > ' + rescue_dir+'pseudo_1.pslx.txt')
		print '                         '+src_dir+'ssaha2 -solexa -best 5 -output pslx -save ' +rescue_dir+'pseudo ' + rescue_dir+'rescue_all.fa > ' + rescue_dir+'pseudo_1.pslx.txt'
		time.sleep(0.1)
		
		
		loc_hash={}
		gplist_=[]
		biplist=open(rescue_dir+'pseudo_list.txt','r')
		while 1:
			linee=biplist.readline()
			if not linee: break
			lines2=linee.split('\t')
			gp=lines2[0]
			if gp not in gplist_:
				gplist_.append(gp)

			loc_hash[lines2[0]]=lines2[1] +	'_'+ lines2[2]	+ '_' +	lines2[3]

		lockeys=loc_hash.keys()

		list_file_name = rescue_dir+ 'pseudo_list.txt'
		pslx_file_name = rescue_dir+ 'pseudo_1.pslx.txt'
		out_file2 =open(dir_name+'/'+project_name +'/'+project_name + '_FusionScan_candidate_result.txt','w')
	
		do_rescue_pr(pslx_file_name, align_dir,	loc_hash, lockeys)
		
		gp_ct_dict=make_align_image(sequ, dir_name,src_dir, temp_dir, align_dir, info_loc_hash, project_name, host, port, out_file2, loc_hash,	gplist_, 60)
		out_file2.close()
		
		ofile =	open(dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed1.txt','w')
		ofile.write("Hgene\tHchr\tHstrand\tHbpt\tTgene\tTchr\tTstrand\tTbp\tHg-Tg\t# of	seed read\tBPseq\n")
		flist=glob.glob(dir_name+'/'+project_name+'/*FusionScan_candidate_result.txt')
	
		ct_dict={}
		t_dict={}
		for f in range(len(flist)):
			fns= flist[f].split('/')
			fn=fns[len(fns)-1]
			#fn1 = fn.split('/')[1]
			fn2= fn.split('_FusionScan')[0]
			ifile =	open(flist[f],'r')
			while 1:
				linei =	ifile.readline()
				if not linei: break
		
				lines =	linei.split('\t')
				gp = lines[0]#+'_'+lines[4]
				if (lines[0]!='Hgene'):
					ct = int(lines[10])
					#os.system('mv '+dir_name +'/'+project_name+'/'+project_name+'_*/alignment_dir/' +gp + '*.RNAalign ' + dir_name+'/'+project_name+'/alignment_dir/')
					if gp_ct_dict.has_key(gp):
						ctt= gp_ct_dict[gp]
						lls=linei.split('\t')
						ttt=lls[0]+'\t'+lls[1]+'\t'+lls[2]+'\t'+lls[3]+'\t'+lls[4]+'\t'+lls[5]+'\t'+lls[6]+'\t'+lls[7]+'\t'+lls[8]+'\t'+lls[9]+'\t'+str(ctt)+'\t'+lls[11]
						t_dict[gp] = ttt
					else:
						t_dict[gp] = linei

			ifile.close()		

		gplist = gp_ct_dict.keys()
		ctlist = gp_ct_dict.values()
		ctlist2=[]
	
		for gp in gplist:
			ct = gp_ct_dict[gp]
			if ct_dict.has_key(ct):
				temp = ct_dict[ct]
				temp.append(gp)
				ct_dict[ct]= temp
			else:
				temp=[]
				temp.append(gp)
				ct_dict[ct]= temp
		
		for ct in ctlist:
			if ct not in ctlist2:
				ctlist2.append(ct)
		
		ctlist2.sort(reverse=True)
		for ct in ctlist2:
			if ct_dict.has_key(ct):
				temp = ct_dict[ct]
				for t in temp:
					if t_dict.has_key(t):
						ofile.write(t_dict[t])
		
		ofile.close()
	
		candi_n= os.popen('wc -l ' + dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed1.txt')
		candi_num=0
		for linec in candi_n.xreadlines():
			candi_num =int(linec.split(' ')[0])
		if candi_num!=0:
			candi_num-=1
		seed_criteria_c=0
		ofile5=open(dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed'+str(min_seed)+'.txt','w')
		resulti=open(dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed1.txt','r')
		ofile5.write("Hgene\tHchr\tHstrand\tHbpt\tTgene\tTchr\tTstrand\tTbp\tHg-Tg\t# of seed read\tBPseq\n")

		while 1:
			liner=resulti.readline()
			if not liner: break
			liners=liner.split('\t')
			if liners[0]!='Hgene':
				seedrc2=liners[10]
				if seedrc2[0]!='#':
					seedrc=int(seedrc2)
					if seedrc>=min_seed:
						seed_criteria_c+=1
						ofile5.write(liner)
		total=0
		filtered=0
		kalbok=0
		notsimilar=0
		notrepeat=0
		ifile3 = open(dir_name+'/'+project_name	+ '/log_dir/log_count_'+project_name,'r')
		while 1:
			lined=ifile3.readline()
			if not lined: break
			lines=lined.split('\t')
			total+=int(lines[0])
			filtered+=int(lines[1])
			kalbok+=int(lines[2])
			notsimilar+=int(lines[3])
			notrepeat+=int(lines[4])
	
		ofile4 = open(dir_name+'/'+project_name + '/count_per_each_step_'+project_name,'w')
		ofile4.write("")
		ofile4.write(" # total reads : "+str(total)+"\n")
		ofile4.write("--> # filtered reads by alignment information : "+str(filtered)+"\n")
		ofile4.write("    --> # fusion genes having read aligned exon boundary totally : "+str(kalbok)+"\n")
		ofile4.write("        --> # fusion genes having no sequence homology : "+str(notsimilar)+"\n")
		ofile4.write("            --> # fusion genes not aligned repeat regions	: "+str(notrepeat)+"\n")
		ofile4.write("                --> # fusion genes having at least 1 seed : "+ str(candi_num)+"\n")
		ofile4.write("                   --> # fusion genes having at least "+str(min_seed)+" seeds : "+ str(seed_criteria_c)+"\n")
		print '-----------------------------------------------------------------------------------------------------'
		print "#<FusionScan Result for "+project_name+">"
		print "## total reads : "+str(total)
		print "#--> # filtered reads by alignment information : "+str(filtered)
		print "#    --> # fusion genes having read aligned exon boundary totally : "+str(kalbok)
		print "#        --> # fusion genes having no sequence homology : "+str(notsimilar)
		print "#            --> # fusion genes not aligned repeat regions : "+str(notrepeat)
		print "#                --> # fusion genes having at least 1 seed : "+ str(candi_num)
		print "#                   --> # fusion genes having at least "+str(min_seed)+" seeds : "+ str(seed_criteria_c)
		time11=datetime.datetime.now()
		timea='['+str(time11).split('.')[0]+']'
		print timea + '                             '+dir_name+'/'+project_name + '/count_per_each_step_'+project_name

	elif (thread_num==1):
		print '==================================================================================================='
		
		os.system('mkdir -p '+dir_name+'/'+project_name)
		os.system('mkdir -p '+dir_name+'/'+project_name+'/fa_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/pslx_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/log_dir')
                os.system('mkdir -p ' + dir_name+'/'+project_name+'/genefa_dir')
                os.system('mkdir -p ' + dir_name+'/'+project_name+'/rescue_dir')
                os.system('mkdir -p ' + dir_name+'/'+project_name+'/alignment_dir')
                os.system('mkdir -p ' + dir_name+'/'+project_name+'/temp_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'/rescue_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'/genefa_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'/genefa_dir2')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'/readfa_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'/info_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'/alignment_dir')
		os.system('mkdir -p '+dir_name+'/'+project_name+'/'+project_name+'/temp_dir')
		os.system('mv ' +  dir_name+'/'+project_name + '/'+project_name + '_trimmed_refMrnaUn_hg19Un_artifactRemoved_collapsed.fa '+ dir_name+'/'+project_name + '/fa_dir/'+project_name+'.fa')
		time11=datetime.datetime.now()
		timea='['+str(time11).split('.')[0]+']'
		print timea + ' 5. Align to hg19'
		os.system(dir_name+'/src_dir/ssaha2 -solexa -skip 6 -cmatch 20 -best 5 -output pslx -save ' +dir_name+'/src_dir/ssaha2_index_hg19/hg19 ' + dir_name+'/'+project_name + '/fa_dir/'+project_name+'.fa > ' + dir_name+'/'+project_name + '/pslx_dir/'+project_name+'.pslx.txt')
		print '                         '+dir_name+'/'+project_name + '/pslx_dir/'+project_name+'.pslx.txt'

		time11=datetime.datetime.now()
		timea='['+str(time11).split('.')[0]+']'
		print timea + ' 6. Run FusionRNAscan program'
		ofile3 = open(dir_name+'/'+project_name	+ '/log_dir/log_count_'+project_name,'w')
		print '                         java -jar ' + dir_name + '/src_dir/FusionScan.jar '+ dir_name + ' '+ dir_name+'/'+project_name + '/pslx_dir/'+project_name+'.pslx.txt ' +dir_name+'/'+project_name + '/'+project_name + '_trimmed_artifactRemoved.fq' + ' ' +project_name + ' ' +project_name + ' ' + str(read_length) + ' ' + host + ' ' + port + ' > ' +dir_name+'/'+project_name + '/log_dir/log_'+project_name
		os.system('java -jar ' + dir_name + '/src_dir/FusionScan.jar '+ dir_name + ' '+ dir_name+'/'+project_name + '/pslx_dir/'+project_name+'.pslx.txt ' +dir_name+'/'+project_name + '/'+project_name + '_trimmed_artifactRemoved.fq' + ' ' +project_name + ' ' +project_name + ' ' + str(read_length) + ' ' + host + ' ' + port + ' > ' +dir_name+'/'+project_name + '/log_dir/log_'+project_name)
		#os.system('java -cp ' + dir_name + ' FusionScan.Main '+ dir_name + ' '+ dir_name+'/'+project_name + '/pslx_dir/'+project_name+'.pslx.txt ' +dir_name+'/'+project_name + '/'+project_name + '_trimmed_artifactRemoved.fq' + ' ' +project_name + ' ' +project_name + ' ' + str(read_length) + ' ' + host + ' ' + port + ' > ' +dir_name+'/'+project_name + '/log_dir/log_'+project_name)
		time.sleep(0.1)
		os.system('mkdir -p ' +	dir_name+'/'+project_name+'/alignment_dir')
		
                rescue_dir=dir_name+'/'+project_name+'/rescue_dir/'
                align_dir=dir_name+'/'+project_name+'/alignment_dir/'
                temp_dir=dir_name+'/'+project_name+'/temp_dir/'
                src_dir=dir_name+'/src_dir/'

		rescue_dic={}
		cout=open(dir_name+'/'+project_name+'/rescue_dir/rescue_all.fa','w')
		catlist=glob.glob(dir_name+'/'+project_name+'/'+project_name+'/rescue_dir/rescue_all.fa')
		for cfile in catlist:
			cin=open(cfile,'r')
			ful=''
			while 1:
				cline=cin.readline()
				if not cline: 
					cout.write(ful)
					cout.flush()
					break
				if cline[0]=='>':
					cout.write(ful)
					cout.flush()
					gpct=0
					gp=cline[1:].split('-')[0]
					if rescue_dic.has_key(gp):
						gpct=rescue_dic[gp]
						gpct+=1
						rescue_dic[gp]=gpct
					else:
						rescue_dic[gp]=1
						gpct=1
					ful='>'+gp+'-'+str(gpct)+'\n'
					
				else:
					ful=ful+cline
			cin.close()
		cout.close()

		# rescue step to here.
		gp_chr_hash={}
		info_loc_hash={}
	
		os.system('cp '	+dir_name+'/'+project_name + '/'+project_name +'/genefa_dir/* '+ dir_name+'/'+project_name+'/genefa_dir/')
		infoi=open(dir_name+'/'+project_name + '/'+project_name+'/info_dir/info_'+project_name+'.txt','r')
		while 1:
			readi=infoi.readline()
			if not readi: break
			lines2=readi[:-1].split('\t')
			gp = lines2[0]
			hg = lines2[1]
			hchr = lines2[2]
			hstr = lines2[3]
			hbp= lines2[4]
			tg = lines2[5]
			tchr = lines2[6]
			tstr = lines2[7]
			tbp= lines2[8]
			BPseq =	lines2[9]
			info_loc_hash[gp]=hg + "\t"+ hchr + "\t" +hstr + "\t"+ hbp +"\t" + tg + "\t" + tchr+ "\t" + tstr +"\t" + tbp + "\t" +BPseq

		sequ = ''
		tt = ''
		tlen = 0
		writerg	= open(rescue_dir+'pseudo_list.txt','w')
		writergg = open(rescue_dir+'pseudo.fa','w')

		print info_loc_hash
		
		listofFiles=glob.glob(dir_name+'/'+project_name+'/genefa_dir/*.fa')
		fcount=0
		for h in range(len(listofFiles)):

			ff=listofFiles[h]
			gps2=ff.split('/')[-1]
			gps3=gps2.split('_')
			gps=gps2.split('.fa')
			fcount+=1
			ddd=ff[len(ff)-2:]

			if ddd=="fa" and len(gps3)==3:
				gp=gps[0]
				ingf=open(ff,'r')
				cc=0
				bplen=0
				seqq=''
				while 1:
					linegf=ingf.readline()
					if not linegf: break
					if cc==0:
						linegf
						linegfs=linegf.split('>')[1]
						linegfss=linegfs.split('_')
						bplen=int(linegfss[1])
						gp_chr_hash[gp]=linegfss[2]+'_'+linegfss[3]
						cc+=1
					else:
						seqq=linegf[:-1]
						sequ=sequ+seqq
						tt = tt	+ gp + '\t'+ str(tlen) + '\t' + str(bplen + tlen) + '\t' +str(len(seqq)+tlen)+ '\n'
						tlen=tlen+len(seqq)
			
		writerg.write(tt)
		writerg.flush()
		writerg.close()
		writergg.write('>pseudo\n' + sequ+'\n')
		writergg.flush()
		writergg.close()
		
		os.system(src_dir+'ssaha2Build -solexa -save ' + rescue_dir+ 'pseudo '+ rescue_dir+'pseudo.fa')
		print '                         '+src_dir+'ssaha2Build -solexa -save ' + rescue_dir + 'pseudo '+ rescue_dir+'pseudo.fa'
		time.sleep(0.1)
		os.system(src_dir+'ssaha2 -solexa -best 5 -output pslx -save ' +rescue_dir+'pseudo ' + rescue_dir+'rescue_all.fa > ' + rescue_dir+'pseudo_1.pslx.txt')
		print '                         '+src_dir+'ssaha2 -solexa -best 5 -output pslx -save ' +rescue_dir+'pseudo ' + rescue_dir+'rescue_all.fa > ' + rescue_dir+'pseudo_1.pslx.txt'
		time.sleep(0.1)
		
		
		loc_hash={}
		gplist_=[]
		biplist=open(rescue_dir+'pseudo_list.txt','r')
		while 1:
			linee=biplist.readline()
			if not linee: break
			lines2=linee.split('\t')
			gp=lines2[0]
			if gp not in gplist_:
				gplist_.append(gp)

			loc_hash[lines2[0]]=lines2[1] + '_'+ lines2[2] + '_' + lines2[3]

		lockeys=loc_hash.keys()

		list_file_name = rescue_dir+ 'pseudo_list.txt'
		pslx_file_name = rescue_dir+ 'pseudo_1.pslx.txt'
		out_file2 =open(dir_name+'/'+project_name +'/'+project_name + '_FusionScan_candidate_result.txt','w')
	
		do_rescue_pr(pslx_file_name, align_dir, loc_hash, lockeys)
		
		gp_ct_dict=make_align_image(sequ, dir_name,src_dir, temp_dir, align_dir, info_loc_hash, project_name, host, port, out_file2, loc_hash, gplist_, 60)
		out_file2.close()
		
		ofile = open(dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed1.txt','w')
		ofile.write("Hgene\tHchr\tHstrand\tHbpt\tTgene\tTchr\tTstrand\tTbp\tHg-Tg\t# of	seed read\tBPseq\n")
		flist=glob.glob(dir_name+'/'+project_name+'/*FusionScan_candidate_result.txt')
	
		ct_dict={}
		t_dict={}
		for f in range(len(flist)):
			fns= flist[f].split('/')
			fn=fns[len(fns)-1]
			#fn1 = fn.split('/')[1]
			fn2= fn.split('_FusionScan')[0]
			ifile =	open(flist[f],'r')
			while 1:
				linei =	ifile.readline()
				if not linei: break
		
				lines =	linei.split('\t')
				gp = lines[0]#+'_'+lines[4]
				if (lines[0]!='Hgene'):
					ct = int(lines[10])
					#os.system('mv '+dir_name +'/'+project_name+'/'+project_name+'_*/alignment_dir/' +gp + '*.RNAalign ' + dir_name+'/'+project_name+'/alignment_dir/')
					if gp_ct_dict.has_key(gp):
						ctt= gp_ct_dict[gp]
						lls=linei.split('\t')
						ttt=lls[0]+'\t'+lls[1]+'\t'+lls[2]+'\t'+lls[3]+'\t'+lls[4]+'\t'+lls[5]+'\t'+lls[6]+'\t'+lls[7]+'\t'+lls[8]+'\t'+lls[9]+'\t'+str(ctt)+'\t'+lls[11]
						t_dict[gp] = ttt
					else:
						t_dict[gp] = linei

			ifile.close()

		gplist = gp_ct_dict.keys()
		ctlist = gp_ct_dict.values()
		ctlist2=[]
	
		for gp in gplist:
			ct = gp_ct_dict[gp]
			if ct_dict.has_key(ct):
				temp = ct_dict[ct]
				temp.append(gp)
				ct_dict[ct]= temp
			else:
				temp=[]
				temp.append(gp)
				ct_dict[ct]= temp
		
		for ct in ctlist:
			if ct not in ctlist2:
				ctlist2.append(ct)
		
		ctlist2.sort(reverse=True)
		for ct in ctlist2:
			if ct_dict.has_key(ct):
				temp = ct_dict[ct]
				for t in temp:
					if t_dict.has_key(t):
						ofile.write(t_dict[t])
		
		ofile.close()
	
		candi_n= os.popen('wc -l ' + dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed1.txt')
		candi_num=0
		for linec in candi_n.xreadlines():
			candi_num =int(linec.split(' ')[0])
		if candi_num!=0:
			candi_num-=1
		seed_criteria_c=0
		ofile5=open(dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed'+str(min_seed)+'.txt','w')
		resulti=open(dir_name+'/'+project_name+'/all_candidates_for_'+project_name+'_morethan_seed1.txt','r')
		ofile5.write("Hgene\tHchr\tHstrand\tHbpt\tTgene\tTchr\tTstrand\tTbp\tHg-Tg\t# of seed read\tBPseq\n")

		while 1:
			liner=resulti.readline()
			if not liner: break
			liners=liner.split('\t')
			if liners[0]!='Hgene':
				seedrc2=liners[10]
				if seedrc2[0]!='#':
					seedrc=int(seedrc2)
					if seedrc>=min_seed:
						seed_criteria_c+=1
						ofile5.write(liner)
		total=0
		filtered=0
		kalbok=0
		notsimilar=0
		notrepeat=0
		ifile3 = open(dir_name+'/'+project_name + '/log_dir/log_count_'+project_name,'r')
		while 1:
			lined=ifile3.readline()
			if not lined: break
			lines=lined.split('\t')
			total+=int(lines[0])
			filtered+=int(lines[1])
			kalbok+=int(lines[2])
			notsimilar+=int(lines[3])
			notrepeat+=int(lines[4])
	
		ofile4 = open(dir_name+'/'+project_name	+ '/count_per_each_step_'+project_name,'w')
		ofile4.write("")
		ofile4.write(" # total reads : "+str(total)+"\n")
		ofile4.write("--> # filtered reads by alignment information : "+str(filtered)+"\n")
		ofile4.write("    --> # fusion genes having read aligned exon boundary totally : "+str(kalbok)+"\n")
		ofile4.write("        --> # fusion genes having no sequence homology : "+str(notsimilar)+"\n")
		ofile4.write("            --> # fusion genes not aligned repeat regions : "+str(notrepeat)+"\n")
		ofile4.write("                --> # fusion genes having at least 1 seed : "+ str(candi_num)+"\n")
		ofile4.write("                   --> # fusion genes having at least "+str(min_seed)+" seeds : "+ str(seed_criteria_c)+"\n")
		print '---------------------------------------------------------------------------------------------------------------------'
		print "#<FusionScan Result for "+project_name+">"
		print "## total reads : "+str(total)
		print "#--> # filtered reads by alignment information : "+str(filtered)
		print "#    --> # fusion genes having read aligned exon boundary totally : "+str(kalbok)
		print "#        --> # fusion genes having no sequence homology : "+str(notsimilar)
		print "#            --> # fusion genes not aligned repeat regions : "+str(notrepeat)
		print "#                --> # fusion genes having at least 1 seed : "+ str(candi_num)
		print "#                   --> # fusion genes having at least "+str(min_seed)+" seeds : "+ str(seed_criteria_c)
		time11=datetime.datetime.now()
		timea='['+str(time11).split('.')[0]+']'
		print timea + '                             '+dir_name+'/'+project_name + '/count_per_each_step_'+project_name

	if known_fusion_search=='go':
		
		print '======================================================================================================'
		print 'GO to search known kinase related fusion	genes !'
		os.system('python find_known_novel_partner_of_fusion_kinase.py ' + dir_name+'/'+project_name + '/pslx_dir/')


else:
		print usage
		sys.exit(1)

