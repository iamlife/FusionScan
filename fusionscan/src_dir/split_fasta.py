import sys, os, string
import re
iname = sys.argv[1]
project_name = sys.argv[2]
thread_num = int(sys.argv[3])
per_num = int(sys.argv[4])
inf =open(iname,'r')

FullSeqlist=[]
FullSeq =''

c=0
while 1:
	line = inf.readline()
	if not line: 
		FullSeqlist.append(FullSeq)
		break

	if line[0] == '>':
		FullSeqlist.append(FullSeq)
		FullSeq = line
		
	if line[0] !='>':
		FullSeq = FullSeq + line

inf.close()

pn = project_name.split('_')[0]
for a in range(thread_num): # file number count
        ouf = open(project_name+'/fa_dir/'+pn+'_'+str(a+1)+'.fa','w')
        for i in range(a*per_num,(a+1)*per_num): # seq count per file
                ouf.write(FullSeqlist[i])




