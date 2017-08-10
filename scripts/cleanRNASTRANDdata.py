import sys,os,getopt,re
import rcluster,Numeric,dircache
from rcluster import *
from RNAstats import *



## read RNAstrand file


print "Read dataset"

f=open(sys.argv[1])
content=f.readlines()
f.close()

id_re = re.compile("# File (?P<id>.*)\n")
ssa_re = re.compile("(?P<ssa>[\(\)\.]+)\s*\n")

data={}

myid='none'
for line in content:
    cell=id_re.match(line)
    if cell:
        myid=cell.group('id').replace(' ','').replace('.dp','')
        data[myid]=''
    cell=ssa_re.match(line)
    if cell:
        data[myid]+=cell.group('ssa').replace('\n','')

## remove bad structures

key2rm=[]
for id,ssa in data.iteritems():
    contiguity,number_isolate_bp,number_of_bps = ss_stats(ssa)
    min_hairpin_length = hairpincheck(ssa)
    if min_hairpin_length<3 or number_isolate_bp>=1:
        #print 'remove',id,ssa,min_hairpin_length,number_isolate_bp
        key2rm.append(id)

print 'remove ',len(key2rm),' structures'
for id in key2rm:
    del data[id]

outfile= 'clean_' + sys.argv[1]
f=open(outfile,'w')
for idc,content in sorted(data.iteritems(),key=lambda (k,v): (len(v),k)):
    f.write(idc + " " + content+'\n')
f.close()













