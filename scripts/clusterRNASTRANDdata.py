import sys,os,getopt,re
import rcluster,Numeric,dircache
from rcluster import *
from RNAstats import *


## read RNAstrand file


print "Read dataset"

f=open(sys.argv[1])
content=f.readlines()
f.close()

data={}
index2name=[]
name2index={}
count=0
for line in content:
    fields=line.replace('\n','').split(' ')
    data[fields[0]]=fields[1]
    index2name.append(fields[0])
    name2index[fields[0]]=count
    count+=1

## compute distance matrix

print "Compute distance matrix"

size_array = len(data)
dist_array = Numeric.zeros((size_array,size_array),Numeric.Float)

RNAdistanceinfile='rnadistance.in'
f=open(RNAdistanceinfile,'w')
for id in index2name:
    ssa=data[id]
    f.write(ssa+'\n')
f.close()

command="RNAdistance -Xm -Df < %s" % (RNAdistanceinfile)
flux = os.popen(command)
buffer = flux.readlines()
flux.close()

maxdist=0
counter=0
for index1 in range(1,len(buffer)-1,1):
    i1=index1
    row=buffer[index1].rstrip('\n').rstrip().lstrip()
    fields=row.split(' ')
    for index2 in range(len(fields)):
        i2=index2
        dist=float(fields[index2])
        dist_array[i1,i2]=dist
        dist_array[i2,i1]=dist
        counter+=1
        if dist>maxdist:    
            maxdist=dist
        if False and dist==0.0 and i1!=i2:
            print i1,' ',i2
            print data[index2name[i1]],index2name[i1]
            print data[index2name[i2]],index2name[i2]

print 'max distance between two targets: ',maxdist
print 'size: ',(size_array*(size_array-1.0))/2.0,';total assigned: ',counter

method='single'

print "Generate R script";
numberofclusters=100
cutheight=10
if method=='kmeans':
    (Rfilename,datafilename) = makeKMEANSfile(dist_array,numberofclusters);
elif method=='ward' or method=='average' or method=='centroid' or method=='single' or method=='median' or method=='complete':
    (Rfilename,datafilename) = makeHCLUSTERfile(dist_array,method,numberofclusters,cutheight);
    
print "Run R script";
if method=='kmeans':
    Rresults = runKMEANSfile(Rfilename);
elif method=='ward' or method=='average' or method=='centroid' or method=='single' or method=='median' or method=='complete':
    Rresults = runHCLUSTERfile(Rfilename);
print "Extract medoids";

medioddata={}
maxdiameter=0
minseparation=maxdist
maxclustsize=0

for iclust,clustdata in Rresults.iteritems():
    
    if method=='ward' or method=='average' or method=='centroid' or method=='single' or method=='median' or method=='complete':
        clustdata['medoid']=computeCentroid(clustdata['content'],dist_array)
    
    medoidname=index2name[clustdata['medoid']-1]
    if method=='kmeans':
        medioddiameter=clustdata['diameter']
        separation=clustdata['separation']
        clustersize=clustdata['size']
        print iclust,medoidname,medioddiameter,separation,clustersize
        if medioddiameter>maxdiameter:
            maxdiameter=medioddiameter
        if separation<minseparation:
            minseparation=separation
    else:
        clustersize=len(clustdata['content'])
        print iclust,medoidname,clustersize
    medioddata[iclust]=(medoidname,data[medoidname])
    if clustersize>maxclustsize:
        maxclustsize=clustersize

if method=='kmeans':
    print "Max diameter: ",maxdiameter
    print "Min separation: ",minseparation
print "Max cluster size: ",maxclustsize

pair=(0,0)
avesumdistc=0.0
mindistc=maxdist
for iclust1,clustdata1 in Rresults.iteritems():
    sumdistc=0.0
    for iclust2,clustdata2 in Rresults.iteritems():
        id1=clustdata1['medoid']-1
        id2=clustdata2['medoid']-1
        if iclust1!=iclust2:
            mydist=dist_array[id1,id2]
            sumdistc+=mydist
            if mydist<mindistc:
                mindistc=mydist
                pair=(id1,id2)
            if mydist==0:
                print 'SAME? ',id1,' ',id2
                print data[index2name[id1]],index2name[id1]
                print data[index2name[id2]],index2name[id2]
    avesumdistc+=(sumdistc)/(len(Rresults)-1)
print "Min distance between two medoids: ",mindistc,'; ',pair[0],'(',index2name[pair[0]],') ',pair[1],'(',index2name[pair[1]],')'
print "Average distance between medoids: ",avesumdistc/len(Rresults)
#print data[pair[0]]
#print data[pair[1]]

## write dataset

print "Build datafile"
outfile='rnastrand_dataset_filtered.txt'
f=open(outfile,'w')
for idc,content in sorted(medioddata.iteritems(),key=lambda (k,v): (len(v[1]),k)):
    f.write(content[1]+'\n')
f.close()
print 'Results stored in ',outfile













