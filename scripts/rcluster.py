import os, sys, getopt, re;
import time,Numeric,dircache,Bio.Cluster;

displaymode=True

###############################################################################

def makeHCLUSTERfile(matdist,method,nbclusters,cutheight):
    
    Rfilename="temp.R";  
    datafilename="sample_distances.dat";  
    fulldatafilename="%s/%s" % (os.getcwd(),datafilename);
    
    f=open(Rfilename,'w');
    fdata=open(datafilename,'w');
    
    ## header
    f.write("library(cluster)\n");
    
    ## build row vector (store in a file to avoid overload)
    dima=len(matdist);
    currentsize=0;
    for irow in range(dima):
        for icol in range(dima):
            currentsize+=1;
            fdata.write(str(matdist[irow][icol]) + ' ');
        fdata.write("\n");
    fdata.write("\n");
    
    ## build matrix
    tablename="tabledist";
    distmatrixname="matdist";
    buffer = "%s <- read.table(\"%s\")\n" % (tablename,fulldatafilename);
    f.write(buffer);
    buffer = "%s <- as.dist(%s, diag = FALSE, upper = FALSE)\n" % (distmatrixname,tablename) 
    f.write(buffer);
    
    ## cluster
    treename = "mytree";
    clustername = "myclusters";
    buffer = "%s <- hclust(%s,\"%s\")\n" % (treename,distmatrixname,method);
    f.write(buffer);
    if method=='ward':
        buffer = "%s <- cutree(%s,h=%d)\n" % (clustername,treename,cutheight);
    else:
        buffer = "%s <- cutree(%s,k=%d)\n" % (clustername,treename,nbclusters);
    f.write(buffer);
    
    f.write("print(\"clusters\")\n");
    buffer="%s\n" % (clustername);
    f.write(buffer);
        
    f.close();
    fdata.close();
    
    return (Rfilename,datafilename);

###############################################################################

def runHCLUSTERfile(rfilename):
    
    results = {};
    
    command="R -f %s --slave" % (rfilename);
    flux = os.popen(command);
    buffer = flux.readlines();
    flux.close();

    # read medoid ids
    if not buffer[0].startswith("[1] \"clusters\""):
        print "WARNING: Corrupted R output. Missing cluster info.";
        sys.exit(1);
    
    clustertable={}
    for k in range(1,len(buffer)):
        line=buffer[k].lstrip('\ \t\n').rstrip('\ \t\n')
        if k%2==1:
            fields=line.split(' ')
            ids=[]
            for fc in fields:
                if fc!='':
                    ids.append(fc)
        else:
            fields=line.split(' ')
            iclust=[]
            for fc in fields:
                if fc!='':
                    iclust.append(int(fc))
            if len(iclust)!=len(ids):
                print ids,iclust
                print 'ABORT'
                sys.exit(1)
            for i in range(len(ids)):
                clustertable[ids[i]]=int(iclust[i])

    # build cluster
    for rid,number in clustertable.iteritems():
        id=int(rid[1:])
        if results.has_key(number):
            results[number]['content'].append(id)
        else:   
            results[number]={}
            results[number]['medoid']=''
            results[number]['content']=[id]

    # return cluster infos
    return results;

###############################################################################

def computeCentroid(listofmembers,matdist):
    centroid=''
    mindist=10000000
    totaldist={}
    for elem1 in listofmembers:
        totaldist[elem1]=0
        for elem2 in listofmembers:
            if elem1!=elem2:
                totaldist[elem1]+=matdist[elem1,elem2]
        if totaldist[elem1]<=mindist:
            centroid=elem1

    return centroid

###############################################################################

def makeKMEANSfile(matdist,nbclusters,displaymode=False):
    
    Rfilename="temp.R";  
    datafilename="sample_distances.dat";  
    fulldatafilename="%s/%s" % (os.getcwd(),datafilename);
    
    f=open(Rfilename,'w');
    fdata=open(datafilename,'w');
    
    ## header
    f.write("library(cluster)\n");
    
    if displaymode:
        addonfilename="arrayPlot.R";
        path1="%s" % (addonfilename);
        path2="./scripts/%s" % (addonfilename);
        path3="../%s" % (addonfilename);
        foundplot=False;
        if os.path.exists(path1):
            fullpath=os.path.abspath(path1);
            foundplot=True;
        elif os.path.exists(path2):
            fullpath=os.path.abspath(path2);
            foundplot=True;
        elif os.path.exists(path3):
            fullpath=os.path.abspath(path3);
            foundplot=True;
        else:
            print "WARNING: Cannot find \"%s\"" % (addonfilename);
        if foundplot:
            buffer="source(\"%s\")\n" % (fullpath);
            f.write(buffer);
    
    ## build row vector
    dima=len(matdist);
    currentsize=0;
    for irow in range(dima):
        for icol in range(dima):
            currentsize+=1;
            fdata.write(str(matdist[irow][icol]) + ' ');
        fdata.write("\n");
    fdata.write("\n");
    
    ## build matrix
    distmatrixname="matdist";
    buffer = "%s <- read.table(\"%s\")\n" % (distmatrixname,fulldatafilename);
    f.write(buffer);
    
    ## cluster
    clustername = "myclusters";
    buffer = "%s <- pam(%s,%d)\n" % (clustername,distmatrixname,nbclusters);
    f.write(buffer);
    
    if displaymode and foundplot:
        f.write("arrayPlot(matdist)\n");
    
    f.write("print(\"medoids\")\n");
    buffer="%s$id.med\n" % (clustername);
    f.write(buffer);
    
    f.write("print(\"infos\")\n");
    buffer="%s$clusinfo\n" % (clustername);
    f.write(buffer);
    
    f.write("print(\"clustering\")\n");
    buffer="%s$clustering\n" % (clustername);
    f.write(buffer);
    
    f.close();
    fdata.close();
    
    return (Rfilename,datafilename);

###############################################################################

def runKMEANSfile(rfilename):
    
    results = {};
    
    command="R -f %s --slave" % (rfilename);
    flux = os.popen(command);
    buffer = flux.readlines();
    flux.close();
    
    currentline = buffer[0];
    cmptline=0;
    
    # read medoid ids
    if currentline.startswith("[1] \"medoids\""):
        cmptline+=1;
        currentline = buffer[cmptline];
        nclus=1;
        while not currentline.startswith("[1] \""):
            starti=currentline.find('[');
            decomp=currentline[starti:].replace('\n','').split(' ');
            # skip first arg is non relevant
            for k in range(1,len(decomp)):
                if decomp[k] != '':
                    results[nclus] = {};
                    results[nclus]['medoid'] = int(decomp[k]);
                    nclus+=1;
            cmptline+=1;
            currentline = buffer[cmptline];
    else:
        print "WARNING: Corrupted R output. Missing medoid ids.";
        sys.exit(1);
    
    ## read cluster infos
    if currentline.startswith("[1] \"infos\""):
        cmptline+=2; # skip header
        currentline = buffer[cmptline];
        nofclusters=0;
        
        while not currentline.startswith("[1] \"") and cmptline < len(buffer):
            nofclusters+=1;
            starti=currentline.find('[');
            closei=currentline.find(']');
            nclus=int(currentline[starti+1:closei-1]);
            if nclus!=nofclusters:
                print "WARNING: Corrupted cluster numbering (At %d while %d expected)." % (nclus,nofclusters);
                sys.exit(1);
            decomp=currentline[starti:].split(' ');
            kinfo=0;
            for k in range(1,len(decomp)):
                if decomp[k] != '':
                    if kinfo==0:
                        results[nclus]['size']=int(decomp[k]);
                    elif kinfo==1:
                        results[nclus]['max_diss']=float(decomp[k]);
                    elif kinfo==2:
                        results[nclus]['av_diss']=float(decomp[k]);
                    elif kinfo==3:
                        results[nclus]['diameter']=float(decomp[k]);
                    elif kinfo==4:
                        results[nclus]['separation']=float(decomp[k]);
                    else:
                        print "WARNING: corrupted line\"%s\"" % (currentline.replace('\n',''));
                        sys.exit(1);
                    kinfo+=1;
            # read done check and goto the next line
            if kinfo != 5:
                print "WARNING: corrupted line\"%s\"" % (currentline.replace('\n',''));
                sys.exit(1);
            cmptline+=1;
            currentline = buffer[cmptline];
        # end while
        
        if nofclusters != len(results):
            print "WARNING: Missing clusters (%d found while %d expected)." % (nofclusters,len(results));
            sys.exit(1);
    else:
        print "WARNING: Corrupted R output. Missing clusters info.";
        sys.exit(1);
    
    # read clustering
    if currentline.startswith("[1] \"clustering\""):
        cmptline+=1;
        currentline = buffer[cmptline];
        nsample=0;
        while not currentline.startswith("[1] \"") and cmptline < len(buffer):
            starti=currentline.find('[');
            decomp=currentline[starti:].split(' ');
            # skip first arg is non relevant
            for k in range(1,len(decomp)):
                if decomp[k] != '':
                    nsample+=1;
                    nclus=int(decomp[k]);
                    if not results[nclus].has_key('content'):
                        results[nclus]['content'] = [];
                    results[nclus]['content'].append(nsample);
                    nclus+=1;
            cmptline+=1;
            if cmptline < len(buffer):
                currentline = buffer[cmptline];
            else:
                currentline="";
    else:
        print "WARNING: Corrupted R output. Missing cluster content.";
        sys.exit(1);
    
    # return cluster infos
    return results;

