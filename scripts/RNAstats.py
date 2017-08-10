import os,sys;
from  math import log;

###############################################################################

# count different loop types

def loop_counter(structure):

    stack=[];
    prev_openk  = -1;
    prev_closek = -1;
    sse_counter = {'hairpin':0, 'stack':0, 'bulge':0, 'internal': 0, 'multi': 0, 'root':0, 'bp':0, 'freein':0, 'freeout':0};
    level_open = 0;
    bp_counter=0

    lindex = range(len(structure));
    for k in lindex:
        carac = structure[k];
        if carac=='(':
            sse_counter['bp']+=1
            if len(stack)==0:
                sse_counter['root']+=1;

            if prev_openk>=0 and prev_closek>=0:
                level_open = len(stack);

            stack.append(k);
            prev_openk  = -1;
            prev_closek = -1;
        elif carac==')':
            closek=k;
            openk=stack.pop();

            if prev_openk<0 and prev_closek<0:
                sse_counter['hairpin']+=1;
            #sse_counter['freein']-=closek-openk-1
            else:
                leftgap  = prev_openk-openk-1;
                rightgap = closek-prev_closek-1;

                if leftgap==0 and rightgap==0:
                    sse_counter['stack']+=1;
                else:
                    if level_open==0:
                        if leftgap==0 or rightgap==0:
                            sse_counter['bulge']+=1;
                        else:
                            sse_counter['internal']+=1;
                    else:
                        sse_counter['multi']+=1;
            prev_openk = openk;
            prev_closek = closek;
        elif carac=='.':
            if len(stack)>0:
                sse_counter['freein']+=1;
            else:
                sse_counter['freeout']+=1;
        
    return sse_counter;

###############################################################################

# compute contiguity, number of isolate base pair and total number of base pairs.

def hairpincheck(structure):
    
    previi=-1;
    min_hairpin_length=len(structure);
    
    lindex = range(len(structure));
    for k in lindex:
        carac = structure[k];        
        if carac=='(':
            previi=k;
        elif carac==')':
            if previi>=0:
                hairpin_length=k-previi-1;
                if hairpin_length<min_hairpin_length:
                    min_hairpin_length=hairpin_length;
            previi=-1;
    
    return min_hairpin_length;

###############################################################################

# compute contiguity, number of isolate base pair and total number of base pairs.

def ss_stats(structure):

    stack=[];
    total_stem_length = 0;
    current_stem_length = 0;
    number_bp   = 0;
    number_contiguous_stacks = 0;
    prev_openk  = -1;
    prev_closek = -1;
    stack_len=0;
    number_isolate_bp = 0;
    min_hairpin_length=len(structure);
    
    lindex = range(len(structure));
    for k in lindex:
        carac = structure[k];
        if carac=='(':
            number_bp += 1;
            stack.append(k);
            last_bp = carac;
        elif carac==')':
            number_bp += 1;
            closek=k;
            openk=stack.pop();

            if prev_openk!=openk+1 or prev_closek!=closek-1:
                number_contiguous_stacks += 1;

            if len(stack)==0:
                total_stem_length += k - openk + 1;
                if total_stem_length<min_hairpin_length:
                    min_hairpin_length=total_stem_length;
                prev_openk  = -1;
                prev_closek = -1;
 
            if openk>0:
                if structure[openk-1]!='(' and structure[openk+1]!='(':
                    number_isolate_bp += 1;
            else:
                if structure[openk+1]!='(':
                    number_isolate_bp += 1;

            if closek+1>=len(structure):
                if structure[closek-1]!=')':
                    number_isolate_bp += 1;
            else:
                if structure[closek-1]!=')' and structure[closek+1]!=')':
                    number_isolate_bp += 1;

            prev_openk = openk;
            prev_closek = closek;
            last_bp = carac;

    if number_contiguous_stacks > 0:
        contiguity = float(number_bp + total_stem_length)/number_contiguous_stacks;
    else:
        contiguity = 0.0;

    number_of_bps = int(number_bp/2);

    return contiguity,number_isolate_bp,number_of_bps;

#############################################################################

# compute entropy of sequences

def compute_entropy(listseqs):

    lenseq=0;
    outdic={}
    ntfrac={};
    for i in range(len(listseqs[0])):
        ntfrac[i]={};
        ntfrac[i]['A']=0.0;
        ntfrac[i]['C']=0.0;
        ntfrac[i]['G']=0.0;
        ntfrac[i]['U']=0.0;

    if len(listseqs) > 0:
        for seq in listseqs:
            if lenseq==0:
                lenseq=len(seq);
            else:
                if lenseq!=len(seq):
                    print "sequencelist corrupted";
                    sys.exit(1);
            for i in range(len(seq)):
                ntfrac[i][seq[i]]+=1.0;
    # compute entropy
            
    outdic['nseqs'] = len(listseqs);
    ntlist=['A','C','G','U'];
    outdic['avg']=0.0;
    for i in range(lenseq):
        outdic[i]=0.0;
        for nt in ntlist:
            if ntfrac[i][nt] > 0:
                myfrac = float(ntfrac[i][nt])/outdic['nseqs'];
                val = -1.0 * myfrac * log(myfrac,4);
            else:
                myfrac = 0.0;
                val = 0.0
            outdic[i] += val;
            outdic['avg'] += val;
    outdic['avg'] /= lenseq;

    return outdic['avg'];

