import pandas as pd

# Functions
# Get location
def GetLoc(x, refpoint):
    if refpoint == "TSS":
        if x["probe.strand"]=="+":
            if x["str.start"]<x["gene.start"]: return "UPSTREAM"
            else: return "DOWNSTREAM"
        else:
            if x["str.start"]>x["gene.start"]: return "UPSTREAM"
            else: return "DOWNSTREAM"
    else:
        if x["probe.strand"]=="+":
            if x["str.start"]<x["gene.stop"]: return "UPSTREAM"
            else: return "DOWNSTREAM"
        else:
            if x["str.start"]>x["gene.stop"]: return "UPSTREAM"
            else: return "DOWNSTREAM"

# Get perc of STRs that are significant
def GetEnrich(data, p_cutoff, norm=True):
    if norm:
        num_strs = data[["chrom","str.start"]].drop_duplicates().shape[0]
    else: num_strs = 1 # return raw number
    if p_cutoff == "signif":
        num_enrich = data[(data["best_str"]) & (data["signif"])][["chrom","str.start"]].drop_duplicates().shape[0]
    else:
        num_enrich = data[(data["best_str"]) & (data["p.wald"]<=p_cutoff)][["chrom","str.start"]].drop_duplicates().shape[0]
    return num_enrich*1.0/num_strs

# Get number of STRs in each bin
def GetNumStrsProfile(data, location, distcol, distbins):
    upstream_num = []
    downstream_num = []
    for i in range(len(distbins)-1):
        l = distbins[i]
        u = distbins[i+1]
        ups = data[(location=="UPSTREAM") & (data[distcol]>l) & (data[distcol]<u)]
        downs = data[(location=="DOWNSTREAM") & (data[distcol]>l) & (data[distcol]<u)]
        upstream_num.append(ups[["chrom","str.start"]].drop_duplicates().shape[0])
        downstream_num.append(downs[["chrom","str.start"]].drop_duplicates().shape[0])
    num = upstream_num[::-1] + downstream_num
    bins = map(lambda x: -1*x, list(distbins[::-1])[:-1]) + list(distbins)
    return bins, num
    
# Function to get profile around TSS/TES using different p.value as cutoff
def GetEnrichmentProfile(data, location, distcol, p_cutoff, distbins):
    upstream_enrich = []
    downstream_enrich = []
    for i in range(len(distbins)-1):
        l = distbins[i]
        u = distbins[i+1]
        ups = data[(location=="UPSTREAM") & (data[distcol]>l) & (data[distcol]<u)]
        downs = data[(location=="DOWNSTREAM") & (data[distcol]>l) & (data[distcol]<u)]
        upstream_enrich.append(GetEnrich(ups, p_cutoff))
        downstream_enrich.append(GetEnrich(downs, p_cutoff))
    enrich = upstream_enrich[::-1] + downstream_enrich
    bins = map(lambda x: -1*x, list(distbins[::-1])[:-1]) + list(distbins)
    return bins, enrich

def Smooth(x, extend=1):
    newx = []
    for i in range(len(x)):
        left = i-extend
        if left < 0: left = 0
        right = i+extend
        if right >= len(x): right = len(x)-1
        num = right-left + 1
        newx.append(sum(x[left:right+1])*1.0/num)
    return newx
