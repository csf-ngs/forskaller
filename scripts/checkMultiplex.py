#!/usr/bin/env python

#check fsk3 for multiplexes with identical index combinations as in input file

import sys

def readIndices(indfile):
    indices = []
    with open(indfile) as inf:
         for line in inf:
             if line.lower().startswith("ind") or line.lower().startswith("bc"):
                continue 
             items = line.split()
             indices.append(items[0:2])
    return indices


def ind2q(indices, pos):
    li = ["\'"+i[pos]+"\'" for i in indices]
    return ",".join(li)

def createQuery(indices):
    first = " adaptor_tag in ("+ind2q(indices,0)+") "
    second = " adaptor_secondary_tag in ("+ind2q(indices, 1)+") "
    count = str(len(indices))
    dual = first+" and "+second
    single = first
    indp = dual
    q = "select multi_id, count(sample_id) from multiplex_samples join samples on sample_id = samples.id where multi_id in (select multi_id from multiplex_samples join samples on sample_id = samples.id where multi_id in (select multi_id from multiplex_samples join samples on sample_id = samples.id where "+indp+") and "+indp+" group by multi_id having count(sample_id) = "+count+") group by multi_id having count(sample_id) = "+count+";"
    return q



indfile = sys.argv[1]
indices = readIndices(indfile)
#print(indices)
istr = [i[0]+"\t"+i[1] for i in indices]
print("\n".join(istr))
q = createQuery(indices)
print(q)




