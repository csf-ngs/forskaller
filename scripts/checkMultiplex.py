#!/usr/bin/env python

#check fsk3 for multiplexes with identical index combinations as in input file

import sys

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def readIndices(indfile, rc):
    indices = []
    with open(indfile) as inf:
         for line in inf:
             if line.lower().startswith("ind") or line.lower().startswith("bc") or line.lower().startswith("tag"):
                continue 
             items = line.split()
             bcs = items[0:2]
             if rc == "true":
               bcs[1] = reverse_complement(bcs[1])
             indices.append(bcs)
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
    #indp = dual
    indp = single
    q = "select multi_id, count(sample_id) from multiplex_samples join samples on sample_id = samples.id where multi_id in (select multi_id from multiplex_samples join samples on sample_id = samples.id where multi_id in (select multi_id from multiplex_samples join samples on sample_id = samples.id where "+indp+") and "+indp+" group by multi_id having count(sample_id) = "+count+") group by multi_id having count(sample_id) = "+count+";"
    return q



indfile = sys.argv[1]
rc = sys.argv[2]
indices = readIndices(indfile, rc)
#print(indices)
istr = [i[0]+"\t"+i[1] for i in indices]
print("\n".join(istr))
q = createQuery(indices)
print(q)




