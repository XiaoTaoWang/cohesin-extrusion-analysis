import bisect

def parse_peaks(fil):

    D = {}
    with open(fil, 'r') as source:
        for line in source:
            c, s, e = line.rstrip().split()[:3]
            s, e = int(s), int(e)
            if not c in D:
                D[c] = []
            D[c].append((s, e))
    
    for c in D:
        D[c].sort()
    
    return D

def bisect_search(p, ref):

    cache = set()
    if not p[0] in ref:
        return cache
    
    List = ref[p[0]]
    idx = max(0, bisect.bisect(List, p[1:])-1)
    for q in List[idx:]:
        if q[1] <= p[1]:
            continue
        if q[0] >= p[2]:
            break
        cache.add((p[0],)+q)
    
    return cache

def load_convergent_loops(infil):

    loops = []
    with open(infil, 'r') as source:
        for line in source:
            c1, s1, e1, c2, s2, e2, label = line.rstrip().split()
            if label == 'convergent':
                s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
                loops.append((c1, s1, e1, c2, s2, e2))
    
    return loops

def calculate_distance(start, end, candi):

    min_dis = 0
    for t in candi:
        d1 = abs(t[1] - start)
        d2 = abs(t[1] - end)
        d3 = abs(t[2] - start)
        d4 = abs(t[2] - end)
        if min([d1, d2, d3, d4]) > min_dis:
            min_dis = min([d1, d2, d3, d4])
    
    return min_dis

loops = load_convergent_loops('K562_CTCF_merge_distance_Convergent_Loops_all_7col.bedpe')
anchors = parse_peaks('K562.CTCF-loop-anchors.bed')
loading_sites = parse_peaks('K562.cohesin-loading-sites.new.bed')
idx1 = 0
idx2 = 0
collect = []
for c1, s1, e1, c2, s2, e2 in loops:
    cache = bisect_search((c1, e1, s2), anchors)
    min_dis = calculate_distance(e1, s2, cache)
    if min_dis > 10000:
        continue # remove loops with nested loops
    
    tmp_collect = []
    cache = bisect_search((c1, e1, s2), loading_sites)
    for t in cache:
        d1 = abs(t[1] - e1)
        d2 = abs(t[1] - s2)
        d3 = abs(t[2] - e1)
        d4 = abs(t[2] - s2)
        if min([d1, d2, d3, d4]) > 0.25*(s2-e1):
            if not len(tmp_collect):
                idx1 += 1
                tmp_collect.append([c1, s1, e1, 'cr{0}_S'.format(idx1)])
                tmp_collect.append([c2, s2, e2, 'cr{0}_E'.format(idx1)])
            idx2 += 1
            tmp_collect.append([c1, t[1], t[2], 'cr{0}_M-{1}'.format(idx1, idx2)])
    
    collect.extend(tmp_collect)

with open('K562-loading-sites-within-convergent-loops.apa.bed', 'w') as out:
    for t in collect:
        out.write('\t'.join(list(map(str, t)))+'\n')