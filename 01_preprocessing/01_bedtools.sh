awk '$7 > 5' K562_CTCF_0h_NIPBL_CUT_TAG005_peaks.narrowPeak | \
    sort -k8,8nr | \
    awk '!seen[$1,$2,$3]++' | \
    head -n 20000 > NIPBL_sig5_top20k.narrowPeak

sort -k8,8g -k7,7nr K562_WT_CTCF_ChIA-PET_distance_50k3MB_Convergent_Loops_PET9.bedpe | \
    head -n 6000 > distance_Convergent_Loops_top10k.bedpe

nipblp="NIPBL_sig5_top20k.narrowPeak"
rad21p="RAD21_top50k.narrowPeak"
rnapiip="K562_CTCF_0h_Pol2_CUT_TAG004_peaks.narrowPeak"
h3k27acp="K562_H3K27ac_chip-seq_ENCFF544LXB.bed"
h3k4me1p="K562_H3K4me1_chip-seq_ENCFF540NGG.bed"

bedtools intersect -a ${nipblp} -b ${rad21p} -u | \
bedtools intersect -a - -b ${rnapiip} -u | \
bedtools intersect -a - -b ${h3k27acp} -u | \
bedtools intersect -a - -b ${h3k4me1p} -u > k562-chip-nipbl-overlap-rad21-rnapii-h3k27ac-h3k4me1.bed