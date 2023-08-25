for vars1 in allouts_1/*_peaks_id.bed_summit_freeblackhot_sorted_merged_1.bed
do
	BEDNAME_1=$(echo "${vars1}" | cut -d '/' -f 2 | cut -d '_' -f 1)
	for vars2 in allouts_2/*_peaks_id.bed_summit_freeblackhot_sorted_merged_2.bed 
	do
		BEDNAME_2=$(echo "${vars2}" | cut -d '/' -f 2 | cut -d '_' -f 1)
		closestdistance=$(bedtools closest -a ${vars1} -b ${vars2} -d | awk '{if($9 >= 0 && $9 <= 150) print $0}' | wc -l)
		totpeak1=$(cat ${vars1} | wc -l)
		totpeak2=$(cat ${vars2} | wc -l)
		echo "${BEDNAME_1}	${BEDNAME_2}	${totpeak1}	${totpeak2}	${closestdistance}"
	done
done
