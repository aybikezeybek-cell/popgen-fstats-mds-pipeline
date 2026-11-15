#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# End-to-End Ancient DNA f-Statistics Pipeline (EIGENSTRAT -> AdmixTools)
# Usage:  bash scripts/01_run_admixtools.sh <dataset_prefix> <outgroup>
# ==============================================================================

PREFIX="${1:-v66.p1_1240K.aadr.patch.PUB}"
OUTGROUP="${2:-Mbuti}"

WORK_DIR="outputs/intermediate"
LOG_DIR="outputs/logs"
SUBSET="${WORK_DIR}/anatolia_subset"

mkdir -p "$WORK_DIR" "$LOG_DIR"

POPS=(
  "Turkey_PPN"
  "Turkey_N"
  "Turkey_Catalhoyuk_EN"
  "Turkey_Catalhoyuk_MN"
  "Turkey_Catalhoyuk_LN"
  "Jordan_PPNB"
  "Iran_GanjDareh_N"
)

# 1. Create population list including outgroup
printf "%s\n" "${POPS[@]}" "$OUTGROUP" > "${WORK_DIR}/poplist.txt"

# 2. Extract lightweight EIGENSTRAT subset via convertf
cat << EOF > "${WORK_DIR}/par_extract.txt"
genotypename:    ${PREFIX}.geno
snpname:         ${PREFIX}.snp
indivname:       ${PREFIX}.ind
outputformat:    EIGENSTRAT
genotypeoutname: ${SUBSET}.geno
snpoutname:      ${SUBSET}.snp
indivoutname:    ${SUBSET}.ind
poplistname:     ${WORK_DIR}/poplist.txt
EOF

convertf -p "${WORK_DIR}/par_extract.txt"

# 3. Generate unique unordered pairs for Outgroup f3(P1, P2; Outgroup)
F3_POP="${WORK_DIR}/f3_popfile.txt"
> "$F3_POP"
N=${#POPS[@]}
for ((i=0; i<N; i++)); do
  for ((j=i+1; j<N; j++)); do
    printf "%s\t%s\t%s\n" "${POPS[$i]}" "${POPS[$j]}" "$OUTGROUP" >> "$F3_POP"
  done
done

# 4. Generate non-redundant quartets for f4(W, X; Y, Outgroup)
F4_POP="${WORK_DIR}/f4_popfile.txt"
> "$F4_POP"
for ((i=0; i<N; i++)); do
  for ((j=i+1; j<N; j++)); do
    for ((k=0; k<N; k++)); do
      if [[ $k -ne $i && $k -ne $j ]]; then
        printf "%s\t%s\t%s\t%s\n" "${POPS[$i]}" "${POPS[$j]}" "${POPS[$k]}" "$OUTGROUP" >> "$F4_POP"
      fi
    done
  done
done

# 5. Run qp3Pop & qpDstat
cat << EOF > "${WORK_DIR}/par_f3.txt"
genotypename: ${SUBSET}.geno
snpname:      ${SUBSET}.snp
indivname:    ${SUBSET}.ind
popfilename:  ${F3_POP}
inbreed:      YES
EOF

cat << EOF > "${WORK_DIR}/par_f4.txt"
genotypename: ${SUBSET}.geno
snpname:      ${SUBSET}.snp
indivname:    ${SUBSET}.ind
popfilename:  ${F4_POP}
f4mode:       YES
EOF

qp3Pop -p "${WORK_DIR}/par_f3.txt" > "${LOG_DIR}/f3_results.log"
qpDstat -p "${WORK_DIR}/par_f4.txt" > "${LOG_DIR}/f4_results.log"
