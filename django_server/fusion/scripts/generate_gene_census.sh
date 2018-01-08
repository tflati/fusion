cat "$1" | tail -n +2 | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,"http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln="$1}'
