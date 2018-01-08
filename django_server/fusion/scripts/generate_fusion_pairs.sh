cat $file | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$4,$2,$5,"http://cancer.sanger.ac.uk/cosmic/fusion/overview?fid="$2"&gid="$5}'
