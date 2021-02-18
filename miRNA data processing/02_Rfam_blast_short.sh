PATH=${PATH}

blastn=/glusterfs/home/song_yunjie/software/ncbi-blast-2.5.0+/bin/blastn
db=/glusterfs/home/wu_liping0101/database/rfam/rfam

$blastn -db $db -query dir_mapper_seq_*/reads_nr.fa -task blastn-short -outfmt 6 -out ${prefix}.rfam.out

sed '1i qseid\tsseid\tpident\tlength\tmismatch\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' ${prefix}.rfam.out > ${prefix}.rfam.xls
