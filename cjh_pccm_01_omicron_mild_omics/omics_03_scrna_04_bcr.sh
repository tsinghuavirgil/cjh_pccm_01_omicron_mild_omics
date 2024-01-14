##Using IgBLAST
contig_fasta=$1
contig_annotation=$2
sample=$3

AssignGenes.py igblast --outdir ./  -s $contig_fasta \
 -b /PROJ2/development/maxingyong/software/igblast_changeo/share/igblast \
 --organism human --exec /PROJ2/development/maxingyong/software/micromamba/envs/genomics/bin/igblastn \
  --loci ig --nproc 8 --format blast

MakeDb.py igblast -i $sample\_igblast.fmt7 \
 -s $contig_fasta \
 -r /PROJ2/development/maxingyong/software/igblast_changeo/share/germlines/imgt/human/vdj/imgt_human_*.fasta \
 --10x $contig_annotation --extended


#####Parsing IMGT output
sh sed_work.sh

for i in `cat file.list `; do tail -n  +2 $i >>merge.txt; done

ls add_sample_name/* | head -n 1 | xargs head -n 1 >header

cat header merge.txt >all_pass.txt

rm merge.txt

ParseDb.py select -d all_pass.txt -f locus -u "IGH" \
        --logic all --regex --outname heavy
ParseDb.py select -d all_pass.txt -f locus -u "IG[LK]" \
        --logic all --regex --outname light

awk -F "_" '{ if (++count[$1"_"$2"_"$3] == 1) print }' $input_file heavy_parse-select.txt >heavy_parse-select_filter.txt
awk -F "_" '{ if (++count[$1"_"$2"_"$3] == 1) print }' $input_file light_parse-select.txt >light_parse-select_filter.txt

###Parsing V(D)J data
dist_cutoff=0.1

DefineClones.py -d ./heavy_parse-select_filter.txt --act set --model ham --norm len --dist $dist_cutoff

python ./light_cluster.py -d ./heavy_parse-select_filter_clone-pass.tsv -e ./light_parse-select_filter.txt  -o all_clone-pass.tsv

DefineClones.py -d all_clone-pass.tsv --act set --model ham --norm len --dist $dist_cutoff

