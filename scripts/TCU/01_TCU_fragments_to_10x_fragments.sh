# Software versions
# bedtools v2.30.0
# bgzip (htslib) 1.12
# wd: liam_manuscript_reproducibility
# Fragment files from GEO are not 10X formatted fragments files!
# ArchR expects 10X formatted fragment files, in practice this means that the gzipped BED files (that is what you can download from GEO) need to be sorted, and then compressed with bgzip and and index created with tabix.
# Useful threads that pointed me to this solution: 
# My file was first not sorted, this led tabix to fail: https://www.biostars.org/p/386768/
# Here 10x fragments format is described: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
# Mention that ArchR expects 10x formatted fragment files: https://github.com/GreenleafLab/ArchR/issues/727

echo "Convert ATAC-seq data from Mimitou et al., 2021 for ASAP-seq and DOGMA-seq T cell stimulation experiments from fragments.tsv.gz to 10-like fragments.tsv.gz"

for datatype in data/original/Mimitou2021/*/ ; do
    echo "$datatype"
    for directory in $datatype*/ ; do
        sample_cond=`basename ${directory}`
        echo "$sample_cond"
        if [ "$sample_cond" = "from_asap_large_data_files" ];
        then continue
        else
        fragments_zip=`echo "${directory}ATAC/"*fragments*`
        echo "$fragments_zip"
        file_name=`basename ${fragments_zip}`
        echo "$file_name"
        echo "unzip, sort and bgzip"
        echo "zcat $fragments_zip | bedtools sort | bgzip -c > data/derived/Mimitou2021/$file_name"
        zcat $fragments_zip | bedtools sort | bgzip -c > data/derived/Mimitou2021/$file_name
        echo "create index"
        echo "tabix -p bed data/derived/Mimitou2021/$file_name"
        tabix -p bed data/derived/Mimitou2021/$file_name
        echo ""
        fi
    done
done