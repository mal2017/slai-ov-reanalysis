echo 'sample_name,individual,Experiment,Instrument,sample_preparation_protocol' > ../config/sample_table.csv 
cut -d',' -f 1,7,17,20,28 SraRunTable.txt | tail -n+2 >> ../config/sample_table.csv 

echo 'sample_name,fq1_uri,fq2_uri,library_type,cdna_read,umi_read,barcode_read' > tmp_subsample_table.csv
tail -n+2 ../config/sample_table.csv | \
    cut -f 1 -d ',' | awk '{FS = ",";OFS = "," ;print $1,"data/"$1"_1.fastq","data/"$1"_2.fastq","10x","read 2","read 1"}' >> \
    tmp_subsample_table.csv

Rscript make_full_sample_table.R

rm tmp_subsample_table.csv