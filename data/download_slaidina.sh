vdb-config --prefetch-to-cwd

fasterq-dump -e 4 -p --split-files --include-technical $(cat SRR_Acc_List.txt | tr '\n' ' ')
#  --split-spot

ALLOWED=$(ls *_2.fastq | cut -f 1 -d "_" | tr " " "|")

echo 'sample_name,individual,Experiment,Instrument,sample_preparation_protocol' > ../config/sample_table.csv 
cut -d',' -f 1,7,17,20,28 SraRunTable.txt | tail -n+2 | grep -e "$ALLOWED" >> ../config/sample_table.csv 

echo 'sample_name,fq1_uri,fq2_uri,library_type,cdna_read,umi_read,barcode_read' > ../config/subsample_table.csv
tail -n+2 ../config/sample_table.csv | \
    cut -f 1 -d ',' | awk '{FS = ",";OFS = "," ;print $1,"data/"$1"_1.fastq","data/"$1"_2.fastq","10xV2","read 2","read 1"}' >> \
    ../config/subsample_table.csv
