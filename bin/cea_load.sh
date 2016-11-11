#!/usr/bin/env bash

#./../cassbase/bin/cea_load.sh /home/vagrant/Downloads/cass_phenotype.xls /home/vagrant/cxgn/cassbase/bin/ true true cass_index_CASS_6Genotypes_Sampling_2015
#./../cassbase/bin/cea_load.sh /home/vagrant/Downloads/t100_phenotypes.xls /home/vagrant/cxgn/cassbase/bin/ false true cass_index_T100

input_file_name=$1
processing_file_dir=$2 #/home/vagrant/cxgn/cassbase/bin
clear_old_database=$3
clear_old_index_files=$4
directory_name=$5 #cass_index_<Project_name>

perl /home/vagrant/cxgn/cassbase/bin/transform_phenotypes_to_expression_atlas_lucy_index.pl -i ${input_file_name} -o ${processing_file_dir}/lucy.tsv -c ${processing_file_dir}/pre_corr.tsv -f ${processing_file_dir}/corr.tsv -p ${processing_file_dir}/project.txt -d ${processing_file_dir}/desc.tsv

#DEBUG=false
if ${clear_old_database}; then
dropdb -U postgres cea
createdb -U postgres cea
psql -U postgres -d cea -h localhost -a -f /home/vagrant/cxgn/Tea/import_project/create_tea_schema.sql
fi

if ${clear_old_index_files}; then  
rm -R /home/vagrant/cxgn/cea_indexed_files/correlation_indexes/${directory_name}/*
rm -R /home/vagrant/cxgn/cea_indexed_files/expression_indexes/${directory_name}/*
rm -R /home/vagrant/cxgn/cea_indexed_files/loci_and_description_index/${directory_name}/*
fi

perl /home/vagrant/cxgn/Tea/import_project/TEA_import_project_metadata_treatments.pl -d cea -H localhost -u postgres -p postgres -t ${processing_file_dir}/project.txt -n

perl /home/vagrant/cxgn/Tea/import_project/index_correlation_file.pl ${processing_file_dir}/corr.tsv /home/vagrant/cxgn/cea_indexed_files/correlation_indexes/${directory_name}/
perl /home/vagrant/cxgn/Tea/import_project/index_expression_file.pl ${processing_file_dir}/lucy.tsv /home/vagrant/cxgn/cea_indexed_files/expression_indexes/${directory_name}/
perl /home/vagrant/cxgn/Tea/import_project/index_description_file.pl ${processing_file_dir}/desc.tsv /home/vagrant/cxgn/cea_indexed_files/loci_and_description_index/${directory_name}/

