#!/usr/bin/env bash

#./../cassbase/bin/cea_load.sh /home/vagrant/Downloads/cass_phenotype.xls /home/vagrant/cxgn /home/vagrant/cxgn/cassbase/bin true true cass_index_CASS_6Genotypes_Sampling_2015
#./../cassbase/bin/cea_load.sh /home/vagrant/Downloads/t100_phenotypes.xls /home/vagrant/cxgn /home/vagrant/cxgn/cassbase/bin false true cass_index_T100


trial_id=$1 #/home/vagrant/Downloads/cass_phenotype.xls
root_dir=$2 #/home/vagrant/cxgn
processing_file_dir=$3 #/home/vagrant/cxgn/cassbase/bin
clear_old_database=$4 #true or false
clear_old_index_files=$5 #true or false
directory_name=$6 #cass_index_<Project_name>

wget --output-document=${processing_file_dir}/phenotype_download.xls "http://0:3000/brapi/v1/studies/${trial_id}/table?format=xls&observationLevel=plant"
sleep 10

perl ${root_dir}/cassbase/bin/transform_phenotypes_to_expression_atlas_lucy_index.pl -i ${processing_file_dir}/phenotype_download.xls -o ${processing_file_dir}/lucy.tsv -c ${processing_file_dir}/pre_corr.tsv -f ${processing_file_dir}/corr.tsv -p ${processing_file_dir}/project.txt -d ${processing_file_dir}/desc.tsv

if ${clear_old_database}; then
dropdb -U postgres cea
createdb -U postgres cea
#DROP SCHEMA public CASCADE;
#CREATE SCHEMA public;
psql -U postgres -d cea -h localhost -a -f ${root_dir}/Tea/import_project/create_tea_schema.sql
fi

if ${clear_old_index_files}; then
rm -R ${root_dir}/indexed_files/correlation_indexes/${directory_name}/*
rm -R ${root_dir}/indexed_files/expression_indexes/${directory_name}/*
rm -R ${root_dir}/indexed_files/loci_and_description_index/${directory_name}/*
fi

perl ${root_dir}/Tea/import_project/TEA_import_project_metadata_treatments.pl -d cea -H localhost -u postgres -p postgres -t ${processing_file_dir}/project.txt -n

perl ${root_dir}/Tea/import_project/index_correlation_file.pl ${processing_file_dir}/corr.tsv ${root_dir}/indexed_files/correlation_indexes/${directory_name}/
perl ${root_dir}/Tea/import_project/index_expression_file.pl ${processing_file_dir}/lucy.tsv ${root_dir}/indexed_files/expression_indexes/${directory_name}/
perl ${root_dir}/Tea/import_project/index_description_file.pl ${processing_file_dir}/desc.tsv ${root_dir}/indexed_files/loci_and_description_index/${directory_name}/

