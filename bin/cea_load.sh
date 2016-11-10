#!/usr/bin/env bash

perl /home/vagrant/cxgn/cassbase/bin/transform_phenotypes_to_expression_atlas_lucy_index.pl -i /home/vagrant/Downloads/cass_phenotype.xls -o /home/vagrant/cxgn/cassbase/bin/lucy.tsv -c /home/vagrant/cxgn/cassbase/bin/pre_corr.tsv -f /home/vagrant/cxgn/cassbase/bin/corr.tsv -p /home/vagrant/cxgn/cassbase/bin/project.txt -d /home/vagrant/cxgn/cassbase/bin/desc.tsv

dropdb -U postgres cea
createdb -U postgres cea
psql -U postgres -d cea -h localhost -a -f /home/vagrant/cxgn/Tea/import_project/create_tea_schema.sql

perl /home/vagrant/cxgn/Tea/import_project/TEA_import_project_metadata_treatments.pl -d cea -H localhost -u postgres -p postgres -t /home/vagrant/cxgn/cassbase/bin/project.txt -n

rm -R /home/vagrant/cxgn/cea_indexed_files/correlation_indexes/cass_index/*
rm -R /home/vagrant/cxgn/cea_indexed_files/expression_indexes/cass_index/*
rm -R /home/vagrant/cxgn/cea_indexed_files/loci_and_description_index/cass_index/*

perl /home/vagrant/cxgn/Tea/import_project/index_correlation_file.pl /home/vagrant/cxgn/cassbase/bin/corr.tsv /home/vagrant/cxgn/cea_indexed_files/correlation_indexes/cass_index/
perl /home/vagrant/cxgn/Tea/import_project/index_expression_file.pl /home/vagrant/cxgn/cassbase/bin/lucy.tsv /home/vagrant/cxgn/cea_indexed_files/expression_indexes/cass_index/
perl /home/vagrant/cxgn/Tea/import_project/index_description_file.pl /home/vagrant/cxgn/cassbase/bin/desc.tsv /home/vagrant/cxgn/cea_indexed_files/loci_and_description_index/cass_index/

