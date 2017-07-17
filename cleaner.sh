#!/bin/sh

rm dump.rdb
echo 'removed dump.rdb'

rm output/rj/graph_dereplicated.fasta
echo 'removed output/rj/graph_dereplicated.fasta'

rm output/rj/dereplicated_database/dump.rdb
echo 'removed output/rj/dereplicated_database/dump.rdb'

rm output/rj/graph.des
echo 'removed output/rj/graph.des'

rm output/rj/graph.esq
echo 'removed output/rj/graph.esq'

rm output/rj/graph.rlt
echo 'removed output/rj/graph.rlt'

rm output/rj/graph.sds
echo 'removed output/rj/graph.sds'

rm output/rj/subgraph_temp.edge.list
echo 'removed output/rj/subgraph_temp.edge.list'

rm output/rj/subgraph_temp_dereplicated.fasta
echo 'removed output/rj/subgraph_temp_dereplicated.fasta'

rm output/rj/subgraph_temp.0.spm
echo 'removed output/rj/subgraph_temp.0.spm'

rm output/rj/subgraph_temp.des
echo 'removed output/rj/subgraph_temp.des'

rm output/rj/subgraph_temp.esq
echo 'removed output/rj/subgraph_temp.esq'

rm output/rj/subgraph_temp.rlt
echo 'removed output/rj/subgraph_temp.rlt'

rm output/rj/subgraph_temp.sds
echo 'removed output/rj/subgraph_temp.sds'

echo 'Database flushed'
redis-cli FLUSHALL