#!/bin/sh

rm dump.rdb
echo 'removed dump.rdb'

rm output/rj/dereplicated_database/dereplicated.fasta
echo 'removed output/rj/dereplicated_database/dereplicated.fasta'

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

echo 'Database flushed'
redis-cli FLUSHALL