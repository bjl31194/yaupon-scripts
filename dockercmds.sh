#! /bin/bash

docker pull ajing/chemtreemap
docker run -t -i -p 8000:8000 --name chemtreemap ajing/chemtreemap /bin/bash
docker cp ~/yaupon/metabolomics/smiles.txt chemtreemap:examples/smiles.txt

docker exec -it chemtreemap bin/bash
cd examples
