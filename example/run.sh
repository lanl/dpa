export DPA_HOME=${PWD}/..
export PDB_ID=1jef

curl -o ${DPA_HOME}/structuredata/${PDB_ID}_raw.pdb https://files.rcsb.org/download/${PDB_ID}.pdb

grep -e HETATM ${DPA_HOME}/structuredata/${PDB_ID}_raw.pdb | grep -v HOH > ${DPA_HOME}/structuredata/ligand_${PDB_ID}.pdb

awk '/^ATOM/' ${DPA_HOME}/structuredata/${PDB_ID}_raw.pdb > ${DPA_HOME}/structuredata/${PDB_ID}.pdb
python ../src/generate_surface_points.py ${DPA_HOME}/structuredata/${PDB_ID}.pdb  ${DPA_HOME}/structuredata/${PDB_ID}.surf 
echo "${PDB_ID}" > pdb_list.txt
perl ../src/ggdpa.pl -f pdb_list.txt
perl ../src/dpa.pl -f pdb_list.txt 1 1 -topp 0.98 -cutoff 6 6 -walldpa -wcpdb
perl ../src/extract_dpa_results.pl ${PDB_ID}
