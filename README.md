# datasets_ncbi
Program recieves inforamtion for GSE or PRJNA identifiers from ncbi.

Usage:
python datasets_ncbi.py [email] [identifier]

Input: email address and valid GSE or PRJNA identifier

Output:
1. If PRJNA was given, file named "[PRJNA_id].csv"
2. If GSE was given:
  a. File named "[GSE_id]_microarray.csv"
  b. Optinal file "[GSE_id]_rnaseq.csv" if srp identifier was found in GSE datasets.
 
