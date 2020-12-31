import sys,os
from urllib.error import HTTPError
from Bio import Entrez
import csv


#microarry object for storing data
class microarray():
    def __init__(self,
                 experiment_id= "",
                 platform_id= "",
                 suppfile= "",
                 ftplink = "",
                 affy = False,
                 rnaseq= False):
        self.experiment_id= experiment_id
        self.platform_id= platform_id
        self.suppfile = suppfile
        self.ftplink = ftplink
        self.affy = affy
        self.rnaseq = rnaseq


# input: term_if, db
# output: list of all uid that connected to the [term_id]
def get_uid_list_from_esearch(term_id, db):
    try:
        handle = Entrez.esearch(db=db, term=term_id, rettype='json')
    except IOError as err:
        print(f'Other error occurred: {err}')
        handle.close()
        exit(2)
    esearch_res = Entrez.read(handle)
    handle.close()
    uid_list = esearch_res["IdList"]

    # [term_id] didn't return any uids- not a valid [term_id]
    if len(uid_list) < 1:
        print("{} is not valid".format(term_id))
        exit(2)
    return uid_list


# input: uid list as a string, db
# output: esummary_res json as object
def get_summary_obj_from_esummary(uid, db):
    try:
        handle = Entrez.esummary(db=db, id=uid, rettype='json')
    except IOError as err:
        print(f'Other error occurred: {err}')
        handle.close()
        exit(2)
    esummary_res = Entrez.read(handle)
    handle.close()
    return esummary_res


# input: uid list
# output: microarray object with relevant information
def create_microarray_obj(uid_list, db):
    # creating empty object
    microarray_obj = microarray()
    fields_filling_counter = 0
    for uid in uid_list:
        # all relevant information was obtained
        if fields_filling_counter >= 2:
            break
        experiment_summary = get_summary_obj_from_esummary(uid, db)
        uid_obj = experiment_summary[0]
        accession = uid_obj["Accession"]

        # gaining information of GSE
        if accession[:3] == "GSE":
            microarray_obj.experiment_id = accession
            microarray_obj.platform_id = "GPL"+uid_obj["GPL"]
            microarray_obj.suppfile = uid_obj["suppFile"]
            microarray_obj.ftplink = uid_obj["FTPLink"]
            extrelations_list = uid_obj["ExtRelations"]

            # searching for SRP and create runinfo file once found
            if len(extrelations_list) > 0:
                for item in extrelations_list:
                    target_object = item["TargetObject"]
                    if target_object[:3] == "SRP":
                        microarray_obj.rnaseq = True
                        file_name = accession+"_rnaseq.csv"
                        create_runinfo_file(file_name, target_object)

        # gaining information of GPL
        if accession[:3] == "GPL":
            if "Affymetrix" in uid_obj["title"]:
                microarray_obj.affy = True
            fields_filling_counter += 1

    return microarray_obj


# input: file_name, target_object_id
# output: [file_name] file
def create_runinfo_file(file_name, target_object_id):
    try:
        os.system("wget -O ./{file} 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={target}'".format(
            file=file_name, target=target_object_id
        ))
    except HTTPError as http_err:
        print(f'HTTP error occurred: {http_err}')
    except Exception as err:
        print(f'Other error occurred: {err}')


# input: microarray object
# output: [GSE_experiment_id_microarray].csv file
def print_microarray_to_file(experiment_summary):
    file_name = experiment_summary.experiment_id+"_microarray.csv"
    with open(file_name, "w") as csv_file:
        fieldnames = ['experiment_id', 'platform_id', 'suppfile', 'ftplink', 'affy', 'rnaseq']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(
            {
                'experiment_id': experiment_summary.experiment_id,
                'platform_id': experiment_summary.platform_id,
                'suppfile': experiment_summary.suppfile,
                'ftplink': experiment_summary.ftplink,
                'affy': experiment_summary.affy,
                'rnaseq': experiment_summary.rnaseq
            }
        )
        print("{file} was created successfully".format(file=file_name))


def handle_gse_identifier(identifier_id):
    uid_list = get_uid_list_from_esearch(identifier_id, "gds")
    microarray_obj = create_microarray_obj(uid_list, "gds")
    print_microarray_to_file(microarray_obj)


def handle_prjna_identifier(identifier_id):
    file_name = identifier_id+".csv"
    create_runinfo_file(file_name, identifier_id)


def main(email, identifier_id):
    Entrez.email = email
    if len(identifier_id) > 3 and identifier_id[:3] == "GSE":
        print("Receiving datasets for GSE: {}".format(identifier_id))
        handle_gse_identifier(identifier_id)
    elif len(identifier_id) > 5 and identifier_id[:5] == "PRJNA":
        print("Receiving datasets for PRJNA: {}".format(identifier_id))
        handle_prjna_identifier(identifier_id)
    else:
        print("identifier id: '{}' is not supported".format(identifier_id))


if __name__ == "__main__":
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usgae: datasets_ncbi [email] [identifier_id]")