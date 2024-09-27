# under dev for test

from Bio import Entrez
import requests

Entrez.email = "seung1.yoo@gmail.com"

def _get_ensembl_id(refseq_id):
    handle = Entrez.esearch(db="nucleotide", term=refseq_id)
    record = Entrez.read(handle)
    print(record)

    if not record['IdList']:
        print(f"No results found for {refseq_id}")
        return None

    nucleotide_id = record['IdList'][0]
    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="xml")
    records = Entrez.read(handle)

    for feature in records[0]['GBSeq_feature-table']:
        if feature['GBFeature_key'] == 'CDS':
            for qualifier in feature['GBFeature_quals']:
                if qualifier['GBQualifier_name'] == 'db_xref' and 'Ensembl' in qualifier['GBQualifier_value']:
                    ensembl_id = qualifier['GBQualifier_value'].split(":")[1]
                    return ensembl_id

    return None



def get_ensembl_id(refseq_id):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/homo_sapiens/{refseq_id}?"

    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)

    if not response.ok:
        response.raise_for_status()

    data = response.json()
    print(data)

    ensembl_id_s = list()
    for entry in data:
        if entry['type'] == 'gene':
            ensembl_id_s.append(entry['id'])

    return ensembl_id_s


def get_mane_info(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?expand=1"

    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)

    if not response.ok:
        response.raise_for_status()

    data = response.json()
    print(data)

    # MANE 정보 확인
    mane_info = {}
    if 'Transcript' in data:
        for transcript in data['Transcript']:
            if 'is_mane_select' in transcript and transcript['is_mane_select']:
                mane_info['MANE Select'] = transcript['id']
            if 'is_mane_plus_clinical' in transcript and transcript['is_mane_plus_clinical']:
                mane_info['MANE Plus Clinical'] = transcript['id']

    return mane_info


if __name__ == "__main__":
    outfh = open("hg38.HouseKeepingGenes.Ensembl.MANE", "w")
    for line in open("hg38.HouseKeepingGenes.bed"):
        items = line.rstrip().split("\t")
        refseq_id = items[3]
        ensembl_id_s = get_ensembl_id(refseq_id)

        if ensembl_id_s:
            for ensembl_id in ensembl_id_s:
                mane_info = get_mane_info(ensembl_id)
                print(f"Ensembl ID for {refseq_id}: {ensembl_id}")
                outfh.write(f"{refseq_id}\t{ensembl_id}\n")
        else:
            print(f"Ensembl ID not found for {refseq_id}")
            outfh.write(f"{refseq_id}\tN/A\n")
    outfh.close()

