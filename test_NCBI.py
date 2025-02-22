# !/usr/bin/ python3
#-*- coding : utf-8 -*- 

from Bio import Entrez 
import json 

Entrez.email = "lin.christine0@gmail.com"


def GeneSymbol(filename) : 
    dico={}
    fichier = open(filename,'r')
    lines = fichier.readlines()
    for line in lines : 
        line= line.strip()
        if line : 
            li = line.split(',')
            gene =li[0]
            species = li[1]
            if '_gca' in species : 
                species = species.split('_gca')[0]  #retirer le "_gca" et ce qui suit 
            dico[species] = gene  

    return dico


def get_gene_id (gene_name, organism) : 

    # extraction des ID pour chaque gene ID 
    handle=Entrez.esearch(db="gene", term = f'{organism} [Orgn] AND {gene_name} [gene]')
    record = Entrez.read(handle)
    ids=record["IdList"]
    gene_id = ids[0]
    handle.close()

    # extraction du nom officiel 
    handle = Entrez.esummary(db="gene", id= gene_id, retmode="json")
    result = json.loads(handle.read())
    handle.close()

    official_name = result["result"].get(gene_id, {}).get("description")

    return ids, official_name


def get_linked_ids (gene_ids, target_db ) : 

    # Requete pour avoir les IDs dans les autres db à partir du gene ID 
    handle = Entrez.elink(dbfrom = "gene", db= target_db, id = gene_ids, retmode= 'json')
    result = json.loads(handle.read())
    handle.close()

    # extraction des IDs pour les différents db 
    linked_ids = []
    for linkset in result.get("linksets", []) : 
        for linksetdb in linkset.get("linksetdbs", []) :
            linked_ids.extend(linksetdb.get("links", []))
    
    return linked_ids


def get_seq_info (seq_ids, db) : 
    
    sequences = []
    
    # Requete pour avoir les RNA number acces ou transcrit number access selon la db 
    handle = Entrez.esummary(db = db, id = ",".join(seq_ids), retmode = 'json')
    result = json.loads(handle.read())
    handle.close()

    # extraction des informatiosn RNA number access ou transcrit number access
    for seq_id in seq_ids : 
        seq_info = result.get("result", {}).get(str(seq_id), {})
        if seq_info : 
            acc_version = seq_info.get("accessionversion","")  
            if acc_version.startswith(("NP_", "NM_", "XP_", "XM_")) : 
                if acc_version not in sequences : 
                    sequences.append(acc_version)
    
    return sequences


def extract_info (gene_symbol, organism) : 

    gene_id, name = get_gene_id(gene_symbol, organism)

    if gene_id :
        # extraction des ids pour les db prot et nucléotide
        protein_ids = get_linked_ids(gene_id, "protein")  # db protein
        nucleotide_ids = get_linked_ids(gene_id, 'nuccore') # db nucleotide

        # extraction des numéros d'access prot et transcrit
        protein_info = get_seq_info(protein_ids, "protein")
        nucleotide_info = get_seq_info(nucleotide_ids, 'nuccore')

    return   {"Gene Symbol " : gene_symbol,
              "Name" : name,
              "Gene ID " : gene_id, 
              "Protein" : protein_info, 
              "Transcript" : nucleotide_info
              }

def main (dico_sp) :
    dico_NCBI={}
    # faire un dico avec tout les infos pour tout les gene/espèce
    for species, gene in dico_sp.items() : 
        res = extract_info(gene_symbol= gene, organism= species)
        dico_NCBI[species] = res
    return dico_NCBI


#fichier = input("Entrez une fichier avec pour chaque ligne gene,espece : ")
#dico = GeneSymbol(fichier)
#print(main(dico_sp= dico))
