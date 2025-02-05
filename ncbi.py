# !/usr/bin/ python3
#-*- coding : utf-8 -*- 

from Bio import Entrez, SeqIO
import re
import sys
import os.path

def Request (species_info, db, requete, retention) :
    """
    Effectue les requêtes API pour récupérer les informations sur le NCBI.
    """


    Entrez.email = "Chevreau.julien21@gmail.com@univ-rouen.fr"

    # Va chercher les IDs dans la base de donnée demandée
    handle = Entrez.esearch(db=db, term=requete, retmax=retention)
    records = Entrez.read(handle)
    identifiers = records["IdList"]

    # Cherche les données à partir des IDs. Résultat en format GenBank
    handle = Entrez.efetch(db=db, id=identifiers, retmax="100",
                            rettype="gbk", retmode="text")
    text = handle.read()

    # Construit un fichier temporaire pour être lu ensuite
    filename = "tmp_%s.gbk" % (db)
    with open(filename,"w") as gbk:
        gbk.write(text)

    return(filename)

def Gene (species_info) :
    """
    Récupère les informations relatives au gène d'intérêt
    """

    # Requête pour la base de donnée
    requete = "%s[Orgn] AND %s[gene Name] AND Refseq [Keyword]" % (species_info["species"], species_info["gene_symbol"])
    gene_file = Request(species_info,"gene", requete, 1)

    # Ouverture et lecture du fichier temporaire contenant les informations du gène
    with open(gene_file,"r") as file:
        content = "".join(file.readlines())

        # Expressions régulières pour trouver ID et nom officiel
        pattern_id = re.search(r"ID:\s(\d*)",content)
        pattern_official = re.search(r"Name:\s(.*)\s\[",content)
        pattern_official_alt = re.search(r"(.*)\s\[.*\]",content)

        # Population du dictionnaire en cas de match
        if pattern_official != None :
            species_info["official_name"] = pattern_official.group(1)
        elif pattern_official_alt != None : # Pattern alternatif dans le cas où le premier n'est pas trouvé
            species_info["official_name"] = pattern_official_alt.group(1)

        if pattern_id != None :
            species_info["gene_id"] = pattern_id.group(1)
    return(species_info)

def ProTranscript (species_info):
    """
    Récupère les informations relatives aux protéines et transcrits traduits
    """

    # Requête qui récupère les entrées officielles RefSeq
    requete = "%s[Orgn] AND %s[gene Name] AND srcdb refseq [properties]" % (species_info["species"], species_info["gene_symbol"])

    # Passage de la requête à la fonction qui va chercher les informations
    prot_file = Request(species_info, "protein",requete ,60)

    # Parcours du fichier temporaire créé pour en extraire les informations
    for seq_record in SeqIO.parse(prot_file, "genbank"):
        try : # Si db_source n'existe pas pour l'entrée en cours, passe à la suivante
            db_source = seq_record.annotations["db_source"].split(" ")

            if db_source[0] == "REFSEQ:" :
                if "transcript_id" not in species_info :
                    species_info["transcript_id"] = [db_source[-1]]
                    species_info["prot_id"] = [seq_record.id]
                else :
                    species_info["transcript_id"].append(db_source[-1])
                    species_info["prot_id"].append(seq_record.id)
        except :
            continue

    # Indique si aucun transcrit/protéine n'est trouvé            
    if "transcript_id" not in species_info :
                species_info["transcript_id"] = ["No data found"]
                species_info["prot_id"] = ["No data found"]
    return(species_info)

def Info (species_info):
    """
    Aggrège les différentes fonctions utilisées pour aller chercher les informations 
    """
    # Formate le nom de l'espèce pour le NCBI. Utile pour les bactéries.
    species_info["species"] = "_".join(species_info["species"].split("_")[0:2])

    # Appelle les fonctions pour récupérer les informations
    species_info = Gene(species_info)
    species_info = ProTranscript(species_info)
    return(species_info)

# En cas de lancement par ligne de commande : lit le fichier en entrée
if __name__ == '__main__':
    gene_symbols = sys.argv[1]
    if os.path.isfile(gene_symbols): # Vérifie que l'argument soit un fichier
        with open(gene_symbols, "r") as infos:
            for line in infos:
                symbol = line.split(",")[0]
                current_species = line[:-1].split(",")[1]

                # Dictionnaire initial à passer en argument aux sous-scripts
                species_info = {}
                species_info["species"] = current_species
                species_info["gene_symbol"] = symbol

                print("\t",Info(species_info))
    else :
        print("Erreur : Veuillez donnez un nom de fichier accessible comme seul argument \n"\
        "Format de chaque ligne : [Symbole de gène],[Espèce] \n" \
        "Ex: RAD51,homo_sapiens")