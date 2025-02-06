import requests, sys, json

fichier = "LC_GeneSymbols_45.txt"

def liste_gene_uniprot(fichier):
    dico_espece_gene = {}

    with open(fichier, "r") as file:
        for ligne in file:
            li = ligne.strip().split(",")  # Suppression des espaces inutiles et séparation

            if len(li) == 2:  # Vérification que la ligne est correcte
                gene = li[0].strip()
                espece = li[1].strip().replace("_", " ").lower()  # Tout en minuscule

                # Nettoyage des noms d'espèces s'ils contiennent "gca"
                if "gca" and "coli" in espece:  
                    espece = "Escherichia coli (strain K12)"  # Séparation et nettoyage
                
                espece = espece.capitalize()  # Mettre seulement la première lettre en majuscule
                
                dico_espece_gene[gene] = espece  # Stockage dans le dictionnaire

    return dico_espece_gene 

dico_espece_gene = liste_gene_uniprot(fichier)
print(dico_espece_gene)


def extraire_info_uniprot(dico_espece_gene):
    dico_uniprot_id = {}
    dico_uniprot_nom = {}
    dico_uniprot_pdb = {}
    dico_uniprot = {'gene_symbol':'','uniprot_id':'','protein_name':'','pdb_id':''}
    uniprot =[]
    pdb_entries = []
    for gene, espece in dico_espece_gene.items():
        query = f"gene:{gene} AND reviewed:true AND {espece}"
        params = {
        "query": query,
        "fields": "accession,protein_name,xref_pdb",
        "sort": "accession desc",
        "size": 3  # Limite à 3 résultats
    }

        headers = {
        "accept": "application/json"
    }
        base_url = "https://rest.uniprot.org/uniprotkb/search"

        response = requests.get(base_url, headers=headers, params=params)
    
        if not response.ok:
            print(f"Erreur pour {gene} ({espece}) : {response.status_code}")
            continue  # Passe au gène suivant au lieu d'arrêter

    # Extraction des données JSON
        data = response.json()
        results = data.get("results", [])

        if not results:
            query = f"gene:{gene} AND {espece}"  # Suppression de "reviewed:true"
            params = {
        "query": query,
        "fields": "accession,protein_name,xref_pdb",
        "sort": "accession desc",
        "size": 1  # Limite à 3 résultats
    }
            response = requests.get(base_url, headers=headers, params=params)
        
            if not response.ok:
                continue
        
            data = response.json()
            results = data.get("results", [])


    # Affichage des résultats pour ce gène

        if not results:
            continue
        else:
            for entry in results:
                accession = entry.get('primaryAccession', 'N/A')
                protein_name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', 'N/A')
                for xref in entry.get("uniProtKBCrossReferences"):
                    if xref.get("database") == "PDB":
                        pdb_entries.append(xref.get("id"))

            # Stockage dans les dictionnaires
                dico_uniprot_id[gene] = accession
                dico_uniprot_nom[gene] = protein_name['value']
                dico_uniprot_pdb[gene] = pdb_entries
                dico_uniprot['gene_symbol']=gene
                dico_uniprot['uniprot_id']=accession
                dico_uniprot['protein_name']=protein_name['value']
                dico_uniprot['pdb_id']=pdb_entries
                uniprot.append(dico_uniprot)
    return dico_uniprot_id,dico_uniprot_nom,dico_uniprot_pdb, uniprot


dico_uniprot_id, dico_uniprot_nom, dico_uniprot_pdb, uniprot = extraire_info_uniprot(dico_espece_gene)

print(dico_uniprot_id)
print(dico_uniprot_nom)
print(dico_uniprot_pdb)
print(uniprot)
#fusionner les 3 dictionnaire en 1
# nom de gene = gene_symbol