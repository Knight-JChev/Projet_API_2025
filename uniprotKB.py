import requests, sys, json

def extraire_info_uniprot(dico_espece_gene):

    dico_uniprot = {'gene_symbol':'','uniprot_id':'','protein_name':'','pdb_id':''}
    uniprot =[]
    pdb_entries = []
    for espece, gene in dico_espece_gene.items():
        print(espece)
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
                dico_uniprot['gene_symbol']=gene
                dico_uniprot['uniprot_id']=accession
                dico_uniprot['protein_name']=protein_name['value']
                dico_uniprot['pdb_id']=pdb_entries
                uniprot.append(dico_uniprot)
    return  uniprot


uniprot = extraire_info_uniprot(dico)

print(uniprot)
# nom de gene = gene_symbol
