import requests, json

# creation dico 
liste = input("Entrez une fichier avec pour chaque ligne gene,espece : ")

# Création d'un dico avec {"espece" : "gene"}
def GeneSymbol(filename) : 
    dico={}
    fichier = open(filename,'r')
    lines = fichier.readlines()
    for line in lines : 
        line= line.strip()
        if line : 
            li = line.split(',')
            species = li[1]
            gene =li[0]
            dico[species]=gene
    return dico 




# Extraire les uniprot_ID
def uniprot_ID (dico_espece) : 

    dico_uniprot={}

    for species, gene in dico_espece.items() :  
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+{species}&format=json&fields=accession"
        response = requests.get(url)

        if response.ok : 
            data= response.json()
            results =data.get("results", [])

            if results : 
                dico_uniprot[species] = results[0].get("primaryAccession") 
            else : 
                dico_uniprot[species] = None
            
        else : 
            print(f"Aucune information trouver pour {species} : {gene} ")
 
    return dico_uniprot



def QuickGO (uniprot_id) :

    dico_aspect={}
    # requete pour avoir le goID
    url = f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId={uniprot_id}&limit=200&page=1"
    response = requests.get(url, headers={"Accept": "application/json"})
    if response.ok:
        data = response.json()
        results = data.get("results", [])
        for res in results :
            aspect = res.get("goAspect")                
            goID = res.get("goId")

            # lien go term 
            

            #requete API pour avoir Go term 
            url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{goID}"
            response = requests.get(url, headers={"Accept": "application/json"})
            if response.ok:
                data = response.json()
                results = data.get("results", [])
                if results : 
                    go_name=results[0].get("name")
                    go_link = f"https://amigo.geneontology.org/amigo/term/{goID}"
                    go= f"name: {go_name}, link: {go_link}"
                if aspect not in dico_aspect:
                    dico_aspect[aspect]={}   # Créer un nouveau sous-dictionnaire si l'aspect n'existe pas encore
                dico_aspect[aspect][goID] = set()
                dico_aspect[aspect][goID].add(go)
                                     
    return  dico_aspect

def main_GO (dico_uniprot_id) : 

    dico ={}
    for species, uniprot in dico_uniprot_id.items() :
        if uniprot : 
            res = QuickGO(uniprot)
            dico[species]=res
            
    return dico 




#dico_GeneSymbol = GeneSymbol(liste)
#dico_uniprot_ID = uniprot_ID(dico_GeneSymbol)
#print(main_GO(dico_uniprot_ID))