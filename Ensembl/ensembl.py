# !/usr/bin/ python3
#-*- coding : utf-8 -*-

import requests
import sys
import os.path

## Fonction pour requête API et sortie en Json
def RequestJson (server, ext):
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    #Retour en cas d'erreur
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return(r.json())

## Fonction pour récupérer la division Ensembl en cas d'échec de la méthode principale
def RetrieveDivAlt(species) :

    # Retrouve les divisions ensembl 
    server = "https://rest.ensembl.org"
    ext = "/info/divisions?" 
 
    divisions = RequestJson(server,ext)[::-1] # Inverse pour chercher dans les espèce bactériennes d'abord

    # Retrouve la division associé à l'espèce
    j = 0
    div = None
    while not div and j<len(divisions) :
        
        ext = f"/info/species?division={divisions[j]}"     
        decoded = RequestJson(server,ext)

        for i in range(len(decoded["species"])) :
            if species in decoded["species"][i]["name"]:
                div = decoded["species"][i]["division"]
        j+=1
    return(div)

## Fonction pour récupérer la division Ensembl 
def RetrieveDiv (species) :
    # Requête des infos concernant l'espèce
    server = "https://rest.ensembl.org" # Serveur Ensembl pour requête API
    ext = f"/info/assembly/{species}?" # Requête d'un assemblage de l'espèce
    decoded = RequestJson(server, ext) # Résultat requête

    accession = decoded['assembly_accession'] # Récupération d'une accession d'un génome    

    # Requête des infos concernant l'accession
    ext = f"/info/genomes/assembly/{accession}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}) 

    # En cas d'erreur utiliser une méthode alternative
    if not r.ok:
        division = RetrieveDivAlt(species)
    else : 
        division = r.json()["division"] # Résultat API

    
    return(division[7:].lower())

## Fonction pour récupérer les informations relatives au symbole dans l'espèce spécifiée
def InfoGene (species_info) :
    
    current_species = species_info["species"] # Raccourci   

    # Serveur Ensembl pour requête API
    server = "https://rest.ensembl.org"

    
    # Requête vers l'espèce et symbole de gène courants
    ext = f"/lookup/symbol/{current_species}/{species_info["gene_symbol"]}?expand=1"     
    decoded = RequestJson(server,ext)

    # Construction du dictionnaire contenant les infos
    species_info["gene_id"] = decoded["id"]
    species_info["division"] = RetrieveDiv(current_species)
    species_info["transcript_id"] = list()
    species_info["prot_id"] = list()

    # Ajout des id de protéines, transcrits traduits et non traduits
    for i in decoded["Transcript"] :
            if "Translation" in i :
                species_info["prot_id"].append(i["Translation"]["id"])
                species_info["transcript_id"].append(i["id"])
            else : # En cas de non traduction ajout d'une mention au transcrit
                species_info["transcript_id"].append(i["id"] + "(Not translated)")

    div = species_info["division"]
    if div == "vertebrates":
        div = "www"

    species_info["gene_browser"] = f"https://{div}.ensembl.org/{current_species.capitalize()}/Location/View?db=core;g={species_info["gene_id"]}"
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

                print(InfoGene(species_info))
    else :
        print("Erreur : Veuillez donnez un nom de fichier accessible comme seul argument")