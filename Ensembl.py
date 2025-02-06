import requests
import sys
import json
import os

def extract_division():
  """
Récupère les divisions taxonomiques d'Ensembl et leurs espèces associées.
    
Returns:
  dict: Un dictionnaire associant chaque espèce à sa division taxonomique.
  """

  print("Extraction des divisions d'Ensembl en cours...")

  server = "https://rest.ensembl.org"
  ext = "/info/divisions?"
 
  r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
  if not r.ok:
    r.raise_for_status()
    sys.exit()
 
  Divisions = r.json()

  dic_div = {}

  for division in Divisions:
    ext_species = "/info/genomes/division/" + division + "?"
    
    r = requests.get(server+ext_species, headers={ "Content-Type" : "application/json"})
  
    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    for x in decoded:
      name = x["name"]
      division = x["division"]
      if not name in dic_div:
        dic_div[name] = division

  return dic_div

def ensembl_extract(dico):
  """
Extrait les informations d'Ensembl pour chaque espèce et gène du dictionnaire fourni.
  """

  dic_div = extract_division()
  print("Divisions d'Ensembl extraites")
  results = []

  server = "https://rest.ensembl.org"

  for espece in dico : 
    ext = "/lookup/symbol/" + espece + "/" + dico[espece] + "?expand=1"
    r = requests.get(server + ext, headers={ "Content-Type" : "application/json"})
  
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    else:
      print(f"Requête effectuée pour l'espèce {espece}")
      
    decoded = r.json()

    gene_symbol = decoded["display_name"]
    species = decoded["species"]

    gene_access_number = decoded["id"]

    data = decoded["Transcript"]

    rna_access_number = []
    protein_access_number = []

    for transcript in data:
      rna_access_number.append(transcript["id"])
      if "Translation" in transcript:  # Vérifier si une clé "Translation" existe
        ensp = transcript["Translation"].get("id")      # Récupérer l'ID (ENSP)
        if ensp:  # Vérifier que les deux valeurs existent
          protein_access_number.append(ensp)
    
    prot_link = "/Transcript/ProteinSummary?p="
    if dic_div[espece.lower()] == "EnsemblVertebrates" :
      div = "www."
    elif dic_div[espece.lower()] == "EnsemblPlants":
      div = "plants."
    elif dic_div[espece.lower()] == "EnsemblProtists":
      div = "protists."
    elif dic_div[espece.lower()] == "EnsemblMetazoa":
      div = "metazoa."
    elif dic_div[espece.lower()] == "EnsemblFungi":
      div = "fungi."
    else:
      div = "bacteria."
      for x in protein_access_number : 
        prot_link = "/Transcript/ProteinSummary_" + x + "?db=core;p="
    
    genome_browser_link = "https://" + div + "ensembl.org/" + species + "/Location/View?db=core;g=" + gene_access_number
    gene_access_number_link = "https://" + div + "ensembl.org/" + species + "/Gene/Summary?db=core;g=" + gene_access_number
    
    if type(rna_access_number) == list:
      rna_access_number_link = []
      for x in rna_access_number :
        rna_access_number_link.append("https://" + div + "ensembl.org/" + species + "/Transcript/Summary?t=" + x)
    else:
      rna_access_number_link = "https://" + div + "ensembl.org/" + species + "/Transcript/Summary?t=" + rna_access_number
    
    if type(protein_access_number) == list:
      protein_access_number_link = []
      for x in protein_access_number :
        protein_access_number_link.append("https://" + div + "ensembl.org/" + species + prot_link + x)
    else:
      protein_access_number_link = "https://" + div + "ensembl.org/" + species + prot_link + protein_access_number

    results.append({
      "gene_symbol": gene_symbol,
      "species": species,
      "gene_access_number": gene_access_number,
      "rna_access_number": rna_access_number,
      "protein_access_number": protein_access_number,
      "genome_browser_link": genome_browser_link,
      "gene_access_number_link": gene_access_number_link,
      "rna_access_number_link": rna_access_number_link,
      "protein_access_number_link": protein_access_number_link
    })
      
  return results

if __name__ == '__main__':
  print("Ce fichier est destiné à être importé comme un module.")