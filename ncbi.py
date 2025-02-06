from Bio import Entrez
from Bio import SeqIO
import re
import os

def get_nested_value(data:dict, keys:list, default=None):
  """
Accède aux valeurs imbriquées dans un dictionnaire en suivant une liste de clés.

Args:
  data (dict): Dictionnaire à parcourir.
  keys (list): Liste des clés à suivre pour trouver la valeur.
  default: Valeur par défaut à retourner si une clé est absente.

Returns:
  any: Valeur trouvée ou valeur par défaut.
  """
  for key in keys:
    if isinstance(data, dict) and key in data:
      data = data[key]
    else:
      return default
  return data

def get_gene_name(record:dict):
  """
Permet de récupérer le nom du gène depuis différents champs.

Args:
  record (dict): Enregistrement NCBI (résultat d'une requête).

Returns:
  str: Nom du gène ou "Not available" si aucun nom n'est trouvé.
  """
  # Chemin 1 : Nom dans Gene-ref_formal-name
  gene_name = get_nested_value(record, [
      "Entrezgene_gene", "Gene-ref", "Gene-ref_formal-name", "Gene-nomenclature", "Gene-nomenclature_name"
  ])
  if gene_name:
    return gene_name

  # Chemin 2 : Description dans Gene-ref_desc
  gene_desc = get_nested_value(record, ["Entrezgene_gene", "Gene-ref", "Gene-ref_desc"])
  if gene_desc:
    return gene_desc

  # Chemin 3 : Protéine associée
  prot_name = get_nested_value(record, ["Entrezgene_prot", "Prot-ref", "Prot-ref_name"])
  if prot_name:
    return prot_name[0] if isinstance(prot_name, list) else prot_name

  # Chemin 4 : RNA associé
  rna_name = get_nested_value(record, ["Entrezgene_rna", "RNA-ref_ext", "RNA-ref_ext_name"])
  if rna_name:
    return rna_name

  return "Not available"

def ncbi_extract(dico:dict, mail:str):
  """
Extrait les informations génomiques à partir de NCBI pour une liste d'espèces et de gènes.
    
Args:
  dico (dict): Dictionnaire contenant les espèces et gènes.
  mail (str): Adresse email pour l'accès à l'API NCBI.
    
Returns:
  list: Liste de dictionnaires contenant les résultats pour chaque espèce et gène.
  """
  Entrez.email = mail
  results = []

  # Vérification préalable des dossiers nécessaires
  if not os.path.exists("gbk") :
    os.mkdir("gbk")

  # Créer un fichier de résumé des résultats
  for espece in dico :
    gene = dico[espece]

    # Extraction du nom d'espèce
    regex = r"^[A-Za-z]*_{1}[a-z]*"
    match = re.search(regex, espece, re.MULTILINE)
    if match:
      espece2 = match.group()
    
    # Requête à la base de données NCBI Gene
    myterm = gene+"[Gene] AND " + espece2 + "[Orgn]"

    handle = Entrez.esearch(db="gene", term=myterm, retmax=100)
    records = Entrez.read(handle)
    identifiers = records["IdList"]

    if not identifiers:
      continue
    else:
      print(f"Requête effectuée pour l'espèce {espece2}")
    
      handle = Entrez.efetch(db="gene", id=identifiers[0], rettype="xml", retmode="text", retmax=100)
      records = Entrez.read(handle)

      # Liste pour stocker les résultats
      rna_access_number = []
      rna_pred_access_number = []
      prot_access_number = []
      prot_pred_access_number = []
      rna_link = []
      rna_pred_link = []
      prot_link = []
      prot_pred_link = []

      # Analyse des informations du gène
      for record in records:
        gene_symbol = get_nested_value(record, [
          "Entrezgene_gene", "Gene-ref", "Gene-ref_locus"
          ], default="Not available")

        gene_name = get_gene_name(record)

        gene_id = get_nested_value(record, [
              "Entrezgene_track-info", "Gene-track", "Gene-track_geneid"
          ], default="Not available")
      
        gene_link = "https://www.ncbi.nlm.nih.gov/gene/"+gene_id

        # Recherche des numéros d'accesssions des protéines pour la Bactérie (pas de transcrits trouvés dans NCBI Nucleotide)
        orga = record['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_orgname']['OrgName']['OrgName_div']
        if orga == "BCT":
          prot_access_number.append(record["Entrezgene_locus"][0]["Gene-commentary_products"][0]["Gene-commentary_accession"])

      # Requête à la base de données NCBI Nucleotide

      handle = Entrez.esearch(db="nucleotide", term=myterm, retmax=100)
      records = Entrez.read(handle)
      identifiers = records["IdList"]

      handle = Entrez.efetch(db="nucleotide", id=identifiers, rettype="gb", retmode="text", retmax=100)
      text = handle.read()
      genbankfile = "./gbk/" + gene + "_" + espece2 + "_RefSeq.gbk"
      with open(genbankfile, "w") as genbank_file:
        genbank_file.write(text)

      # Dictionnaire qui permettra de vérifier la présence de séquences codantes
      feat_type = {"CDS":0}

      for seq_record in SeqIO.parse(genbankfile, "genbank"):
        # Récupération des séquences d'ARN (numéro d'accession commençant par NM)
        if seq_record.name.startswith("NM_"):
          ncbit = seq_record.name
          rna_access_number.append(ncbit)

          for feat in seq_record.features:

            if feat.type in feat_type:
              feat_type[feat.type] += 1
            else:
              feat_type[feat.type] = 1

            if feat.type == "CDS":
              for qual in feat.qualifiers:
                if qual == "gene":
                  if gene in feat.qualifiers[qual]:
                    continue
                  else:
                    break
                elif qual == "protein_id":
                  ncbip = feat.qualifiers[qual]
                  prot_access_number.append(ncbip[0])

          if feat_type["CDS"] == 0:
            prot_access_number.append("No protein translated from this transcript")
        
          feat_type = {"CDS":0}

        elif seq_record.name.startswith("XM_"):
          ncbitp = seq_record.name
          rna_pred_access_number.append(ncbitp)

          for feat in seq_record.features:

            if feat.type in feat_type:
              feat_type[feat.type] += 1
            else:
              feat_type[feat.type] = 1

            if feat.type == "CDS":
              for qual in feat.qualifiers:
                if qual == "gene":
                  if gene in feat.qualifiers[qual]:
                    continue
                  else:
                    break
                elif qual == "protein_id":
                  ncbipp = feat.qualifiers[qual]
                  prot_pred_access_number.append(ncbipp[0])

          if feat_type["CDS"] == 0:
            prot_pred_access_number.append("No protein translated from this transcript")
        
          feat_type = {"CDS":0}

      if rna_access_number:
        for nm in rna_access_number:
          rna_link.append("https://www.ncbi.nlm.nih.gov/nuccore/"+nm)

      if rna_pred_access_number:
        for xm in rna_pred_access_number:
          rna_pred_link.append("https://www.ncbi.nlm.nih.gov/nuccore/"+xm)

      if prot_access_number:
        for np in prot_access_number:
          prot_link.append("https://www.ncbi.nlm.nih.gov/protein/" + np)

      if prot_pred_access_number:
        for xp in prot_pred_access_number:
          prot_pred_link.append("https://www.ncbi.nlm.nih.gov/protein/" + xp)

      results.append({
        "species": espece,
        "gene_symbol": gene_symbol,
        "gene_name": gene_name,
        "gene_id": gene_id,
        "gene_link": gene_link,
        "rna_id": rna_access_number if rna_access_number else "No data found",
        "prot_id": prot_access_number if prot_access_number else "No data found",
        "rna_pred_id": rna_pred_access_number if rna_pred_access_number else "No data found",
        "prot_pred_id": prot_pred_access_number if prot_pred_access_number else "No data found",
        "rna_link": rna_link if rna_link else None,
        "prot_link": prot_link if prot_link else None,
        "rna_pred_link": rna_pred_link if rna_pred_link else None,
        "prot_pred_link": prot_pred_link if prot_pred_link else None
      })

  return results

if __name__ == '__main__':
  print("Ce fichier est destiné à être importé comme un module.")