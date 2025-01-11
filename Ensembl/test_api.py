# !/usr/bin/ python3
#-*- coding : utf-8 -*- 

import requests, sys, json
 
server = "https://rest.ensembl.org"
ext = "/info/species?Diviision=EnsemblProtist"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
test_file = open("test_species.txt","a")
test_file.write(json.dumps(decoded,indent=2))
test_file.close()