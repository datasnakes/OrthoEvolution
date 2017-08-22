#!usr/bin/Python
def _parse_json(json_obj):
    """Edited from Source:  https://github.com/okeeffdp/json_to_newick"""
    try:
        # Test is the key 'name' in the current level of the dictionary.
        newick = json_obj['name']
    except KeyError:
        # Catch no 'name' trees and start the newick string with empty qoutes
        newick = ''

    if 'branch_length' in json_obj:
        newick = newick + ':' + str(json_obj['branch_length'])
    # If there are 'children'
    if 'children' in json_obj:
        # Initialise a list to contain the daughter info
        info = []
        # For each child, treat it as a new dictionary object
        for child in json_obj['children']:
            # parse the newick string straight into it the list
            info.append(_parse_json(child))
        # join all the daughter info together with a comma
        info = ','.join(info)
        # Concatenate all the children together at the start of the parent
        # newick string
        newick = '(' + info + ')' + newick
    newick = newick + ';'
    return newick

# if __name__ == '__main__':
# 	n = '(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F'
#
# 	json = {
# 		"name" : "F",
# 		"children": [
# 			{"name": "A", "branch_length": 0.1},
# 			{"name": "B", "branch_length": 0.2},
# 			{"name": "E","branch_length": 0.5,
# 			"children": [
# 				{"name": "C", "branch_length": 0.3},
# 				{"name": "D", "branch_length": 0.4}
# 				]
# 			}
# 		]
# 	}
# 	print(json_to_newick(json))
#
# 	json2 = {
# 		# "name" : "F",
# 		"children": [
# 			{#"name": "A",
# 			"branch_length": 0.1
# 			},
# 			{#"name": "B",
# 			"branch_length": 0.2
# 			},
# 			{#"name": "E",
# 			"branch_length": 0.5,
# 			"children": [
# 				{#"name": "C",
# 				"branch_length": 0.3
# 				},
# 				{#"name": "D",
# 				"branch_length": 0.4
# 				}
# 				]
# 			}
# 		]
# 	}
# 	print(json_to_newick(json2))
#
# 	json3 = {
# 		# "name" : "F",
# 		"children": [
# 			{#"name": "A",
# 			# "branch_length": 0.1
# 			},
# 			{#"name": "B",
# 			# "branch_length": 0.2
# 			},
# 			{#"name": "E",
# 			# "branch_length": 0.5,
# 			"children": [
# 				{#"name": "C",
# 				# "branch_length": 0.3
# 				},
# 				{#"name": "D",
# 				# "branch_length": 0.4
# 				}
# 				]
# 			}
# 		]
# 	}
# 	print(json_to_newick(json3))
#
# 	expected = '(((Drosophila_melanogaster:1,Caenorhabditis_elegans:1):1,((Ciona_intestinalis:1,Ciona_savignyi:1):1,(((((((((Taeniopygia_guttata:1,Ficedula_albicollis:1):1,((Meleagris_gallopavo:1,Gallus_gallus:1):1,Anas_platyrhynchos:1):1):1,Pelodiscus_sinensis:1):1,Anolis_carolinensis:1):1,((((((Procavia_capensis:1,Loxodonta_africana:1):1,Echinops_telfairi:1):1,(Choloepus_hoffmanni:1,Dasypus_novemcinctus:1):1):1,((((Oryctolagus_cuniculus:1,Ochotona_princeps:1):1,((((Mus_musculus:1,Rattus_norvegicus:1):1,Dipodomys_ordii:1):1,Ictidomys_tridecemlineatus:1):1,Cavia_porcellus:1):1):1,(((Microcebus_murinus:1,Otolemur_garnettii:1):1,(((((Papio_anubis:1,Macaca_mulatta:1):1,Chlorocebus_sabaeus:1):1,((((Pan_troglodytes:1,Homo_sapiens:1):1,Gorilla_gorilla:1):1,Pongo_abelii:1):1,Nomascus_leucogenys:1):1):1,Callithrix_jacchus:1):1,Tarsius_syrichta:1):1):1,Tupaia_belangeri:1):1):1,((Sorex_araneus:1,Erinaceus_europaeus:1):1,(((Pteropus_vampyrus:1,Myotis_lucifugus:1):1,((((Mustela_putorius_furo:1,Ailuropoda_melanoleuca:1):1,Canis_familiaris:1):1,Felis_catus:1):1,Equus_caballus:1):1):1,((((Bos_taurus:1,Ovis_aries:1):1,Tursiops_truncatus:1):1,Vicugna_pacos:1):1,Sus_scrofa:1):1):1):1):1):1,((Macropus_eugenii:1,Sarcophilus_harrisii:1):1,Monodelphis_domestica:1):1):1,Ornithorhynchus_anatinus:1):1):1,Xenopus_tropicalis:1):1,Latimeria_chalumnae:1):1,(((Danio_rerio:1,Astyanax_mexicanus:1):1,(((Tetraodon_nigroviridis:1,Takifugu_rubripes:1):1,((((Poecilia_formosa:1,Xiphophorus_maculatus:1):1,Oryzias_latipes:1):1,Gasterosteus_aculeatus:1):1,Oreochromis_niloticus:1):1):1,Gadus_morhua:1):1):1,Lepisosteus_oculatus:1):1):1,Petromyzon_marinus:1):1):1):1,Saccharomyces_cerevisiae:1);'
# 	mega_json = {"children":[{"name":"","children":[{"name":"","children":[{"name":"Drosophila_melanogaster","branch_length":1},{"name":"Caenorhabditis_elegans","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"Ciona_intestinalis","branch_length":1},{"name":"Ciona_savignyi","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Taeniopygia_guttata","branch_length":1},{"name":"Ficedula_albicollis","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"Meleagris_gallopavo","branch_length":1},{"name":"Gallus_gallus","branch_length":1}],"branch_length":1},{"name":"Anas_platyrhynchos","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Pelodiscus_sinensis","branch_length":1}],"branch_length":1},{"name":"Anolis_carolinensis","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Procavia_capensis","branch_length":1},{"name":"Loxodonta_africana","branch_length":1}],"branch_length":1},{"name":"Echinops_telfairi","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"Choloepus_hoffmanni","branch_length":1},{"name":"Dasypus_novemcinctus","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Oryctolagus_cuniculus","branch_length":1},{"name":"Ochotona_princeps","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Mus_musculus","branch_length":1},{"name":"Rattus_norvegicus","branch_length":1}],"branch_length":1},{"name":"Dipodomys_ordii","branch_length":1}],"branch_length":1},{"name":"Ictidomys_tridecemlineatus","branch_length":1}],"branch_length":1},{"name":"Cavia_porcellus","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Microcebus_murinus","branch_length":1},{"name":"Otolemur_garnettii","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Papio_anubis","branch_length":1},{"name":"Macaca_mulatta","branch_length":1}],"branch_length":1},{"name":"Chlorocebus_sabaeus","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Pan_troglodytes","branch_length":1},{"name":"Homo_sapiens","branch_length":1}],"branch_length":1},{"name":"Gorilla_gorilla","branch_length":1}],"branch_length":1},{"name":"Pongo_abelii","branch_length":1}],"branch_length":1},{"name":"Nomascus_leucogenys","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Callithrix_jacchus","branch_length":1}],"branch_length":1},{"name":"Tarsius_syrichta","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Tupaia_belangeri","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"Sorex_araneus","branch_length":1},{"name":"Erinaceus_europaeus","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Pteropus_vampyrus","branch_length":1},{"name":"Myotis_lucifugus","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Mustela_putorius_furo","branch_length":1},{"name":"Ailuropoda_melanoleuca","branch_length":1}],"branch_length":1},{"name":"Canis_familiaris","branch_length":1}],"branch_length":1},{"name":"Felis_catus","branch_length":1}],"branch_length":1},{"name":"Equus_caballus","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Bos_taurus","branch_length":1},{"name":"Ovis_aries","branch_length":1}],"branch_length":1},{"name":"Tursiops_truncatus","branch_length":1}],"branch_length":1},{"name":"Vicugna_pacos","branch_length":1}],"branch_length":1},{"name":"Sus_scrofa","branch_length":1}],"branch_length":1}],"branch_length":1}],"branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"Macropus_eugenii","branch_length":1},{"name":"Sarcophilus_harrisii","branch_length":1}],"branch_length":1},{"name":"Monodelphis_domestica","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Ornithorhynchus_anatinus","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Xenopus_tropicalis","branch_length":1}],"branch_length":1},{"name":"Latimeria_chalumnae","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Danio_rerio","branch_length":1},{"name":"Astyanax_mexicanus","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Tetraodon_nigroviridis","branch_length":1},{"name":"Takifugu_rubripes","branch_length":1}],"branch_length":1},{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"","children":[{"name":"Poecilia_formosa","branch_length":1},{"name":"Xiphophorus_maculatus","branch_length":1}],"branch_length":1},{"name":"Oryzias_latipes","branch_length":1}],"branch_length":1},{"name":"Gasterosteus_aculeatus","branch_length":1}],"branch_length":1},{"name":"Oreochromis_niloticus","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Gadus_morhua","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Lepisosteus_oculatus","branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Petromyzon_marinus","branch_length":1}],"branch_length":1}],"branch_length":1}],"branch_length":1},{"name":"Saccharomyces_cerevisiae","branch_length":1}],"name":""}
#
# 	observed = json_to_newick(mega_json)
#
# 	print(observed)
# 	print(observed == expected)
