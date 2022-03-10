import re
import csv
import sys
from pathlib import Path
import requests
import json
import os

global exactMatchDic
exactMatchDic = dict()

headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}

def falcon_UMLS_CUI_function():
    value = global_dic["value"]
    output = ""
    url = 'https://labs.tib.eu/sdm/biofalcon/api?mode=short'
    text = str(value).replace("_"," ")
    payload = '{"text":"'+text+'"}'
    r = requests.post(url, data=payload.encode('utf-8'), headers=headers)
    if r.status_code == 200:
        response=r.json()
        if len(response['entities'][1])>0:
            return response['entities'][1][0]
        else:
            return ""
    else:
        return ""

def replaceExactMatch(value):                       
    if value != "":
        replacedValue = exactMatchDic[value]
    else:
        replacedValue = "" 
    return(replacedValue)

def dictionaryCreation():
    directory = Path(os.path.abspath(os.path.join(os.getcwd(), os.path.dirname(__file__)))).parent.absolute()
    with open(str(directory)+"/Sources/label_cui_dictionary.csv",'r') as data:
        for row in csv.DictReader(data):
            exactMatchDic.update({row['SampleOriginLabel']:row['CUI']}) 
dictionaryCreation()

# returns a string in lower case
def tolower(value):
    return value.lower()

# returns the concatenation of two strings
def concat2(value1,value2):
    if bool(value1) and bool(value2):
        result = str(str(value1)+str(value2))
    else:
        result = ""  
    return(result)

# return a string in upper case
def toupper(value):
    return value.upper()


# return a string in title case
def totitle(value):
    return value.title()


# return a string after removing leading and trailing whitespaces
def trim(value):
    return value.strip()


# return a string without s2
def chomp(value, toremove):
    return value.replace(toremove, '')


#return the substring (index2 can be null, index2 can be negative value)
def substring(value, index1, index2):
    if index2 is None:
        return value[int(index1):]
    else:
        return value[int(index1):int(index2)]


#replace value2 by value3
def replaceValue(value, value2, value3):
    return value.replace(value2, value3)


#returns the first appearance of the regex in value
def match(value, regex):
    return re.match(regex, value)[0]

def variantIdentifier(column1, column2,prefix):
    value = ""
    if (str(column1) != "nan"):
        value = re.sub('_.*','',str(column2))+"_"+str(column1).replace("c.","").replace(">", "~")
        value = prefix+value
    return value

def execute_function(row,dic):
    if "tolower" in dic["function"].lower():
        return tolower(row[dic["func_par"]["value"]])
    elif "toupper" in dic["function"]:
        return toupper(row[dic["func_par"]["value"]])
    elif "totitle" in dic["function"]:
        return totitle(row[dic["func_par"]["value"]])
    elif "trim" in dic["function"]:
        return trim(row[dic["func_par"]["value"]])
    elif "chomp" in dic["function"]:
        return chomp(row[dic["func_par"]["value"]],dic["func_par"]["toremove"])
    elif "substring" in dic["function"]:
        if "index2" in dic["func_par"].keys():
            return substring(row[dic["func_par"]["value"]],dic["func_par"]["index1"],dic["func_par"]["index2"])
        else:
            return substring(row[dic["func_par"]["value"]],dic["func_par"]["index1"],None)
    elif "replaceValue" in dic["function"]:
        return replaceValue(row[dic["func_par"]["value"]],dic["func_par"]["value2"],dic["func_par"]["value3"])
    elif "match" in dic["function"]:
        return match(dic["func_par"]["regex"],row[dic["func_par"]["value"]])
    elif "variantIdentifier" in dic["function"]:
        return variantIdentifier(row[dic["func_par"]["column1"]],row[dic["func_par"]["column2"]],dic["func_par"]["prefix"])
    elif "falcon_UMLS_CUI_function" in dic["function"]:
        return falcon_UMLS_CUI_function(row[dic["func_par"]["value"]])
    elif "replaceExactMatch" in dic["function"]:
        return replaceExactMatch(row[dic["func_par"]["value"]])
    else:
        print("Invalid function")
        print("Aborting...")
        sys.exit(1)

def execute_function_mysql(row,header,dic):
    if "tolower" in dic["function"]:
        return tolower(row[header.index(dic["func_par"]["value"])])
    elif "toupper" in dic["function"]:
        return toupper(row[header.index(dic["func_par"]["value"])])
    elif "totitle" in dic["function"]:
        return totitle(row[header.index(dic["func_par"]["value"])])
    elif "trim" in dic["function"]:
        return trim(row[header.index(dic["func_par"]["value"])])
    elif "chomp" in dic["function"]:
        return chomp(row[header.index(dic["func_par"]["value"])],dic["func_par"]["toremove"])
    elif "substring" in dic["function"]:
        if "index2" in dic["func_par"].keys():
            return substring(row[header.index(dic["func_par"]["value"])],dic["func_par"]["index1"],dic["func_par"]["index2"])
        else:
            return substring(row[header.index(dic["func_par"]["value"])],dic["func_par"]["index1"],None)
    elif "replaceValue" in dic["function"]:
        return replaceValue(row[header.index(dic["func_par"]["value"])],dic["func_par"]["value2"],dic["func_par"]["value3"])
    elif "match" in dic["function"]:
        return match(dic["func_par"]["regex"],row[header.index(dic["func_par"]["value"])])
    elif "variantIdentifier" in dic["function"]:
        return variantIdentifier(row[header.index(dic["func_par"]["column1"])],row[header.index(dic["func_par"]["column2"])],dic["func_par"]["prefix"])
    elif "falcon_UMLS_CUI_function" in dic["function"]:
        return falcon_UMLS_CUI_function(row[header.index(dic["func_par"]["value"])])
    elif "replaceExactMatch" in dic["function"]:
        return replaceExactMatch(row[header.index(dic["func_par"]["value"])])
    else:
        print("Invalid function")
        print("Aborting...")
        sys.exit(1)

def create_dictionary(triple_map):
    dic = {}
    inputs = []
    for tp in triple_map.predicate_object_maps_list:
        if "#" in tp.predicate_map.value:
            key = tp.predicate_map.value.split("#")[1]
            tp_type = tp.predicate_map.mapping_type
        elif "/" in tp.predicate_map.value:
            key = tp.predicate_map.value.split("/")[len(tp.predicate_map.value.split("/"))-1]
            tp_type = tp.predicate_map.mapping_type
        if "constant" in tp.object_map.mapping_type:
            value = tp.object_map.value
            tp_type = tp.object_map.mapping_type
        elif "#" in tp.object_map.value:
            value = tp.object_map.value.split("#")[1]
            tp_type = tp.object_map.mapping_type
        elif "/" in tp.object_map.value:
            value = tp.object_map.value.split("/")[len(tp.object_map.value.split("/"))-1]
            tp_type = tp.object_map.mapping_type
        else:
            value = tp.object_map.value
            tp_type = tp.object_map.mapping_type

        dic.update({key : value})
        if (key != "executes") and ([value,tp_type] not in inputs):
            inputs.append([value,tp_type])

    dic["inputs"] = inputs
    return dic
